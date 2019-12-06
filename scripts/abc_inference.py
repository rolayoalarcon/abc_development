import argparse
import json
import timeit
import pandas as pd
import pybedtools
from pyliftover import LiftOver
from joblib import Parallel, delayed
from genomic_functions import activity_of_element, coordinate_conversion, format_bedstyle, rescale_bin, get_contacts, promoter_coord, close_element, selfBin_contact, abc_score, chromosome_conversion


def main():
# READING ARGUMENTS
    parser = argparse.ArgumentParser(description='Arguments for Activity from Contacts')

    parser.add_argument('--enhancers', type=str, help='Enhancer BedFile with the activity reported')
    parser.add_argument('--tss', type=str, help='TSS BedFile')
    parser.add_argument('--promoters', type=str, help='activity for promoter regions')

    parser.add_argument('--hic', type=str, help='Hi-C regularized counts', 
                        default='../external_dataset/K562_filteredRegularized_contactCount.tsv')
    parser.add_argument('--bincoord', type=str, help='Coordinates for bins',
                        default='../external_dataset/K562_filteredBins.bed')
    parser.add_argument('--chain', type=str, help='Chain file for coordinate liftover',
                        default='../external_dataset/hg38ToHg19.over.chain')
    
    parser.add_argument('--chromap', type=str, help='Chromosome mappping file',
                        default='../external_dataset/GRCh37_UCSC2ensembl.txt')
    
    parser.add_argument('-p', type=int, help='Cores to use during processing', default=1)
    parser.add_argument('--scaler', type=int, help='Values to multiply for re-scaling', default=100)
    parser.add_argument('--closer', type=int, help='Cutoff for enhancer vecinity', default=5_000_000)
    parser.add_argument('--gamma', type=int, help='Gamma powerlaw parameter', default=-0.7764771175681618)
    parser.add_argument('--scaleparam', type=int, help='Scale powerlaw parameter', default=10.787505424121239)
    parser.add_argument('--mindist', type=int, help='Minimum distance for powerlaw', default=1_000_000)
    parser.add_argument('--promlength', type=int, help='Promoter length', default=500)
    parser.add_argument('--cutoff', type=int, help='Cutoff probability', default=0)
    parser.add_argument('--outfile', type=str, help='Output file name')

    args = parser.parse_args()

    # ASSIGNING ARGUMENTS TO VARIABLES
    enhancer_bedfile = args.enhancers
    tss_bedfile = args.tss
    promoters_bedfile = args.promoters
    hic_file = args.hic

    num_cores = args.p
    chromosome_mapping = args.chromap
    coord_conversion = args.chain
    filtered_coords = args.bincoord
    SCALER_VALUE = args.scaler
    CLOSER_CUTOFF = args.closer
    SCALE_PL = args.scaleparam
    GAMMA_PL = args.gamma
    DISTMIN_PL = args.mindist
    PROMOTER_LEN = args.promlength
    CUTOFF_POSITIVE = args.cutoff
    output = args.outfile

    '''
    Reading file from input
    '''
    print('Reading files...')
    tss_df = pd.read_csv(tss_bedfile, sep='\t', header=None, names=['chr', 'start', 'end', 'gene_id', 'score', 'strand'])
    promoters_df = pd.read_csv(promoters_bedfile, sep='\t')
    enhancer_df = pd.read_csv(enhancer_bedfile, sep='\t')

    # For some reason, the indexing is faster when read from csv.
    hic_counts = pd.read_csv(hic_file, sep='\t')

    filtered_bins = pybedtools.BedTool(filtered_coords)

    lift_file = LiftOver(coord_conversion)
    chromap_file = pd.read_csv(chromosome_mapping, sep='\t', header=None, names=['chromosome', 'ensembl_chr'],
                               index_col='chromosome')
    try:
        enhancer_df = enhancer_process(enhancer_info=enhancer_df, filtered_bincoord=filtered_bins, 
                                       coordinate_liftover=lift_file, chromosome_relation=chromap_file)

        tss_df = tss_process(tss_info=tss_df, filtered_bincoord=filtered_bins, coordinate_liftover=lift_file,
                             chromosome_relation=chromap_file)

        tss_enhancer_intersected = intersect_elements(tss_intersect=tss_df, enhancer_intersect=enhancer_df,
                                                      closer_value=CLOSER_CUTOFF)

        rescaled_data = rescale_rows(processed_df=tss_df, s_value=SCALER_VALUE, regularized_counts=hic_counts,
                                     num_p=num_cores)

        denom_dict, tss_enhancer_newinf = calculate_denominator(enhancer_tss_info=tss_enhancer_intersected,
                                                                promoter_info=promoters_df,
                                                                scaled_counts=rescaled_data, gamma_powerlaw=GAMMA_PL,
                                                                scale_powerlaw=SCALE_PL, s_value=SCALER_VALUE,
                                                                distance_min=DISTMIN_PL, promoter_length=PROMOTER_LEN,
                                                                num_p=num_cores)

        calculate_abc(enhancer_tss_info=tss_enhancer_newinf, denominator_values=denom_dict,
                      num_p=num_cores, positive_cutoff=CUTOFF_POSITIVE, output_name=output)

    finally:
        pybedtools.cleanup()

def enhancer_process(enhancer_info, filtered_bincoord, coordinate_liftover, chromosome_relation):

    '''

    PROCESSING ENHANCERS

    1. Calculate activity
    2. Convert chromosomes and coordinates
    3. Instersect with bins

    '''

    score_cols = ['CPM_atac_0', 'CPM_atac_1', 'CPM_h3k27ac_0']
    print('Processing enhancers...')
    enhancer_info['activity'] = enhancer_info.apply(activity_of_element, axis=1,
                                                    args=(score_cols),
                                                    result_type='expand')

    # The coordinates for the enhancers are in HG38 but the coordinates for the bin
    # are in b37. Also, we want to find the midpoint of each enhancer to guarantee
    # that only one bin intersects with each enhancer.
    enhancer_info['midpoint_start'] = enhancer_info['start'] + (abs(enhancer_info['end'] - enhancer_info['start']) // 2)
    enhancer_info['midpoint_end'] = enhancer_info['midpoint_start'] + 1
    enhancer_info.rename(columns={'start': 'enhancer_start', 'end': 'enhancer_end'}, inplace=True)

    enhancer_info[['start_conv', 'end_conv']] = enhancer_info.apply(coordinate_conversion,
                                                                    args=('chr', 'midpoint_start',
                                                                          'midpoint_end',
                                                                          coordinate_liftover.convert_coordinate),
                                                                    axis=1, result_type='expand')

    enhancer_info['chromosome_conv'] = chromosome_conversion(enhancer_info.copy(), chromosome_relation)

    enhancer_imp = format_bedstyle(enhancer_info, name_col='enhancer_id', score_col='activity',
                                   start_col='start_conv', end_col='end_conv')

    enhancer_bed = pybedtools.BedTool.from_dataframe(enhancer_imp)

    enhancer_bin_intersect = enhancer_bed.intersect(filtered_bincoord, wo=True)
    enh_bins_intersect_df = enhancer_bin_intersect.to_dataframe(names=['enh_chrom', 'enhConv_start', 'enhConv_end',
                                                                       'enhancer_id', 'activity', 'strand',
                                                                       'enhBin_chrom', 'enhBin_start', 'enhBin_end',
                                                                       'enhBin_id', 'dot', 'strand2',
                                                                       'overlap'], header=None)

    enhancer_set = enhancer_info.loc[enhancer_info['enhancer_id'].isin(enh_bins_intersect_df['enhancer_id'].tolist())]
    enhancer_set.set_index('enhancer_id', inplace=True)
    enh_bins_intersect_df.set_index('enhancer_id', inplace=True)
    enhancer_processed = enhancer_set.join(enh_bins_intersect_df[['enhBin_start', 'enhBin_end', 'enhBin_id', 'enhConv_start', 'enhConv_end']])
    enhancer_processed.reset_index(inplace=True)

    print(enhancer_processed.head())
    return enhancer_processed


def tss_process(tss_info, filtered_bincoord, coordinate_liftover, chromosome_relation):
    '''

    PROCESSING PROMOTERS

    1. Convert coordinates
    2. Intersect with bins

    '''

    print('Processing Transcription Start Sites...')
    tss_info[['start_conv', 'end_conv']] = tss_info.apply(coordinate_conversion,
                                                           args=('chr', 'start', 'end',
                                                                 coordinate_liftover.convert_coordinate),
                                                           axis=1, result_type='expand')
    # print(tss_info.loc[tss_info['gene_id'] == 'ENSG00000100889:U.U.U', ['start_conv', 'end_conv']])

    tss_info['chromosome_conv'] = chromosome_conversion(tss_info.copy(), chromosome_relation)
    print(tss_info.columns)
    tss_imp = format_bedstyle(tss_info, name_col='gene_id', score_col='score',
                              start_col='start_conv', end_col='end_conv')

    tss_bed = pybedtools.BedTool.from_dataframe(tss_imp)

    tss_bin_intersect = tss_bed.intersect(filtered_bincoord, wo=True)
    tss_bins_intersect_df = tss_bin_intersect.to_dataframe(names=['tss_chrom', 'tssConv_start', 'tssConv_end',
                                                                   'gene_id', 'something', 'strand',
                                                                   'tssBin_chrom', 'tssBin_start', 'tssBin_end',
                                                                   'tssBin_id', 'something_2', 'tssBin_strand',
                                                                   'overlap'], header=None)

    tss_set = tss_info.loc[tss_info['gene_id'].isin(tss_bins_intersect_df['gene_id'].tolist())]
    tss_set.set_index('gene_id', inplace=True)
    tss_bins_intersect_df.set_index('gene_id', inplace=True)
    tss_processed = tss_set.join(tss_bins_intersect_df[['tssBin_start', 'tssBin_end', 'tssBin_id', 'tssConv_start', 'tssConv_end']])
    tss_processed.reset_index(inplace=True)

    print(tss_processed.head())
    #print('Number of bins: {}\nExpected minutes {}'.format(num_bins, (num_bins*6)/60))
    return tss_processed


def rescale_rows(processed_df, regularized_counts, s_value, num_p):
    print('Re-scaling rows...')
    print(str(num_p))
    intersected_bins = processed_df['tssBin_id'].drop_duplicates().tolist()
    start = timeit.default_timer()
    dict_list = Parallel(n_jobs=num_p)(delayed(rescale_bin)(bin_id=bin_number,
                                                            regularized_matrix=regularized_counts,
                                                            base_number=s_value) for bin_number in intersected_bins)
    elapsed = timeit.default_timer() - start
    print('Elapsed: {}'.format(elapsed))

    bins_rescaled = {key: value for d in dict_list for key, value in d.items()}

    with pd.HDFStore('scaled_rows.h5') as data_storage:
        for key, value in bins_rescaled.items():
            data_storage['bin_{}'.format(key)] = value

    return bins_rescaled


def intersect_elements(tss_intersect, enhancer_intersect, closer_value):
    '''

    INTERSECTING TSSs WITH ENHANCERS

    '''
    print('Intersecting TSSs with enhancers...')

    tss_bed_df = tss_intersect[['chr', 'start', 'end', 'gene_id', 'tssBin_id', 'strand']]

    enhancer_intersect['binID_activity'] = enhancer_intersect['enhBin_id'].astype(str) + ':' + enhancer_intersect['activity'].astype(str)
    enhancer_intersect['enhancer_midpoint'] = enhancer_intersect['enhancer_id'].astype(str) + ':' + enhancer_intersect['midpoint_start'].astype(str)
    enhancer_bed_df = enhancer_intersect[['chr', 'enhancer_start', 'enhancer_end',
                                          'enhancer_midpoint', 'binID_activity', 'strand']]

    tss_bed = pybedtools.BedTool.from_dataframe(tss_bed_df)
    enhancer_bed = pybedtools.BedTool.from_dataframe(enhancer_bed_df)

    tss_enh_window = tss_bed.window(enhancer_bed, w=closer_value)
    tss_enh_window_df = tss_enh_window.to_dataframe(names=['tss_chrom', 'tss_start', 'tss_end',
                                                           'gene_id', 'tssBin_id', 'strand',
                                                           'enh_chrom', 'enh_start', 'enh_end',
                                                           'enhancer_midpoint', 'binID_activity',
                                                           'enh_strand'], header=None)

    window_expanded = tss_enh_window_df.assign(midpoint=tss_enh_window_df['enhancer_midpoint'].apply(lambda x: int(x.split(':')[1])),
                                               enhancer_id=tss_enh_window_df['enhancer_midpoint'].apply(lambda x: x.split(':')[0]),
                                               enhBin_id=tss_enh_window_df['binID_activity'].apply(lambda x: int(x.split(':')[0])),
                                               activity=tss_enh_window_df['binID_activity'].apply(lambda x: float(x.split(':')[1])))

    results_window = window_expanded[['tss_chrom', 'tss_start', 'tss_end', 'gene_id', 'tssBin_id', 'strand',
                                      'enh_chrom', 'enh_start', 'enh_end', 'enhancer_id', 'midpoint', 'activity',
                                      'enhBin_id']]

    results_window['distance'] = abs(results_window['tss_start'] - results_window['midpoint'])

    print(results_window.head())
    return results_window


def calculate_denominator(enhancer_tss_info, promoter_info, scaled_counts, gamma_powerlaw, scale_powerlaw, s_value, distance_min,
                          promoter_length, num_p):

    '''

    CALCULATING THE DENOMINATOR FOR EACH TSS

    '''
    score_cols = ['CPM_atac_0', 'CPM_atac_1', 'CPM_h3k27ac_0']
    print('Calculating denominator for TSSs...')
    enhancer_tss_info['corrected_contacts'] = enhancer_tss_info.apply(get_contacts,
                                                                      args=(gamma_powerlaw, scale_powerlaw,
                                                                            scaled_counts, s_value, distance_min),
                                                                      axis=1)

    promoter_info['promoter_activity'] = promoter_info.apply(activity_of_element, axis=1,
                                                             args=(score_cols),
                                                             result_type='expand')


    enhancer_tss_info.set_index('gene_id', inplace=True)
    promoter_info.set_index('gene_id', inplace=True)

    enhancer_tss_info = enhancer_tss_info.join(promoter_info[['promoter_activity']])
    enhancer_tss_info.reset_index(inplace=True)

    gene_list = enhancer_tss_info['gene_id'].drop_duplicates().tolist()
    dict_list = Parallel(n_jobs=num_p)(delayed(close_element)(gene_id=gene, genomic_info=enhancer_tss_info,
                                                              gamma_self=gamma_powerlaw, scale_self=scale_powerlaw,
                                                              min_self=distance_min, scale_mult=s_value,
                                                              rescaled_name=scaled_counts) for gene in gene_list)

    correction_values = {key: value for d in dict_list for key, value in d.items()}

    enhancer_tss_info.to_csv('ehancer_tss_info.tsv', sep='\t')

    with open('correction_values.json', 'w') as data_storage:
        json.dump(correction_values, data_storage)

    print(enhancer_tss_info.head())
    return correction_values, enhancer_tss_info


def calculate_abc(enhancer_tss_info, denominator_values, num_p, positive_cutoff, output_name):

    print('Calculating ABC scores...')
    print('Using cutoff: {}'.format(positive_cutoff))
    gene_list = enhancer_tss_info['gene_id'].drop_duplicates().tolist()

    gene_enhancer_scores = Parallel(n_jobs=num_p)(delayed(abc_score)(gene_id=gene,
                                                                     genomic_info=enhancer_tss_info,
                                                                     correction_terms=denominator_values) for gene in gene_list)
    abc_concatenation = pd.concat(gene_enhancer_scores)
    abc_positives = abc_concatenation[abc_concatenation['ABC_score'] > positive_cutoff].copy()

    abc_positives['midpoint_end'] = abc_positives['midpoint'] + 1
    abc_arcs = abc_positives[['tss_chrom', 'tss_start', 'tss_end',
                              'enh_chrom', 'midpoint', 'midpoint_end', 'ABC_score']]

    abc_positives.to_csv(output_name, sep='\t', index=False)
    #abc_arcs.to_csv('abc_arcdata.tsv', sep='\t', index=False)

    print('Done.')


if __name__ == '__main__':
    main()
