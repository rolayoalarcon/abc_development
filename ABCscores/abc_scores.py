import argparse
import math as mt
import numpy as np
import pandas as pd
import pyBigWig as pbw
from pybedtools import BedTool
from pyliftover import LiftOver
from joblib import Parallel, delayed

# CONSTANTS FOR SCRIPT
CLOSER_CUTOFF = 5_000_000
POWER_CUTOFF = 1_000_000
SCALER_VALUE = 100
GAMMA_K562 = -0.7764771175681618
SCALE_K562 = 10.787505424121239
PROMOTER_LENGTH = 500
SCORE_CUTOFF = 0.5
CHROMOSOME_MAPPING = '../data/external/GRCh37_UCSC2ensembl.txt'
COORDINATE_LIFTOVER = '../data/external/hg38ToHg19.over.chain'
FILTERED_BINS_BEDFILE = '../data/external/K562_filteredBins.bed'


'''

FUNCTIONS START

'''


# This function calculates the activity of each enhancer by returning the
# the geometric mean of H3K27 and DHS reads at a given enhancer region.
def activity_of_element(row, chrom_name, start_name, end_name):
    dhs_reads = atac_bw.stats(row[chrom_name], row[start_name], row[end_name], type='sum')[0]

    h3k_reads = h3k27ac_bw.stats(row[chrom_name],
                                 row[start_name], row[end_name], type='sum')[0]

    # Geometric mean
    g_mean = mt.sqrt(dhs_reads * h3k_reads)

    return(g_mean)


# This function uses liftover to convert the coordinates. We take into account
# the name of the columns where the coordinates are.
def coordinate_conversion(row, start_name, end_name):

    start_conversion = hg38_hg19.convert_coordinate(row['chromosome'], row[start_name])
    end_conversion = hg38_hg19.convert_coordinate(row['chromosome'], row[end_name])

    if len(start_conversion) > 0:
        ch_hg19, start_hg19, _, _ = start_conversion[0]
    else:
        ch_hg19, start_hg19 = -1, -1

    if len(end_conversion) > 0:
        _, end_hg19, _, _ = end_conversion[0]
    else:
        end_hg19 = -1

    if ch_hg19 != -1:
        ch_b37 = chromosome_relation.loc[ch_hg19].values[0]
    else:
        ch_b37 = -1

    return(ch_b37, min(start_hg19, end_hg19), max(start_hg19, end_hg19))


# This funciton re-scales the rows of each bin and returns a df with the new
# values
def rescale_row(bin_id):
    bin_row = regularized_counts.loc[(regularized_counts['bin1_id'] == bin_id)
                                     | (regularized_counts['bin2_id'] == bin_id)].copy()

    scaled_row = rescaler(bin_row['regularized_count'].values)

    bin_row['scaled_count'] = scaled_row
    return({bin_id: bin_row})


# Resturn vector so that the max is 100.
def rescaler(numbers_list):
    maximum = max(numbers_list)
    rescaled_numbers = (numbers_list / maximum) * SCALER_VALUE
    return(rescaled_numbers)


# This function returns the re-scaled number of contacts between two bins
def get_contacts(row):
    involved_bins = [row['tssBin_id'], row['enhBin_id']]

    pseudocount = power_law(distance=row['distance'])
    bin_info = bins_rescaled[row['tssBin_id']]
    max_contacts = max(bin_info['regularized_count'].values)

    count = bin_info.loc[(bin_info['bin1_id'] == min(involved_bins))
                         & (bin_info['bin2_id'] == max(involved_bins)), 'scaled_count'].values
    if len(count) > 0:
        count_pseudocount = count[0] + ((pseudocount/max_contacts) * 100)
        if count_pseudocount > 100:
            count_pseudocount = 100
    else:
        count_pseudocount = (pseudocount/max_contacts) * 100

    return(count_pseudocount)


# This function calculates a pseudocount of the expected number of contacts
# given a certain distance.
def power_law(distance, gamma=GAMMA_K562, scale=SCALE_K562):
    if distance > POWER_CUTOFF:
        correction = scale*(distance**gamma)
    elif distance <= POWER_CUTOFF:
        correction = scale*(POWER_CUTOFF**gamma)

    return(correction)


# This function will calculate the denominator for each gene
def close_element(gene_id):
    gene_enhancer_info = tss_enh_close_df[tss_enh_close_df['Ensembl gene id'] == gene_id].copy()
    denominator = sum(gene_enhancer_info['activity'].values * gene_enhancer_info['corrected_conacts'].values)
    promoter_act = gene_enhancer_info['promoter_activity'].values[0]

    promoter_contacts = selfBin_contact(gene_enhancer_info['tssBin_id'].values[0])
    promoter_value = promoter_act * promoter_contacts

    denominator += promoter_value
    return({gene_id: denominator})


def promoter_coord(row):
    if row['strand'] == '+':
        second_coord = row['tss_start'] - PROMOTER_LENGTH
    elif row['strand'] == '-':
        second_coord = row['tss_start'] + PROMOTER_LENGTH
    return(min(row['tss_start'], second_coord), max(row['tss_start'], second_coord))


def selfBin_contact(bin_id):
    pseudocount = power_law(distance=0)
    bin_info = bins_rescaled[bin_id]
    max_contacts = max(bin_info['regularized_count'].values)
    count = bin_info.loc[(bin_info['bin1_id'] == bin_id)
                         & (bin_info['bin2_id'] == bin_id), 'scaled_count'].values

    if len(count) > 0:
        count_pseudocount = count[0] + (pseudocount/max_contacts) * 100
        if count_pseudocount > 100:
            count_pseudocount = 100
    else:
        count_pseudocount = (pseudocount/max_contacts) * 100
    return(count_pseudocount)



def abc_score(gene_id):
    gene_info = tss_enh_close_df[tss_enh_close_df['Ensembl gene id'] == gene_id].copy()
    activity = gene_info['activity'].values * gene_info['corrected_conacts'].values
    abc_scores = activity/correction_terms[gene_info['Ensembl gene id'].values[0]]

    gene_info['ABC_score'] = abc_scores
    return(gene_info)


'''

FUNCTIONS END

'''

# READING ARGUMENTS
parser = argparse.ArgumentParser(description='Arguments for Activity from Contacts')
parser.add_argument('--atac', type=str, help='ATAC BigWig')
parser.add_argument('--h3k27ac', type=str, help='H3K27ac BigWig')
parser.add_argument('--enhancers', type=str, help='Enhancer BedFile')
parser.add_argument('--tss', type=str, help='TSS BedFile')
parser.add_argument('--hic', type=str, help='Hi-C regularized counts')
parser.add_argument('-c', type=int, help='Cores to use during processing', default=1)

args = parser.parse_args()

# ASSIGNING ARGUMENTS TO VARIABLES
atac_file = args.atac
h3k27ac_file = args.h3k27ac
enhancer_bedfile = args.enhancers
tss_bedfile = args.tss
hic_file = args.hic
num_cores = args.c


'''

FIRST STEP WILL BE READING FILES

'''
print('Reading files...')

atac_bw = pbw.open(atac_file)
h3k27ac_bw = pbw.open(h3k27ac_file)
tss_info = pd.read_csv(tss_bedfile, sep='\t')
enhancer_info = pd.read_csv(enhancer_bedfile, sep='\t')
regularized_counts = pd.read_csv(hic_file, sep='\t')
filtered_bins = BedTool(FILTERED_BINS_BEDFILE)
hg38_hg19 = LiftOver(COORDINATE_LIFTOVER)
chromosome_relation = pd.read_csv(CHROMOSOME_MAPPING, sep='\t', header=None,
                                  names=['chromosome', 'ensembl_chr'],
                                  index_col='chromosome')

'''

PROCESSING ENHANCERS

1. Calculate activity
2. Convert chromosomes and coordinates
3. Instersect with bins

'''
print('Processing enhancers...')
enhancer_info['activity'] = enhancer_info.apply(activity_of_element, axis=1,
                                                args=('chromosome', 'start', 'end'),
                                                result_type='expand')

# The coordinates for the enhancers are in HG38 but the coordinates for the bin
# are in b37. Also, we want to find the midpoint of each enhancer to guarantee
# that only one bin intersects with each enhancer.

enhancer_info['midpoint_start'] = enhancer_info['start'] + (abs(enhancer_info['end'] - enhancer_info['start'])//2)
enhancer_info['midpoint_end'] = enhancer_info['midpoint_start'] + 1

enhancer_info.rename(columns={'start': 'enhancer_start',
                              'end': 'enhancer_end'}, inplace=True)

enhancer_info[['chromosome_ensembl', 'start_hg19', 'end_hg19']] = enhancer_info.apply(coordinate_conversion,
                                                                                      args=('midpoint_start', 'midpoint_end'),
                                                                                      axis=1, result_type='expand')

# Preparing dataframe for BedTool

enhancer_imp = enhancer_info[['chromosome_ensembl', 'start_hg19', 'end_hg19',
                              'Enhancer id', 'activity', 'strand']].copy()
enhancer_imp.replace(-1, np.nan, inplace=True)
enhancer_imp.dropna(inplace=True)
enhancer_imp['start_hg19'] = enhancer_imp['start_hg19'].astype(int)
enhancer_imp['end_hg19'] = enhancer_imp['end_hg19'].astype(int)

# Now we intersect enhancers with filtered bins
enhancers_bed = BedTool.from_dataframe(enhancer_imp)

enh_bins_intersect = enhancers_bed.intersect(filtered_bins, wo=True)
enh_bins_intersect_df = enh_bins_intersect.to_dataframe(names=['enh_chrom', 'enh_start', 'enh_end',
                                                               'Enhancer id', 'activity', 'strand',
                                                               'enhBin_chrom', 'enhBin_start', 'enhBin_end',
                                                               'enhBin_id', 'dot', 'strand2',
                                                               'overlap'], header=None)

# No enhancers should intersect with more than one bin
enhancer_list = enh_bins_intersect_df['Enhancer id'].tolist()
enhancer_nodup = enh_bins_intersect_df['Enhancer id'].drop_duplicates().tolist()
assert len(enhancer_list) == len(enhancer_nodup)

'''

PROCESSING PROMOTERS

1. Convert coordinates
2. Intersect with bins
'''
print('Processing promoters...')


tss_info[['chromosome_ensembl', 'start_hg19', 'end_hg19']] = tss_info.apply(coordinate_conversion,
                                                                            args=('start', 'end'),
                                                                            axis=1, result_type='expand')

tss_imp = tss_info[['chromosome_ensembl', 'start_hg19', 'end_hg19',
                    'Ensembl gene id', 'dot', 'strand']].copy()
tss_imp.replace(-1, np.nan, inplace=True)
tss_imp.dropna(inplace=True)
tss_imp['start_hg19'] = tss_imp['start_hg19'].astype(int)
tss_imp['end_hg19'] = tss_imp['end_hg19'].astype(int)

# Intersecting

tss_bed = BedTool.from_dataframe(tss_imp)
tss_bins_intersect = tss_bed.intersect(filtered_bins, wo=True)
tss_bins_intersect_df = tss_bins_intersect.to_dataframe(names=['tss_chrom', 'tss_start', 'tss_end',
                                                               'Ensembl gene id', 'something', 'strand',
                                                               'tssBin_chrom', 'tssBin_start', 'tssBin_end',
                                                               'tssBin_id', 'something_2', 'tssBin_strand',
                                                               'overlap'], header=None)

tss_list = tss_bins_intersect_df['Ensembl gene id'].tolist()
tss_nodup = tss_bins_intersect_df['Ensembl gene id'].drop_duplicates().tolist()
assert len(tss_list) == len(tss_nodup)


'''

PREPARE DATA ROWS
I propose that this section be paralelized, since it takes some time, even with
only one chromosome. We would only have to iterate over the bin IDs

'''
print('Re-scaling rows...')
bins_intersected = tss_bins_intersect_df['tssBin_id'].drop_duplicates().tolist()
dict_list = Parallel(n_jobs=num_cores)(delayed(rescale_row)(bin_number) for bin_number in bins_intersected)
bins_rescaled = {key: value for d in dict_list for key, value in d.items()}


'''

INTERSECTING TSSs WITH ENHANCERS

'''
print('Intersecting TSSs with enhancers...')
# First we will get the coordinates in hg38 back into the information. We also
# get information for the complete span of the enhancer.
enh_bins_intersect_df.set_index('Enhancer id', inplace=True)
enhancer_info.set_index('Enhancer id', inplace=True)

enhancer_bins_info = enh_bins_intersect_df.join(enhancer_info[['chromosome', 'enhancer_start', 'enhancer_end']])
enhancer_bins_info.reset_index(inplace=True)

# We do the same with TSSs
tss_info.set_index('Ensembl gene id', inplace=True)
tss_bins_intersect_df.set_index('Ensembl gene id', inplace=True)

tss_bins_info = tss_bins_intersect_df.join(tss_info[['chromosome', 'start', 'end']])
tss_bins_info.reset_index(inplace=True)

# Look for enhancers within 5Mb of each TSS.
tss_bed_df = tss_bins_info[['chromosome', 'start', 'end',
                            'Ensembl gene id', 'tssBin_id', 'strand']]

enh_bed_df = enhancer_bins_info[['chromosome', 'enhancer_start', 'enhancer_end',
                                 'enhBin_id', 'activity', 'strand']]

# Searching
tss_bed = BedTool.from_dataframe(tss_bed_df)
enh_bed = BedTool.from_dataframe(enh_bed_df)

tss_enh_close = tss_bed.window(enh_bed, w=CLOSER_CUTOFF)
tss_enh_close_df = tss_enh_close.to_dataframe(names=['tss_chrom', 'tss_start', 'tss_end',
                                                     'Ensembl gene id', 'tssBin_id', 'strand',
                                                     'enh_chrom', 'enh_start', 'enh_end',
                                                     'enhBin_id', 'activity', 'enh_strand'], header=None)

tss_enh_close_df['distance'] = abs(tss_enh_close_df['tss_start'] - (tss_enh_close_df['enh_end'] - 250))

'''

CALCULATING THE DENOMINATOR FOR EACH TSS

'''
print('Calculating denominator for TSSs')
# Number of contacts
tss_enh_close_df['corrected_conacts'] = tss_enh_close_df.apply(get_contacts, axis=1)
tss_enh_close_df[['promoter_start', 'promoter_end']] = tss_enh_close_df.apply(promoter_coord, axis=1,
                                                                              result_type='expand')

tss_enh_close_df['promoter_activity'] = tss_enh_close_df.apply(activity_of_element, axis=1,
                                                            args=('tss_chrom', 'tss_start', 'tss_end'),
                                                            result_type='expand')

# Calculate denominator with funciton
# This can also be Parallel
gene_list = tss_enh_close_df['Ensembl gene id'].drop_duplicates().tolist()
correction_list = Parallel(n_jobs=num_cores)(delayed(close_element)(gene_id) for gene_id in gene_list)
correction_terms = {key: value for d in correction_list for key, value in d.items()}


'''

CALCULATING ABC SCORE FOR EACH GENE-ENHANCER PAIR

'''
print('Calculating ABC scores...')
gene_enhancer_scores = Parallel(n_jobs=num_cores)(delayed(abc_score)(gene_id) for gene_id in gene_list)
abc_concatenation = pd.concat(gene_enhancer_scores)
abc_positives = abc_concatenation[abc_concatenation['ABC_score'] > SCORE_CUTOFF].copy()

abc_positives.to_csv('abc_test_results.txt', sep='\t')
print('Done.')
