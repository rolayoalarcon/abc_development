import os
import yaml
import pandas as pd
import argparse

NORMALIZATION_FACTOR = 1_000_000

# Chromosomes whose reads will be counted for CPM calculation
READING_CHROMOSOMES = ['chr{}'.format(i) for i in range(1, 23)]
READING_CHROMOSOMES.extend(['chrY', 'chrX'])

# columns to be set as index when joining
INDEX_COLS = ['chr', 'start', 'end', 'id', 'length', 'strand']


def norm_reads(cov_dir, stat_dir, prefix, assay, region, selected_chrom, norm_factor):
    cov_file = '{assay}_{pre}_{gr}.tsv'.format(assay=assay, pre=prefix, gr=region)
    stat_file = '{prefix}.stats'.format(prefix=prefix)

    cov_df = pd.read_csv(os.path.join(cov_dir, cov_file), sep='\t', header=None, names=['chr', 'start', 'end',
                                                                                        'id', 'score', 'strand',
																		                'number_reads',
                                                                                        'intersected_bases', 'length',
                                                                                        'percentage_mapped'])
    coverage_df = cov_df.loc[cov_df['chr'].isin(selected_chrom)]

    stat_df = pd.read_csv(os.path.join(stat_dir, stat_file), sep='\t', header=None,
                          names=['chr', 'length', 'mapped_reads', 'unmapped_reads'])
    alnstats_df = stat_df.loc[stat_df['chr'].isin(selected_chrom)]

    cpm_df = count_cpm(coverage_df, alnstats_df, norm_factor)

    return cpm_df


def count_cpm(coverage, statistics, n_fact):
    total_reads = statistics['mapped_reads'].sum()
    coverage['CPM'] = (coverage['number_reads']/total_reads) * n_fact

    return coverage


def concatenate_cpm(region_dictionary, index_cols):
    # Renaming the columns so we can join them
    for key in list(region_dictionary.keys()):
        for i, df in enumerate(region_dictionary[key]):
            df.rename(columns={'CPM': 'CPM_{assay}_{num}'.format(assay=key,
                                                                 num=i)}, inplace=True)

    # index columns and only select those columns with cpm data
    all_df = [dsub.set_index(index_cols)[select_cpm(dsub.columns, 'CPM')] for d in list(region_dictionary.values()) for dsub in d]

    # there must always be at least two dataframes
    template = all_df[0].join(all_df[1])

    if len(all_df):
        for i in range(2, len(all_df)):
            template = template.join(all_df[i])

    return template


def select_cpm(column_list, selected):
    cols = [c for c in column_list if c.startswith(selected)]
    return cols

def main():
    parser = argparse.ArgumentParser(description='Arguments for coverage calculation')
    parser.add_argument('--coverage_dir', type=str, help='directory with coverage information')
    parser.add_argument('--stats_dir', type=str, help='directory with alignment stats information')
    parser.add_argument('--info', type=str, help='yaml file with sample relations')
    parser.add_argument('-o', type=str, help='output directory')
    args = parser.parse_args()

    with open(args.info, 'r') as metainfo:
        experiment_info = yaml.safe_load(metainfo)

    promoters_dictionary = {}
    enhancers_dictionary = {}

    print(READING_CHROMOSOMES)

    try:
        print('Calculating CPM')
        for epi_mark in list(experiment_info.keys()):
            print('\t\tpromoters...')
            promoters_cpm = {epi_mark: [norm_reads(cov_dir=args.coverage_dir, stat_dir=args.stats_dir,
                                                   prefix=pre, assay=epi_mark, region='promoters',
                                                   selected_chrom=READING_CHROMOSOMES,
                                                   norm_factor=NORMALIZATION_FACTOR) for pre in experiment_info[epi_mark]]}
            promoters_dictionary.update(promoters_cpm)

            print('\t\tenhancers...')
            enhancers_cpm = {epi_mark: [norm_reads(cov_dir=args.coverage_dir, stat_dir=args.stats_dir,
                                                   prefix=pre, assay=epi_mark, region='enhancers',
                                                   selected_chrom=READING_CHROMOSOMES,
                                                   norm_factor=NORMALIZATION_FACTOR) for pre in experiment_info[epi_mark]]}

            enhancers_dictionary.update(enhancers_cpm)
        print('\t\tconcatenating...')
        concatenated_prom = concatenate_cpm(promoters_dictionary, INDEX_COLS)
        concatenated_enha = concatenate_cpm(enhancers_dictionary, INDEX_COLS)

        promoters_out = os.path.join(args.o, 'promoters_cpm.tsv')
        concatenated_prom.to_csv(promoters_out, sep='\t')

        enhancers_out = os.path.join(args.o, 'enhancers_cpm.tsv')
        concatenated_enha.to_csv(enhancers_out, sep='\t')

    finally:
        print('done')


if __name__ == '__main__':
    main()