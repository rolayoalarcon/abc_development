import numpy as np
from scipy import stats


# This function calculates the activity of each enhancer by returning the
# the geometric mean of H3K27 and DHS reads at a given enhancer region.
def activity_of_element(row, activity_cols_1, activity_cols_2, activity_cols_3):
    activity_cols = [activity_cols_1, activity_cols_2, activity_cols_3]
    activity_reads = row[activity_cols].tolist()
    geom_mean = stats.gmean(activity_reads)
    return geom_mean


# This function uses liftover to convert the coordinates. We take into account
# the name of the columns where the coordinates are.
def coordinate_conversion(row, chrom_name, start_name, end_name, conversion_function):

    start_conversion = conversion_function(*row[[chrom_name, start_name]])
    end_conversion = conversion_function(*row[[chrom_name, end_name]])

    if len(start_conversion) > 0:
        ch_new, start_new, _, _ = start_conversion[0]
    else:
        ch_new, start_new = -1, -1

    if len(end_conversion) > 0:
        _, end_new, _, _ = end_conversion[0]
    else:
        end_new = -1

    return (min(start_new, end_new), max(start_new, end_new))


def chromosome_conversion(dataframe, chromosome_relation):
    dataframe.set_index('chr', inplace=True)
    dataframe_joined = dataframe.join(chromosome_relation)

    return dataframe_joined['ensembl_chr'].values


def format_bedstyle(dataframe, name_col, score_col, start_col, end_col):

    dataframe_imp = dataframe[['chromosome_conv', start_col, end_col,
                               name_col, score_col, 'strand']].copy()

    dataframe_imp.replace(-1, np.nan, inplace=True)
    dataframe_imp.dropna(inplace=True)
    dataframe_imp[start_col] = dataframe_imp[start_col].astype(int)
    dataframe_imp[end_col] = dataframe_imp[end_col].astype(int)
    return dataframe_imp


def rescale_bin(bin_id, regularized_matrix, base_number):
    binrow = regularized_matrix.loc[(regularized_matrix['bin1_id'] == bin_id)
                                    | (regularized_matrix['bin2_id'] == bin_id)]

    scaled_df = binrow.assign(scaled_count=(binrow['regularized_count']/binrow['regularized_count'].max())*base_number)
    return {bin_id: scaled_df}


def get_contacts(row, gamma_param, scale_param, rescaled_info, base_number, min_dist):
    tss_bin = row['tssBin_id']
    enhancer_bin = row['enhBin_id']
    separation = row['distance']

    #print(tss_bin)

    involved_bins = [tss_bin, enhancer_bin]

    pseudocount = power_law(distance=separation, gamma=gamma_param, scale=scale_param, min_cutoff=min_dist)


    bin_info = rescaled_info[tss_bin]

    max_contacts = bin_info['regularized_count'].max()

    count = bin_info.loc[(bin_info['bin1_id'] == min(involved_bins))
                         & (bin_info['bin2_id'] == max(involved_bins)), 'scaled_count'].values
    if len(count) > 0:
        count_pseudocount = count[0] + ((pseudocount/max_contacts) * base_number)
        if count_pseudocount > base_number:
            count_pseudocount = base_number
    else:
        count_pseudocount = (pseudocount/max_contacts) * base_number

    return count_pseudocount


# This function calculates a pseudocount of the expected number of contacts
# given a certain distance.
def power_law(distance, gamma, scale, min_cutoff):
    if distance > min_cutoff:
        correction = scale*(distance**gamma)
    elif distance <= min_cutoff:
        correction = scale*(min_cutoff**gamma)

    return correction


def promoter_coord(row, promoter_l, start_name):
    if row['strand'] == '+':
        second_coord = row[start_name] - promoter_l
    elif row['strand'] == '-':
        second_coord = row[start_name] + promoter_l
    return (min(row[start_name], second_coord), max(row[start_name], second_coord))


def close_element(gene_id, genomic_info, gamma_self, scale_self, min_self, rescaled_name, scale_mult):
    #print(gene_id)
    gene_info = genomic_info[genomic_info['gene_id'] == gene_id].copy()

    denominator = sum(gene_info['activity'].values * gene_info['corrected_contacts'].values)
    promoter_act = gene_info['promoter_activity'].values[0]

    promoter_contacts = selfBin_contact(gene_info['tssBin_id'].values[0], gamma_self,
                                        scale_self, min_self, rescaled_name, scale_mult)

    promoter_value = promoter_act * promoter_contacts
    denominator += promoter_value

    return {gene_id: denominator}


def selfBin_contact(bin_id, gamma_param, scale_param, min_dist, rescaled_info, base_number):
    pseudocount = power_law(distance=0, gamma=gamma_param, scale=scale_param, min_cutoff=min_dist)

    bin_info = rescaled_info[bin_id]
    max_contacts = bin_info['regularized_count'].max()

    count = bin_info.loc[(bin_info['bin1_id'] == bin_id) & (bin_info['bin2_id'] == bin_id), 'scaled_count'].values

    if len(count) > 0:
        count_pseudocount = count[0] + (pseudocount/max_contacts) * base_number
        if count_pseudocount > base_number:
            count_pseudocount = base_number
    else:
        count_pseudocount = (pseudocount/max_contacts) * base_number
    return count_pseudocount


def abc_score(gene_id, genomic_info, correction_terms):
    gene_info = genomic_info[genomic_info['gene_id'] == gene_id].copy()
    activity = gene_info['activity'].values * gene_info['corrected_contacts'].values

    abc_scores = activity/correction_terms[gene_info['gene_id'].values[0]]

    gene_info['ABC_score'] = abc_scores
    return gene_info
