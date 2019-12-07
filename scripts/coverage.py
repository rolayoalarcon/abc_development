import os
import yaml
import argparse
import subprocess

COVERAGE_COMMAND = 'bedtools coverage -a {bedfile} -b {bamfile} > {output}'
IDXSTATS_COMMAND = 'samtools idxstats {bamfile} > {statout}'
SEPARATOR = ' && '


def create_commands_coverage(key, bamdirectory, peak_file, promoter_file, command_str, yaml_info, output_directory, separate):

    promoter_commands = [command_str.format(bedfile=promoter_file,
                                            bamfile='{}.bam'.format(os.path.join(bamdirectory, prefix)),
                                            output=format_coverage(output_directory, prefix, key, 'promoters')) for prefix in yaml_info[key]]
    peaks_commands = [command_str.format(bedfile=peak_file,
                                         bamfile='{}.bam'.format(os.path.join(bamdirectory, prefix)),
                                         output=format_coverage(output_directory, prefix, key, 'enhancers')) for prefix in yaml_info[key]]

    promoter_oneline = separate.join(promoter_commands)
    peaks_oneline = separate.join(peaks_commands)

    assay_line = separate.join([promoter_oneline, peaks_oneline])

    return assay_line


def create_commands_stat(key, bamdirectory, statdirectory,  command_str, yaml_info, separate):
    stats_commands = [command_str.format(bamfile='{}.bam'.format(os.path.join(bamdirectory, prefix)),
                                         statout='{}.stats'.format(os.path.join(statdirectory, prefix))) for prefix in yaml_info[key]]
    assay_stats = separate.join(stats_commands)

    return assay_stats


def format_coverage(outdir, bam_pre, astype, genomic_region):
    outfile = '{assay}_{pre}_{gr}.tsv'.format(assay=astype,
                                          pre=bam_pre,
                                          gr=genomic_region)
    output = os.path.join(outdir, outfile)

    return output


def execute_instruction(complete_command):
    print('Executing:\n {}'.format(complete_command))
    p = subprocess.Popen(complete_command, shell=True)

    (stdoutdata, stderrdata) = p.communicate()

    if stderrdata != None:
        raise RuntimeError("Command failed.")

def main():
    parser = argparse.ArgumentParser(description='Arguments for coverage calculation')
    parser.add_argument('--peaks', type=str, help='ATAC peak bedfile')
    parser.add_argument('--promoters', type=str, help='Promoter bedfile')
    parser.add_argument('--bam_directory', type=str, help='directory with bam files')
    parser.add_argument('--info', type=str, help='yaml file with sample relations')
    parser.add_argument('--o_coverage', type=str, help='coverage output directory')
    parser.add_argument('--o_stats', type=str, help='alignment stats output directory')
    args = parser.parse_args()

    # This list will contain all the commands for each individual assay type
    coverage_list = []
    stats_list = []

    # Open yaml file
    with open(args.info, 'r') as metainfo:
        experiment_info = yaml.safe_load(metainfo)

    try:
        print('Creating command for each epigenetic mark...')
        for epi_mark in list(experiment_info.keys()):

            # Coverage command
            mark_command = create_commands_coverage(key=epi_mark, bamdirectory=args.bam_directory, peak_file=args.peaks,
                                                    promoter_file=args.promoters, command_str=COVERAGE_COMMAND,
                                                    yaml_info=experiment_info, output_directory=args.o_coverage,
                                                    separate=SEPARATOR)

            coverage_list.append(mark_command)

            # stats command
            mark_stats = create_commands_stat(key=epi_mark, bamdirectory=args.bam_directory,
                                              statdirectory=args.o_stats, command_str=IDXSTATS_COMMAND,
                                              yaml_info=experiment_info, separate=SEPARATOR)
            stats_list.append(mark_stats)

        print('Creating complete commands...')
        complete_coverage = SEPARATOR.join(coverage_list)
        complete_stats = SEPARATOR.join(stats_list)

        execute_instruction(complete_coverage)
        execute_instruction(complete_stats)

    finally:
        print('done')


if __name__ == '__main__':
    main()