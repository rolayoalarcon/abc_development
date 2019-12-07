# abc_development
This is my implementation of the Activity By Contact model for predicting regulatory interactions between enhancers and transcribed genes. This method is described [here](https://www.biorxiv.org/content/10.1101/529990v1). At the time that I started developing this, their original code was not available yet.
  
  
# How to use

So we start with: 
- bedfile for promoter regions
- bedfile for ATAC or DHS peaks
- Bam files for ATAC (or DHS) and H3K27ac


### Getting coverage
So we start by calculating the number of reads present in our promoter and enhancer regions. To do this we use the the **coverage.py** script. At the same time, we determine the total number of reads aligned with idxstats.

python coverage.py --peaks ../data/atac_peaks/atac_peaks_chr22.txt --promoters ../data/genomic_information/gene_promoters_all.bed --bam_directory ../data/bam_files/ --info ../data/metainformation/config.yaml --o_coverage ../data/coverage --o_stats ../data/aln_stats/

### Determining CPM
After this, we normalize the number of mapped reads by calculating the Counts per Million. For this we use the **cpm_columns.py**. 

python cpm_columns.py --coverage_dir ../data/coverage/ --stats_dir ../data/aln_stats/ --info ../data/metainformation/config.yaml -o ../data/cpm_files/


### Getting predictions
Now we use the **abc_inference.py**. For this script we need:
- cpm fileswe generated in the previous step
- Hi-C data (unfortunatley this data is too heavy to upload)

python abc_inference.py --enhancers ../data/cpm_files/enhancers_cpm.tsv --promoters ../data/cpm_files/promoters_cpm.tsv --tss ../data/genomic_information/gene_tss_all.txt --hic ../data/external/K562_filteredRegularized_contactCount.tsv --bincoord ../data/external/K562_filteredBins.bed -o results.tsv -p 4

### Done :)
