# HiFiTargetEnrichmentQC
QC scripts for the Twist Alliance Dark Genes Panel

This is a Snakemake scripts to benchmark data output per gene.
The  data required to run the benchmark is 
1- Bam file
2- Phased SNVs and Indels file
3- SVs file

For Snakemake to run naming conventiones needs to be follwod 

1. Create a link to the bam file
2. Create a link to the SV vcf file e.g., HG002.SVs.vcf NOTE: the first part of the name must equals the bam file name so in this case HG002
3. Create a SNVs directory e.g., mkdir -p HG002_SNVs NOTE: the first part of the name must equals the bam file name so in this case HG002
4. in the previous direcotry link the SNVs file with the following name: phased_merge_output_AF_filtered.vcf.gz also link the phased_merge_output_AF_filtered.vcf.gz.tbi