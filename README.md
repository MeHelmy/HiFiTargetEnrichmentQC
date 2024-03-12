# HiFiTargetEnrichmentQC

QC scripts for the Twist Alliance Dark Genes Panel

This is a Snakemake script to benchmark data output per gene. The data required to run the benchmark are:

1. Bam file
2. Phased SNVs and Indels file
3. SVs file

For Snakemake to run, naming conventions need to be followed:

1. Create a link to the bam file.
2. Create a link to the SV VCF file, e.g., HG002.SVs.vcf. NOTE: The first part of the name must match the bam file name, so in this case, HG002.
3. Create a SNVs directory, e.g., mkdir -p HG002_SNVs. NOTE: The first part of the name must match the bam file name, so in this case, HG002.
4. In the previous directory, link the SNVs file with the following name: phased_merge_output_AF_filtered.vcf.gz. Also link the phased_merge_output_AF_filtered.vcf.gz.tbi.

An example of how to run the script:
The `HG002.stat` is an example of the output needed. Replace HG002 with your sample name if needed.
