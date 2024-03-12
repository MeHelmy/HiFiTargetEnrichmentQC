# import Lib
############
import os, glob
from snakemake.utils import min_version

############################

# Snake Version
###############
min_version("6.15.5")
##############

# Config File
#############
if os.path.isfile("config.yaml"):
    configfile: "config.yaml"
else:
    sys.exit("Looks like there is no config.yaml file in " + os.getcwd() + " make sure there is one or at least specify one with the --configfile commandline parameter.")
#############


# GET WORKING DIRECTORY DEFAULT IS CURRENT DIRECTORY
####################################################
script_location = os.getcwd()
data_dir =  config["sample_directory"] if config['sample_directory'] else os.getcwd()
REFERENCES = config["reference"]
SNVs_dir = config['SNVs_dir_name']
af = config['af']
#############

# Preparing conda environments.
###############################
STAT_ENV=script_location+"/envs/stat.yaml"
#############

## ------------------------------------------------------------------------------------ ##
## Main Rules
## ------------------------------------------------------------------------------------ ##

rule all:
    # input: "{myfile}.mosdepth.done",
    # input: "{myfile}_coverage.bed",
    # input: "{myfile}_zero_covered.bed",
    # input: "{myfile}_average_zero_cover.bed",
    # input: "{myfile}_average_zero_min_max_cover.bed", "{myfile}_SNVs_benchamrk.done"
    # input: "{myfile}_average_zero_min_max_cover.bed", "{myfile}_SNVs_benchamrk_values.tsv"
    # input: "{myfile}_coverage_and_benchamrk.tsv"
    input:  "{myfile}_average_zero_min_max_cover_exon.bed",
     "{myfile}_base_stat.txt",
     # "{{myfile}}_{plots}".format(plots=config['plot_dir']),
     # "{myfile}_variant_coverage.bed",
     "{{myfile}}_{SV_bench_dir}".format(SV_bench_dir=config['SV_OUT_DIR']), "{{myfile}}_{SNVs_bench_split}.tsv".format(SNVs_bench_split=config['SNVs_split_file_name'])
    output: touch("{myfile}.stat")

rule Mosdepth:
    """
    Input for Rule MergeCoverageWithBam
    """
    input: "{myfile}.bam"
    output: "{myfile}.mosdepth.done"
    params:
        bed_file=config['bed_file'],
        coverage_prefix="{myfile}",
        mq = config['MAPQ']
    threads: config['mosdepth_threads']
    log: "{myfile}.mosdepth.log"
    conda: STAT_ENV
    shell:"""
        mosdepth -t {threads}  --by {params.bed_file} --mapq {params.mq} {params.coverage_prefix} {input} > {log} 2>&1 && touch "{output}"
    """

rule MergeCoverageWithBam:
    """
    Input to Rule ZeroCoverage
    """
    input: "{myfile}.mosdepth.done"
    output: "{myfile}_coverage.bed"
    message: "merge average coverage per gene with gene statistics"
    params:
        bed_file=config['bed_region']
    conda: STAT_ENV
    shell:"""
        tail -n +2 {params.bed_file} | bedtools intersect -wo -a - -b  {wildcards.myfile}.regions.bed.gz -f 1 -r | cut -f1-12,17 > {output}
    """

rule ZeroCoverage:
    """
    Input to Rule MergeZeroCoverageWithStat
    """
    input: "{myfile}_coverage.bed"
    output: "{myfile}_zero_covered.bed"
    message: "Calculating zero covered percentage"
    log: "{myfile}_zero_covered.log"
    conda: STAT_ENV
    params:
        covmin = config['covmin'],
        cov0 = config['cov0'],
        covmin_value = config['covmin_value1'],
    shell:"""
        bedtools intersect -wo -a {wildcards.myfile}.per-base.bed.gz -b {input} | awk 'BEGIN{{OFS="\\t";}}$4==0 {{print $5,$6,$7,$9, $(NF)}}' | datamash -s -g 1,2,3,4 sum 5 1> {params.cov0}.bed 2>>{log} &&\
        bedtools intersect -wo -a {wildcards.myfile}.per-base.bed.gz -b {input} | awk -v value={params.covmin_value} 'BEGIN{{OFS="\\t";}}$4<value {{print $5,$6,$7,$9, $(NF)}}' | datamash -s -g 1,2,3,4 sum 5 1> {params.covmin}.bed 2>>{log} &&\
        bedtools intersect -wao -a {params.covmin}.bed -b {params.cov0}.bed -r  -f 1 | bedtools sort | cut -f 1-5,10 | sed 's/\s-1/\t0/' 1> {output} 2>> {log} &&\
        rm {params.covmin}.bed {params.cov0}.bed
    """

rule MergeZeroCoverageWithStat:
    """
    Input to Rule MinAndMaxCover
    Count both bases with zero coverage and bases less than X coverage, X value from params
    """
    input:
        average_cover = "{myfile}_coverage.bed",
        zero_cover = "{myfile}_zero_covered.bed",
    output: "{myfile}_average_zero_cover.bed"
    message: "Running {rule} to merge coverage with statistics genes file"
    log: "{myfile}_average_zero_cover.log"
    conda: STAT_ENV
    shell:"""
        bedtools intersect -wao -a {input.average_cover} -b {input.zero_cover} -f 1 -r |  cut -f1-13,18-19 | sed 's/\t-1/\t0/;s/\t\./\t0/g' 1> {output} 2> {log}
    """

rule MinAndMaxCover:
    """
    Input to Rule UpdateSNVsStatWithBenchmark
    """
    input: "{myfile}_average_zero_cover.bed"
    output: "{myfile}_average_zero_min_max_cover.bed"
    message: "Calculating all coverages for genes"
    log: "{myfile}_average_zero_min_max_cover.log"
    conda: STAT_ENV
    shell:"""
        bedtools intersect -wo -a {wildcards.myfile}.per-base.bed.gz -b {input} | datamash -s -g 5,6,7,9 min 4 max 4 | bedtools intersect -wo -a {input} -b - -f 1 -r |  cut -f1-15,20,21 1>{output} 2>{log}
    """

rule SNVsBenchmark:
    """
    Input Rule to CollectSNVsBenchmark
    """
    input:
        bed_file = "{genes_dir}/{{gene}}.bed".format(genes_dir=config['genes_dir']),
        snv = "{{myfile}}_SNVs/phased_merge_output_{af}.vcf.gz".format(af="AF_filtered"),
        # truth_file = "{myfile}_gold_filtered.vcf.gz",
        # snv = "{myfile}_SNVs/phased_merge_output.vcf.gz"
    output: directory("{SNV_bench_dir}/{{gene}}__{{myfile}}".format(SNV_bench_dir=config['output_dir_name']))
    threads: config['SNVs_bench_threads']
    params:
        GOLD_SNP=config['gold_vcf'],
        ref=config['REF_38'],
        extra=config['extra_rtg'],
        samples_to_compare=config['samples_to_compare']
    conda: STAT_ENV
    shell:"""
      rtg RTG_MEM=16G vcfeval --baseline {params.GOLD_SNP}  --bed-regions {input.bed_file} -c {input.snv} -o {output} -t {params.ref} -T {threads} {params.extra} {params.samples_to_compare} 
    """

rule CollectSNVsBenchmark:
    """
    Input Rule to WriteSNVsBenchmark
    """
    input: lambda wildcards: expand(["{{SNV_bench_dir}}/{gene}__{{myfile}}".format(gene=os.path.basename(i).split('.')[0]) for i in glob.glob(config['genes_dir']+"/*.bed")], SNV_bench_dir=config['output_dir_name'], myfile=wildcards.myfile)
    output: touch("{myfile}_SNVs_benchamrk.done")

rule WriteSNVsBenchmark:
    """
    Input Rule to UpdateSNVsStatWithBenchmark
    """
    input:
        Bench_done = "{myfile}_SNVs_benchamrk.done",
        snv = "{{myfile}}_SNVs/phased_merge_output_{af}.vcf.gz".format(af="AF_filtered"),
        SNVs_benchmark_dir_list = glob.glob(config['output_dir_name']+"/*")
    output: "{myfile}_SNVs_benchamrk_values.tsv"
    message: "Writing SNVs Benchamrk for {rule}"
    conda: STAT_ENV
    params:
        genes=config['genes_dir']
    shell:"""
        for i in {input.SNVs_benchmark_dir_list};
            do
            filename=$(basename $i)
            gene_name=$(echo $filename | sed 's/__{wildcards.myfile}.*//g')
            variants=$(bcftools view -H {input.snv} -R {params.genes}/${{gene_name}}.bed | wc -l)
            snps=$(bcftools view -H -v snps {input.snv} -R {params.genes}/${{gene_name}}.bed | wc -l)
            indels=$(bcftools view -H -v indels {input.snv} -R {params.genes}/${{gene_name}}.bed | wc -l)
            Bench_values=$(tail -n +3  ${{i}}/summary.txt | head -n 1 | awk 'BEGIN{{OFS="\t"}}{{print $(NF-2), $(NF-1), $(NF)}}')
            echo "${{gene_name}}\t${{variants}}\t${{snps}}\t${{indels}}\t${{Bench_values}}";
            done  >> {output}
    """

            # for i in {input.SNVs_benchamrk_dir_list};
            #     do
            #     filename=$(basename $i)
            #     gene_name=$(echo $filename | sed 's/__{wildcards.myfile}.*//g')
            #     variants=$(bcftools view -H {wildcards.myfile}_SNVs/phased_merge_output.vcf.gz -R {params.genes}/${{gene_name}}.bed | wc -l)
            #     snps=$(bcftools view -H -v snps {wildcards.myfile}_SNVs/phased_merge_output.vcf.gz -R {params.genes}/${{gene_name}}.bed | wc -l)
            #     indels=$(bcftools view -H -v indels {wildcards.myfile}_SNVs/phased_merge_output.vcf.gz -R {params.genes}/${{gene_name}}.bed | wc -l)
            #     Bench_values=$(tail -n +3  ${{i}}/summary.txt | head -n 1 | awk 'BEGIN{{OFS="\t"}}{{print $(NF-2), $(NF-1), $(NF)}}')
            #     echo "${{gene_name}}\t${{variants}}\t${{snps}}\t${{indels}}\t${{Bench_values}}";
            #     done  >> {output}

rule UpdateSNVsStatWithBenchmark:
    """
    Input Rule for PhaseBenchmark
    """
    input:
        benchmark="{myfile}_SNVs_benchamrk_values.tsv",
        stat="{myfile}_average_zero_min_max_cover.bed"
    output: "{myfile}_coverage_and_benchamrk.tsv"
    message: "merging benchamrk with coverage"
    conda: STAT_ENV
    params:
        maxcov = config['covmin_value1'],
        genes_without_benchmark = config['genes_without_benchmark'],
    shell:"""
        temp_out=UpdateSNVsStatWithBenchmark_temp_$(date +%s).tsv
        python {data_dir}/merge_coverage_bench.py {input.benchmark} {input.stat} {params.maxcov} > $temp_out
        python {data_dir}/scripts/update_pandas_df.py $temp_out {output} {params.genes_without_benchmark} && rm -rf $temp_out
    """

rule PhaseBenchmark:
    """
    Input Rule for BaseStat
    """
    input:
        benchmark="{myfile}_coverage_and_benchamrk.tsv",
        phased_snvs="{{myfile}}_{SNVs_dir}/phased_merge_output_{af}.vcf.gz".format(SNVs_dir=SNVs_dir, af="AF_filtered"),
    output: "{myfile}_coverage_and_benchamrk_phase.tsv"
    message: "Add phasing benchamrk to statistics matrix"
    conda: STAT_ENV
    params:
        genes_region = config['bed_file'],
    shell:"""
        tmp="PhaseBenchmark_phase_stat_"$(date +%s).tsv
        bcftools view -e 'GT="mis" || GT="RR"' {input.phased_snvs} | bedtools intersect -wo -a {params.genes_region} -b - | bedtools sort | cut -f 4,13,14 | python {data_dir}/scripts/phase_analysis.py >  $tmp &&\
        python {data_dir}/scripts/update_phase_info.py {input.benchmark} $tmp {output} && rm -rf $tmp
    """

rule ExoneAverageCoverage:
    """
    Input Rule to MergeExonCoverageWithBam
    """
    input: "{myfile}.bam"
    output: "{myfile}.mosdepth.exon.done"
    params:
        bed_file=config['exon_file'],
        coverage_prefix="{myfile}.exon",
        mq = config['MAPQ'],
    threads: config['mosdepth_threads']
    log: "{myfile}.mosdepth.exon.log"
    message: "Running {rule} to calculate exon coverage"
    conda: STAT_ENV
    shell:"""
        mosdepth -t {threads}  --by {params.bed_file} --mapq {params.mq} {params.coverage_prefix} {input} > {log} 2>&1 && touch "{output}"
    """

rule MergeExonCoverageWithBam:
    """
    Input Rule to ZeroExonCoverage
    """
    input: "{myfile}.mosdepth.exon.done"
    output: "{myfile}_coverage_exon.bed"
    message: "Running {rule} merge average coverage per exon with exon file"
    params:
        bed_file=config['exon_file']
    conda: STAT_ENV
    shell:"""
        bedtools intersect -wo -a {params.bed_file} -b  {wildcards.myfile}.exon.regions.bed.gz -f 1 -r | cut -f 1-7,12 | uniq > {output}
    """

rule ZeroExonCoverage:
    """
    Input Rule to MergeZeroExonCoverageWithStat
    """
    input: "{myfile}_coverage_exon.bed"
    output: "{myfile}_zero_covered_exon.bed"
    message: "Running {rule} calculating zero covered percentage for exon"
    log: "{myfile}_zero_covered_exon.log"
    conda: STAT_ENV
    shell:"""
        bedtools intersect -wo -a {wildcards.myfile}.exon.per-base.bed.gz -b {input} | awk 'BEGIN{{OFS="\\t";}}$4==0 {{print $5,$6,$7,$11, $(NF)}}' | datamash -g 1,2,3,4 sum 5 1> {output} 2> {log}
    """

rule MergeZeroExonCoverageWithStat:
    """
    Input Rule to MinAndMaxExonCover
    """
    input:
        average_cover = "{myfile}_coverage_exon.bed",
        zero_cover = "{myfile}_zero_covered_exon.bed",
    output: "{myfile}_average_zero_cover_exon.bed"
    message: "Running {rule} to merge exon coverage with statistics genes file"
    log: "{myfile}_average_zero_cover_exon.log"
    params:
        covmin_value = config["covmin_value1"]
    conda: STAT_ENV
    shell:"""
        tmp="MergeZeroExonCoverageWithStat_average_and_zero_"$(date +%s).bed
        tmp2="MergeZeroExonCoverageWithStat_lt_minval_"$(date +%s).bed
        bedtools intersect -wao -a {input.average_cover} -b {input.zero_cover} -f 1 -r |  cut -f1-8,13 | sed 's/\s-1/\t0/' 1> $tmp 2> {log} &&\
        bedtools intersect -wo -a {wildcards.myfile}.exon.per-base.bed.gz -b $tmp |  awk -v value={params.covmin_value} 'BEGIN{{OFS="\\t";}}$4<value {{print $5,$6,$7, $8,$9, $10, $(NF)}}' | datamash -s -g 1,2,3,4,5,6 sum 7 > $tmp2
        bedtools intersect -wao -a $tmp -b $tmp2 -r  -f 1 | bedtools sort | awk 'BEGIN{{OFS="\\t";}}{{if($16=="."){{$16=0}} print}}' | datamash -g 1,2,3,4,5,6,7,8,9 sum 16 >        {output} && rm -rf $tmp $tmp2
    """

rule MinAndMaxExonCover:
    """
    Input Rule to all
    """
    input: "{myfile}_average_zero_cover_exon.bed"
    output: "{myfile}_average_zero_min_max_cover_exon.bed"
    message: "Running {rule} calculating all coverages for exon"
    log: "{myfile}_average_zero_min_max_cover_exon.log"
    conda: STAT_ENV
    shell:"""
        bedtools intersect -wo -a {wildcards.myfile}.exon.per-base.bed.gz -b {input} |  datamash -s -g   5,6,7,8,9,10,11 min 4 max 4 | bedtools intersect -wo -a {input} -b - -f 1 -r |  cut -f1-10,18,19 | uniq | sed '1i Chr\\tStart\\tEnd\\tgene_id\\ttranscript_id\\tgene_name\\texon_id\\tAverage_cov\\tZero_cover\\tlt8\\tMin_bp_cover\\tMax_bp_cover' 1>{output} 2>{log}
    """

rule bam_intersect_gene:
    """
    Input Rule to BaseStat
    """
    input:
        full_bam = "{myfile}.bam",
        full_bam_bai = "{myfile}.bam.bai",
    output:
        bam_genes = "{{myfile}}_{buffer}_.bam".format(buffer=config['buffer']),
    message: "Intersecting {input.full_bam} with CMRGs"
    log: "{{myfile}}_{buffer}_.log".format(buffer=config['buffer']),
    conda: STAT_ENV
    params:
        genes_region = config['bed_file'],
        ref = config['reference'],
        buffer = config['buffer'],
    shell:"""
        (bedtools slop -b {params.buffer} -i {params.genes_region} -g {params.ref}.fai | bedtools intersect -wa -ubam  -abam {input.full_bam} -b - > {output} && samtools index {output}) 2> {log}
    """

rule BaseStat:
    """
    Input Rule to all
        - Calcaute on targte reads
        - Percntage of bases > 8, 10, amd 20x
        - Plot all graphes
        - Print genes with no coverage (average coverage is < 1)
    """
    input:
        full_bam = "{myfile}.bam",
        full_bam_bai = "{myfile}.bam.bai",
        reads_on_target_bam = "{{myfile}}_{buffer}_.bam".format(buffer=config['buffer']),
        mosdepth = "{myfile}.mosdepth.done",
        full_stat_file = "{myfile}_coverage_and_benchamrk_phase.tsv",
        variant_uncovered = "{myfile}_variant_uncovered.tsv",
    output:
        stat = "{myfile}_base_stat.txt",
        plots = directory("{{myfile}}_{plots}".format(plots=config['plot_dir'])),
    message: "Calculating on target reads and bases >= X coverage"
    log: "{myfile}_base_stat.log"
    conda: STAT_ENV
    params:
        read_flag = config['read_flag'],
        covmin_value1 = config['covmin_value1'],
        covmin_value2 = config['covmin_value2'],
        covmin_value3 = config['covmin_value3'],
        percentage_uncovered_genes = config['percentage_uncovered_genes'],
        genes_without_benchmark = config['genes_without_benchmark'],
    shell:"""
        tmp_df_uncovered_genes=BaseStat_uncover_genes_df_"$(date +%s)".tsv &&\
        total_reads=$(samtools view -c -F {params.read_flag} {input.full_bam}) &&\
        reads_on_target=$(samtools view -c -F {params.read_flag} {input.reads_on_target_bam}) &&\
        percentage_on_target=$(printf  "%.2f%%" $(echo "($reads_on_target/$total_reads)*100" | bc -l) ) &&\
        echo -e "##Reads_on_target\\tTotal_reads\\tPercentage_on_target\\n#reads\\t$reads_on_target\\t$total_reads\\t$percentage_on_target\\n" >> {output.stat} &&\
        bases_in_region=$(samtools view  -F  {params.read_flag} {input.reads_on_target_bam} | awk '{{s+=length($10)}}END{{print s}}') &&\
        bases_in_total=$(samtools view  -F  {params.read_flag} {input.full_bam} | awk '{{s+=length($10)}}END{{print s}}') &&\
        percentage_bases_on_target=$(printf  "%.2f%%" $(echo "($bases_in_region/$bases_in_total)*100" | bc -l) ) &&\
        echo -e "##Bases_on_target\\tTotal_bases\\tPercentage_bases_on_target\\n#bases\\t$bases_in_region\\t$bases_in_total\\t$percentage_bases_on_target\\n" >> {output.stat} &&\
        number_of_gene_bases=$( tail -n +2 {input.full_stat_file} | cut -f 4 | datamash sum 1 ) &&\
        bases_gt_first_value=$(bedtools intersect -wo -a {wildcards.myfile}.per-base.bed.gz -b {wildcards.myfile}_coverage.bed | awk -v value={params.covmin_value1} 'BEGIN{{OFS="\\t";}}$4>value {{print $5,$6,$7,$9, $(NF)}}' | datamash sum 5) &&\
        bases_gt_second_value=$(bedtools intersect -wo -a {wildcards.myfile}.per-base.bed.gz -b {wildcards.myfile}_coverage.bed | awk -v value={params.covmin_value2} 'BEGIN{{OFS="\\t";}}$4>value {{print $5,$6,$7,$9, $(NF)}}' | datamash sum 5) &&\
        bases_gt_third_value=$(bedtools intersect -wo -a {wildcards.myfile}.per-base.bed.gz -b {wildcards.myfile}_coverage.bed | awk -v value={params.covmin_value3} 'BEGIN{{OFS="\\t";}}$4>value {{print $5,$6,$7,$9, $(NF)}}' | datamash sum 5) &&\
        echo -e "\\n##Percentage of bases with>x coverage" >> {output.stat} &&\
        printf "#per\\tPercentage of bases with %d coverage or more: %.2f%%\\n" {params.covmin_value1} $(echo "($bases_gt_first_value/$number_of_gene_bases)*100" | bc -l) >> {output.stat} &&\
        printf "#per\\tPercentage of bases with %d coverage or more: %.2f%%\\n" {params.covmin_value2} $(echo "($bases_gt_second_value/$number_of_gene_bases)*100" | bc -l) >> {output.stat} &&\
        printf "#per\\tPercentage of bases with %d coverage or more: %.2f%%\\n" {params.covmin_value3} $(echo "($bases_gt_third_value/$number_of_gene_bases)*100" | bc -l) >> {output.stat} &&\
        echo -e "\\n##Uncovered genes" >> {output.stat} &&\
        awk 'NR==1 {{print "#UN\\t"$0}} NR>1 && $13<1 {{print "#UN\\t"$0}} END{{print "\\n"}}' {input.full_stat_file} >> {output.stat} &&\
        echo -e "\\n##Genes where 50% of gene body < 8x coverage" >> {output.stat} &&\
        awk -v v={params.percentage_uncovered_genes} 'NR==1 {{printf "#50\\t%s\\tGene_percent_lt8\\n", $0}}  NR>1 && $14/$4*100 > v {{printf "#50\\t%s\\t%2.f\\n", $0, $14/$4*100}} END{{print "\\n"}}' {input.full_stat_file} >> {output.stat} &&\
        python {script_location}/scripts/add_uncover_genes_annotation.py {input.full_stat_file} {input.variant_uncovered} $tmp_df_uncovered_genes 2>> {log} &&\
        mv $tmp_df_uncovered_genes {input.full_stat_file} && rm -rf $tmp_df_uncovered_genes 2>> {log} &&\
        python {script_location}/scripts/stat_plot.py {input.full_stat_file} {output.plots} {params.genes_without_benchmark} 2>> {log}
    """

# bases_in_region=$(samtools view -h  -F  {params.read_flag} {input.reads_on_target_bam} | bedtools genomecov -ibam - -d | awk '{{s+=$3}}END{{print s}}') &&\
# bases_in_total=$(samtools view -h -F  {params.read_flag} {input.full_bam} | bedtools genomecov -ibam - -d | awk '{{s+=$3}}END{{print s}}') &&\

rule VariantCoverage:
    """
    Input Rule to BaseStat and CreateTruthSetVCF
    """
    input:
        bam_genes = "{{myfile}}_{buffer}_.bam".format(buffer=config['buffer']),
    output:
        variant_coverage = "{myfile}_variant_coverage.bed",
        variant_uncovered = "{myfile}_variant_uncovered.tsv",
    message: "Calculating coverage for GIAB variants"
    log: "{myfile}_variant_coverage.log"
    params:
        variant_bed = config["gold_variant_bed"],
        genes_bed = config["bed_file"],
    conda: STAT_ENV
    shell:"""
    bedtools coverage -b {input.bam_genes} -a {params.variant_bed} |  awk '{{printf "%s\\t%d\\t%d\\t%s\\t%s\\t%d\\t%.2f\\n",  $1, $2, $3, $6, $7, $(NF - 3), $(NF)}}' | bedtools intersect -wo -a - -b {params.genes_bed} | cut -f 1-7,11 1>{output.variant_coverage} 2>{log} &&\
    awk '$6<1' {output.variant_coverage} | datamash -g 8 count 6 1> {output.variant_uncovered} 2>>{log}
    """

rule CreateTruthSetVCF:
    """
    Input to SNVsBenchmark
    """
    input:
        variant_coverage = "{myfile}_variant_coverage.bed",
    output:
        truth_vcf = "{myfile}_gold_filtered.vcf.gz",
        truth_vcf_tbi = "{myfile}_gold_filtered.vcf.gz.tbi",
    message: "Filter gold vcf to contain covered variants only"
    log: "{myfile}_gold_filtered.log"
    params:
        truth_vcf = config['gold_vcf'],
        min_cover = config['min_snv_coverage_to_benchmark'],
    conda: STAT_ENV
    shell:"""
    awk '$6<{params.min_cover}' {input.variant_coverage} | bedtools intersect -header -v -a {params.truth_vcf} -b - | bgzip > {output.truth_vcf} && tabix {output.truth_vcf} 2>{log}
    """

rule SVsBenchmark:
    """
    Input Rule to all
    """
    input:
        SV_file = "{myfile}.SVs.vcf",
    output:
        Benchamrk_dir = directory("{{myfile}}_{SV_bench_dir}".format(SV_bench_dir=config['SV_OUT_DIR']))
    message: "Benchamarking SVs"
    log: "{myfile}_SV_benchamrk.log"
    params:
        GOLD_SV = config["GOLD_SV"],
        BED_SV = config["BED_SV"],
        extra = config['extra'],
    conda: STAT_ENV
    shell:"""
        tmp="SVsBenchamrk_SV_DEL_INS_"$(date +%s).vcf.gz &&\
        bedtools sort  -header -i {input.SV_file}| bcftools view -i 'INFO/SVTYPE="DEL" | INFO/SVTYPE="INS"' -Oz -o "$tmp" &&\
        tabix "$tmp" &&\
        truvari bench -b {params.GOLD_SV} -c $tmp -f {REFERENCES} -o {output.Benchamrk_dir} --includebed {params.BED_SV} {params.extra} > {log} &&\
        rm -rf $tmp
    """

rule SNVsBenchmarkSeparately:
    """
    Input Rule to all
    Here we identify fp, fn, tp, and tn per variant.
    """
    input:
        # SNV_bech_dir = "{SNV_bench_dir}".format(SNV_bench_dir=config['output_dir_name']),
        SNV_bech_dir = "{myfile}_SNVs_benchamrk.done"
    output:
        Benchamrk_dir = "{{myfile}}_{SNVs_bench_split}.tsv".format(SNVs_bench_split=config['SNVs_split_file_name'])
    message: "Benchamarking SNVs Separately"
    log: "{myfile}_SNVs_bech_split.log"
    params:
        min_quality = config['SNVs_min_quality'],
        benchamrk_dir = config['output_dir_name'],
    conda: STAT_ENV
    shell:"""
        python {script_location}/scripts/short_variant_benchamrk_split.py {params.benchamrk_dir} {params.min_quality} > {output.Benchamrk_dir}
    """
## ------------------------------------------------------------------------------------ ##
## Success and failure messages
## ------------------------------------------------------------------------------------ ##

onsuccess:
	shell("mkdir -p snake_log && find . -maxdepth 1  \( -name 'snakejob*' -or -name 'slurm*' \) -type f -exec mv -t snake_log {{}} \;")

onerror:
	shell("mkdir -p snake_log && find . -maxdepth 1  \( -name 'snakejob*' -or -name 'slurm*' \) -type f -exec mv -t snake_log {{}}  \;")
