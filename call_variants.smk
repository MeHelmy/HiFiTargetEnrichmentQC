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


# Preparing conda environements.
###############################
SNIFFLES_ENV = script_location+"/envs/sniffles2.yaml"
CLAIR_ENV = script_location+"/envs/clair3.yaml"
BCFTOOLS = script_location+"/envs/bcftools.yaml"
#############




rule all:
    # input: "{myfile}.SVs.vcf", "{{myfile}}_{SNVs_dir}".format(SNVs_dir=SNVs_dir), "{{myfile}}_{SNVs_dir}/phased_merge_output_{af}.vcf.gz".format(SNVs_dir=SNVs_dir, af="AF_filtered")
    input: "{myfile}.SVs.vcf",
           "{{myfile}}_{SNVs_dir}/phased_merge_output_{af}.vcf.gz".format(SNVs_dir=SNVs_dir, af="AF_filtered")
    output: touch("{myfile}.done")

rule Sniffles2:
    """
        Calling SVs using Sniffles2
    """
    input:
        bam = "{myfile}.bam",
        bai = "{myfile}.bam.bai",
    output: "{myfile}.SVs.vcf",
    message: "Calling SVs form {input.bam}"
    log: "{myfile}.SVs.log"
    conda: SNIFFLES_ENV
    params:
        reference = REFERENCES,
    threads: config['sniffles_threads']
    shell:"""
        echo "my env is here $CONDA_PREFIX" &&\
        sniffles --input {input.bam} --vcf {output} --reference {params.reference} --threads {threads} > {log} 2>&1
    """

rule Clair3:
    """
        Calling SNVs using Clair3
    """
    input:
        bam = "{myfile}.bam",
        bai = "{myfile}.bam.bai",
    # output: directory("{{myfile}}_{SNV_dir_name}".format(SNV_dir_name=SNVs_dir)),
    output: "{{myfile}}_{SNV_dir_name}/phased_merge_output.vcf.gz".format(SNV_dir_name=SNVs_dir),
    # output: "{{myfile}}_{SNV_dir_name}/phased_merge_output.vcf.gz".format(SNV_dir_name=SNVs_dir),
    log: "{{myfile}}_{SNV_dir_name}.log".format(SNV_dir_name=SNVs_dir),
    message: "Calling SVs form {input.bam}"
    conda: CLAIR_ENV
    params:
        reference = REFERENCES,
        platform = config['reads_technology'],
        bed_region = config['bed_file'],
        snv_out_dir = SNVs_dir
    threads: config['clair_threads']
    shell:"""
        bam={input.bam}
        bam_name=${{bam/.bam/}};
        run_clair3.sh\
         --bam_fn={input.bam}\
         --output=${{bam_name}}_{params.snv_out_dir}\
         --ref_fn={params.reference}\
         --threads {threads}\
         --platform={params.platform}\
         --model_path=$CONDA_PREFIX/bin/models/{params.platform}\
         --bed_fn={params.bed_region}\
         --parallel=$(which parallel)\
         --enable_phasing > {log} 2>&1
    """

rule FilterVariants:
    input:  "{{myfile}}_{SNVs_dir}/phased_merge_output.vcf.gz".format(SNVs_dir=SNVs_dir)
    output: "{{myfile}}_{SNVs_dir}/phased_merge_output_{{af}}.vcf.gz".format(SNVs_dir=SNVs_dir, af="AF_filtered")
    message: "Filtering SNVs self.fail('message')e {input} using allele frequency of {af}"
    conda: BCFTOOLS
    params:
        af = af
    threads: config['clair_threads']
    log: "{{myfile}}_{SNVs_dir}/phased_merge_output_{{af}}.log".format(SNVs_dir=SNVs_dir, af="AF_filtered")
    shell:"""
        bcftools view --threads {threads} -Oz -i 'AF>={params.af}' -o {output} {input} 2>{log} && tabix {output}
    """
## ------------------------------------------------------------------------------------ ##
## Success and failure messages
## ------------------------------------------------------------------------------------ ##

onsuccess:
	shell("mkdir -p snake_log && find . -maxdepth 1  \( -name 'snakejob*' -or -name 'slurm*' \) -type f -exec mv -t snake_log {{}}  \;")

onerror:
	shell("mkdir -p snake_log && find . -maxdepth 1  \( -name 'snakejob*' -or -name 'slurm*' \) -type f -exec mv -t snake_log {{}}  \;")
