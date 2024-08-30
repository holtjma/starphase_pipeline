
import csv

import sys
sys.path.append('../scripts')
from PipelineConfig import *

from SharedParsing import load_batches

#######################################################################
# Functions to pull files from pipeline area
#######################################################################

SAMPLE_METADATA = load_batches()

def get_bam_file(wildcards):
    '''
    Gets the haplotagged BAM file from our metadata
    '''
    sample = wildcards.sample
    return SAMPLE_METADATA[sample]['tagged_bam']

def get_phased_vcf_file(wildcards):
    '''
    Gets the phased VCF file from our metadata
    '''
    sample = wildcards.sample
    return SAMPLE_METADATA[sample]['phased_vcf']

def get_sample_id(wildcards):
    '''
    Returns the VCF sample ID, this should match the sample label typically
    '''
    sample_id = wildcards.sample
    return sample_id

#######################################################################
# Primary root rules
#######################################################################

def get_peddy_cohort():
    ret = []
    for sample in SAMPLE_METADATA:
        # adding the HTML file triggers the rule to run
        ret.append(f'{PIPELINE_FOLDER}/peddy/{sample}.html')
    return ret

rule peddy_cohort:
    input:
        get_peddy_cohort()

def get_starphase_cohort():
    ret = []
    for sample in SAMPLE_METADATA:
        # just add the diplotypes to trigger the rule
        ret.append(f'{PIPELINE_FOLDER}/starphase/{sample}.diplotypes.json')
    return ret

rule starphase_cohort:
    input:
        get_starphase_cohort()

localrules: collation_file
rule collation_file:
    output:
        tsv=f'{PIPELINE_FOLDER}/collation/all_batches.tsv'
    params:
        script=f'{SCRIPTS_FOLDER}/GenerateCollationFile.py'
    shell: '''
        python3 {params.script} \
            -o {output.tsv}
    '''

#######################################################################
# Ancestry checks
#######################################################################

rule peddy:
    input:
        vcf=get_phased_vcf_file
    output:
        ped="{pipeline}/peddy/{sample}.ped",
        html="{pipeline}/peddy/{sample}.html",
        # TODO: there are other outputs, but we don't need to specify them since this is a one-of analysis
    params:
        prefix="{pipeline}/peddy/{sample}",
        sample_id=get_sample_id,
        sites='hg38'
    resources:
        mem_mb=4*1024,
        runtime=5 #minutes
    threads: 1
    conda: f"{ENV_FOLDER}/peddy.yaml"
    log: "{pipeline}/logs/peddy/{sample}.log"
    benchmark: "{pipeline}/benchmark/peddy/{sample}.log"
    shell: '''
        echo -e "{params.sample_id}\t{params.sample_id}\t0\t0\t0\t0" > {output.ped}
        peddy \
            --procs {threads} \
            --sites {params.sites} \
            --prefix {params.prefix} \
            {input.vcf} \
            {output.ped}
        '''

#######################################################################
# PGx tools
#######################################################################

rule starphase_bam:
    input:
        database=STARPHASE_DATABASE,
        reference=REFERENCE_FASTA,
        bam=get_bam_file,
        vcf=get_phased_vcf_file
    output:
        json="{pipeline}/starphase/{sample}.diplotypes.json",
        pharmcat_tsv="{pipeline}/starphase/{sample}.for_pharmcat.tsv",
        debug_folder=directory("{pipeline}/starphase/{sample}_debug")
    params:
        options='' # may need to adjust if someone wants to use targeted
    resources:
        # partition="short",
        mem_mb=16*1024, # megabytes, increase if needed
        runtime=60 # minutes, typically much less than this
    threads: 1
    conda: f"{ENV_FOLDER}/starphase.yaml"
    log: "{pipeline}/logs/starphase/{sample}.log"
    benchmark: "{pipeline}/benchmark/starphase/{sample}.log"
    shell: '''
        pbstarphase \
            diplotype \
            -v \
            {params.options} \
            --threads {threads} \
            --database {input.database} \
            --reference {input.reference} \
            --vcf {input.vcf} \
            --bam {input.bam} \
            --output-calls {output.json} \
            --pharmcat-tsv {output.pharmcat_tsv} \
            --output-debug {output.debug_folder}
        '''