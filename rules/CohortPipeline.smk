
import csv
import glob

# get at the PipelineConfig
import sys
script_folder = os.path.join(os.path.dirname(workflow.basedir), 'scripts')
sys.path.append(script_folder)
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

def get_mosdepth_cohort():
    ret = []
    for sample in SAMPLE_METADATA:
        # we just need the summary file
        ret.append(f'{PIPELINE_FOLDER}/mosdepth/{sample}.mosdepth.summary.txt')
    return ret

rule mosdepth_cohort:
    input:
        get_mosdepth_cohort()

def get_peddy_cohort():
    ret = []
    for sample in SAMPLE_METADATA:
        # add the two files we need, even though they come from the same rule
        ret.append(f'{PIPELINE_FOLDER}/peddy/{sample}.het_check.csv')
        ret.append(f'{PIPELINE_FOLDER}/peddy/{sample}.sex_check.csv')
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

def get_batch_files():
    batch_files = sorted(glob.glob(f'{DATA_FOLDER}/cohort_batches/*.tsv'))
    # template.tsv is in here, but that shouldn't be a problem
    return batch_files

# The next two rules run fast and could be localrules, but we want to capture the logs in the event of failure
rule collation_file:
    input:
        batch_files=get_batch_files()
    output:
        tsv=f'{PIPELINE_FOLDER}/collation/all_batches.tsv'
    params:
        script=f'{SCRIPTS_FOLDER}/GenerateCollationFile.py'
    log: f'{PIPELINE_FOLDER}/logs/collation/all_batches.log'
    shell: '''
        python3 {params.script} \
            -o {output.tsv} \
            > {log} 2>&1
    '''

rule aggregate_summary:
    input:
        peddy_files=get_peddy_cohort(),
        starphase_files=get_starphase_cohort(),
        coverage_files=get_mosdepth_cohort(),
        tsv=f'{PIPELINE_FOLDER}/collation/all_batches.tsv'
    output:
        tsv=f'{PIPELINE_FOLDER}/aggregate/aggregate_summary.tsv'
    params:
        script=f'{SCRIPTS_FOLDER}/AggregateData.py'
    log: f'{PIPELINE_FOLDER}/logs/aggregate/aggregate_summary.log'
    shell: '''
        python3 {params.script} \
            -i {input.tsv} \
            -o {output.tsv} \
            > {log} 2>&1
    '''

#######################################################################
# QC checks
#######################################################################

rule mosdepth:
    input:
        bam=get_bam_file
    output:
        summary="{pipeline}/mosdepth/{sample}.mosdepth.summary.txt",
        dist="{pipeline}/mosdepth/{sample}.mosdepth.global.dist.txt"
    params:
        prefix="{pipeline}/mosdepth/{sample}",
        chrom="chr1" # we just need a rough estimate to rule out really low coverage, chr1 will do
    resources:
        mem_mb=8*1024,
        runtime=2*60 #minutes
    threads: 4
    conda: f"{ENV_FOLDER}/mosdepth.yaml"
    log: "{pipeline}/logs/mosdepth/{sample}.log"
    benchmark: "{pipeline}/benchmark/mosdepth/{sample}.tsv"
    shell: '''
        mosdepth \
            -n \
            --threads {threads} \
            --fast-mode \
            --chrom {params.chrom} \
            {params.prefix} \
            {input.bam} \
            > {log} 2>&1
        '''

#######################################################################
# Ancestry checks
#######################################################################

rule peddy:
    input:
        vcf=get_phased_vcf_file
    output:
        ped="{pipeline}/peddy/{sample}.ped",
        het_check="{pipeline}/peddy/{sample}.het_check.csv",
        sex_check="{pipeline}/peddy/{sample}.sex_check.csv",
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
            {output.ped} \
            > {log} 2>&1
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
            --output-debug {output.debug_folder} \
            > {log} 2>&1
        '''
