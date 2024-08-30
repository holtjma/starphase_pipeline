
import argparse as ap
import os

from PipelineConfig import *

PREPRINT_SNAKEFILE = f'{SCRIPTS_FOLDER}/preprint/PreprintExtras.smk'

if __name__ == '__main__':
    description = 'Runs the extra tools for the pre-print'
    p = ap.ArgumentParser(description=description, formatter_class=ap.RawTextHelpFormatter)

    #standard options we have for all snakemake pipelines
    p.add_argument('-x', '--execute', dest='execute', action='store_true', default=False, help='execute the commands (default: False)')
    p.add_argument('-C', '--clean', dest='clean', action='store_true', default=False, help='clean all output files (default: False)')
    p.add_argument('-A', '--all', dest='target_all', action='store_true', default=False, help='target all outputs files (default: False)')

    #things we want to generate
    p.add_argument('-s', '--starphase-cohort', dest='starphase_cohort', action='store_true', default=False, help='generate the StarPhase cohort results (default: False)')
    p.add_argument('-p', '--peddy-cohort', dest='peddy_cohort', action='store_true', default=False, help='generate the Peddy cohort results (default: False)')
    
    args = p.parse_args()

    #generate all the targets from the above commands
    targets = []
    
    if args.target_all or args.starphase_cohort:
        targets.append('starphase_cohort')

    if args.target_all or args.peddy_cohort:
        targets.append('peddy_cohort')
    
    snakemake_frags = [
        'snakemake',
        #standard options
        '--printshellcmds',
        '--keep-going',
        '--use-conda',
        #'--use-singularity', '--singularity-args "-B /pbi -B /home --cleanenv"',
        '--jobs', '1000',
        #things we need to pass/configure
        '--profile', SNAKEMAKE_PROFILE,
        '--snakefile', SNAKEFILE
    ]

    if not args.execute:
        print('Dry run only: enabled')
        snakemake_frags.append('--dry-run')
    
    if args.clean:
        print('Clean-up: enabled')
        snakemake_frags.append('--delete-all-output')

    if len(targets) > 0:
        snakemake_frags += targets
        cmd = ' '.join(snakemake_frags)
        print(cmd)
        ret_code = os.system(cmd)
        if ret_code != 0:
            raise Exception(f'Snakemake error {ret_code}, see above')
    else:
        print('WARNING: No targets specified, exiting without doing any work.')
