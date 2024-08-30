
# Root config items
SNAKEMAKE_PROFILE = '/home/UNIXHOME/mholt/snakemake_profiles/slurm'
REFERENCE_FASTA = '/pbi/collections/appslabht/reference/human_GRCh38_no_alt_analysis_set.fasta'

# can be derived from file location
import os
PROJECT_DIR = os.path.dirname(os.path.dirname(__file__))

#anything derived from the root config
RULES_FOLDER = f'{PROJECT_DIR}/rules'
ENV_FOLDER = f'{PROJECT_DIR}/conda_envs'
PIPELINE_FOLDER = f'{PROJECT_DIR}/pipeline'
DATA_FOLDER = f'{PROJECT_DIR}/data'
RESULTS_FOLDER = f'{PROJECT_DIR}/results'
SCRIPTS_FOLDER = f'{PROJECT_DIR}/scripts'

# sub-paths
STARPHASE_DATABASE = f'{DATA_FOLDER}/starphase_db/v0.14.1/pbstarphase_20240826.json.gz'
SNAKEFILE = f'{RULES_FOLDER}/CohortPipeline.smk'

# put all batch in these configs
BATCH_FILES = [
    f'{DATA_FOLDER}/cohort_batches/preprint_cohort.tsv'
]
