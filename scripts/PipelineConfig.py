
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
SNAKEMAKE_PROFILE = f'{DATA_FOLDER}/snakemake_profile'
REFERENCE_FASTA = f'{DATA_FOLDER}/human_GRCh38_no_alt_analysis_set.fasta'
STARPHASE_DATABASE = f'{DATA_FOLDER}/starphase_db/v0.14.1/pbstarphase_20240826.json.gz'
SNAKEFILE = f'{RULES_FOLDER}/CohortPipeline.smk'
