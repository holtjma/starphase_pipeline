# StarPhase cohort pipeline
This repo contains a pipeline and collection of scripts to generate aggregate summary statistics the StarPhase publication.
Contributing sites each run this pipeline on their local data files and then the aggregate summaries are collected to generate final figures.

## Pipeline overview
At a high level, the pipeline performs the following core steps:

1. `mosdepth` - Run on individual datasets. This generates summary coverage statistics that are used to filter out low coverage datasets. This tends to be the most computationally expensive step of the whole pipeline.
2. `peddy` - Run on individual datasets. This generates summary ancestry and sex predictions for the datasets.
3. `pbstarphase` - Run on individual datasets. This generates the PGx diplotype calls of interest.
4. Aggregation script - Run on the full cohort. This loads all the other output files, filters out low coverage datasets, and then aggregates PGx diplotype calls based on predicted ancestry and sex. The final output will be located in `./pipeline/aggregate/aggregate_summary.tsv`.

## User instructions
This cohort pipeline assumes you are operating in a cluster environment and have conda installed.
To install conda, we recommend the [Bioconda setup instructions](https://bioconda.github.io).
Mamba is also a viable alternative.

### Running the pipeline
This is the main process for generating anonymized, ancestry-level aggregate statistics from a cohort.

1. Clone this repository to your local compute infrastructure:
```
git clone https://github.com/holtjma/starphase_pipeline.git
cd starphase_pipeline
```
2. Create the cohort file following the template provided in `./data/cohort_batches/template.tsv`. Note this must be saved to a different file. This requires populating the unique sample ID as well as the locations of the haplotagged BAM and phased VCF file. The provided sample ID is expected to match the values inside the provided BAM and VCF.
3. The pipeline expects a reference FASTA file to be located in the `./data` sub-folder. It can be soft-linked for ease-of-use:
```
ln -s {path_to_reference} ./data/human_GRCh38_no_alt_analysis_set.fasta
```
4. Create and activate the provided conda environment:
```
conda env create --name starphase_pipeline --file ./conda_envs/starphase_pipeline.yaml
# follow prompts to install the environment
conda activate starphase_pipeline
```
5. Set up the run conditions for snakemake. There are two options:
  * **Recommended**: Set up a snakemake profile for your compute environment. Several templates are provided on [snakemake submission profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles). This is provided via the `--profile {SNAKEMAKE_PROFILE}` option, and passed through to the snakemake invocation.
  * Alternate: This pipeline can be run on a single machine as well by specifying `--cores {CORES}` instead of a snakemake profile. While this is easier to set up, it forces all work onto a single machine which may be slow for larger cohorts. This is option is passed through to the snakemake invocation.
6. Run the full pipeline, which will produce the final cohort summary at `./pipeline/aggregate/aggregate_summary.tsv`:
```
# dry-run of the pipeline
python3 ./scripts/RunCohortPipeline.py --profile {SNAKEMAKE_PROFILE} -A
# executes the commands
python3 ./scripts/RunCohortPipeline.py --profile {SNAKEMAKE_PROFILE} -Ax 
```
7. Review `./pipeline/aggregate/aggregate_summary.tsv` which will have the following format:
```
gene	ancestry	sex	starphase_diplotype	count
ABCG2	AFR	female	rs2231142 reference (G)/rs2231142 reference (G)	33
...
```

### Generating aggregate figures
This step is primarily for aggregating and visualization results from multiple cohorts.
However, it can be run on a single cohort file.

1. Copy all aggregate summary files (the output from [generating a cohort summary](#generating-a-cohort-summary)) into the `./data/cohort_aggregate_files` folder with distinct filenames.
2. Run `python3 ./scripts/AnalyzeData.py`, which will generate a collection of visualizations in the `./results` sub-folder.
