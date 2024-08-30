# StarPhase cohort pipeline
This repo contains a pipeline and collection of scripts to generate aggregate summary statistics for cohorts from different contributing sites.

## User instructions
This cohort pipeline assumes you are operating in a cluster environment and have conda installed.
To install conda, we recommend the [Bioconda setup instructions](https://bioconda.github.io).

### Generating a cohort summary
This is the main process for generating anonymized, ancestry-level aggregate statistics from a cohort.

1. Clone this repository to your local compute infrastructure:
```
git clone https://github.com/holtjma/starphase_pipeline.git
cd starphase_pipeline
```
2. Create a cohort file following the template provided in `./data/cohort_batches/template.tsv`. This requires populating the unique sample ID as well as the locations of the haplotagged BAM and phased VCF file. Note: the provided sample ID is expected to match the values inside the provided BAM and VCF.
3. The pipeline expects a reference FASTA file and [snakemake submission profile](https://snakemake.readthedocs.io/en/stable/executing/cli.html#profiles) to be located in the `./data` sub-folder. They can be soft-linked for ease-of-use:
```
ln -s {path_to_profile} ./data/snakemake_profile
ln -s {path_to_reference} ./data/human_GRCh38_no_alt_analysis_set.fasta
```
  * This pipeline can be run on a local machine as well. If this is desired, the `snakemake` invocation inside `./scripts/RunCohortPipeline.py` can be modified for your local use case.
4. Create and activate the provided conda environment:
```
conda env create --name starphase_pipeline --file ./conda_envs/starphase_pipeline.yaml
# follow prompts to install the environment
conda activate starphase_pipeline
```
5. Run the full pipeline, which will produce the final cohort summary at `./pipeline/aggregate/aggregate_summary.tsv`:
```
python3 ./scripts/RunCohortPipeline -Ax
```
6. Review `./pipeline/aggregate/aggregate_summary.tsv` which will have the following format:
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
