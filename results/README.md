# Results
This folder contains the supporting population results for the StarPhase publication.
These results are automatically generated using the `scripts/AnalyzeData.py` script.
The list of files generated:

* `ancestry_images` - Contains the automatically generated images for ancestry-specific distribution of haplotypes.
  * `{gene}_distribution.png` - For the given `{gene}`, the distribution of haplotypes for each ancestry. Only the top 20 are shown. If more than 20 are present, then the top 19 are shown and the rest are labeled as "Other".
  * `{gene}_{seq_type}_delta_distribution.png` - This is just for the HLA genes. For the given `{seq_type}` (DNA or cDNA), this will show the edit distance between the consensus sequence and the database sequence. These are grouped by exact match, off-by-one, minor and major delta values.
  * `cpic_stack.png` - A custom joint image created for the publication.
  * `CYP2D6_cn_distribution.png` - The distribution of CYP2D6 haplotype copy number values. Full CYP2D6 and hybrids are counted separately.
  * `CYP2D6_dip_func_distribution.png` - The distribution of metabolizer categories for the _diplotypes_ of each CYP2D6 population. These were generated using the [lookup tables provided by PharmGKB](https://www.pharmgkb.org/page/cyp2d6RefMaterials).
  * `CYP2D6_impact_distribution.png` - The distribution of functional categories for the _haplotypes_ of each CYP2D6 population. These were generated using the [lookup tables provided by PharmGKB](https://www.pharmgkb.org/page/cyp2d6RefMaterials).
  * `HLA_joint.png` - A joint image created for the publication. It contains the four datasets from `{gene}_{seq_type}_delta_distribution.png`.
* `combined_aggregate.tsv` - A file containing the combined metrics from all input cohorts. These are the "raw" data values used to create all results. If a single input cohort TSV is used, this file should be identical to that single input.
* `haplotype_aggregate.tsv` - A file containing the aggregate haplotype metrics from all input cohorts. These are metrics that have been processed from the combined cohort, and primarily represent the raw data values going into the generated images.
* `observed_haps.png` - The number of haplotypes observed at least once in the whole cohort. The percentage of the database is above each bar and used to color the bar (red-yellow-green color scale).
