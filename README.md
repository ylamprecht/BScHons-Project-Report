# BScHons-Project-Report

## Bioinformatics pipeline for whole genome sequencing (WGS) data generated to investigate genetic variants associated with COVID-19 vaccination status
This project aims to identify host genetic variants associated with severe COVID-19 in unvaccinated individuals from African populations, to determine how these may be influenced by vaccination. Through genetic screening, the findings could help prioritise certain individuals for vaccination based on their genetic risk factors. The specific objectives are:

### 1. Comparison of Variant Prioritisation Tools
to evaluate the performance of eVai and Franklin classifying relevant genetic variants and thereafter select the most appropriate tool.

#### Scripts:
- **Data Cleaning**: `clean_eVai.py`, `clean_franklin.py`, and `clean_VIPR.py` scripts refine each tool's output to retain only relevant columns, facilitating subsequent analysis.
- **Merging**: `merge_with_VIPR.py` merges eVai or Franklin outputs with VIPR pathogenicity scores, preparing data for classification comparisons.
- **Comparisons**:
  - **General Comparison**: `compare_all_variants.R` generates Venn diagrams comparing the variants and genes identified by eVai and Franklin.
  - **Classification-Specific Comparison**: `compare_common_variants.R` compares specific variant classifications (Pathogenic/Likely Pathogenic, VUS, Benign/Likely Benign) between eVai and Franklin, with results presented as Venn diagrams.
  - **Interrater Agreement**: `Kappa_comparison.R` calculates Cohen’s weighted Kappa to assess interrater agreement between eVai and Franklin on variant classifications.
  - **VIPR Comparison**: `compare_VIPR.R` produces boxplots to show the distribution of VIPR pathogenicity scores across eVai’s and Franklin’s variant classifications.

#### Supplementary Material:
- Weighted Kappa values for 20 severe COVID-19 samples.
- Venn diagrams showing gene, variant, and classification overlaps between eVai and Franklin for each of the 20 severe COVID-19 samples.

### 2. Design a bioinformatics pipeline
that identifies candidate host genetic variants associated with an increased risk of severe/critical COVID-19, using WGS.

#### Pipeline Steps:

1. **Variant Calling**:
   - **Genome Preparation**: Run `prep_genome.pbs` once on the Centre for High Performance Computing (CHPC) with the reference genome file `hg38.fa`.
   - **Calling Variants**: Submit `variant_calling.pbs`, adjusting sample number and specifying the directory containing the raw FASTQ files. Transfer the resulting VCF file to the local computer for further analysis.

2. **Variant Classification**:
   - Upload VCF files to [Franklin](https://franklin.genoox.com/) as "Inherited Disease Single Cases". Franklin applies default filters, excluding synonymous, low-confidence, and common variants.
   - Export the results and organise them into cohort-specific folders for further comparison.

3. **Cohort Comparison**:
   - Run `processing.py` in the directory containing the cohort-specific CSV files downloaded from Franklin. The script:
     - Merges and refines CSV files.
     - Counts variant occurrences by classification.
     - Combines cohort data to create eight output files summarising variant distribution across cohorts.

These reseults allow the identification of variants of specific classifications present in unvaccinated individuals with severe COVID-19 but absent in other cohorts.

## Requirements
- **Python** for data refinement and cohort comparison.
- **R** for variant comparison and statistical analysis.
- **PBS Script Execution**: Access to the Centre for High Performance Computing (CHPC) for running `prep_genome.pbs` and `variant_calling.pbs`.
