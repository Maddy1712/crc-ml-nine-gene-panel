# Input data

The processed inputs are **not stored in git** (size). Download them from Zenodo
(**[add DOI/link]**) and place them in this folder:

| File | Content |
|---|---|
| `ENA_counts1_cleaned.tsv` | gene x sample raw counts, ENA cohort 1 |
| `ENA_counts2_cleaned.tsv` | gene x sample raw counts, ENA cohort 2 |
| `ENA_counts3.tsv` | gene x sample raw counts, ENA cohort 3 |
| `TCGA_COAD_READ_Merged_Count_Matrix.tsv` | GDC STAR gene counts, COAD + READ merged (701 samples) |
| `All_Metadata_Combined.csv` | sample IDs and tumor/normal labels for all cohorts |

**Cohort order:** ENA1 = PRJNA596460 (40 samples), ENA2 = PRJNA658001 (131 samples), ENA3 = PRJNA218851 (50 samples).

## Expected formats

- Count files are tab-separated. Column 1 holds the gene ID (Ensembl IDs; version suffixes such as `.12` are stripped
  automatically). All other columns are integer raw counts, one per sample. The TCGA file may have a
  `gene_name` column in position 2; it is used only to map IDs to symbols.
- The metadata CSV needs one column with the sample ID (`sample`, `Sample_ID`, `sample_name` or `id`)
  and one with the label (`condition`, `group`, `status`, `sample_type` or `type`).
  Accepted labels: Normal / normal tissue / control / healthy / non-tumor / 11 -> **Normal**;
  Tumor / tumour / cancer / tumor tissue / 01 -> **Tumor**.
- Sample IDs in the count files must match the metadata IDs (non-alphanumeric characters are ignored when matching).

## How the inputs were generated

- ENA cohorts: `scripts/preprocessing/01_download_fastq.sh` -> `02_qc_trim_align.sh` -> `03_featurecounts.R`
  (then describe here any "cleaning" applied to produce the `*_cleaned.tsv` files).
- TCGA: STAR gene-level counts downloaded from the GDC (https://portal.gdc.cancer.gov), COAD and READ merged.
