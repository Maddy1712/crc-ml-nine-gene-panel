# Leakage-controlled ML identification of a nine-gene transcriptomic panel for colorectal cancer

Code accompanying the manuscript:

> Tiwari M, Mani A. *Machine Learning-Based Identification and Cross-Cohort Validation of a
> Stability-Selected Nine-Gene Transcriptomic Panel for Colorectal Cancer.* (submitted)
> Archived release: Zenodo DOI **[add after release]**

## What the pipeline does

1. **Development data:** three public ENA RNA-seq cohorts (PRJNA596460, PRJNA658001, PRJNA218851; 221 samples).
   **External test:** TCGA-COAD/READ (701 samples), never used for filtering, selection, tuning or thresholds.
2. **Nested cross-validation** (5 outer x 3 inner folds, stratified by condition and cohort). Inside each
   training partition only: low-expression filtering, cohort-adjusted DESeq2, Elastic Net alpha tuning, feature selection.
3. **Stability selection:** genes chosen in >= 3 of 5 outer folds form the final panel.
4. **Five classifiers** on the panel: KNN, ridge logistic regression, random forest, SVM, XGBoost.
5. **Leave-one-cohort-out (LOCO)** validation across the three ENA cohorts.
6. **Locked external validation** in TCGA: panel, scaling and trained models fixed, no refitting.

Seed: `2026`. Key settings are at the top of `R/crc_ml_pipeline.R`.

## Repository layout

```
R/crc_ml_pipeline.R              main analysis (nested CV, LOCO, TCGA test, figures, tables)
scripts/00_install_packages.R    installs R dependencies
scripts/preprocessing/           FASTQ download, QC/trimming/alignment, featureCounts (templates)
data/README.md                   input file formats and where to download them
docs/REVISION_NOTES.md           changes made for release + points to disclose
CITATION.cff, LICENSE            how to cite / MIT license
```

## Quick start

```bash
git clone https://github.com/maddy1712/crc-ml-nine-gene-panel.git
cd crc-ml-nine-gene-panel
Rscript scripts/00_install_packages.R
# download the processed inputs from Zenodo into ./data (see data/README.md)
Rscript R/crc_ml_pipeline.R            # or: Rscript R/crc_ml_pipeline.R <data_dir> <out_dir>
```

Outputs go to `results/` (CSV tables, PDF/PNG figures, `LOCKED_MODEL_ENA_panel.rds`, `R_SessionInfo.txt`).
The run fits DESeq2 dozens of times, so allow substantial time.

## Mapping of outputs to the manuscript

| Manuscript item | Output file(s) in `results/` |
|---|---|
| Table I (final panel) | `FINAL_BIOMARKER_PANEL.csv` |
| Table II (nested CV, LOCO) | `Nested_CV_Model_Summary.csv`, `LOCO_Metrics_By_Cohort.csv` |
| Table III (TCGA) | `TCGA_EXTERNAL_TEST_METRICS.csv` |
| Figures | `Figure_*.pdf/png` (fill in the figure-number mapping here) |

## Reproducibility notes

- Exact numbers can differ slightly with different R, package (notably `xgboost`, `glmnet`, `DESeq2`) or OS versions.
  The full session of the reported run is in `docs/sessionInfo_original_run.txt`. Key versions: R 4.3.2 (Windows 11),
  DESeq2 1.42.1, glmnet 5.0, xgboost 3.2.1.1, randomForest 4.7-1.2, e1071 1.7-17, class 7.3-22, pROC 1.19.0.1,
  caret 7.0-1, dplyr 1.1.4, ggplot2 4.0.3.
- Input data are public (ENA, GDC). Processed matrices and metadata are on Zenodo: **[add link]**.

## Contact

Corresponding author: Ashutosh Mani (amani@mnnit.ac.in), Motilal Nehru National Institute of Technology, Prayagraj.
