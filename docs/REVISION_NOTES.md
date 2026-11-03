# Revision notes (release version vs. the script used for the manuscript)

## Changes (none alters a numerical computation)
1. Header documents usage and inputs/outputs.
2. Input folder and output folder are now arguments (`DATA_DIR`, `OUTDIR`); defaults `data/` and `results/`.
   The old output folder name `CRC_ML_FINAL_RESULTS_with_geoVAlidation` was replaced by `results`.
3. Section 8 (single-fit Elastic Net, comparison only): scaling now uses ENA only. TCGA was previously passed
   as a dummy "test" matrix; the training output is identical, but TCGA no longer appears anywhere in feature selection.
4. Section 13 (extra-cohort validation, GSE292858) is switched off by default (`RUN_EXTRA_VALIDATION <- FALSE`).
   It is not part of the reported results.
5. Missing-input error message points to `data/README.md`.

## Verify before release
Run the original and the revised script on your machine and compare the key tables:
```bash
diff old_results/Nested_CV_Metrics_All_Folds.csv  results/Nested_CV_Metrics_All_Folds.csv
diff old_results/TCGA_EXTERNAL_TEST_METRICS.csv   results/TCGA_EXTERNAL_TEST_METRICS.csv
diff old_results/FINAL_BIOMARKER_PANEL.csv        results/FINAL_BIOMARKER_PANEL.csv
```
They should be identical.

## Points to disclose in the manuscript (Methods)
- Before any differential-expression step, the gene universe was restricted to genes present in both the ENA and the TCGA
  count matrices (`common_genes`). This step uses gene identifiers only (no expression values or labels), so it does not leak
  outcome information, but state it so the "TCGA not used before validation" claim is exact.
- Per-sample log2-CPM normalisation uses each sample's own library size and no labels.
- Panel-level TCGA statistics in `FINAL_BIOMARKER_PANEL.csv` are descriptive only and were not used to choose genes.
