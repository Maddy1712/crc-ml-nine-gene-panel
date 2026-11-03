# Installs the packages used by the pipeline. Run once:  Rscript scripts/00_install_packages.R
if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
BiocManager::install(c("DESeq2", "Rsubread"), ask = FALSE, update = FALSE)
install.packages(c("glmnet", "randomForest", "e1071", "xgboost", "class",
                   "pROC", "ggplot2", "dplyr", "caret"))
# For exact reproducibility, record versions with renv:
#   install.packages("renv"); renv::init()   # then commit renv.lock
# Versions of the reported run: R 4.3.2, DESeq2 1.42.1 (Bioconductor 3.18), glmnet 5.0, xgboost 3.2.1.1,
# randomForest 4.7-1.2, e1071 1.7-17, pROC 1.19.0.1, caret 7.0-1, dplyr 1.1.4, ggplot2 4.0.3
# (see docs/sessionInfo_original_run.txt). Use R 4.3.2 for the closest match.
