############################################################
# CRC ML BIOMARKER DISCOVERY - LEAKAGE-CONTROLLED PIPELINE
# Development: ENA1+ENA2+ENA3   External test: TCGA COAD+READ
# Nested CV: training-only DEG -> Elastic Net -> 5 ML models
#
# USAGE (fresh R session, from the repository root):
#   Rscript R/crc_ml_pipeline.R                      # uses ./data and writes ./results
#   Rscript R/crc_ml_pipeline.R <data_dir> <out_dir> # custom locations
# or, interactively:  source("R/crc_ml_pipeline.R")
#
# INPUT : gene-level count matrices + metadata (see data/README.md)
# OUTPUT: tables, figures, locked model and sessionInfo in <out_dir>
# Upstream processing (FastQC/Trimmomatic/HISAT2/featureCounts) is in
# scripts/preprocessing/.
############################################################

############################ 1. PACKAGES
required_packages <- c("DESeq2", "glmnet", "randomForest", "e1071", "xgboost",
                       "class", "pROC", "ggplot2", "dplyr", "caret")
missing_pk <- required_packages[!vapply(required_packages, requireNamespace,
                                        logical(1), quietly = TRUE)]
if (length(missing_pk) > 0) {
  stop("Install missing packages first: ", paste(missing_pk, collapse = ", "),
       "\n  BiocManager::install(c('", paste(missing_pk, collapse = "','"), "'))")
}
suppressPackageStartupMessages({
  library(DESeq2); library(glmnet); library(randomForest); library(e1071)
  library(xgboost); library(pROC); library(ggplot2); library(dplyr); library(caret)
})
cat("xgboost version:", as.character(packageVersion("xgboost")),
    "| DESeq2 version:", as.character(packageVersion("DESeq2")), "\n")

############################ 2. SETTINGS
options(stringsAsFactors = FALSE)
SEED <- 2026
set.seed(SEED)
cli_args <- commandArgs(trailingOnly = TRUE)
DATA_DIR <- if (length(cli_args) >= 1) cli_args[1] else "data"
OUTDIR   <- if (length(cli_args) >= 2) cli_args[2] else "results"
dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
# Optional locked validation on an extra cohort (section 13). Off by default; the
# manuscript results do not depend on it.
RUN_EXTRA_VALIDATION <- FALSE

N_OUTER <- 5
N_INNER <- 3
PADJ_CUTOFF <- 0.05
LOG2FC_CUTOFF <- 1
MIN_DEG <- 10
MAX_DEG_FOR_ML <- 500
ALPHA_GRID <- c(0.25, 0.50, 0.75)
MIN_ENET_GENES <- 5
MAX_ENET_GENES <- 50
N_RF <- 500
N_XGB_ROUNDS <- 200
# TRUE = L2-penalised (ridge) logistic regression. Plain glm() on a near-separable ENA training set gives
# extreme coefficients that do not transfer across platforms (it collapsed on TCGA: AUC 0.62, sens 0.04).
LR_PENALISED <- TRUE

ena_files <- file.path(DATA_DIR, c("ENA_counts1_cleaned.tsv", "ENA_counts2_cleaned.tsv", "ENA_counts3.tsv"))
tcga_file <- file.path(DATA_DIR, "TCGA_COAD_READ_Merged_Count_Matrix.tsv")
metadata_file <- file.path(DATA_DIR, "All_Metadata_Combined.csv")

miss <- c(ena_files, tcga_file, metadata_file)
miss <- miss[!file.exists(miss)]
if (length(miss) > 0) stop("Missing input files (see data/README.md). Looked in: ", normalizePath(DATA_DIR, mustWork = FALSE),
                         "\n", paste(miss, collapse = "\n"))

############################ 3. HELPERS
clean_gene_names <- function(x) {
  x <- toupper(trimws(sub("\\.[0-9]+$", "", as.character(x))))
  x[x == ""] <- NA
  x
}

aggregate_duplicate_genes <- function(mat) {
  mat <- as.matrix(mat)
  genes <- clean_gene_names(rownames(mat))
  keep <- !is.na(genes)
  mat <- mat[keep, , drop = FALSE]
  genes <- genes[keep]
  if (anyDuplicated(genes)) {
    mat <- rowsum(mat, group = genes, reorder = FALSE)
  } else {
    rownames(mat) <- genes
  }
  mat
}

read_count_file <- function(file) {
  cat("Reading:", file, "\n")
  x <- read.delim(file, header = TRUE, sep = "\t", check.names = FALSE,
                  stringsAsFactors = FALSE)
  if (ncol(x) < 2) stop("Invalid count file: ", file)
  genes <- x[[1]]
  # Drop annotation columns (e.g. TCGA 'gene_name' in column 2) - keep numeric sample columns only
  is_num <- vapply(x[-1], is.numeric, logical(1))
  if (any(!is_num)) {
    cat("  Dropping non-numeric annotation column(s):",
        paste(colnames(x)[-1][!is_num], collapse = ", "), "\n")
  }
  mat <- as.matrix(x[, -1, drop = FALSE][, is_num, drop = FALSE])
  suppressWarnings(mode(mat) <- "numeric")
  n_na <- sum(is.na(mat))
  if (n_na > 0) {
    warning(file, ": ", n_na, " non-numeric/NA cells set to 0. ",
            "Check the file has ONLY gene IDs in column 1 (no extra annotation columns).")
  }
  rownames(mat) <- genes
  mat <- aggregate_duplicate_genes(mat)
  mat[is.na(mat)] <- 0
  mat[mat < 0] <- 0
  round(mat)
}

find_metadata_column <- function(df, candidates) {
  idx <- match(tolower(candidates), tolower(colnames(df)))
  idx <- idx[!is.na(idx)]
  if (length(idx) == 0) return(NULL)
  colnames(df)[idx[1]]
}

normalize_condition <- function(x) {
  x <- tolower(trimws(as.character(x)))
  out <- rep(NA_character_, length(x))
  out[x %in% c("normal", "normal tissue", "normal sample", "control", "healthy",
               "non-tumor", "non tumor", "non_tumor", "11")] <- "Normal"
  out[x %in% c("tumor", "tumour", "cancer", "tumor tissue", "tumour tissue", "01")] <- "Tumor"
  out
}

match_samples <- function(matrix_samples, metadata_samples) {
  res <- match(matrix_samples, metadata_samples)
  c1 <- gsub("[^A-Za-z0-9]", "", matrix_samples)
  c2 <- gsub("[^A-Za-z0-9]", "", metadata_samples)
  miss <- is.na(res)
  res[miss] <- match(c1[miss], c2)
  res
}

safe_auc <- function(y, pred) {
  y <- factor(y, levels = c("Normal", "Tumor"))
  if (length(unique(y)) < 2 || length(unique(pred)) < 2 || anyNA(pred)) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(response = y, predictor = pred,
                                 levels = c("Normal", "Tumor"),
                                 direction = "<", quiet = TRUE)))
}

calculate_metrics <- function(y, prob, threshold = 0.5) {
  y <- factor(y, levels = c("Normal", "Tumor"))
  pred <- factor(ifelse(prob >= threshold, "Tumor", "Normal"),
                 levels = c("Normal", "Tumor"))
  cm <- table(Actual = y, Predicted = pred)
  TN <- cm["Normal", "Normal"]; TP <- cm["Tumor", "Tumor"]
  FP <- cm["Normal", "Tumor"];  FN <- cm["Tumor", "Normal"]
  sens <- if (TP + FN > 0) TP / (TP + FN) else NA_real_
  spec <- if (TN + FP > 0) TN / (TN + FP) else NA_real_
  prec <- if (TP + FP > 0) TP / (TP + FP) else NA_real_
  f1 <- if (is.na(prec) || is.na(sens) || prec + sens == 0) NA_real_ else
    2 * prec * sens / (prec + sens)
  data.frame(Accuracy = (TP + TN) / sum(cm), Balanced_Accuracy = mean(c(sens, spec)), Sensitivity = sens, Specificity = spec,
             Precision = prec, F1 = f1, AUC = safe_auc(y, prob))
}

get_class_weights <- function(y) {
  y <- factor(y, levels = c("Normal", "Tumor"))
  tab <- table(y)
  w <- rep(1, length(y))
  w[y == "Normal"] <- length(y) / (2 * tab[["Normal"]])
  w[y == "Tumor"]  <- length(y) / (2 * tab[["Tumor"]])
  w
}

scale_train_test <- function(train_x, test_x) {
  train_x <- as.matrix(train_x); test_x <- as.matrix(test_x)
  means <- colMeans(train_x, na.rm = TRUE)
  sds <- apply(train_x, 2, sd, na.rm = TRUE)
  sds[is.na(sds) | sds == 0] <- 1
  list(train = sweep(sweep(train_x, 2, means, "-"), 2, sds, "/"),
       test  = sweep(sweep(test_x,  2, means, "-"), 2, sds, "/"),
       means = means, sds = sds)
}

# FIX: library size must come from the FULL gene set, not from the DEG subset.
# Per-sample normalisation uses no labels, so computing it once is leakage-free.
log_cpm_full <- function(counts) {
  counts <- as.matrix(counts)
  lib <- colSums(counts)
  lib[lib <= 0] <- 1
  log2(sweep(counts, 2, lib, "/") * 1e6 + 1)
}

save_gg <- function(p, name, w = 8, h = 6) {
  ggsave(file.path(OUTDIR, paste0(name, ".pdf")), p, width = w, height = h)
  ggsave(file.path(OUTDIR, paste0(name, ".png")), p, width = w, height = h, dpi = 300)
}

############################ 4. LOAD DATA
ena1 <- read_count_file(ena_files[1])
ena2 <- read_count_file(ena_files[2])
ena3 <- read_count_file(ena_files[3])
cat("ENA1/2/3 samples:", ncol(ena1), ncol(ena2), ncol(ena3), "\n")

common_ena <- Reduce(intersect, list(rownames(ena1), rownames(ena2), rownames(ena3)))
if (length(common_ena) < 1000) stop("Only ", length(common_ena), " genes common across ENA files.")

ENA <- cbind(ena1[common_ena, , drop = FALSE],
             ena2[common_ena, , drop = FALSE],
             ena3[common_ena, , drop = FALSE])
ena_cohort <- rep(c("ENA1", "ENA2", "ENA3"), c(ncol(ena1), ncol(ena2), ncol(ena3)))
colnames(ENA) <- make.unique(trimws(colnames(ENA)))
names(ena_cohort) <- colnames(ENA)

TCGA <- read_count_file(tcga_file)
cat("TCGA samples:", ncol(TCGA), "\n")

# Ensembl -> gene symbol map (TCGA file has a gene_name column) used for labels in tables/figures
tcga_hdr <- names(read.delim(tcga_file, nrows = 1, check.names = FALSE))
symbol_map <- character(0)
if (length(tcga_hdr) > 2 && tolower(tcga_hdr[2]) %in% c("gene_name", "symbol", "gene_symbol")) {
  ann <- read.delim(tcga_file, colClasses = c("character", "character", rep("NULL", length(tcga_hdr) - 2)),
                    check.names = FALSE)
  ids <- clean_gene_names(ann[[1]])
  ok <- !is.na(ids) & !duplicated(ids) & nzchar(ann[[2]])
  symbol_map <- setNames(ann[[2]][ok], ids[ok])
}
sym <- function(g) { out <- unname(symbol_map[g]); out[is.na(out)] <- g[is.na(out)]; out }

metadata <- read.csv(metadata_file, header = TRUE, check.names = FALSE,
                     stringsAsFactors = FALSE)
cat("Metadata columns:", paste(colnames(metadata), collapse = ", "), "\n")

sample_col <- find_metadata_column(metadata, c("sample", "Sample_ID", "sampleid", "sample.id",
                                               "sample_name", "sample_name_id", "id"))
condition_col <- find_metadata_column(metadata, c("condition", "group", "status",
                                                  "sample_type", "sample.type", "type"))
if (is.null(sample_col)) stop("Cannot identify sample-ID column in metadata.")
if (is.null(condition_col)) stop("Cannot identify condition column in metadata.")

meta <- metadata
meta$SampleID <- make.unique(trimws(as.character(meta[[sample_col]])))
meta$Condition <- normalize_condition(meta[[condition_col]])

# ---- ENA metadata
ena_match <- match_samples(colnames(ENA), meta$SampleID)
ena_keep <- !is.na(ena_match)
if (sum(ena_keep) < 20) stop("Only ", sum(ena_keep), " ENA samples matched metadata.")
ENA <- ENA[, ena_keep, drop = FALSE]
ena_meta <- meta[ena_match[ena_keep], , drop = FALSE]
rownames(ena_meta) <- colnames(ENA)
ena_meta$Cohort <- ena_cohort[colnames(ENA)]
if (anyNA(ena_meta$Condition)) {
  stop("Unrecognised ENA conditions. Raw values: ",
       paste(unique(meta[[condition_col]][ena_match[ena_keep]][is.na(ena_meta$Condition)]),
             collapse = ", "))
}
ena_meta$Condition <- factor(ena_meta$Condition, levels = c("Normal", "Tumor"))

# ---- TCGA metadata
tcga_match <- match_samples(colnames(TCGA), meta$SampleID)
tcga_keep <- !is.na(tcga_match)
if (sum(tcga_keep) < 20) stop("Only ", sum(tcga_keep), " TCGA samples matched metadata.")
TCGA <- TCGA[, tcga_keep, drop = FALSE]
tcga_meta <- meta[tcga_match[tcga_keep], , drop = FALSE]
rownames(tcga_meta) <- colnames(TCGA)
if (anyNA(tcga_meta$Condition)) {
  stop("Unrecognised TCGA conditions. Raw values: ",
       paste(unique(meta[[condition_col]][tcga_match[tcga_keep]][is.na(tcga_meta$Condition)]),
             collapse = ", "))
}
tcga_meta$Condition <- factor(tcga_meta$Condition, levels = c("Normal", "Tumor"))

cat("\nENA:\n");  print(table(ena_meta$Cohort, ena_meta$Condition))
cat("\nTCGA:\n"); print(table(tcga_meta$Condition))

write.csv(data.frame(Dataset = c("ENA", "TCGA"),
                     Samples = c(ncol(ENA), ncol(TCGA)),
                     Normal = c(sum(ena_meta$Condition == "Normal"), sum(tcga_meta$Condition == "Normal")),
                     Tumor = c(sum(ena_meta$Condition == "Tumor"), sum(tcga_meta$Condition == "Tumor"))),
          file.path(OUTDIR, "Dataset_Summary.csv"), row.names = FALSE)

# ---- common genes
common_genes <- intersect(rownames(ENA), rownames(TCGA))
if (length(common_genes) < 1000) stop("Only ", length(common_genes), " genes shared by ENA and TCGA.")
ENA <- ENA[common_genes, , drop = FALSE]
TCGA <- TCGA[common_genes, , drop = FALSE]
cat("Common ENA-TCGA genes:", length(common_genes), "\n")

y <- ena_meta$Condition
cohort <- ena_meta$Cohort

# Per-sample normalisation (label-free) computed once on the full gene set
ENA_lcpm <- log_cpm_full(ENA)
TCGA_lcpm <- log_cpm_full(TCGA)

############################ 5. CORE FUNCTIONS (defined ONCE, before use)

filter_low_expression <- function(counts, labels, min_cpm = 1) {
  cpm <- sweep(counts, 2, pmax(colSums(counts), 1), "/") * 1e6
  min_n <- max(2, min(table(labels)))   # expressed in at least the smaller class size
  rownames(counts)[rowSums(cpm >= min_cpm) >= min_n]
}

run_deseq <- function(dds) {
  out <- tryCatch(suppressMessages(DESeq(dds, quiet = TRUE)), error = function(e) NULL)
  if (is.null(out)) {
    out <- tryCatch(suppressMessages(DESeq(dds, fitType = "mean", sfType = "poscounts",
                                           quiet = TRUE)),
                    error = function(e) {
                      message("  DESeq2 failed: ", conditionMessage(e)); NULL })
  }
  out
}

get_training_degs <- function(counts, labels, cohort, verbose = FALSE) {
  counts <- as.matrix(counts)
  labels <- factor(labels, levels = c("Normal", "Tumor"))
  ok <- !is.na(labels)
  counts <- counts[, ok, drop = FALSE]
  labels <- droplevels(labels[ok])
  cohort <- as.character(cohort[ok])
  if (nlevels(labels) < 2 || any(table(labels) < 2)) {
    message("  DEG: need >=2 samples in both classes."); return(character(0))
  }
  keep <- filter_low_expression(counts, labels)
  if (length(keep) < 20) { message("  DEG: too few genes after filtering."); return(character(0)) }
  counts <- round(counts[keep, , drop = FALSE])
  counts <- counts[rowSums(counts) > 0, , drop = FALSE]

  coldata <- data.frame(condition = labels, cohort = factor(cohort), row.names = colnames(counts))
  design_formula <- if (nlevels(coldata$cohort) >= 2) ~ cohort + condition else ~ condition
  # FIX: fall back if cohort is confounded with condition (not full rank)
  mm <- model.matrix(design_formula, coldata)
  if (qr(mm)$rank < ncol(mm)) design_formula <- ~ condition
  if (verbose) cat("  DESeq2 design:", deparse(design_formula), "| genes:", nrow(counts), "\n")

  dds <- tryCatch(DESeqDataSetFromMatrix(counts, coldata, design_formula),
                  error = function(e) { message("  DESeq2 object error: ", conditionMessage(e)); NULL })
  if (is.null(dds)) return(character(0))
  dds <- run_deseq(dds)
  if (is.null(dds)) return(character(0))

  res <- tryCatch(results(dds, contrast = c("condition", "Tumor", "Normal"), alpha = PADJ_CUTOFF),
                  error = function(e) { message("  results() error: ", conditionMessage(e)); NULL })
  if (is.null(res)) return(character(0))
  res <- as.data.frame(res)
  res$gene <- rownames(res)
  res <- res[!is.na(res$padj) & !is.na(res$log2FoldChange), , drop = FALSE]
  res <- res[order(res$padj, -abs(res$log2FoldChange)), , drop = FALSE]  # FIX: was ascending |LFC|

  degs <- res$gene[res$padj < PADJ_CUTOFF & abs(res$log2FoldChange) > LOG2FC_CUTOFF]
  if (length(degs) < MIN_DEG) degs <- head(res$gene, 100)
  if (length(degs) > MAX_DEG_FOR_ML) degs <- head(degs, MAX_DEG_FOR_ML)  # already ranked
  if (verbose) cat("  DEGs returned:", length(degs), "\n")
  out <- unique(degs)
  attr(out, "res") <- res      # full training-only DESeq2 table (used for volcano plot)
  out
}

fit_elastic_net <- function(x, y, alpha_value) {
  x <- as.matrix(x)
  y <- factor(y, levels = c("Normal", "Tumor"))
  if (nrow(x) != length(y)) stop("x/y row mismatch.")
  if (anyNA(x) || anyNA(y)) stop("NA in Elastic Net input.")
  if (nlevels(droplevels(y)) < 2) stop("Elastic Net needs both classes.")
  x <- x[, apply(x, 2, sd) > 1e-8, drop = FALSE]
  if (ncol(x) == 0) stop("No variable genes remain.")
  if (is.na(alpha_value) || alpha_value < 0 || alpha_value > 1) alpha_value <- 0.5
  nfolds_use <- min(5, max(3, min(table(y))))
  set.seed(1000)
  fit <- cv.glmnet(x = x, y = y, family = "binomial", alpha = alpha_value,
                   type.measure = "auc", nfolds = nfolds_use,
                   weights = get_class_weights(y), standardize = FALSE, maxit = 1e5)
  fit$training_columns <- colnames(x)
  fit
}

# FIX: one consistent selector (same lambda for selection AND ranking, no zero-coef padding)
select_enet_genes <- function(fit) {
  get_nz <- function(s) {
    cf <- as.matrix(coef(fit, s = s))
    cf <- cf[rownames(cf) != "(Intercept)", , drop = FALSE]
    d <- data.frame(Gene = rownames(cf), Coefficient = cf[, 1])
    d <- d[d$Coefficient != 0, , drop = FALSE]
    d$AbsCoefficient <- abs(d$Coefficient)
    d[order(d$AbsCoefficient, decreasing = TRUE), , drop = FALSE]
  }
  lam <- "lambda.1se"
  d <- get_nz(lam)
  if (nrow(d) < MIN_ENET_GENES) { lam <- "lambda.min"; d <- get_nz(lam) }
  if (nrow(d) == 0) stop("Elastic Net selected no genes.")
  d <- head(d, MAX_ENET_GENES)
  d$Lambda <- lam
  rownames(d) <- NULL
  d
}

tune_elastic_net <- function(counts, lcpm, labels, cohort, inner_folds = 3) {
  labels <- factor(labels, levels = c("Normal", "Tumor"))
  set.seed(SEED)
  fold_id <- caret::createFolds(labels, k = inner_folds, list = FALSE)
  auc_mat <- matrix(NA_real_, nrow = length(ALPHA_GRID), ncol = inner_folds)

  for (k in seq_len(inner_folds)) {
    val_idx <- which(fold_id == k)
    tr_idx <- setdiff(seq_along(labels), val_idx)
    if (nlevels(droplevels(labels[tr_idx])) < 2 || nlevels(droplevels(labels[val_idx])) < 2) next

    # DEGs depend only on the inner-train fold -> compute ONCE per fold (not per alpha)
    genes <- get_training_degs(counts[, tr_idx, drop = FALSE], labels[tr_idx], cohort[tr_idx])
    if (length(genes) < 2) next
    sc <- scale_train_test(t(lcpm[genes, tr_idx, drop = FALSE]),
                           t(lcpm[genes, val_idx, drop = FALSE]))

    for (a in seq_along(ALPHA_GRID)) {
      fit <- tryCatch(fit_elastic_net(sc$train, labels[tr_idx], ALPHA_GRID[a]),
                      error = function(e) NULL)
      if (is.null(fit)) next
      prob <- tryCatch(as.numeric(predict(fit, newx = sc$test[, fit$training_columns, drop = FALSE],
                                          s = "lambda.1se", type = "response")),
                       error = function(e) rep(NA_real_, length(val_idx)))
      auc_mat[a, k] <- safe_auc(labels[val_idx], prob)
    }
  }
  results <- data.frame(Alpha = ALPHA_GRID,
                        Mean_AUC = rowMeans(auc_mat, na.rm = TRUE),
                        SD_AUC = apply(auc_mat, 1, sd, na.rm = TRUE))
  if (all(is.nan(results$Mean_AUC) | is.na(results$Mean_AUC))) {
    warning("Alpha tuning produced no valid AUC; using alpha = 0.5.")
    return(list(alpha = 0.5, results = results))
  }
  list(alpha = results$Alpha[which.max(results$Mean_AUC)], results = results)
}

train_models <- function(x, y) {
  x <- as.matrix(x)
  y <- factor(y, levels = c("Normal", "Tumor"))
  w <- get_class_weights(y)
  models <- list()

  # ---- Logistic regression (glm, ridge fallback if rank-deficient / failed)
  safe_names <- make.names(colnames(x), unique = TRUE)
  df <- as.data.frame(x); colnames(df) <- safe_names; df$Class <- y
  lr <- if (LR_PENALISED) NULL else
    tryCatch(suppressWarnings(glm(Class ~ ., data = df, family = binomial(), weights = w)),
             error = function(e) NULL)
  if (!is.null(lr) && !anyNA(coef(lr))) {
    models$Logistic_Regression <- list(type = "glm", model = lr, names = safe_names)
  } else {
    if (!LR_PENALISED) message("  glm failed/rank-deficient -> using ridge logistic regression.")
    ridge <- cv.glmnet(x, y, family = "binomial", alpha = 0, weights = w,
                       standardize = FALSE, nfolds = 5)
    models$Logistic_Regression <- list(type = "ridge", model = ridge)
  }

  # ---- Random forest
  set.seed(SEED)
  models$Random_Forest <- list(type = "rf",
                               model = randomForest(x = x, y = y, ntree = N_RF, importance = TRUE))

  # ---- SVM
  set.seed(SEED)
  models$SVM <- list(type = "svm",
                     model = e1071::svm(x = x, y = y, kernel = "radial", probability = TRUE,
                                        class.weights = c(Normal = 1, Tumor = 1)))

  # ---- XGBoost: classic xgb.train API (works on xgboost 1.x, 2.x and 3.x;
  #      avoids the changed xgboost() interface that returns classes/different args)
  y_xgb <- as.numeric(y == "Tumor")
  spw <- if (sum(y_xgb == 1) > 0) sum(y_xgb == 0) / sum(y_xgb == 1) else 1
  dtrain <- xgb.DMatrix(data = x, label = y_xgb)
  params <- list(objective = "binary:logistic", eval_metric = "auc", max_depth = 3,
                 eta = 0.05, subsample = 0.8, colsample_bytree = 0.8,
                 min_child_weight = 1, scale_pos_weight = spw)
  set.seed(SEED)
  models$XGBoost <- list(type = "xgb",
                         model = xgb.train(params = params, data = dtrain,
                                           nrounds = N_XGB_ROUNDS, verbose = 0))

  # ---- KNN (lazy learner: prediction needs training data, passed in predict_model)
  models$KNN <- list(type = "knn", model = NULL)
  models
}

predict_model <- function(obj, x_train, x_test, y_train) {
  x_test <- as.matrix(x_test)
  switch(obj$type,
    glm = {
      nd <- as.data.frame(x_test); colnames(nd) <- obj$names
      as.numeric(suppressWarnings(predict(obj$model, newdata = nd, type = "response")))
    },
    ridge = as.numeric(predict(obj$model, newx = x_test, s = "lambda.min", type = "response")),
    rf = as.numeric(predict(obj$model, newdata = x_test, type = "prob")[, "Tumor"]),
    svm = {
      p <- attr(predict(obj$model, newdata = x_test, probability = TRUE), "probabilities")
      as.numeric(p[, "Tumor"])
    },
    xgb = as.numeric(predict(obj$model, xgb.DMatrix(data = x_test))),
    knn = {
      k <- max(3, min(11, floor(sqrt(nrow(x_train)))))
      if (k %% 2 == 0) k <- k + 1
      pred <- class::knn(train = x_train, test = x_test,
                         cl = factor(y_train, levels = c("Normal", "Tumor")), k = k, prob = TRUE)
      p <- attr(pred, "prob")
      as.numeric(ifelse(pred == "Tumor", p, 1 - p))
    },
    stop("Unknown model type: ", obj$type))
}

############################ 6. NESTED CROSS-VALIDATION
cat("\n==== NESTED", N_OUTER, "-FOLD CV ====\n")

set.seed(SEED)
strata <- interaction(y, cohort, drop = TRUE)       # balance condition AND cohort across folds
outer_folds <- caret::createFolds(strata, k = N_OUTER, list = TRUE, returnTrain = FALSE)

outer_metrics <- list(); outer_predictions <- list(); selected_gene_list <- list()

for (k in seq_len(N_OUTER)) {
  cat("\n---- OUTER FOLD", k, "of", N_OUTER, "----\n")
  val_idx <- outer_folds[[k]]
  tr_idx <- setdiff(seq_len(ncol(ENA)), val_idx)
  tr_counts <- ENA[, tr_idx, drop = FALSE]
  tr_y <- y[tr_idx]; va_y <- y[val_idx]; tr_coh <- cohort[tr_idx]

  if (nlevels(droplevels(va_y)) < 2) stop("Outer fold ", k, " validation set has one class only.")
  cat("Train N/T:", sum(tr_y == "Normal"), sum(tr_y == "Tumor"),
      "| Val N/T:", sum(va_y == "Normal"), sum(va_y == "Tumor"), "\n")

  cat("Tuning alpha (inner CV)...\n")
  tuning <- tune_elastic_net(tr_counts, ENA_lcpm[, tr_idx, drop = FALSE], tr_y, tr_coh, N_INNER)
  write.csv(tuning$results, file.path(OUTDIR, paste0("OuterFold_", k, "_ElasticNet_Tuning.csv")),
            row.names = FALSE)
  cat("Selected alpha:", tuning$alpha, "\n")

  cat("DEG selection (outer-train only)...\n")
  degs <- get_training_degs(tr_counts, tr_y, tr_coh, verbose = TRUE)
  if (length(degs) < 2) stop("Too few DEGs in outer fold ", k, ". See messages above.")

  sc <- scale_train_test(t(ENA_lcpm[degs, tr_idx, drop = FALSE]),
                         t(ENA_lcpm[degs, val_idx, drop = FALSE]))
  enet <- fit_elastic_net(sc$train, tr_y, tuning$alpha)
  sel <- select_enet_genes(enet)
  genes <- sel$Gene
  cat("Selected genes:", length(genes), "(", sel$Lambda[1], ")\n")
  write.csv(sel, file.path(OUTDIR, paste0("OuterFold_", k, "_Selected_Genes.csv")), row.names = FALSE)
  selected_gene_list[[k]] <- genes

  train_x <- sc$train[, genes, drop = FALSE]
  val_x <- sc$test[, genes, drop = FALSE]

  cat("Training 5 models...\n")
  models <- train_models(train_x, tr_y)

  for (mn in names(models)) {
    prob <- predict_model(models[[mn]], train_x, val_x, tr_y)
    m <- calculate_metrics(va_y, prob)
    m$Model <- mn; m$Outer_Fold <- k; m$N_Features <- length(genes)
    outer_metrics[[length(outer_metrics) + 1]] <- m
    outer_predictions[[length(outer_predictions) + 1]] <- data.frame(
      Sample = colnames(ENA)[val_idx], Actual = as.character(va_y), Probability = prob,
      Prediction = ifelse(prob >= 0.5, "Tumor", "Normal"), Model = mn, Outer_Fold = k)
    cat(sprintf("  %-20s AUC=%.3f\n", mn, m$AUC))
  }
}

cv_metrics <- do.call(rbind, outer_metrics)
cv_predictions <- do.call(rbind, outer_predictions)
write.csv(cv_metrics, file.path(OUTDIR, "Nested_CV_Metrics_All_Folds.csv"), row.names = FALSE)
write.csv(cv_predictions, file.path(OUTDIR, "Nested_CV_Predictions.csv"), row.names = FALSE)

cv_summary <- cv_metrics %>%
  group_by(Model) %>%
  summarise(Mean_Accuracy = mean(Accuracy, na.rm = TRUE), Mean_Balanced_Accuracy = mean(Balanced_Accuracy, na.rm = TRUE), SD_Accuracy = sd(Accuracy, na.rm = TRUE),
            Mean_Sensitivity = mean(Sensitivity, na.rm = TRUE), SD_Sensitivity = sd(Sensitivity, na.rm = TRUE),
            Mean_Specificity = mean(Specificity, na.rm = TRUE), SD_Specificity = sd(Specificity, na.rm = TRUE),
            Mean_Precision = mean(Precision, na.rm = TRUE), Mean_F1 = mean(F1, na.rm = TRUE),
            Mean_AUC = mean(AUC, na.rm = TRUE), SD_AUC = sd(AUC, na.rm = TRUE), .groups = "drop")
write.csv(cv_summary, file.path(OUTDIR, "Nested_CV_Model_Summary.csv"), row.names = FALSE)
print(as.data.frame(cv_summary))

# ---- gene stability
gf <- sort(table(unlist(selected_gene_list)), decreasing = TRUE)
gene_frequency_df <- data.frame(Gene = names(gf), Frequency = as.integer(gf),
                                Percentage = as.integer(gf) / N_OUTER * 100)
gene_frequency_df$Symbol <- sym(gene_frequency_df$Gene)
write.csv(gene_frequency_df, file.path(OUTDIR, "ElasticNet_Gene_Selection_Frequency.csv"), row.names = FALSE)
stable_genes <- gene_frequency_df$Gene[gene_frequency_df$Frequency >= ceiling(N_OUTER * 0.6)]
if (length(stable_genes) == 0) stable_genes <- head(gene_frequency_df$Gene, 20)
write.csv(gene_frequency_df[gene_frequency_df$Gene %in% stable_genes, ],
          file.path(OUTDIR, "Stable_Final_Candidate_Genes.csv"), row.names = FALSE)
cat("\nStable genes:", paste(stable_genes, collapse = ", "), "\n")

############################ 6b. LEAVE-ONE-COHORT-OUT (LOCO) VALIDATION
# Cross-cohort generalisation inside the development data: train on two ENA cohorts
# (DEG selection, alpha tuning, Elastic Net, 5 models - all training-only), test on the third.
cat("\n==== LEAVE-ONE-COHORT-OUT VALIDATION ====\n")
run_fold <- function(tr_idx, val_idx, tag) {
  tr_counts <- ENA[, tr_idx, drop = FALSE]
  tr_y <- y[tr_idx]; va_y <- y[val_idx]; tr_coh <- cohort[tr_idx]
  tuning <- tune_elastic_net(tr_counts, ENA_lcpm[, tr_idx, drop = FALSE], tr_y, tr_coh, N_INNER)
  degs <- as.character(get_training_degs(tr_counts, tr_y, tr_coh, verbose = TRUE))
  if (length(degs) < 2) stop("too few DEGs")
  sc <- scale_train_test(t(ENA_lcpm[degs, tr_idx, drop = FALSE]), t(ENA_lcpm[degs, val_idx, drop = FALSE]))
  enet <- fit_elastic_net(sc$train, tr_y, tuning$alpha)
  genes <- select_enet_genes(enet)$Gene
  tx <- sc$train[, genes, drop = FALSE]; vx <- sc$test[, genes, drop = FALSE]
  models <- train_models(tx, tr_y)
  ms <- list(); ps <- list()
  for (mn in names(models)) {
    prob <- predict_model(models[[mn]], tx, vx, tr_y)
    m <- calculate_metrics(va_y, prob)
    m$Model <- mn; m$Held_Out <- tag; m$N_Features <- length(genes); m$N_Test <- length(val_idx)
    ms[[mn]] <- m
    ps[[mn]] <- data.frame(Sample = colnames(ENA)[val_idx], Actual = as.character(va_y), Probability = prob,
                           Prediction = ifelse(prob >= 0.5, "Tumor", "Normal"), Model = mn, Held_Out = tag)
    cat(sprintf("  %-20s AUC=%.3f\n", mn, m$AUC))
  }
  list(metrics = do.call(rbind, ms), preds = do.call(rbind, ps), genes = genes)
}

loco_ok <- list()
for (cc in sort(unique(cohort))) {
  cat("\n---- Held-out cohort:", cc, "(train on the other cohorts) ----\n")
  r <- tryCatch(run_fold(which(cohort != cc), which(cohort == cc), cc),
                error = function(e) { message("LOCO failed for ", cc, ": ", conditionMessage(e)); NULL })
  if (!is.null(r)) loco_ok[[cc]] <- r
}
loco_freq <- table(character(0))
if (length(loco_ok) > 0) {
  loco_metrics <- do.call(rbind, lapply(loco_ok, `[[`, "metrics")); rownames(loco_metrics) <- NULL
  loco_preds <- do.call(rbind, lapply(loco_ok, `[[`, "preds")); rownames(loco_preds) <- NULL
  write.csv(loco_metrics, file.path(OUTDIR, "LOCO_Metrics_By_Cohort.csv"), row.names = FALSE)
  write.csv(loco_preds, file.path(OUTDIR, "LOCO_Predictions.csv"), row.names = FALSE)
  loco_summary <- loco_metrics %>% group_by(Model) %>%
    summarise(Mean_AUC = mean(AUC, na.rm = TRUE), Min_AUC = min(AUC, na.rm = TRUE),
              Mean_Sensitivity = mean(Sensitivity, na.rm = TRUE), Mean_Specificity = mean(Specificity, na.rm = TRUE),
              Mean_Balanced_Accuracy = mean(Balanced_Accuracy, na.rm = TRUE), .groups = "drop")
  write.csv(loco_summary, file.path(OUTDIR, "LOCO_Summary.csv"), row.names = FALSE)
  print(as.data.frame(loco_summary))
  loco_freq <- table(unlist(lapply(loco_ok, `[[`, "genes")))
  loco_genes_df <- data.frame(Gene = names(loco_freq), Symbol = sym(names(loco_freq)),
                              LOCO_Folds = as.integer(loco_freq), stringsAsFactors = FALSE)
  loco_genes_df <- loco_genes_df[order(-loco_genes_df$LOCO_Folds), ]
  write.csv(loco_genes_df, file.path(OUTDIR, "LOCO_Selected_Genes.csv"), row.names = FALSE)
} else {
  warning("All LOCO folds failed; see messages above.")
}

############################ 7. FINAL DEG ANALYSIS (ALL ENA, cohort-adjusted DESeq2)
# This is the reportable tumor-vs-normal DEG result. TCGA is not used here.
cat("\n==== FINAL DEG ANALYSIS ON ALL ENA ====\n")
final_degs <- get_training_degs(ENA, y, cohort, verbose = TRUE)
if (length(final_degs) < 2) stop("Final DEG analysis failed. See messages above.")
deg_all <- attr(final_degs, "res")
names(deg_all)[names(deg_all) == "gene"] <- "Gene"
deg_all$Symbol <- sym(deg_all$Gene)
deg_all$Status <- ifelse(deg_all$padj < PADJ_CUTOFF & deg_all$log2FoldChange > LOG2FC_CUTOFF, "Up in tumor",
                  ifelse(deg_all$padj < PADJ_CUTOFF & deg_all$log2FoldChange < -LOG2FC_CUTOFF,
                         "Down in tumor", "Not significant"))
deg_all <- deg_all[, c("Gene", "Symbol", "baseMean", "log2FoldChange", "lfcSE", "stat",
                       "pvalue", "padj", "Status")]
rownames(deg_all) <- NULL
deg_sig <- deg_all[deg_all$Status != "Not significant", , drop = FALSE]
write.csv(deg_all, file.path(OUTDIR, "DEG_Full_Results_ENA.csv"), row.names = FALSE)
write.csv(deg_sig, file.path(OUTDIR, "DEG_Significant_ENA.csv"), row.names = FALSE)
cat("Genes tested:", nrow(deg_all), "| significant DEGs:", nrow(deg_sig),
    "| up:", sum(deg_all$Status == "Up in tumor"),
    "| down:", sum(deg_all$Status == "Down in tumor"), "\n")

############################ 8. SINGLE-FIT ELASTIC NET (comparison only)
cat("\n==== SINGLE-FIT ELASTIC NET ON ALL ENA (for comparison) ====\n")
degs_vec <- as.character(final_degs)
final_tuning <- tune_elastic_net(ENA, ENA_lcpm, y, cohort, N_INNER)
final_alpha <- final_tuning$alpha
write.csv(final_tuning$results, file.path(OUTDIR, "Final_ElasticNet_Tuning.csv"), row.names = FALSE)
en_sc <- scale_train_test(t(ENA_lcpm[degs_vec, , drop = FALSE]), t(ENA_lcpm[degs_vec, , drop = FALSE]))  # ENA only; TCGA not touched here
final_enet <- fit_elastic_net(en_sc$train, y, final_alpha)
enet_single <- select_enet_genes(final_enet)
enet_single$Symbol <- sym(enet_single$Gene)
write.csv(enet_single, file.path(OUTDIR, "ElasticNet_SingleFit_All_ENA_Genes.csv"), row.names = FALSE)
cat("Single-fit alpha:", final_alpha, "| genes:", paste(enet_single$Symbol, collapse = ", "), "\n")

############################ 9. FINAL BIOMARKER PANEL
# ONE definition used everywhere (models, TCGA test, tables, figures):
# genes selected in >= 60% of the nested-CV outer folds (stability selection, ENA only).
panel_genes <- head(stable_genes, MAX_ENET_GENES)
panel_genes <- panel_genes[panel_genes %in% rownames(ENA_lcpm) & panel_genes %in% rownames(TCGA_lcpm)]
if (length(panel_genes) < 2) {
  warning("Fewer than 2 stable genes; using single-fit Elastic Net genes as the panel.")
  panel_genes <- enet_single$Gene
}

gene_stats <- function(lcpm, cond, genes) {
  cond <- factor(as.character(cond), levels = c("Normal", "Tumor"))
  do.call(rbind, lapply(genes, function(g) {
    x <- lcpm[g, ]
    xt <- x[cond == "Tumor"]; xn <- x[cond == "Normal"]
    w <- suppressWarnings(wilcox.test(xt, xn))
    auc_raw <- as.numeric(w$statistic) / (length(xt) * length(xn))
    data.frame(Gene = g, Mean_Normal = mean(xn), Mean_Tumor = mean(xt),
               Delta_log2CPM = mean(xt) - mean(xn), Wilcoxon_P = w$p.value,
               AUC = max(auc_raw, 1 - auc_raw))
  }))
}
st_ena <- gene_stats(ENA_lcpm, y, panel_genes)
st_tcga <- gene_stats(TCGA_lcpm, tcga_meta$Condition, panel_genes)
fr <- gene_frequency_df[match(panel_genes, gene_frequency_df$Gene), ]

panel_tbl <- data.frame(
  Gene = panel_genes, Symbol = sym(panel_genes),
  Selected_in_Outer_Folds = fr$Frequency, Percent_of_Folds = fr$Percentage,
  In_SingleFit_ElasticNet = panel_genes %in% enet_single$Gene,
  SingleFit_EN_Coefficient = enet_single$Coefficient[match(panel_genes, enet_single$Gene)],
  DESeq2_log2FC = deg_all$log2FoldChange[match(panel_genes, deg_all$Gene)],
  DESeq2_padj = deg_all$padj[match(panel_genes, deg_all$Gene)],
  ENA_Mean_Normal = st_ena$Mean_Normal, ENA_Mean_Tumor = st_ena$Mean_Tumor,
  ENA_Delta_log2CPM = st_ena$Delta_log2CPM, ENA_Wilcoxon_P = st_ena$Wilcoxon_P, ENA_AUC = st_ena$AUC,
  TCGA_Mean_Normal = st_tcga$Mean_Normal, TCGA_Mean_Tumor = st_tcga$Mean_Tumor,
  TCGA_Delta_log2CPM = st_tcga$Delta_log2CPM, TCGA_Wilcoxon_P = st_tcga$Wilcoxon_P, TCGA_AUC = st_tcga$AUC,
  Direction = ifelse(st_ena$Delta_log2CPM > 0, "Up in tumor", "Down in tumor"),
  Concordant_ENA_TCGA = sign(st_ena$Delta_log2CPM) == sign(st_tcga$Delta_log2CPM),
  stringsAsFactors = FALSE)
panel_tbl <- panel_tbl[order(-panel_tbl$Selected_in_Outer_Folds, -panel_tbl$ENA_AUC), , drop = FALSE]
rownames(panel_tbl) <- NULL
panel_tbl$Selected_in_LOCO_Folds <- as.integer(ifelse(panel_tbl$Gene %in% names(loco_freq), loco_freq[panel_tbl$Gene], 0L))
panel_genes <- panel_tbl$Gene
panel_tbl$Label <- make.unique(panel_tbl$Symbol)
lab_map <- setNames(panel_tbl$Label, panel_tbl$Gene)
write.csv(panel_tbl, file.path(OUTDIR, "FINAL_BIOMARKER_PANEL.csv"), row.names = FALSE)
writeLines(paste0(seq_len(nrow(panel_tbl)), ". ", panel_tbl$Symbol, " (", panel_tbl$Gene, ") - ",
                  panel_tbl$Direction, " - selected in ", panel_tbl$Selected_in_Outer_Folds, "/", N_OUTER,
                  " folds - ENA AUC ", sprintf("%.3f", panel_tbl$ENA_AUC),
                  " - TCGA AUC ", sprintf("%.3f", panel_tbl$TCGA_AUC)),
           file.path(OUTDIR, "FINAL_BIOMARKER_PANEL.txt"))
cat("\nFINAL BIOMARKER PANEL (", length(panel_genes), "genes):\n")
print(panel_tbl[, c("Symbol", "Gene", "Direction", "Selected_in_Outer_Folds", "ENA_AUC", "TCGA_AUC")])

############################ 10. FINAL MODELS + LOCKED TCGA TEST
cat("\n==== TRAIN ON ALL ENA (PANEL) -> LOCKED TCGA TEST ====\n")
panel_sc <- scale_train_test(t(ENA_lcpm[panel_genes, , drop = FALSE]),
                             t(TCGA_lcpm[panel_genes, , drop = FALSE]))
panel_train_x <- panel_sc$train
panel_tcga_x <- panel_sc$test

# Diagnostic: how far does each panel gene move between ENA and TCGA after ENA-based scaling?
# Large |ENA - TCGA| gaps within the same class flag genes that are platform-sensitive.
shift_tbl <- do.call(rbind, lapply(c("Normal", "Tumor"), function(cl) {
  data.frame(Gene = panel_genes, Symbol = sym(panel_genes), Class = cl,
             ENA_mean_z = colMeans(panel_train_x[y == cl, , drop = FALSE]),
             TCGA_mean_z = colMeans(panel_tcga_x[tcga_meta$Condition == cl, , drop = FALSE]),
             ENA_sd_z = apply(panel_train_x[y == cl, , drop = FALSE], 2, sd),
             TCGA_sd_z = apply(panel_tcga_x[tcga_meta$Condition == cl, , drop = FALSE], 2, sd),
             row.names = NULL)
}))
shift_tbl$Mean_shift_z <- shift_tbl$TCGA_mean_z - shift_tbl$ENA_mean_z
write.csv(shift_tbl, file.path(OUTDIR, "Panel_ENA_vs_TCGA_Shift_Diagnostic.csv"), row.names = FALSE)
final_models <- train_models(panel_train_x, y)
# Locked artefact: panel + ENA scaling + trained models. Nothing below may modify it.
tryCatch(saveRDS(list(panel_genes = panel_genes, symbols = sym(panel_genes), means = panel_sc$means,
                      sds = panel_sc$sds, models = final_models, train_x = panel_train_x, train_y = y),
                 file.path(OUTDIR, "LOCKED_MODEL_ENA_panel.rds")),
         error = function(e) message("Could not save locked model: ", conditionMessage(e)))

tcga_metrics_l <- list(); tcga_pred_l <- list()
for (mn in names(final_models)) {
  prob <- predict_model(final_models[[mn]], panel_train_x, panel_tcga_x, y)
  m <- calculate_metrics(tcga_meta$Condition, prob)
  m$Model <- mn; m$N_Features <- length(panel_genes)
  tcga_metrics_l[[mn]] <- m
  tcga_pred_l[[mn]] <- data.frame(Sample = colnames(TCGA), Actual = as.character(tcga_meta$Condition),
                                  Probability = prob,
                                  Prediction = ifelse(prob >= 0.5, "Tumor", "Normal"), Model = mn)
}
tcga_metrics <- do.call(rbind, tcga_metrics_l); rownames(tcga_metrics) <- NULL
tcga_predictions <- do.call(rbind, tcga_pred_l); rownames(tcga_predictions) <- NULL
write.csv(tcga_metrics, file.path(OUTDIR, "TCGA_EXTERNAL_TEST_METRICS.csv"), row.names = FALSE)
write.csv(tcga_predictions, file.path(OUTDIR, "TCGA_EXTERNAL_TEST_PREDICTIONS.csv"), row.names = FALSE)
print(tcga_metrics)

############################ 11. FIGURES (each wrapped so one failure does not stop the rest)
fig <- function(name, expr) {
  tryCatch({ force(expr); invisible(TRUE) },
           error = function(e) { message("FIGURE FAILED [", name, "]: ", conditionMessage(e)); invisible(FALSE) })
}
col_cond <- c(Normal = "#2166AC", Tumor = "#B2182B")
col_stat <- c("Up in tumor" = "#B2182B", "Down in tumor" = "#2166AC", "Not significant" = "grey75")

make_heatmap <- function(lcpm, cond, genes, labels, title, file, h = 6) {
  mat <- t(scale(t(lcpm[genes, , drop = FALSE])))
  mat[!is.finite(mat)] <- 0
  mat <- pmax(pmin(mat, 3), -3)
  cond <- as.character(cond)
  ord_g <- if (nrow(mat) > 1) hclust(dist(mat))$order else 1
  sgn <- sign(rowMeans(mat[, cond == "Tumor", drop = FALSE]) - rowMeans(mat[, cond == "Normal", drop = FALSE]))
  ord_s <- order(cond, colSums(mat * sgn))
  m2 <- mat[ord_g, ord_s, drop = FALSE]
  df <- as.data.frame(as.table(m2))
  df$Condition <- setNames(cond, colnames(lcpm))[as.character(df$Var2)]
  df$Gene <- factor(unname(labels[as.character(df$Var1)]), levels = unname(labels[rownames(m2)]))
  p <- ggplot(df, aes(Var2, Gene, fill = Freq)) + geom_tile() +
    facet_grid(. ~ Condition, scales = "free_x", space = "free_x") +
    scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-3, 3), name = "z-score") +
    labs(title = title, x = "Samples", y = NULL) + theme_minimal() +
    theme(axis.text.x = element_blank(), axis.ticks.x = element_blank(), panel.grid = element_blank())
  save_gg(p, file, 10, h)
}

make_boxplots <- function(lcpm, cond, genes, title, file) {
  m <- t(lcpm[genes, , drop = FALSE])
  long <- data.frame(Expression = as.vector(m),
                     Gene = factor(rep(unname(lab_map[genes]), each = nrow(m)), levels = unname(lab_map[genes])),
                     Condition = factor(rep(as.character(cond), times = length(genes)), levels = c("Normal", "Tumor")))
  nc <- min(4, length(genes)); nr <- ceiling(length(genes) / nc)
  p <- ggplot(long, aes(Condition, Expression, fill = Condition)) +
    geom_boxplot(outlier.shape = NA, alpha = 0.6) + geom_jitter(width = 0.15, size = 0.5, alpha = 0.3) +
    facet_wrap(~Gene, scales = "free_y", ncol = nc) + scale_fill_manual(values = col_cond) +
    labs(title = title, y = "log2 CPM", x = NULL) + theme_classic() + theme(legend.position = "none")
  save_gg(p, file, 3 * nc + 1, 2.6 * nr + 1)
}

roc_df <- function(actual, prob, model) {
  r <- pROC::roc(factor(actual, levels = c("Normal", "Tumor")), prob,
                 levels = c("Normal", "Tumor"), direction = "<", quiet = TRUE)
  data.frame(FPR = 1 - r$specificities, TPR = r$sensitivities, Model = model)
}

cm_table <- function(d, model) {
  cm <- as.data.frame(table(Actual = factor(d$Actual, levels = c("Tumor", "Normal")),
                            Predicted = factor(d$Prediction, levels = c("Normal", "Tumor"))))
  cm$Pct <- round(100 * cm$Freq / ave(cm$Freq, cm$Actual, FUN = sum), 1)
  cm$Model <- model
  cm
}
cm_plot <- function(cm, title) {
  ggplot(cm, aes(Predicted, Actual, fill = Pct)) + geom_tile(colour = "white") +
    geom_text(aes(label = paste0(Freq, "\n(", Pct, "%)")), size = 4.5) +
    scale_fill_gradient(low = "white", high = "#4A90D9", limits = c(0, 100), name = "% of actual") +
    labs(title = title) + theme_classic()
}

# ---- A. DATA OVERVIEW
fig("A1 sample composition", {
  cd <- data.frame(Cohort = c(as.character(cohort), rep("TCGA", ncol(TCGA))),
                   Condition = c(as.character(y), as.character(tcga_meta$Condition)))
  cd <- as.data.frame(table(cd$Cohort, cd$Condition)); names(cd) <- c("Cohort", "Condition", "N")
  save_gg(ggplot(cd, aes(Cohort, N, fill = Condition)) + geom_col() +
            geom_text(aes(label = N), position = position_stack(vjust = 0.5), colour = "white") +
            scale_fill_manual(values = col_cond) + labs(title = "Samples per cohort and condition", y = "Samples") +
            theme_classic(), "Figure_A1_Sample_Composition", 7, 5)
})
fig("A2 PCA ENA batch", {
  vv <- apply(ENA_lcpm, 1, var)
  top_var <- names(sort(vv, decreasing = TRUE))[seq_len(min(1000, length(vv)))]
  pc1 <- prcomp(t(ENA_lcpm[top_var, ]), scale. = TRUE)
  pv1 <- round(100 * pc1$sdev^2 / sum(pc1$sdev^2), 1)
  save_gg(ggplot(data.frame(PC1 = pc1$x[, 1], PC2 = pc1$x[, 2], Cohort = cohort, Condition = y),
                 aes(PC1, PC2, colour = Cohort, shape = Condition)) + geom_point(size = 2.5, alpha = 0.8) +
            labs(title = "ENA PCA (1000 most variable genes): cohort vs condition",
                 x = paste0("PC1 (", pv1[1], "%)"), y = paste0("PC2 (", pv1[2], "%)")) + theme_classic(),
          "Figure_A2_PCA_ENA_Cohort_Batch", 8, 6)
})

# ---- B. DIFFERENTIAL EXPRESSION
fig("B1 DEG counts", {
  dc <- as.data.frame(table(deg_all$Status)); names(dc) <- c("Status", "N")
  dc <- dc[dc$Status != "Not significant", ]
  save_gg(ggplot(dc, aes(Status, N, fill = Status)) + geom_col() + geom_text(aes(label = N), vjust = -0.3) +
            scale_fill_manual(values = col_stat) +
            labs(title = paste0("DEGs: padj < ", PADJ_CUTOFF, ", |log2FC| > ", LOG2FC_CUTOFF), x = NULL, y = "Genes") +
            theme_classic() + theme(legend.position = "none"), "Figure_B1_DEG_Counts", 5, 5)
})
fig("B2 volcano", {
  lab_top <- function(st, n) { d <- deg_all[deg_all$Status == st, ]; head(d[order(d$padj, -abs(d$log2FoldChange)), ], n) }
  lab_df <- rbind(lab_top("Up in tumor", 10), lab_top("Down in tumor", 10))
  lab_df <- lab_df[!lab_df$Gene %in% panel_genes, ]
  pan_df <- deg_all[deg_all$Gene %in% panel_genes, ]; pan_df$Label <- unname(lab_map[pan_df$Gene])
  deg_all$nlp <- -log10(pmax(deg_all$padj, 1e-300)); lab_df$nlp <- -log10(pmax(lab_df$padj, 1e-300))
  pan_df$nlp <- -log10(pmax(pan_df$padj, 1e-300))
  save_gg(ggplot(deg_all, aes(log2FoldChange, nlp, colour = Status)) + geom_point(alpha = 0.45, size = 0.8) +
            scale_colour_manual(values = col_stat) +
            geom_vline(xintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed") +
            geom_hline(yintercept = -log10(PADJ_CUTOFF), linetype = "dashed") +
            geom_text(data = lab_df, aes(label = Symbol), colour = "grey20", size = 2.8, vjust = -0.6, check_overlap = TRUE) +
            geom_point(data = pan_df, shape = 21, colour = "black", fill = NA, size = 3, stroke = 1) +
            geom_text(data = pan_df, aes(label = Label), colour = "black", fontface = "bold", size = 3.4, vjust = -1) +
            labs(title = "Volcano plot: Tumor vs Normal (ENA, DESeq2 adjusted for cohort)",
                 subtitle = "Circled/bold = final biomarker panel; grey labels = top 10 up and down",
                 x = "log2 fold change (Tumor / Normal)", y = "-log10 adjusted p") + theme_classic(),
          "Figure_B2_Volcano_DEGs", 9, 7)
})
fig("B3 MA plot", {
  pan_df <- deg_all[deg_all$Gene %in% panel_genes, ]; pan_df$Label <- unname(lab_map[pan_df$Gene])
  save_gg(ggplot(deg_all, aes(log10(baseMean + 1), log2FoldChange, colour = Status)) +
            geom_point(alpha = 0.4, size = 0.8) + scale_colour_manual(values = col_stat) +
            geom_hline(yintercept = c(-LOG2FC_CUTOFF, LOG2FC_CUTOFF), linetype = "dashed") +
            geom_point(data = pan_df, shape = 21, colour = "black", fill = NA, size = 3, stroke = 1) +
            geom_text(data = pan_df, aes(label = Label), colour = "black", size = 3, vjust = -1, check_overlap = TRUE) +
            labs(title = "MA plot (ENA)", x = "log10(mean normalised count + 1)", y = "log2 fold change") +
            theme_classic(), "Figure_B3_MA_Plot", 8, 6)
})
fig("B4 top DEGs", {
  sig <- deg_sig[deg_sig$baseMean >= 10, ]
  up <- sig[sig$log2FoldChange > 0, ]; up <- head(up[order(-up$log2FoldChange), ], 15)
  dn <- sig[sig$log2FoldChange < 0, ]; dn <- head(dn[order(dn$log2FoldChange), ], 15)
  td <- rbind(up, dn); td$Label <- make.unique(td$Symbol)
  save_gg(ggplot(td, aes(reorder(Label, log2FoldChange), log2FoldChange, fill = Status)) + geom_col() +
            coord_flip() + scale_fill_manual(values = col_stat) +
            labs(title = "Top 15 up- and down-regulated DEGs (by |log2FC|, baseMean >= 10)", x = NULL,
                 y = "log2 fold change") + theme_classic(), "Figure_B4_Top_DEGs", 8, 8)
})
fig("B5 DEG heatmap", {
  top50 <- head(deg_sig$Gene[order(deg_sig$padj, -abs(deg_sig$log2FoldChange))], 50)
  lab50 <- setNames(make.unique(sym(top50)), top50)
  make_heatmap(ENA_lcpm, y, top50, lab50, "Top 50 DEGs by adjusted p - ENA", "Figure_B5_Heatmap_Top50_DEGs", 10)
})

# ---- C. FEATURE SELECTION / STABILITY
fig("C1 tuning", {
  save_gg(ggplot(final_tuning$results, aes(Alpha, Mean_AUC)) +
            geom_errorbar(aes(ymin = Mean_AUC - SD_AUC, ymax = Mean_AUC + SD_AUC), width = 0.03) +
            geom_line() + geom_point(size = 3) +
            labs(title = "Elastic Net alpha tuning (inner CV, all ENA)", y = "Mean AUC") + theme_classic(),
          "Figure_C1_ElasticNet_Alpha_Tuning", 6, 5)
})
fig("C2 features per fold", {
  nf <- unique(cv_metrics[, c("Outer_Fold", "N_Features")])
  save_gg(ggplot(nf, aes(factor(Outer_Fold), N_Features)) + geom_col() +
            geom_text(aes(label = N_Features), vjust = -0.3) +
            labs(title = "Elastic Net genes selected per outer fold", x = "Outer fold", y = "Genes") +
            theme_classic(), "Figure_C2_Features_Per_Fold", 6, 5)
})
fig("C3 selection frequency", {
  tf <- head(gene_frequency_df, 30); tf$Label <- make.unique(tf$Symbol)
  tf$InPanel <- ifelse(tf$Gene %in% panel_genes, "In final panel", "Not in panel")
  save_gg(ggplot(tf, aes(reorder(Label, Frequency), Frequency, fill = InPanel)) + geom_col() + coord_flip() +
            scale_fill_manual(values = c("In final panel" = "#B2182B", "Not in panel" = "grey60")) +
            labs(title = "Elastic Net selection frequency across outer folds", x = NULL,
                 y = paste0("Folds selected (of ", N_OUTER, ")"), fill = NULL) + theme_classic(),
          "Figure_C3_Gene_Selection_Frequency", 8, 8)
})

# ---- D. FINAL BIOMARKER PANEL
fig("D1 panel heatmap ENA", make_heatmap(ENA_lcpm, y, panel_genes, lab_map,
    "Final biomarker panel - ENA (development)", "Figure_D1_Heatmap_Panel_ENA", max(3.5, 0.45 * length(panel_genes) + 2)))
fig("D2 panel heatmap TCGA", make_heatmap(TCGA_lcpm, tcga_meta$Condition, panel_genes, lab_map,
    "Final biomarker panel - TCGA (external)", "Figure_D2_Heatmap_Panel_TCGA", max(3.5, 0.45 * length(panel_genes) + 2)))
fig("D3 boxplots ENA", make_boxplots(ENA_lcpm, y, panel_genes, "Final biomarker panel - ENA", "Figure_D3_Boxplots_Panel_ENA"))
fig("D4 boxplots TCGA", make_boxplots(TCGA_lcpm, tcga_meta$Condition, panel_genes, "Final biomarker panel - TCGA",
                                      "Figure_D4_Boxplots_Panel_TCGA"))
fig("D5 concordance", {
  cc <- data.frame(Label = panel_tbl$Label, ENA = panel_tbl$ENA_Delta_log2CPM, TCGA = panel_tbl$TCGA_Delta_log2CPM)
  save_gg(ggplot(cc, aes(ENA, TCGA, label = Label)) + geom_hline(yintercept = 0, colour = "grey60") +
            geom_vline(xintercept = 0, colour = "grey60") + geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
            geom_point(size = 3) + geom_text(vjust = -0.8, size = 3.2, check_overlap = TRUE) +
            labs(title = "Biomarker effect size: ENA vs TCGA", x = "ENA: mean log2CPM (Tumor - Normal)",
                 y = "TCGA: mean log2CPM (Tumor - Normal)") + theme_classic(), "Figure_D5_Panel_ENA_vs_TCGA_Concordance", 7, 6)
})
fig("D6 single gene AUC", {
  ad <- rbind(data.frame(Label = panel_tbl$Label, Dataset = "ENA", AUC = panel_tbl$ENA_AUC),
              data.frame(Label = panel_tbl$Label, Dataset = "TCGA", AUC = panel_tbl$TCGA_AUC))
  save_gg(ggplot(ad, aes(AUC, reorder(Label, AUC), colour = Dataset)) + geom_point(size = 3.5, alpha = 0.85) +
            xlim(0.5, 1) + labs(title = "Single-gene discrimination (Tumor vs Normal)", x = "AUC", y = NULL) +
            theme_classic(), "Figure_D6_Panel_SingleGene_AUC", 7, 5)
})
fig("D7 correlation", {
  cmat <- cor(t(ENA_lcpm[panel_genes, , drop = FALSE]))
  dimnames(cmat) <- list(unname(lab_map[panel_genes]), unname(lab_map[panel_genes]))
  cdf <- as.data.frame(as.table(cmat))
  save_gg(ggplot(cdf, aes(Var1, Var2, fill = Freq)) + geom_tile() + geom_text(aes(label = sprintf("%.2f", Freq)), size = 3) +
            scale_fill_gradient2(low = "#2166AC", mid = "white", high = "#B2182B", limits = c(-1, 1), name = "r") +
            labs(title = "Panel gene correlation (ENA)", x = NULL, y = NULL) + theme_minimal() +
            theme(axis.text.x = element_text(angle = 45, hjust = 1)), "Figure_D7_Panel_Correlation_ENA", 7, 6)
})
fig("D8 PCA panel", {
  comb <- cbind(ENA_lcpm[panel_genes, , drop = FALSE], TCGA_lcpm[panel_genes, , drop = FALSE])
  pc2 <- prcomp(t(comb), scale. = TRUE); pv2 <- round(100 * pc2$sdev^2 / sum(pc2$sdev^2), 1)
  save_gg(ggplot(data.frame(PC1 = pc2$x[, 1], PC2 = pc2$x[, 2],
                            Dataset = rep(c("ENA", "TCGA"), c(ncol(ENA), ncol(TCGA))),
                            Condition = c(as.character(y), as.character(tcga_meta$Condition))),
                 aes(PC1, PC2, colour = Condition, shape = Dataset)) + geom_point(size = 2, alpha = 0.7) +
            scale_colour_manual(values = col_cond) +
            labs(title = "PCA of final panel: ENA + TCGA", x = paste0("PC1 (", pv2[1], "%)"),
                 y = paste0("PC2 (", pv2[2], "%)")) + theme_classic(), "Figure_D8_PCA_Panel_ENA_TCGA", 8, 6)
})

# ---- E. MODEL PERFORMANCE
fig("E1 nested CV boxplot", {
  save_gg(ggplot(cv_metrics, aes(Model, AUC, fill = Model)) + geom_boxplot(alpha = 0.5, outlier.shape = NA) +
            geom_jitter(width = 0.1, size = 2) + coord_cartesian(ylim = c(min(0.5, cv_metrics$AUC, na.rm = TRUE), 1)) +
            labs(title = "Nested CV AUC across outer folds (ENA)", x = NULL) + theme_classic() +
            theme(legend.position = "none", axis.text.x = element_text(angle = 30, hjust = 1)),
          "Figure_E1_NestedCV_AUC_by_Model", 8, 6)
})
fig("E2 nested CV ROC", {
  rd <- do.call(rbind, lapply(unique(cv_predictions$Model), function(mn) {
    d <- cv_predictions[cv_predictions$Model == mn, ]; roc_df(d$Actual, d$Probability, mn) }))
  save_gg(ggplot(rd, aes(FPR, TPR, colour = Model)) + geom_path(linewidth = 0.9) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
            labs(title = "Nested CV ROC (pooled out-of-fold predictions, ENA)", x = "1 - Specificity",
                 y = "Sensitivity") + theme_classic(), "Figure_E2_NestedCV_ROC", 8, 6)
})
fig("E3 internal vs external", {
  comp <- rbind(data.frame(Dataset = "ENA nested CV", Model = cv_summary$Model, AUC = cv_summary$Mean_AUC),
                data.frame(Dataset = "TCGA external", Model = tcga_metrics$Model, AUC = tcga_metrics$AUC))
  save_gg(ggplot(comp, aes(Model, AUC, fill = Dataset)) + geom_col(position = position_dodge(0.8), width = 0.7) +
            geom_text(aes(label = sprintf("%.3f", AUC)), position = position_dodge(0.8), vjust = -0.3, size = 3) +
            coord_cartesian(ylim = c(0, 1.05)) + labs(title = "Internal vs external AUC", x = NULL) +
            theme_classic() + theme(axis.text.x = element_text(angle = 30, hjust = 1)),
          "Figure_E3_Internal_vs_External_AUC", 9, 6)
})
fig("E4 TCGA ROC", {
  rd <- do.call(rbind, lapply(unique(tcga_predictions$Model), function(mn) {
    d <- tcga_predictions[tcga_predictions$Model == mn, ]; roc_df(d$Actual, d$Probability, mn) }))
  write.csv(rd, file.path(OUTDIR, "TCGA_ROC_Curve_Data.csv"), row.names = FALSE)
  save_gg(ggplot(rd, aes(FPR, TPR, colour = Model)) + geom_path(linewidth = 0.9) +
            geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
            labs(title = "External validation ROC (TCGA)", x = "1 - Specificity", y = "Sensitivity") + theme_classic(),
          "Figure_E4_TCGA_ROC_Curves", 8, 6)
})
fig("E5 TCGA AUC", {
  save_gg(ggplot(tcga_metrics, aes(reorder(Model, AUC), AUC)) + geom_col(fill = "grey40") +
            geom_text(aes(label = sprintf("%.3f", AUC)), vjust = -0.3) + ylim(0, 1.05) +
            labs(title = "External TCGA AUC", x = NULL) + theme_classic() +
            theme(axis.text.x = element_text(angle = 30, hjust = 1)), "Figure_E5_TCGA_AUC_Comparison", 8, 6)
})
fig("E6 TCGA metrics", {
  ml <- do.call(rbind, lapply(c("Sensitivity", "Specificity", "Balanced_Accuracy", "AUC"), function(m)
    data.frame(Model = tcga_metrics$Model, Metric = m, Value = tcga_metrics[[m]])))
  save_gg(ggplot(ml, aes(Metric, Value, fill = Model)) + geom_col(position = position_dodge(0.85), width = 0.8) +
            geom_text(aes(label = sprintf("%.2f", Value)), position = position_dodge(0.85), vjust = -0.3, size = 2.6) +
            coord_cartesian(ylim = c(0, 1.08)) + labs(title = "TCGA external performance (class-imbalance-aware metrics)", x = NULL) +
            theme_classic(), "Figure_E6_TCGA_Metrics_By_Model", 10, 6)
})
fig("E7 confusion matrices", {
  all_cm <- do.call(rbind, lapply(unique(tcga_predictions$Model), function(mn)
    cm_table(tcga_predictions[tcga_predictions$Model == mn, ], mn)))
  for (mn in unique(all_cm$Model)) {
    save_gg(cm_plot(all_cm[all_cm$Model == mn, ], paste("TCGA confusion matrix:", mn)),
            paste0("Figure_E7_TCGA_ConfusionMatrix_", gsub("[^A-Za-z0-9]+", "_", mn)), 6, 5)
  }
  save_gg(cm_plot(all_cm, "TCGA confusion matrices (all models)") + facet_wrap(~Model),
          "Figure_E7_TCGA_ConfusionMatrices_All", 11, 7)
})
fig("E8 probability distribution", {
  save_gg(ggplot(tcga_predictions, aes(Actual, Probability, fill = Actual)) +
            geom_violin(alpha = 0.5, scale = "width") + geom_jitter(width = 0.12, size = 0.4, alpha = 0.4) +
            geom_hline(yintercept = 0.5, linetype = "dashed") + facet_wrap(~Model) +
            scale_fill_manual(values = col_cond) +
            labs(title = "TCGA predicted tumor probability by true class", x = NULL, y = "P(Tumor)") +
            theme_classic() + theme(legend.position = "none"), "Figure_E8_TCGA_Probability_Distribution", 10, 7)
})

fig("E9 LOCO", {
  save_gg(ggplot(loco_metrics, aes(Held_Out, AUC, fill = Model)) +
            geom_col(position = position_dodge(0.85), width = 0.8) +
            geom_text(aes(label = sprintf("%.2f", AUC)), position = position_dodge(0.85), vjust = -0.3, size = 2.6) +
            coord_cartesian(ylim = c(0, 1.08)) +
            labs(title = "Leave-one-cohort-out validation (train on two ENA cohorts, test on the third)",
                 x = "Held-out cohort", y = "AUC") + theme_classic(), "Figure_E9_LOCO_AUC_By_Cohort", 9, 6)
})

# ---- F. FEATURE IMPORTANCE (panel genes only)
fig("F1 RF importance", {
  rf_imp <- importance(final_models$Random_Forest$model)
  rf_imp_df <- data.frame(Gene = rownames(rf_imp), rf_imp, row.names = NULL)
  rf_imp_df$Symbol <- unname(lab_map[rf_imp_df$Gene])
  rf_imp_df <- rf_imp_df[order(rf_imp_df$MeanDecreaseAccuracy, decreasing = TRUE), ]
  write.csv(rf_imp_df, file.path(OUTDIR, "RandomForest_Variable_Importance.csv"), row.names = FALSE)
  save_gg(ggplot(rf_imp_df, aes(reorder(Symbol, MeanDecreaseAccuracy), MeanDecreaseAccuracy)) + geom_col(fill = "grey40") +
            coord_flip() + labs(title = "Random forest variable importance", x = NULL) + theme_classic(),
          "Figure_F1_RandomForest_Importance", 7, 5)
})
fig("F2 XGB importance", {
  xi <- xgb.importance(model = final_models$XGBoost$model)
  xi <- as.data.frame(xi); xi$Symbol <- unname(lab_map[xi$Feature])
  write.csv(xi, file.path(OUTDIR, "XGBoost_Feature_Importance.csv"), row.names = FALSE)
  save_gg(ggplot(xi, aes(reorder(Symbol, Gain), Gain)) + geom_col(fill = "grey40") + coord_flip() +
            labs(title = "XGBoost feature importance (gain)", x = NULL) + theme_classic(),
          "Figure_F2_XGBoost_Importance", 7, 5)
})
fig("F3 single-fit coefficients", {
  ec <- enet_single; ec$Label <- make.unique(ec$Symbol)
  ec$Direction <- ifelse(ec$Coefficient > 0, "Higher in tumor", "Higher in normal")
  save_gg(ggplot(ec, aes(reorder(Label, Coefficient), Coefficient, fill = Direction)) + geom_col() + coord_flip() +
            scale_fill_manual(values = c("Higher in tumor" = "#B2182B", "Higher in normal" = "#2166AC")) +
            labs(title = "Single-fit Elastic Net coefficients (all ENA; comparison only)", x = NULL) + theme_classic(),
          "Figure_F3_SingleFit_ElasticNet_Coefficients", 7, max(4, 0.3 * nrow(ec) + 2))
})

############################ 12. TABLES AND SUMMARY
expr_out <- function(lcpm, cond, file) {
  d <- as.data.frame(t(lcpm[panel_genes, , drop = FALSE]))
  colnames(d) <- unname(lab_map[panel_genes])
  d$Sample <- rownames(d); d$Condition <- as.character(cond)
  write.csv(d, file.path(OUTDIR, file), row.names = FALSE)
}
expr_out(ENA_lcpm, y, "Panel_Expression_ENA_log2CPM.csv")
expr_out(TCGA_lcpm, tcga_meta$Condition, "Panel_Expression_TCGA_log2CPM.csv")

model_comparison <- bind_rows(
  cv_summary %>% transmute(Dataset = "ENA Nested CV", Model, Accuracy = Mean_Accuracy,
                           Sensitivity = Mean_Sensitivity, Specificity = Mean_Specificity,
                           Precision = Mean_Precision, F1 = Mean_F1, AUC = Mean_AUC, SD_AUC),
  tcga_metrics %>% transmute(Dataset = "TCGA External Test", Model, Accuracy, Sensitivity,
                             Specificity, Precision, F1, AUC, SD_AUC = NA_real_))
write.csv(model_comparison, file.path(OUTDIR, "FINAL_Model_Performance_Comparison.csv"), row.names = FALSE)

writeLines(c("CRC BIOMARKER DISCOVERY - FINAL SUMMARY", "",
             paste0("ENA samples: ", ncol(ENA), " (Normal ", sum(y == "Normal"), " / Tumor ", sum(y == "Tumor"), ")"),
             paste0("TCGA samples: ", ncol(TCGA), " (Normal ", sum(tcga_meta$Condition == "Normal"),
                    " / Tumor ", sum(tcga_meta$Condition == "Tumor"), ")"),
             paste0("Genes shared ENA-TCGA: ", length(common_genes)),
             paste0("DESeq2 (all ENA, cohort-adjusted): ", nrow(deg_all), " genes tested; ", nrow(deg_sig),
                    " DEGs (", sum(deg_all$Status == "Up in tumor"), " up, ", sum(deg_all$Status == "Down in tumor"), " down)"),
             paste0("FINAL BIOMARKER PANEL (selected in >= 60% of outer folds), n = ", length(panel_genes), ":"),
             paste0("  ", panel_tbl$Symbol, " (", panel_tbl$Gene, ") - ", panel_tbl$Direction, collapse = "\n"),
             paste0("Single-fit Elastic Net genes (comparison only): ", paste(enet_single$Symbol, collapse = ", ")),
             "TCGA was never used for DEG selection, feature selection, tuning or thresholding."),
           file.path(OUTDIR, "FINAL_ANALYSIS_SUMMARY.txt"))
############################ 13. OPTIONAL: LOCKED MODEL-LEVEL VALIDATION ON ADDITIONAL COHORTS
# Applies the LOCKED panel + scaling + trained models (no refitting) to another RNA-seq cohort,
# e.g. GSE292858. Put its raw count matrix and a metadata CSV (sample ID + condition) in the
# working directory and adjust the file names below. Skipped automatically if files are absent.
EXTRA_COHORTS <- list(GSE292858 = file.path(DATA_DIR, c(counts = "GSE292858_counts.tsv", meta = "GSE292858_metadata.csv")))
id_from_symbol <- setNames(names(symbol_map), toupper(symbol_map))

validate_external <- function(name, counts_file, meta_file) {
  cat("\n==== LOCKED VALIDATION:", name, "====\n")
  ext <- read_count_file(counts_file)
  if (mean(rownames(ext) %in% rownames(ENA)) < 0.5 && length(symbol_map) > 0) {   # symbols -> Ensembl
    m <- id_from_symbol[rownames(ext)]; keep <- !is.na(m)
    ext <- ext[keep, , drop = FALSE]; rownames(ext) <- unname(m[keep])
    ext <- aggregate_duplicate_genes(ext)
  }
  em <- read.csv(meta_file, check.names = FALSE, stringsAsFactors = FALSE)
  sc_col <- find_metadata_column(em, c("sample", "Sample_ID", "sampleid", "sample_name", "id", "geo_accession"))
  cd_col <- find_metadata_column(em, c("condition", "group", "status", "sample_type", "type"))
  if (is.null(sc_col) || is.null(cd_col)) stop("Cannot find sample/condition columns in ", meta_file)
  ids <- trimws(as.character(em[[sc_col]])); cnd <- normalize_condition(em[[cd_col]])
  idx <- match_samples(colnames(ext), ids)
  keep <- !is.na(idx) & !is.na(cnd[idx])
  ext <- ext[, keep, drop = FALSE]
  cond <- factor(cnd[idx[keep]], levels = c("Normal", "Tumor"))
  cat(name, "samples:", ncol(ext), "| Normal", sum(cond == "Normal"), "| Tumor", sum(cond == "Tumor"), "\n")

  lc <- log_cpm_full(ext)
  avail <- intersect(panel_genes, rownames(lc))
  gs <- gene_stats(lc, cond, avail); gs$Symbol <- sym(gs$Gene)
  write.csv(gs, file.path(OUTDIR, paste0("Validation_", name, "_GeneLevel.csv")), row.names = FALSE)

  miss <- setdiff(panel_genes, rownames(lc))
  if (length(miss) > 0) {
    message("  Panel genes missing in ", name, ": ", paste(sym(miss), collapse = ", "),
            " -> model-level validation skipped (gene-level table written).")
    return(invisible(NULL))
  }
  x <- sweep(sweep(t(lc[panel_genes, , drop = FALSE]), 2, panel_sc$means, "-"), 2, panel_sc$sds, "/")
  ms <- list(); ps <- list()
  for (mn in names(final_models)) {
    prob <- predict_model(final_models[[mn]], panel_train_x, x, y)
    m <- calculate_metrics(cond, prob); m$Model <- mn; ms[[mn]] <- m
    ps[[mn]] <- data.frame(Sample = colnames(ext), Actual = as.character(cond), Probability = prob,
                           Prediction = ifelse(prob >= 0.5, "Tumor", "Normal"), Model = mn)
  }
  mt <- do.call(rbind, ms); rownames(mt) <- NULL
  pr <- do.call(rbind, ps); rownames(pr) <- NULL
  write.csv(mt, file.path(OUTDIR, paste0("Validation_", name, "_Model_Metrics.csv")), row.names = FALSE)
  write.csv(pr, file.path(OUTDIR, paste0("Validation_", name, "_Predictions.csv")), row.names = FALSE)
  print(mt)
  fig(paste(name, "ROC"), {
    rd <- do.call(rbind, lapply(unique(pr$Model), function(mn) {
      d <- pr[pr$Model == mn, ]; roc_df(d$Actual, d$Probability, mn) }))
    save_gg(ggplot(rd, aes(FPR, TPR, colour = Model)) + geom_path(linewidth = 0.9) +
              geom_abline(slope = 1, intercept = 0, linetype = "dashed") +
              labs(title = paste0("Locked-model validation: ", name), x = "1 - Specificity", y = "Sensitivity") +
              theme_classic(), paste0("Figure_G_", name, "_Model_ROC"), 8, 6)
  })
  fig(paste(name, "gene AUC"), {
    save_gg(ggplot(gs, aes(AUC, reorder(Symbol, AUC))) + geom_point(size = 3.5) + xlim(0.5, 1) +
              labs(title = paste0("Single-gene AUC in ", name), x = "AUC", y = NULL) + theme_classic(),
            paste0("Figure_G_", name, "_Gene_AUC"), 7, 5)
  })
  invisible(mt)
}
for (nm in if (RUN_EXTRA_VALIDATION) names(EXTRA_COHORTS) else character(0)) {
  f <- EXTRA_COHORTS[[nm]]
  if (all(file.exists(f))) {
    tryCatch(validate_external(nm, f[["counts"]], f[["meta"]]),
             error = function(e) message("Validation of ", nm, " failed: ", conditionMessage(e)))
  } else {
    message("Skipping ", nm, ": ", paste(f[!file.exists(f)], collapse = ", "), " not found.")
  }
}

capture.output(sessionInfo(), file = file.path(OUTDIR, "R_SessionInfo.txt"))

cat("\nDONE. Files written to:", normalizePath(OUTDIR), "\n")
cat(paste(sort(list.files(OUTDIR)), collapse = "\n"), "\n")
cat("Any 'FIGURE FAILED' messages above name a figure that did not render.\n")
