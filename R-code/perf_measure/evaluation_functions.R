################################################################################
# EVALUATION METRICS
# 
# Top-level categories included:
#  - Distributional Similarity
#    (univariate, bivariate, multivariate)
#  - Downstream Utility
#    (predictive performance, inference performance (tbd))
#  - Overfitting/Memorization
#  - Privacy (MIA/AIA)
################################################################################

# If needed: install dcortools package patch to fix ncc integer overflow issue
#devtools::install_github("jkapar/dcortools", ref = "ncc-preptoterms_fast")

suppressPackageStartupMessages({
  library(ranger)
  library(data.table)
  library(dcortools)
  library(pROC)
  library(foreach)
  library(caret)
  library(gower)
  library(approxOT)
  library(kernlab)
  library(philentropy)
  library(binom)
  library(HiClimR)
})

# ==============================================================================
# Distributional similarity
# ==============================================================================

# ------------------------------------------------------------------------------
# Multivariate measures
# ------------------------------------------------------------------------------

### --- Wasserstein distance (WD) ----------------------------------------------
WD <- function(syn, real_train, real_test = NULL,
               WD_mixed_data_handling = "mixed_distance", WD_base_normalization = real_train, WD_approx_method = "networkflow",
               NN_n_chunks_x = 1, NN_chunks_parallel = F, ...) {
  
  # convert to data.table (maybe vectors are given in univariate case)
  syn <- as.data.table(copy(syn))
  real_train <- as.data.table(copy(real_train))
  if (!is.null(real_test)) real_test <- as.data.table(copy(real_test))
  
  if (ncol(syn) == 1) {
    # univariate case
    if (!is.numeric(syn[[1]])) {
      # 1D Wasserstein for categorical variables with Hamming distance (used in mixed distance) equals TVD 
      res <- TV_Distance(syn[[1]], real_train[[1]], if(!is.null(real_test)) real_test[[1]])
      res[, metric_info := "WD"]
    } else {
      data_processed <- encode_mixed_dt(syn, real_train, real_test, base_normalization = WD_base_normalization, output = "matrix")
      WD_TrainSyn <- wasserstein(na.omit(data_processed$real_train), na.omit(data_processed$syn), p = 1, ground_p = 1, method = "univariate")
      res <- data.table(metric = "Univariate distance", metric_info = "WD", pair = "Train-Syn", result = WD_TrainSyn)
      if (!is.null(real_test)) {
        WD_TestSyn <- wasserstein(na.omit(data_processed$real_test), na.omit(data_processed$syn), p = 1, ground_p = 1, method = "univariate")
        WD_TrainTest <- wasserstein(na.omit(data_processed$real_train), na.omit(data_processed$real_test), p = 1, ground_p = 1, method = "univariate")
        res <- data.table(metric = "Univariate distance", metric_info = "WD", pair = c("Train-Syn", "Test-Syn", "Train-Test"), result = c(WD_TrainSyn, WD_TestSyn, WD_TrainTest))
      }
    }
  } else {
    # multivariate case
    if (WD_mixed_data_handling == "mixed_distance") {
      d_TrainSyn <- mixed_distance(real_train, syn, base_normalization = WD_base_normalization, n_chunks_x = NN_n_chunks_x, parallel_chunks_x = NN_chunks_parallel)
      cost_matrix_TrainSyn <- matrix(d_TrainSyn, nrow = nrow(real_train), ncol = nrow(syn), byrow = TRUE)
      WD_TrainSyn <- wasserstein(cost = cost_matrix_TrainSyn, a = rep(1/nrow(real_train), nrow(real_train)), b = rep(1/nrow(syn), nrow(syn)),
                                 p = 1, ground_p = 1, method = WD_approx_method)
      res <- data.table(metric = "WD", pair = "Train-Syn", result = WD_TrainSyn)
      if (!is.null(real_test)) {
        d_TestSyn <- mixed_distance(real_test, syn, base_normalization = WD_base_normalization, n_chunks_x = NN_n_chunks_x, parallel_chunks_x = NN_chunks_parallel)
        d_TrainTest <- mixed_distance(real_train, real_test, base_normalization = WD_base_normalization,  n_chunks_x = NN_n_chunks_x, parallel_chunks_x = NN_chunks_parallel)
        cost_matrix_TestSyn <- matrix(d_TestSyn, nrow = nrow(real_test), ncol = nrow(syn), byrow = TRUE)
        cost_matrix_TrainTest <- matrix(d_TrainTest, nrow = nrow(real_train), ncol = nrow(real_test), byrow = TRUE)
        WD_TestSyn <- wasserstein(cost = cost_matrix_TestSyn, a = rep(1/nrow(real_test), nrow(real_test)), b = rep(1/nrow(syn), nrow(syn)),
                                  p = 1, ground_p = 1, method = WD_approx_method)
        WD_TrainTest <- wasserstein(cost = cost_matrix_TrainTest, a = rep(1/nrow(real_train), nrow(real_train)), b = rep(1/nrow(real_test), nrow(real_test)),
                                    p = 1, ground_p = 1, method = WD_approx_method)
        res <- data.table(metric = "WD", pair = c("Train-Syn", "Test-Syn", "Train-Test"), result = c(WD_TrainSyn, WD_TestSyn, WD_TrainTest))
      }
    } else if (WD_mixed_data_handling == "normalize_onehot") {
      data_processed <- encode_mixed_dt(syn, real_train, real_test, base_normalization = WD_base_normalization, output = "matrix")
      WD_TrainSyn <- wasserstein(data_processed$real_train, data_processed$syn, p = 1, ground_p = 1, method = WD_approx_method)
      res <- data.table(metric = "WD", pair = "Train-Syn", result = WD_TrainSyn)
      if (!is.null(real_test)) {
        WD_TestSyn <- wasserstein(as.matrix(data_processed$real_test), as.matrix(data_processed$syn), p = 1, ground_p = 1, method = WD_approx_method)
        WD_TrainTest <- wasserstein(as.matrix(data_processed$real_train), as.matrix(data_processed$real_test), p = 1, ground_p = 1, method = WD_approx_method)
        res <- data.table(metric = "WD", pair = c("Train-Syn", "Test-Syn", "Train-Test"), result = c(WD_TrainSyn, WD_TestSyn, WD_TrainTest))
      }
    }
  }
  
  res[, pair := factor(pair, levels = c("Train-Syn", "Test-Syn", "Train-Test"))][]
}

### --- Maximum mean discrepancy (MMD) mit rbf kernel---------------------------
MMD <- function(syn, real_train, real_test = NULL, seed = NULL, ...) {
  if(!is.null(seed)) set.seed(seed)
  
  # convert to data.table (maybe vectors are given in univariate case)
  syn <- as.data.table(copy(syn))
  real_train <- as.data.table(copy(real_train))
  if (!is.null(real_test)) real_test <- as.data.table(copy(real_test))
  
  data_processed <- encode_mixed_dt(syn, real_train, real_test, scale_OH = sqrt(0.5),  output = "matrix")
  sigma <- sigest(data_processed$real_train, fr = 1, scaled = F)[[2]]
  MMD_TrainSyn <- kmmd(data_processed$real_train, data_processed$syn, kernel = laplacedot(sigma))@mmdstats[1]
  res <- data.table(metric = "MMD", pair = "Train-Syn", result = MMD_TrainSyn)
  if (!is.null(real_test)) {
    MMD_TestSyn <- kmmd(data_processed$real_test, data_processed$syn, kernel = rbfdot(sigma))@mmdstats[1]
    MMD_TrainTest <- kmmd(data_processed$real_train, data_processed$real_test, kernel = rbfdot(sigma))@mmdstats[1]
    res <- data.table(metric = "MMD", pair = c("Train-Syn", "Test-Syn", "Train-Test"), result = c(MMD_TrainSyn, MMD_TestSyn, MMD_TrainTest))
  }
  res[, pair := factor(pair, levels = c("Train-Syn", "Test-Syn", "Train-Test"))]
  if (ncol(syn) > 1) {
    res[]
  } else {
    res[, `:=` (metric = "Univariate distance", metric_info = "MMD")][, .(metric, metric_info, pair, result)]
  }
}


### --- Maximum mean discrepancy (MMD) with RF kernel --------------------------
# (needs ranger branch "completely random forests")
MMD_RF <- function(syn, real_train, rf_train_data = "combined", normalized_kernel = FALSE, unbiased = FALSE,
                   rf_min.bucket = 2, rf_num.trees = 50,  rf_replace = FALSE, rf_sample.fraction = 1, seed = NULL, ...) {
  
  
  if (rf_train_data == "real_train") {
    rf_train_data_ <- real_train
  } else if (rf_train_data == "syn") {
    rf_train_data_ <- syn
  } else if (rf_train_data == "combined") {
    rf_train_data_ <- rbind(real_train, syn)
  } else if (is.data.frame(rf_train_data) | is.data.table(rf_train_data)) {
    rf_train_data_ <- rf_train_data
    rf_train_data <- "custom"
  } else {
    stop("rf_train_data must be 'combined', 'real_train', 'syn', or a custom data.frame/data.table")
  }
  
  n <- as.numeric(nrow(real_train))
  m <- as.numeric(nrow(syn))
  
  # train unsupervised random forest (completely random)
  urf <- ranger(x = rf_train_data_, y = sample(1:1000, nrow(rf_train_data_), replace = T), min.bucket = rf_min.bucket, num.trees = rf_num.trees, classification = T, mtry = ncol(real_train),
                splitrule = "extratrees", num.random.splits = 1, sample.fraction = rf_sample.fraction, replace = rf_replace, seed = seed)
  
  
  predict_data <- rbind(real_train, syn)
  if (rf_train_data == "custom") {
    predict_data <- rbind(predict_data, rf_train_data_)
  }
  
  leafIDs <- stats::predict(urf, predict_data, type = 'terminalNodes')$predictions + 1L
  leafIDs <- apply(leafIDs, 2, function(x) match(x, unique(x)))
  max_IDs <- c(0, cumsum(apply(leafIDs, 2, max))[-rf_num.trees])
  leaf_IDs_unique <- sweep(leafIDs, 2, max_IDs, "+")
  leaf_IDs_real_train_unique <- leaf_IDs_unique[1:n, ]
  leaf_IDs_syn_unique <- leaf_IDs_unique[(n+1):(n+m), ]
  if (rf_train_data == "combined") {
    leaf_IDs_rf_train_unique <- leaf_IDs_unique
  } else if (rf_train_data == "real_train") {
    leaf_IDs_rf_train_unique <- leaf_IDs_real_train_unique
  } else if (rf_train_data == "syn") {
    leaf_IDs_rf_train_unique <- leaf_IDs_syn_unique
  } else if (rf_train_data == "custom") {
    leaf_IDs_rf_train_unique <- leaf_IDs_unique[(n+m+1):nrow(leaf_IDs_unique), ]
  }
  
  max_ID_global <- max(leaf_IDs_rf_train_unique)
  
  sum_phi_real_train <- tabulate(leaf_IDs_real_train_unique, max_ID_global)
  sum_phi_syn <- tabulate(leaf_IDs_syn_unique, max_ID_global)
  
  if (normalized_kernel) {
    normalization <- tabulate(leaf_IDs_rf_train_unique, max_ID_global)
    sum_phi_real_train <- sum_phi_real_train / sqrt(normalization)
    sum_phi_syn <- sum_phi_syn / sqrt(normalization)
  } else {
    normalization <- rep(1, length(max_ID_global))
  }
  
  
  if (unbiased) {
    mmd2 <- ((sum(sum_phi_real_train^2) - sum(sum_phi_real_train/sqrt(normalization)))/(n*(n-1)) + 
               (sum(sum_phi_syn^2) - sum(sum_phi_syn/sqrt(normalization)))/(m*(m-1)) - 
               2 * sum(sum_phi_real_train * sum_phi_syn)/(n*m))/rf_num.trees
  } else {
    mmd2 <- sum((sum_phi_real_train/n - sum_phi_syn/m)^2)/rf_num.trees
  }
  mmd2 <- max(0, mmd2)
  mmd <- sqrt(mmd2)
  mmd
}

### --- Detection / C2ST -------------------------------------------------------
detection <- function(syn, real_train = NULL, real_test = NULL, detection_metric = NULL, seed = NULL, ...) {
  
  run_detection <- function(data1, data2) {
    combined <- rbind(
      cbind(data1, is_data1 = 1),
      cbind(data2, is_data1 = 0))
    
    rf_detection <- ranger(x = combined[, -"is_data1"],
                           y = factor(combined[["is_data1"]]),
                           probability = T,
                           seed = seed)
    oob_preds_detection <- rf_detection$predictions[, "1"]
    Accuracy_detection <- mean((oob_preds_detection > 0.5) == combined[["is_data1"]])
    AUROC_detection <- as.numeric(pROC::auc(roc(combined[["is_data1"]], oob_preds_detection, quiet = T)))
    pRMSE_detection <- sqrt(mean((oob_preds_detection - 0.5)^2))
    data.table(metric = "Detection",
               metric_info = c("Accuracy", "AUROC", "pRMSE"),
               result = c(Accuracy_detection, AUROC_detection, pRMSE_detection))
  }
  
  if (!is.null(real_train)) {
    detection_RealTrnSyn <- run_detection(real_train, syn)
    detection_RealTrnSyn[, pair := "Train-Syn"]
  } else {
    detection_RealTrnSyn <- NULL
  }
  
  if (!is.null(real_test)) {
    detection_RealTestSyn <- run_detection(real_test, syn)
    detection_RealTestSyn[, pair := "Test-Syn"]
  } else {
    detection_RealTestSyn <- NULL
  }
  
  if (!is.null(real_train) & !is.null(real_test)) {
    detection_RealTrnRealTest <- run_detection(real_train, real_test)
    detection_RealTrnRealTest[, pair := "Train-Test"]
  } else {
    detection_RealTrnRealTest <- NULL
  }
  
  res <- rbind(detection_RealTrnSyn, detection_RealTestSyn, detection_RealTrnRealTest)
  
  res[, pair := factor(pair, levels = c("Train-Syn", "Test-Syn", "Train-Test"))]
  res_mean <- res[, .(result = mean(result)), keyby = .(metric, metric_info, pair)]
  
  if(is.null(detection_metric)) {
    # return all detection metrics
    return(res_mean)
  } else if (detection_metric %in% c("Accuracy", "AUROC", "pRMSE")) {
    # return only specified detection metric
    return(res_mean[metric_info == detection_metric, ])
  } else {
    stop("Invalid detection_metric. Choose either 'AUROC', 'pRMSE', or NULL to return both.")
  }
}


# ------------------------------------------------------------------------------
# Univariate measures
# ------------------------------------------------------------------------------

# --- Total variation distance (categorical variables; numeric with binning) -
TV_Distance <- function(syn_col, real_train_col, real_test_col = NULL, n_bins = NULL) {
  tabled_data <- bin_and_table(syn_col, real_train_col, real_test_col, n_bins)
  probs_syn_train <- rbind(as.numeric(tabled_data$p_syn), as.numeric(tabled_data$p_real_train))
  d_syn_train <- as.numeric(0.5 * distance(probs_syn_train, method = "manhattan", mute.message = T))
  
  if (!is.null(real_test_col)) {
    probs_syn_test <- rbind(as.numeric(tabled_data$p_syn), as.numeric(tabled_data$p_real_test))
    probs_train_test <- rbind(as.numeric(tabled_data$p_real_train), as.numeric(tabled_data$p_real_test))
    d_syn_test <- as.numeric(0.5 * distance(probs_syn_test, method = "manhattan", mute.message = T))
    d_train_test <- as.numeric(0.5 * distance(probs_train_test, method = "manhattan", mute.message = T))
    data.table(metric = "Univariate distance", metric_info = "TVD",
               pair = factor(c("Train-Syn", "Test-Syn", "Train-Test"), levels = c("Train-Syn", "Test-Syn", "Train-Test")),
               result = c(d_syn_train, d_syn_test, d_train_test))
  } else {
    data.table(metric = "Univariate distance", metric_info = "TVD", pair = factor(c("Train-Syn"), levels = c("Train-Syn", "Test-Syn", "Train-Test")), result = d_syn_train)
  }
}

### --- Univariate distance (per column and mean over columns) -----------------
univariate_distance <- function(syn, real_train, real_test = NULL,
                                UD_num_method = WD,
                                UD_cat_method = WD,
                                UD_n_bins = NULL,
                                UD_results = "mean", ...) {
  
  stopifnot(ncol(real_train) == ncol(syn))
  stopifnot(all(names(real_train) == names(syn)))
  
  is_named_list <- is.list(UD_n_bins) && !is.null(names(UD_n_bins))
  
  univariate_distance_per_col <- foreach(col = names(real_train), .combine = rbind) %do% {
    
    syn_col <- syn[[col]]
    real_train_col <- real_train[[col]]
    real_test_col <- if (!is.null(real_test)) real_test[[col]] else NULL
    if (is.numeric(real_train_col)) {
      UD_n_bins_col <- if (is_named_list) UD_n_bins[[col]] else UD_n_bins
      cbind(metric = "Univariate distance", variable = col, type = factor("numeric", levels = c("numeric", "categorical")), UD_num_method(syn_col, real_train_col, real_test_col, n_bins = UD_n_bins_col)[, -"metric"])
    } else {
      cbind(metric = "Univariate distance", variable = col, type = factor("categorical", levels = c("numeric", "categorical")), UD_cat_method(syn_col, real_train_col, real_test_col)[, -"metric"])
    }
  }
  
  methods <- univariate_distance_per_col[, .(methods = unique(metric_info)), keyby = type][, paste0(methods, collapse = ",")]
  
  univariate_distance_mean <- rbind(
    univariate_distance_per_col[, .(result = mean(result, na.rm = T)), by = .(metric, type, metric_info, pair)],
    univariate_distance_per_col[, .(result = mean(result, na.rm = T), metric_info = paste0(methods,
                                                                                           "; mean")), by = .(metric, pair)][, type := "mixed"]
  )
  
  if (UD_results == "all") {
    rbind(univariate_distance_per_col, univariate_distance_mean, fill = TRUE)
  } else if (UD_results == "univariate") {
    univariate_distance_per_col
  } else if (UD_results == "mean_per_type") {
    univariate_distance_mean
  } else if (UD_results == "mean") {
    univariate_distance_mean[type == "mixed", -"type"]
  } else {
    stop("Invalid choice of UD_results. Choose either 'all', 'univariate', 'mean_per_type', or 'mean'.")
  }
}


# ------------------------------------------------------------------------------
# Bivariate measures (pair-wise dependencies)
# ------------------------------------------------------------------------------

### ---- Correlation distance (per column-pair and aggregated) -----------------
correlation_distance <- function(syn, real_train, real_test = NULL,
                                 CD_method = "dcor",
                                 CD_agg = "MAE",
                                 CD_rescale_pearson = TRUE,
                                 CD_results = "agg", ...) {
  if(ncol(real_train) != ncol(syn) || !(all(colnames(real_train) == colnames(syn)))) {
    stop("Column names of real_train and syn data must match.")
  }
  cor_table_cols <- c("variable1", "variable2", "type.x", "type.y", "cor_type", "cor")
  if(ncol(real_train) == length(cor_table_cols) && all(colnames(real_train) == cor_table_cols) &
     (is.null(real_test) || (ncol(real_test) == length(cor_table_cols) && all(colnames(real_test) == cor_table_cols))) &
     (ncol(syn) == length(cor_table_cols) && all(colnames(syn) == cor_table_cols))) {
    cor_table_real_train <- copy(real_train)
    cor_table_real_test <- copy(real_test)
    cor_table_syn <- copy(syn)
  } else {
    if (CD_method == "classic") {
      cor_table_real_train <- cor_mixed(real_train)
      cor_table_syn <- cor_mixed(syn)
      if (!is.null(real_test)) {
        cor_table_real_test <- cor_mixed(real_test)
      }
      
      correlations <- cbind(cor_table_real_train[, cbind(.SD[,-"cor"], cor_train = cor)],
                            cor_test = if (!is.null(real_test)) cor_table_real_test[, cor] else NA_real_,
                            cor_syn = cor_table_syn[, cor])
      
      if (CD_rescale_pearson) {
        # Rescale Pearson correlations to [0, 1]
        
        cor_cols <- names(which(sapply(correlations, is.numeric)))
        correlations[cor_type == "Pearson", (cor_cols) := (.SD + 1) / 2, .SDcols = cor_cols]
      }
    } else if (CD_method == "dcor") {
      
      cor_table_full_matrix <- matrix(
        NA_real_,
        nrow = ncol(real_train),
        ncol = ncol(real_train),
        dimnames = list(names(real_train), names(real_train))
      )
      
      cols_w_values_real_train <- names(which(sapply(real_train, \(col) any(!is.na(col)))))
      cor_table_real_train_matrix_part <- dcmatrix(real_train[, .SD, .SDcols = cols_w_values_real_train], calc.dcov = F, return.data = F, bias.corr = T, fc.discrete = T, use = "pairwise.complete.obs")$dcor
      cor_table_real_train_matrix <- copy(cor_table_full_matrix)
      cor_table_real_train_matrix[cols_w_values_real_train, cols_w_values_real_train] <- cor_table_real_train_matrix_part
      
      cols_w_values_syn <- names(which(sapply(syn, \(col) any(!is.na(col)))))
      cor_table_syn_matrix_part <- dcmatrix(syn[, .SD, .SDcols = cols_w_values_syn], calc.dcov = F, return.data = F, bias.corr = T, fc.discrete = T, use = "pairwise.complete.obs")$dcor
      cor_table_syn_matrix <- copy(cor_table_full_matrix)
      cor_table_syn_matrix[cols_w_values_syn, cols_w_values_syn] <- cor_table_syn_matrix_part
      
      cor_table_real_train <- cor_table_real_train_matrix[lower.tri(cor_table_real_train_matrix)]
      cor_table_syn <- cor_table_syn_matrix[lower.tri(cor_table_syn_matrix)]
      if (!is.null(real_test)) {
        cols_w_values_real_test <- names(which(sapply(real_test, \(col) any(!is.na(col)))))
        cor_table_real_test_matrix_part <- dcmatrix(real_test[, .SD, .SDcols = cols_w_values_real_test], calc.dcov = F, return.data = F, bias.corr = T, fc.discrete = T, use = "pairwise.complete.obs")$dcor
        cor_table_real_test_matrix <- copy(cor_table_full_matrix)
        cor_table_real_test_matrix[cols_w_values_real_test, cols_w_values_real_test] <- cor_table_real_test_matrix_part
        cor_table_real_test <- cor_table_real_test_matrix[lower.tri(cor_table_real_test_matrix)]
      }
      correlations <- data.table(variable1 = colnames(cor_table_real_train_matrix)[col(cor_table_real_train_matrix)[lower.tri(cor_table_real_train_matrix)]],
                                 variable2 = rownames(cor_table_real_train_matrix)[row(cor_table_real_train_matrix)[lower.tri(cor_table_real_train_matrix)]],
                                 cor_type = "dcor",
                                 cor_train = cor_table_real_train,
                                 cor_test = if (!is.null(real_test)) cor_table_real_test else NA_real_,
                                 cor_syn = cor_table_syn)
      # clip small negative values due to bias correction to 0 for train, test, syn
      cor_cols <- names(which(sapply(correlations, is.numeric)))
      correlations[, (cor_cols) := lapply(.SD, \(x) pmax(x,0)), .SDcols = cor_cols]
      
    } else {
      stop("Unknown method. Use 'classic' or 'dcor'.")
    }
  }
  
  correlations[, cor_diff_syn_train := cor_syn - cor_train]
  correlations[, cor_diff_syn_test := cor_syn - cor_test]
  correlations[, cor_diff_train_test := cor_train - cor_test]
  correlations[, metric_info := cor_type]
  
  # melt to pair variable
  correlations <- melt(correlations, id.vars = c("metric_info",  "variable1", "variable2"),
                       measure.vars = c("cor_train", "cor_test", "cor_syn", "cor_diff_syn_train", "cor_diff_syn_test", "cor_diff_train_test"), 
                       variable.name = "dataset", value.name = "result")
  correlations[dataset == "cor_train", data_split := "Train"]
  correlations[dataset == "cor_test", data_split := "Test"]
  correlations[dataset == "cor_syn", data_split := "Syn"]
  correlations[dataset == "cor_diff_syn_train", `:=` (data_split = NA_character_, pair = "Train-Syn")]
  correlations[dataset == "cor_diff_syn_test", `:=` (data_split = NA_character_, pair = "Test-Syn")]
  correlations[dataset == "cor_diff_train_test", `:=` (data_split = NA_character_, pair = "Train-Test")]
  correlations[, dataset := NULL]
  
  correlations[, metric := fifelse(is.na(data_split), "Correlation difference", "Correlation")]
  setcolorder(correlations, c("metric", "metric_info", "variable1", "variable2", "data_split", "pair", "result"))
  
  if (CD_agg == "MAE") {
    # mean absolute difference
    cor_distance <- correlations[metric == "Correlation difference", .(result = mean(abs(result), na.rm = T)), by = pair]
  } else if (CD_agg == "Frobenius") {
    # Frobenius norm
    cor_distance <- correlations[metric == "Correlation difference", .(result = sqrt(sum((result^2), na.rm = T))), by = pair]
    # Correlation of correlations
  } else if (CD_agg == "cor") {
    cor_distance <- data.table(
      pair = c("Train-Syn", "Test-Syn", "Train-Test"), 
      result = c(cor(correlations[data_split == "Train", result], correlations[data_split == "Syn", result], use = "pairwise.complete.obs"),
                 cor(correlations[data_split == "Test", result], correlations[data_split == "Syn", result], use = "pairwise.complete.obs"),
                 cor(correlations[data_split == "Train", result], correlations[data_split == "Test", result], use = "pairwise.complete.obs"))
    )
  } else {
    stop("Unknown aggregation method. Use 'MAE', 'Frobenius' or 'cor'.")
  }
  
  cor_distance[, pair := factor(pair, levels = c("Train-Syn", "Test-Syn", "Train-Test"))]
  correlation_distance <- data.table(metric = "Correlation distance",
                                     metric_info = paste0(CD_method, "; ", CD_agg),
                                     cor_distance)
  
  if (is.null(real_test)) {
    correlations <- correlations[(metric == "Correlation" & data_split != "Test") | (metric == "Correlation difference" & pair == "Train-Syn"), ]
    correlation_distance <- correlation_distance[pair == "Train-Syn", ]
  }
  
  if (CD_results == "all") {
    rbind(correlations, correlation_distance, fill = TRUE)
  } else if (CD_results == "agg") {
    correlation_distance
  } else {
    stop("Invalid choice of CD_results. Choose either 'agg' or 'all'.")
  }
}

### --- Correlation distance with fastCor ----------------------------------

fastCor_dist_measure <- function (real_train, syn) {
  # char_vars
  char_vars <- names(which(sapply(real_train, is.character)))
  print(char_vars)
  real_train_num <- real_train[, .SD, .SDcols = -char_vars]
  syn_num <- syn[, .SD, .SDcols = -char_vars]
  real_cor <- HiClimR::fastCor(real_train_num)
  cor_dist <- as.matrix(HiClimR::fastCor(syn_num) - real_cor)
  cor_dist[lower.tri(cor_dist, diag = FALSE)] <- NA
  cor_dist_mean <- mean(abs(cor_dist), na.rm = TRUE)
  return(cor_dist_mean)
}


# ==============================================================================
# Downstream utilty
# ==============================================================================

### --- (Generalised) Machine learning utility (Predictive performance) ---------------------------------
MLU <- function(syn, real_train, real_test, MLU_target = NULL, MLU_results = "mean", seed = NULL, MLU_parallel_target = FALSE, ...) {
  
  if (is.null(MLU_target)) {
    target <- names(real_train)
  } else if (length(MLU_target) == 1) {
    target <- MLU_target
    MLU_parallel_target <- FALSE
  } else if (length(MLU_target) < length(names(real_train))) {
    target <- MLU_target
  }
  
  compute_MLU <- function(target) {
    real_train <- real_train[complete.cases(real_train[, .SD, .SDcols = target]), ]
    real_test <- real_test[complete.cases(real_test[, .SD, .SDcols = target]), ]
    syn <- syn[complete.cases(syn[, .SD, .SDcols = target]), ]
    
    
    is_classification <- fifelse(is.numeric(real_train[[target]]), F, T)
    
    rf_gt <- ranger(x = real_train[, .SD, .SDcols = setdiff(names(real_train), target)],
                    y = (real_train[[target]]), 
                    probability = is_classification,
                    seed = seed,
                    num.threads = if (MLU_parallel_target) 2 else NULL)
    
    rf_syn <- ranger(x = syn[, .SD, .SDcols = setdiff(names(syn), target)],
                     y = syn[[target]],
                     probability = is_classification,
                     seed = seed,
                     num.threads = if (MLU_parallel_target) 2 else NULL)
    
    
    preds_gt <- predict(rf_gt, data = real_test[, .SD, .SDcols = setdiff(names(real_test), target)])$predictions
    preds_syn <- predict(rf_syn, data = real_test[, .SD, .SDcols = setdiff(names(real_test), target)])$predictions
    
    if (is_classification && ncol(preds_gt) > ncol(preds_syn)) {
      missing_cols <- setdiff(colnames(preds_gt), colnames(preds_syn))
      preds_syn <- as.data.table(preds_syn)[, (missing_cols) := 0][]
    }
    
    if (is_classification) {
      AUC_gt <- as.numeric(pROC::auc(multiclass.roc((real_test[[target]]), preds_gt, quiet = T)))
      AUC_syn <- as.numeric(pROC::auc(multiclass.roc((real_test[[target]]), preds_syn, quiet = T)))
      data.table(metric = "MLU", variable = target, metric_info = c("gt; AUROC", "syn; AUROC", "MLU"), result = c(AUC_gt, AUC_syn, AUC_syn/AUC_gt))
    } else {
      MAE_gt <- mean(abs(real_test[[target]] - preds_gt))
      MAE_syn <- mean(abs(real_test[[target]] - preds_syn))
      data.table(metric = "MLU", variable = target, metric_info = c("gt; MAE", "syn; MAE", "MLU"), result = c(MAE_gt, MAE_syn, MAE_gt/MAE_syn))
    }
  }
  
  if (MLU_parallel_target) {
    res <- foreach(target = target, .combine = "rbind") %dopar% compute_MLU(target)
  } else {
    res <- foreach(target = target, .combine = "rbind") %do% compute_MLU(target)
  }
  
  
  if (MLU_results == "mean") {
    res[metric_info == "MLU", .(metric = "MLU",
                                metric_info = if(is.null(MLU_target) | length(MLU_target) > 1) "cross-target (mean)" else paste0("target: ", MLU_target), 
                                result = mean(result))]
  } else if (MLU_results == "per_target") {
    res
  }
}



# ==============================================================================
# Overfitting / memorization
# ==============================================================================

### --- Identical match share --------------------------------------------------

IMS <- function(syn, real_train, real_test, ...) {
  
  syn <- copy(syn)
  real_train <- copy(real_train)
  real_test <- copy(real_test)
  
  setkeyv(syn, names(syn))
  setkeyv(real_train, names(real_train))
  setkeyv(real_test, names(real_test))
  
  IM_target <- syn[real_train, .N, by=.EACHI][, .(identical_matches = N)]
  IM_target[, IM := fifelse(identical_matches > 0, 1/identical_matches, 0)]
  IM_baseline <- real_test[real_train, .N, by=.EACHI][, .(identical_matches = N)]
  IM_baseline[, IM := fifelse(identical_matches > 0, 1/identical_matches, 0)]
  
  data.table(metric = "IMS", metric_info = c("IMS(Train->Syn)", "IMS(Train->Test)",
                                             "Share (Identical_matches(Train->Syn) > Identical_matches(Train->Test))"),
             result = c(mean(IM_target$IM), mean(IM_baseline$IM),
                        mean(IM_target[, identical_matches]>IM_baseline[, identical_matches])
             )
  )
}

DCR <- function(syn, real_train, real_test, anchor = c("syn", "train"), DCR_NN_k = 1, DCR_NN_averaging = c("arithmetic", "harmonic"),
                NN_n_chunks_x = 1, NN_parallel_chunks_x = FALSE, ...) {
  
  anchor <- match.arg(anchor)
  DCR_NN_averaging <- match.arg(DCR_NN_averaging)
  
  if (DCR_NN_averaging == "arithmetic") {
    avg_function <- colmeans
  } else if (DCR_NN_averaging == "harmonic") {
    # here: distance-weighted arithmetic = harmonic
    avg_function <- function(x) {
      DCR_NN_k / colmeans(1 / x)
    }
  }
  
  if (anchor == "train") {
    
    if (nrow(syn) != nrow(real_test)) {
      target_size <- min(nrow(syn), nrow(real_test))
      syn <- syn[sample(seq_len(nrow(syn)), size = target_size, replace = FALSE), ]
      real_test <- real_test[sample(seq_len(nrow(real_test)), size = target_size, replace = FALSE), ]
    }
    
    DCR_train_syn <- avg_function(mixed_distance(real_train, syn, base_normalization = real_train, top_n = DCR_NN_k,
                                                 n_chunks_x = NN_n_chunks_x, parallel_chunks_x = NN_parallel_chunks_x)$distance)
    DCR_train_test <- avg_function(mixed_distance(real_train, real_test, base_normalization = real_train, top_n = DCR_NN_k,
                                                  n_chunks_x = NN_n_chunks_x, parallel_chunks_x = NN_parallel_chunks_x)$distance)
    DCR_test_syn <- avg_function(mixed_distance(real_test, syn, base_normalization = real_train, top_n = DCR_NN_k,
                                                n_chunks_x = NN_n_chunks_x, parallel_chunks_x = NN_parallel_chunks_x)$distance)
    DCR_diff_train_memorization <- DCR_train_syn - DCR_train_test
    DCR_share_memorized <- mean(DCR_diff_train_memorization < 0)
    avg_memorization <- abs(sum(DCR_diff_train_memorization[DCR_diff_train_memorization < 0])/length(DCR_diff_train_memorization))
    
    res <- data.table(metric = "DCR", metric_info = c("d(Train->Syn); q0.05", "d(Train->Test); q0.05", "d(Test->Syn); q0.05",
                                                      "d(Train->Syn); median", "d(Train->Test); median", "d(Test->Syn); median",
                                                      "d(Train->Syn); mean", "d(Train->Test); mean", "d(Test->Syn); mean",
                                                      "Share (d(Train->Syn) < d(Train->Test))", "Memorization; mean"),
                      result = c(quantile(DCR_train_syn, 0.05), quantile(DCR_train_test, 0.05), quantile(DCR_test_syn, 0.05),
                                 median(DCR_train_syn), median(DCR_train_test), median(DCR_test_syn, 0.05),
                                 mean(DCR_train_syn), mean(DCR_train_test), mean(DCR_test_syn, 0.05),
                                 DCR_share_memorized, avg_memorization))
  } else if (anchor == "syn") {
    
    if (nrow(real_train) != nrow(real_test)) {
      target_size <- min(nrow(real_train), nrow(real_test))
      real_train <- real_train[sample(seq_len(nrow(real_train)), size = target_size, replace = FALSE), ]
      real_test <- real_test[sample(seq_len(nrow(real_test)), size = target_size, replace = FALSE), ]
    }
    
    DCR_syn_train <- avg_function(mixed_distance(syn, real_train, base_normalization = real_train, top_n = DCR_NN_k,
                                                 n_chunks_x = NN_n_chunks_x, parallel_chunks_x = NN_parallel_chunks_x)$distance)
    DCR_syn_test <- avg_function(mixed_distance(syn, real_test, base_normalization = real_train, top_n = DCR_NN_k,
                                                n_chunks_x = NN_n_chunks_x, parallel_chunks_x = NN_parallel_chunks_x)$distance)
    DCR_test_train <- avg_function(mixed_distance(real_test, real_train, base_normalization = real_train, top_n = DCR_NN_k,
                                                  n_chunks_x = NN_n_chunks_x, parallel_chunks_x = NN_parallel_chunks_x)$distance)
    DCR_diff_syn_overfitting <- DCR_syn_train - DCR_syn_test
    DCR_share_overfitting <- mean(DCR_diff_syn_overfitting < 0)
    avg_overfitting <- abs(sum(DCR_diff_syn_overfitting[DCR_diff_syn_overfitting < 0])/length(DCR_diff_syn_overfitting))
    res <- data.table(metric = "DCR", metric_info = c("d(Syn->Train); q0.05", "d(Syn->Test); q0.05", "d(Test->Train); q0.05",
                                                      "d(Syn->Train); median", "d(Syn->Test); median", "d(Test->Train); median",
                                                      "d(Syn->Train); mean", "d(Syn->Test); mean", "d(Test->Train); mean",
                                                      "Share (d(Syn->Train) < d(Syn->Test))", "Overfitting; mean"),
                      result = c(quantile(DCR_syn_train, 0.05), quantile(DCR_syn_test,0.05), quantile(DCR_test_train, 0.05),
                                 median(DCR_syn_train), median(DCR_syn_test), median(DCR_test_train),
                                 mean(DCR_syn_train), mean(DCR_syn_test), mean(DCR_test_train),
                                 DCR_share_overfitting, avg_overfitting))
  }
  res
}

# ==============================================================================
# Privacy (Membership and attribute inference)
# ==============================================================================

# DCR-based MIA
MIA_DCR <- function(syn, real_train, real_test, add_DCR_diff = TRUE, NN_k = 1, NN_averaging = c("arithmetic", "harmonic"),
                    NN_n_chunks_x = 1, NN_parallel_chunks_x = FALSE, ...) {
  NN_averaging <- match.arg(NN_averaging)
  if (NN_averaging == "arithmetic") {
    avg_function <- colmeans
  } else if (NN_averaging == "harmonic") {
    # here: distance-weighted arithmetic = harmonic
    avg_function <- function(x) {
      NN_k / colmeans(1 / x)
    }
  }
  scores_DCR <- avg_function(mixed_distance(rbind(real_test, real_train), syn, base_normalization = real_train, top_n = NN_k,
                                            n_chunks_x = NN_n_chunks_x, parallel_chunks_x = NN_parallel_chunks_x)$distance)
  
  if (add_DCR_diff) {
    DCR_nonmembers_nonmembers <- avg_function(mixed_distance(real_test, real_test, base_normalization = real_train, top_n = NN_k+1,
                                                             n_chunks_x = NN_n_chunks_x, parallel_chunks_x = NN_parallel_chunks_x)$distance[-1, , drop = F])
    DCR_members_nonmembers <- avg_function(mixed_distance(real_train, real_test, base_normalization = real_train, top_n = NN_k,
                                                          n_chunks_x = NN_n_chunks_x, parallel_chunks_x = NN_parallel_chunks_x)$distance)
    scores_DCR_diff <- c(scores_DCR[seq_len(nrow(real_test))] - DCR_nonmembers_nonmembers,
                         scores_DCR[(nrow(real_test)+1):length(scores_DCR)] - DCR_members_nonmembers)
    list(DCR = -scores_DCR, DCR_Diff = -scores_DCR_diff)
  } else {
    list(DCR = -scores_DCR)
  }
}

MIA_DPI <- function(syn, real_train, real_test, NN_k = 5, eps = 1e-10, NN_n_chunks_x = 1, NN_parallel_chunks_x = FALSE, ...) {
  NN_nonmembers_rest <- mixed_distance(real_test, rbind(real_test, syn), base_normalization = real_train, top_n = NN_k+1,
                                       n_chunks_x = NN_n_chunks_x, parallel_chunks_x = NN_parallel_chunks_x)$index[-1, ,drop = F]
  DPI_nonmembers <- plogis(colsums(NN_nonmembers_rest > nrow(real_test))/(colsums(NN_nonmembers_rest <= nrow(real_test)) + eps))
  NN_members_rest <- mixed_distance(real_train, rbind(real_test, syn) , base_normalization = real_train, top_n = NN_k,
                                    n_chunks_x = NN_n_chunks_x, parallel_chunks_x = NN_parallel_chunks_x)$index
  DPI_members <- plogis(colsums(NN_members_rest > nrow(real_test))/(colsums(NN_members_rest <= nrow(real_test)) + eps))
  scores_DPI <- c(DPI_nonmembers, DPI_members)
  list(DPI = scores_DPI)
}

MIA_RF <- function(syn, real_train, real_test, seed = NULL, ...) {
  combined <- rbind(
    cbind(real_test, member = 0),
    cbind(syn, member = 1))
  
  rf_MIA <- ranger(x = combined[, -"member"],
                   y = factor(combined[["member"]]),
                   probability = T,
                   seed = seed)
  oob_preds_nonmembers <- rf_MIA$predictions[which(combined$member == 0), 2]
  preds_members <- predict(rf_MIA, data = real_train)$predictions[, 2]
  scores_RF <- c(oob_preds_nonmembers, preds_members)
  list(RF = scores_RF)
}

### --- Membership inference attack (detection- or DCR-based) ------------------
MIA <- function(syn, real_train, real_test, MIA_attacks = c(MIA_DCR, MIA_DPI, MIA_RF),
                MIA_FPR_fixed = c(0, 0.001, 0.01, 0.1), K_folds_eps = 5, times_eps = 5, report_best = TRUE, seed = NULL, ...) {
  
  labels <- as.factor(c(rep(0, nrow(real_test)), rep(1, nrow(real_train))))
  
  scores_attacks <- foreach(attack = MIA_attacks, .combine = c) %do% {
    attack(syn, real_train, real_test, seed = seed, ...)
  }
  
  results <- foreach(attack = names(scores_attacks), .combine = rbind) %do% {
    scores <- scores_attacks[[attack]]
    roc_obj <- roc(labels, scores, direction = "<", quiet = TRUE)
    fpr <- 1 - roc_obj$specificities  # FPR
    tpr <- roc_obj$sensitivities      # TPR
    tpr_at_fpr <- approx(fpr, tpr, xout = MIA_FPR_fixed, method = "constant", f = 0, ties = first)$y
    auroc <- as.numeric(auc(roc_obj))
    
    set.seed(seed, kind = "L'Ecuyer-CMRG")
    folds <- createMultiFolds(labels, k = K_folds_eps, times = times_eps)
    
    eps_folds <- foreach(fold = folds, .combine = "rbind") %do% {
      
      # threshold selection
      labels_selection <- labels[-fold]
      scores_selection <- scores[-fold]
      roc_obj_selection <- roc(labels_selection, scores_selection, direction = "<", quiet = TRUE)
      fpr_selection <- 1 - roc_obj_selection$specificities
      tpr_selection <- roc_obj_selection$sensitivities
      n_train_selection <- sum(labels_selection == 1)
      n_test_selection <- sum(labels_selection == 0)
      tpr_CI_selection <- binom.confint(tpr_selection*n_train_selection, n_train_selection, methods="exact")
      fpr_CI_selection <- binom.confint(fpr_selection*n_test_selection, n_test_selection, methods="exact")
      
      # threshold for positive MIA
      ratio_posMIA <- tpr_CI_selection[, "lower"]/fpr_CI_selection[, "upper"]
      max_ID_posMIA_selection <- which.max(ratio_posMIA)
      threshold_posMIA <- roc_obj_selection$thresholds[max_ID_posMIA_selection]
      
      # threshold for negative MIA
      ratio_negMIA <- (1 - fpr_CI_selection[, "upper"])/(1 - tpr_CI_selection[, "lower"])
      max_ID_negMIA_selection <- which.max(ratio_negMIA)
      threshold_negMIA <- roc_obj_selection$thresholds[max_ID_negMIA_selection]
      
      # evaluate on holdout fold
      labels_holdout <- labels[fold]
      scores_holdout <- scores[fold]
      n_train_holdout <- sum(labels_holdout == 1)
      n_test_holdout <- sum(labels_holdout == 0)
      
      # epsilons for positive MIA
      tp_threshold_posMIA_holdout <- sum((scores_holdout >= threshold_posMIA) & (labels_holdout == 1))
      fp_threshold_posMIA_holdout <- sum((scores_holdout >= threshold_posMIA) & (labels_holdout == 0))
      tpr_CI_posMIA_holdout <- binom.confint(tp_threshold_posMIA_holdout, n_train_holdout, methods="exact")
      fpr_CI_posMIA_holdout <- binom.confint(fp_threshold_posMIA_holdout, n_test_holdout, methods="exact")
      
      eps_CI_low_posMIA <- log(tpr_CI_posMIA_holdout[["lower"]]/fpr_CI_posMIA_holdout[["upper"]])
      eps_CI_high_posMIA <- log(tpr_CI_posMIA_holdout[["upper"]]/fpr_CI_posMIA_holdout[["lower"]])
      eps_estimate_posMIA <- log(tpr_CI_posMIA_holdout[["mean"]]/fpr_CI_posMIA_holdout[["mean"]])
      
      # epsilons for negative MIA
      tp_threshold_negMIA_holdout <- sum((scores_holdout >= threshold_negMIA) & (labels_holdout == 1))
      fp_threshold_negMIA_holdout <- sum((scores_holdout >= threshold_negMIA) & (labels_holdout == 0))
      tpr_CI_negMIA_holdout <- binom.confint(tp_threshold_negMIA_holdout, n_train_holdout, methods="exact")
      fpr_CI_negMIA_holdout <- binom.confint(fp_threshold_negMIA_holdout, n_test_holdout, methods="exact")
      
      eps_CI_low_negMIA <- log((1 - fpr_CI_negMIA_holdout[["upper"]])/(1 - tpr_CI_negMIA_holdout[["lower"]]))
      eps_CI_high_negMIA <- log((1 - fpr_CI_negMIA_holdout[["lower"]])/(1 - tpr_CI_negMIA_holdout[["upper"]]))
      eps_estimate_negMIA <- log((1 - fpr_CI_negMIA_holdout[["mean"]])/(1 - tpr_CI_negMIA_holdout[["mean"]]))
      
      # overall epsilons
      eps_CI_low <- max(eps_CI_low_posMIA, eps_CI_low_negMIA)
      eps_CI_high <- max(eps_CI_high_posMIA, eps_CI_high_negMIA)
      eps_estimate <- max(eps_estimate_posMIA, eps_estimate_negMIA)
      
      data.table(eps_CI_low_posMIA = fifelse(is.na(eps_CI_low_posMIA), 0, max(0, eps_CI_low_posMIA)), 
                 eps_estimate_posMIA = fifelse(is.na(eps_estimate_posMIA), 0, max(0, eps_estimate_posMIA)),
                 eps_CI_high_posMIA = fifelse(is.na(eps_CI_high_posMIA), 0, max(0, eps_CI_high_posMIA)),
                 eps_CI_low_negMIA = fifelse(is.na(eps_CI_low_negMIA), 0, max(0, eps_CI_low_negMIA)),
                 eps_estimate_negMIA = fifelse(is.na(eps_estimate_negMIA), 0, max(0, eps_estimate_negMIA)),
                 eps_CI_high_negMIA = fifelse(is.na(eps_CI_high_negMIA), 0, max(0, eps_CI_high_negMIA)),
                 eps_CI_low = fifelse(is.na(eps_CI_low), 0, max(0, eps_CI_low)),
                 eps_estimate = fifelse(is.na(eps_estimate), 0, max(0, eps_estimate)),
                 eps_CI_high = fifelse(is.na(eps_CI_high), 0, max(0, eps_CI_high))
      )
      
    }
    
    eps_mean <- sapply(eps_folds, mean) # max or mean over folds?
    eps_max <- sapply(eps_folds, max)
    
    data.table(metric = "MIA",
               attack = attack,
               metric_info = c(paste0("TPR at FPR=", MIA_FPR_fixed), "AUROC",
                               paste0(names(eps_mean), "_foldmean"), paste0(names(eps_max), "_foldmax")),
               result = c(tpr_at_fpr, auroc, eps_mean, eps_max))
  }
  
  if (report_best) {
    results[, .SD[which.max(result)], by = metric_info]
  } else {
    results
  }
  
  
}


### --- Attribute inference attack ---------------------------------------------
AIA <- function(syn, real_train, real_test, key_vars, target, seed = NULL, ...) {
  
  is_classification <- fifelse(is.numeric(real_train[[target]]), F, T)
  if (is_classification) {
    target_variable_type <- "Categorical"
    n_categories <- uniqueN(rbind(real_train, real_test)[[target]])
  } else {
    target_variable_type <- "Continuous"
    target_range <- range(rbind(real_train, real_test)[[target]], na.rm = T)
    target_mean <- mean(rbind(real_train, real_test)[[target]], na.rm = T)
    target_sd <- sd(rbind(real_train, real_test)[[target]], na.rm = T)
  }
  
  # omit rows with NA in target
  train <- real_train[complete.cases(real_train[, .SD, .SDcols = target]), ]
  test <- real_test[complete.cases(real_test[, .SD, .SDcols = target]), ]
  syn <- syn[complete.cases(syn[, .SD, .SDcols = target]), ]
  
  # syn -> train inference  
  suppressWarnings({rf_syn <- ranger(x = syn[, .SD, .SDcols = key_vars],
                                     y = syn[[target]],
                                     probability = is_classification,
                                     seed = seed)
  })
  preds_rf_syn <- predict(rf_syn, data = train[, .SD, .SDcols = key_vars])$predictions
  
  # test -> train inference (basline)
  suppressWarnings({rf_test <- ranger(x = test[, .SD, .SDcols = key_vars],
                                      y = test[[target]],
                                      probability = is_classification,
                                      seed = seed)
  })
  preds_rf_test <- predict(rf_test, data = train[, .SD, .SDcols = key_vars])$predictions
  
  if (is_classification) {
    # add columns to preds_rf_syn with zero that are missing levels in preds_rf_syn
    preds_rf_syn <- cbind(preds_rf_syn, matrix(0, nrow = nrow(preds_rf_syn), ncol = length(setdiff(levels(train[[target]]), colnames(preds_rf_syn))), dimnames = list(NULL, setdiff(levels(train[[target]]), colnames(preds_rf_syn)))))
    preds_rf_test <- cbind(preds_rf_test, matrix(0, nrow = nrow(preds_rf_test), ncol = length(setdiff(levels(train[[target]]), colnames(preds_rf_test))), dimnames = list(NULL, setdiff(levels(train[[target]]), colnames(preds_rf_test)))))
    AUC_rf_syn <- as.numeric(pROC::auc(multiclass.roc(train[[target]], preds_rf_syn, quiet = T)))
    AUC_rf_test <- as.numeric(pROC::auc(multiclass.roc(train[[target]], preds_rf_test, quiet = T)))
    advantage <- AUC_rf_syn - AUC_rf_test
    result <- data.table(metric = "AIA", target = target, target_info = paste0(target_variable_type, "; ", n_categories, " categories"),
                         key_vars = paste0(key_vars, collapse = ", "), metric_info = c("Syn->Train AUC", "Test->Train AUC", "Advantage  (AUC_Syn - AUC_Test)"), result = c(AUC_rf_syn, AUC_rf_test, advantage))
  } else {
    MAE_rf_syn <- mean(abs(train[[target]] - preds_rf_syn))
    MAE_rf_test <- mean(abs(train[[target]] - preds_rf_test))
    advantage <- MAE_rf_test - MAE_rf_syn
    result <- data.table(metric = "AIA", target = target, target_info = paste0(target_variable_type, "; range: [", round(target_range[1],2), ", ", round(target_range[2],2), "], mean: ", round(target_mean,2), ", sd: ", round(target_sd,2)),
                         key_vars = paste0(key_vars, collapse = ", "), metric_info = c("Syn->Train MAE", "Test->Train MAE", "Advantage (MAE_Test - MAE_Syn)"), result = c(MAE_rf_syn, MAE_rf_test, advantage))
  }
  
  result
  
}
