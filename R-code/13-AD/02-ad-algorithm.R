# ================================================
# HARF synthesizer function for downstream (ds) prediction
# ================================================
harf_ad_pred <- function(data, job, instance, ...) {
  # -------------------------------
  # Train the HARF model 
  # -------------------------------
  # Split data into training (70%) and testing (30%) sets to assess downstream performance
  metab_data <- as.data.frame(instance$data)
  metab_data <- metab_data[complete.cases(metab_data), ]
  metab_clin_feats = instance$metab_clin_feats
  metab_feats = instance$metab_feats
  evidence = instance$evidence
  train_idx <- instance$train_idx
  print(head(train_idx, n = 30))
  test_idx <- instance$test_idx
  print(dim(metab_data[train_idx, metab_feats]))
  start_time <- Sys.time()
  harf_model <- h_arf(
    omx_data = metab_data[train_idx, metab_feats],
    cli_lab_data = metab_data[train_idx, metab_clin_feats],
    feature_ordering = colnames(metab_data),
    parallel = FALSE,   
    verbose = TRUE,
    target = "Braak_bin3",
    ...
  )
  
  training_time_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  
  # -------------------------------
  # Prepare for data synthesis
  # -------------------------------
  # Extract continuous features for performance measures
  cols_features <- as.data.table(metab_data[, .SD, .SDcols = sapply(metab_data, is.numeric)])
  n_iter <- 1   # use 100 for full run
  results_list <- vector("list", n_iter)
  
  # -------------------------------
  # Synthesize data & compute performance measures
  # -------------------------------
  for (i in 1:n_iter) {
    message("Synthesizing data, iteration ", i, "...")
    iter_start <- Sys.time()
    
    synth_single_cell <- h_forge(
      harf_obj = harf_model,
      n_synth = if (evidence) 1 else nrow(metab_data[train_idx, ]),
      evidence = if (evidence) data.frame(Braak_bin3 = metab_data$Braak_bin3[train_idx]) else NULL,
      parallel = FALSE,   
      verbose = TRUE
    )
    
    iter_end <- Sys.time()
    real_train <- as.data.table(metab_data[train_idx, ])
    synth_classical_data <- as.data.table(synth_single_cell)
    # Compute performance measures
    UVD <- tryCatch(
      univariate_distance(
        real_train = real_train[, ..cols_features],
        syn = synth_classical_data[, ..cols_features]
      ),
      error = function(e) {
        data.table(UVD.metric = NA_real_, 
                   UVD.metric_info = NA_real_, 
                   UVD.pair = NA_real_, 
                   UVD.result = NA_real_)
      }
    )
    CD <- tryCatch(
      fastCor_dist_measure(real_train[, ..cols_features],
                           synth_classical_data[, ..cols_features]),
      error = function(e) {
        print(e)
        data.table(CD.metric = NA_real_, 
                   CD.metric_info = NA_real_, 
                   CD.pair = NA_real_, 
                   CD.result = NA_real_)
      }
    )
    MMD_rbk <- tryCatch(
      MMD(real_train, synth_classical_data),
      error = function(e) {
        print(e)
        data.table(MMD.metric = NA_real_, 
                   MMD.pair = NA_real_, 
                   MMD.result = NA_real_)} 
    )
    
    
    # Retrieving training and testing data
    test_data <- metab_data[-train_idx, ]
    real_train <- as.data.frame(metab_data[train_idx, ])
    # Train RF on original training indices
    rf_org <- ranger(
      x = real_train[, -which(colnames(real_train) == "Braak_bin3")],
      y = as.factor(real_train$Braak_bin3),
      num.trees = 1000, 
      min.node.size = 5,
      probability = TRUE
    )
    # Train RF on synthetic data
    synth_single_cell <- as.data.frame(synth_single_cell)
    rf_synth <- ranger(
      x = synth_single_cell[ , -which(colnames(synth_single_cell) == "Braak_bin3")],
      y = as.factor(synth_single_cell$Braak_bin3),
      num.trees = 5000,
      min.node.size = 5,
      probability = TRUE
    )
    # Predict on test data using both models
    pred_rf_org <- predict(rf_org, data = test_data[, -which(colnames(test_data) == "Braak_bin3")])$predictions[ , 1]
    pred_rf_synth <- predict(rf_synth, data = test_data[, -which(colnames(test_data) == "Braak_bin3")])$predictions[ , 1]
    # Evaluate AUC for RF
    roc_rf_org <- roc(test_data$Braak_bin3, pred_rf_org)
    auc_rf_org <- roc_rf_org$auc
    roc_rf_synth <- roc(test_data$Braak_bin3, pred_rf_synth)
    auc_rf_synth <- roc_rf_synth$auc
    auc_rf_diff <- auc_rf_org - auc_rf_synth
    # Compute MCC for RF
    threshold_rf <- coords(roc_rf_synth, "best", ret = "threshold")
    mcc_rf_org <- mltools::mcc(preds = ifelse(pred_rf_org > threshold_rf$threshold, 1, 0),
                               actuals = test_data$Braak_bin3)
    mcc_rf_synth <- mltools::mcc(preds = ifelse(pred_rf_synth > threshold_rf$threshold, 1, 0),
                                 actuals = test_data$Braak_bin3)
    mcc_rf_diff <- mcc_rf_org - mcc_rf_synth
    
    # Train Lasso on original training indices
    x_org_train <- model.matrix(Braak_bin3 ~ . - 1,
                                data = real_train)
    y_org_train <- real_train$Braak_bin3  # Keep as factor
    y_org_train_numeric <- as.numeric(y_org_train == "Late")  # 1 = Late, 0 = Early
    # Fit Lasso on original data
    lasso_org_model <- cv.glmnet(
      x_org_train,
      y_org_train_numeric,
      alpha = 1,
      family = "binomial"
    )
    # Prepare test set
    test_data <- as.data.frame(test_data)
    y_test <- as.numeric(test_data$Braak_bin3 == "Late")  # 1 = Late, 0 = Early
    x_test <- model.matrix(Braak_bin3 ~ . - 1, data = test_data)
    
    # Predict probabilities
    pred_lasso_org <- predict(lasso_org_model, newx = x_test, s = "lambda.min", type = "response")
    
    # Compute AUC
    roc_lasso_org <- roc(y_test, as.vector(pred_lasso_org))
    auc_lasso_org <- roc_lasso_org$auc
    
    # Compute Matthews Correlation Coefficient (MCC) for Lasso
    threshold_lasso <- coords(lasso_auc , "best", ret = "threshold")
    mcc_lasso_org <- mltools::mcc(preds = ifelse(as.vector(pred_lasso_org) > threshold_lasso$threshold, 1, 0),
                                  actuals = test_data$Braak_bin3)
    
    # Fit lasso on synthetic data
    x_synth_train <- model.matrix(Braak_bin3 ~ . - 1, data = synth_single_cell)
    y_synth_train <- as.numeric(synth_single_cell$Braak_bin3 == "Late")  # 1 = Late, 0 = Early
    lasso_synth_model <- cv.glmnet(
      x_synth_train,
      y_synth_train,
      alpha = 1,
      family = "binomial")
    # Predict probabilities on test set
    pred_lasso_synth <- predict(lasso_synth_model, newx = x_test, s = "lambda.min", type = "response")
    # Compute AUC
    roc_lasso_synth <- roc(y_test, as.vector(pred_lasso_synth))
    auc_lasso_synth <- roc_lasso_synth$auc
    # Compute MCC for Lasso on synthetic data
    threshold_lasso_synth <- coords(roc_lasso_synth , "best", ret = "threshold")
    mcc_lasso_synth <- mltools::mcc(preds = ifelse(as.vector(pred_lasso_synth) > threshold_lasso_synth$threshold, 1, 0),
                                    actuals = test_data$Braak_bin3)
    mcc_lasso_diff <- mcc_lasso_org - mcc_lasso_synth
    
    results_list[[i]] <- data.table(
      Data = instance$data_name,
      iteration = i,
      UVD = UVD,
      CD = CD,
      MMD_rbk = MMD_rbk,
      AUC_RF_ORG = auc_rf_org,
      AUC_RF_SYN = auc_rf_synth,
      AUC_RF_DIFF = auc_rf_diff,
      AUC_Lasso_ORG = auc_lasso_org,
      AUC_Lasso_SYN = auc_lasso_synth,
      AUC_Lasso_DIFF = auc_lasso_org - auc_lasso_synth,
      MCC_RF_ORG = mcc_rf_org,
      MCC_RF_SYN = mcc_rf_synth,
      MCC_RF_DIFF = mcc_rf_diff,
      MCC_Lasso_ORG = mcc_lasso_org,
      time = as.numeric(difftime(iter_end, iter_start, units = "mins")) + training_time_minutes,
      Synthesizer = "HARF"
    )
    # Clean up memory
    rm(synth_single_cell)
    gc()
  }
  
  # Combine results
  return(rbindlist(results_list))
}

# ================================================
# ARF synthesizer function for downstream (ds) prediction
# ================================================
arf_ad_pred <- function(data, job, instance, ...) {
  # --------------------------------------
  # Train the ARF model (single-threaded)
  # --------------------------------------
  metab_data <- as.data.frame(instance$data)
  metab_clin_feats = instance$metab_clin_feats
  metab_feats = instance$metab_feats
  evidence = instance$evidence
  train_idx <- instance$train_idx
  test_idx <- instance$test_idx
  start_time <- Sys.time()
  
  classical_arf <- adversarial_rf(
    x = metab_data[train_idx, ],
    prune = TRUE,
    delta = 0,
    verbose = TRUE,
    parallel = FALSE,   
    ...
  )
  
  classical_forde <- forde(
    classical_arf,
    metab_data[train_idx, ],
    parallel = FALSE
  )
  
  training_time_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  # -------------------------------
  # Prepare for data synthesis
  # -------------------------------
  cols_features <- metab_data[, .SD, .SDcols = sapply(metab_data, is.numeric)]
  n_iter <- 1  # change to 100 for full run
  results_list <- vector("list", n_iter)
  for (i in 1:n_iter) {
    message("Synthesizing data, iteration ", i, "...")
    iter_start <- Sys.time()
    
    synth_classical_data <- forge(
      classical_forde,
      n_synth = if (evidence) 1 else nrow(metab_data[train_idx, ]),
      evidence = if (evidence) data.frame(Braak_bin3 = metab_data$Braak_bin3[train_idx]) else NULL,
      parallel = FALSE   
    )
    
    iter_end <- Sys.time()
    real_train <- as.data.table(metab_data[train_idx, ])
    synth_classical_data <- as.data.table(synth_classical_data)
    # Compute performance measures
    UVD <- tryCatch(
      univariate_distance(
        real_train = real_train[, ..cols_features],
        syn = synth_classical_data[, ..cols_features]
      ),
      error = function(e) {
        print(e)
        data.table(UVD.metric = NA_real_, 
                   UVD.metric_info = NA_real_, 
                   UVD.pair = NA_real_, 
                   UVD.result = NA_real_)
      }
    )
    CD <- tryCatch(
      fastCor_dist_measure(real_train[, ..cols_features],
                           synth_classical_data[, ..cols_features]),
      error = function(e) {
        print(e)
        data.table(CD.metric = NA_real_, 
                   CD.metric_info = NA_real_, 
                   CD.pair = NA_real_, 
                   CD.result = NA_real_)
      }
    )
    MMD_rbk <- tryCatch(
      MMD(real_train[, ..cols_features],
          synth_classical_data[, ..cols_features]),
      error = function(e) {
        print(e)
        data.table(MMD.metric = NA_real_, 
                   MMD.pair = NA_real_, 
                   MMD.result = NA_real_)} 
    )
    
    # Retrieving testing data and clusters
    synth_classical_data <- as.data.frame(synth_classical_data)
    test_data <- metab_data[-train_idx, ]
    
    # Train RF on original training indices
    rf_org <- ranger(
      x = metab_data[train_idx, -which(colnames(metab_data) == "Braak_bin3")],
      y = as.factor(metab_data$Braak_bin3[train_idx]),
      num.trees = 5000,
      min.node.size = 5,
      probability = TRUE
    )
    # Train RF on synthetic data
    rf_synth <- ranger(
      x = synth_classical_data[ , -which(colnames(synth_classical_data) == "Braak_bin3")],
      y = as.factor(synth_classical_data$Braak_bin3),
      num.trees = 5000,
      min.node.size = 5,
      probability = TRUE
    )
    # Predict on test data using both models
    pred_rf_org <- predict(rf_org, data = test_data[, -which(colnames(test_data) == "Braak_bin3")])$predictions[ , 1]
    pred_rf_synth <- predict(rf_synth, data = test_data[, -which(colnames(test_data) == "Braak_bin3")])$predictions[ , 1]
    # Evaluate AUC for RF
    roc_rf_org <- roc(test_data$Braak_bin3, pred_rf_org)
    auc_rf_org <- roc_rf_org$auc
    roc_rf_synth <- roc(test_data$Braak_bin3, pred_rf_synth)
    auc_rf_synth <- roc_rf_synth$auc
    auc_rf_diff <- auc_rf_org - auc_rf_synth
    # Compute MCC for RF
    threshold_rf <- coords(roc_rf_synth , "best", ret = "threshold")
    mcc_rf_org <- mltools::mcc(preds = ifelse(pred_rf_org > threshold_rf$threshold, 1, 0),
                               actuals = test_data$Braak_bin3)
    mcc_rf_synth <- mltools::mcc(preds = ifelse(pred_rf_synth > threshold_rf$threshold, 1, 0),
                                 actuals = test_data$Braak_bin3)
    mcc_rf_diff <- mcc_rf_org - mcc_rf_synth
    
    # Train Lasso on original training indices
    x_org_train <- model.matrix(Braak_bin3 ~ . - 1, data = metab_data[train_idx, ])
    y_org_train <- metab_data$Braak_bin3[train_idx]  # Keep as factor
    y_org_train_numeric <- as.numeric(y_org_train == "Late")  # 1 = Late, 0 = Early
    # Fit Lasso on original data
    lasso_org_model <- cv.glmnet(
      x_org_train,
      y_org_train_numeric,
      alpha = 1,
      family = "binomial"
    )
    # Prepare test set
    x_test <- model.matrix(Braak_bin3 ~ . - 1, data = test_data)
    y_test <- as.numeric(test_data$Braak_bin3 == "Late")  # 1 = Late, 0 = Early
    # Predict probabilities
    pred_lasso_org <- predict(lasso_org_model, newx = x_test, s = "lambda.min", type = "response")
    roc_lasso_org <- roc(y_test, as.vector(pred_lasso_org))
    auc_lasso_org <- roc_lasso_org$auc
    # Fit lasso on synthetic data
    x_synth_train <- model.matrix(Braak_bin3 ~ . - 1, data = synth_classical_data)
    y_synth_train <- as.numeric(synth_classical_data$Braak_bin3 == "Late")  # 1 = Late, 0 = Early
    lasso_synth_model <- cv.glmnet(
      x_synth_train,
      y_synth_train,
      alpha = 1,
      family = "binomial"
    )
    # Predict probabilities on test set
    pred_lasso_synth <- predict(lasso_synth_model, newx = x_test, s = "lambda.min", type = "response")
    roc_lasso_synth <- roc(y_test, as.vector(pred_lasso_synth))
    auc_lasso_synth <- roc_lasso_synth$auc
    # Compute MCC for Lasso on synthetic data
    threshold_lasso_synth <- coords(roc_lasso_synth , "best", ret = "threshold")
    mcc_lasso_org <- mltools::mcc(preds = ifelse(as.vector(pred_lasso_org) > threshold_lasso_synth$threshold, 1, 0),
                                  actuals = test_data$Braak_bin3)
    mcc_lasso_synth <- mltools::mcc(preds = ifelse(as.vector(pred_lasso_synth) > threshold_lasso_synth$threshold, 1, 0),
                                    actuals = test_data$Braak_bin3)
    mcc_lasso_diff <- mcc_lasso_org - mcc_lasso_synth
    
    # Compute AUC
    results_list[[i]] <- data.table(
      Data = instance$data_name,
      iteration = i,
      UVD = UVD,
      CD = CD,
      MMD_rbk = MMD_rbk,
      AUC_RF_ORG = auc_rf_org,
      AUC_RF_SYN = auc_rf_synth,
      AUC_RF_DIFF = auc_rf_diff,
      AUC_Lasso_ORG = auc_lasso_org,
      AUC_Lasso_SYN = auc_lasso_synth,
      AUC_Lasso_DIFF = auc_lasso_org - auc_lasso_synth,
      MCC_RF_ORG = mcc_rf_org,
      MCC_RF_SYN = mcc_rf_synth,
      MCC_RF_DIFF = mcc_rf_diff,
      time = as.numeric(difftime(iter_end, iter_start, units = "mins")) + training_time_minutes,
      Synthesizer = "ARF"
    )
    # Clean up memory
    rm(synth_classical_data)
    gc()
  }
  return(rbindlist(results_list))
}

