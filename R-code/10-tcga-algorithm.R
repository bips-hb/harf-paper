# ================================================
# HARF synthesizer function for downstream (ds) prediction
# ================================================
harf_ds_pred <- function(data, job, instance, ...) {
  # -------------------------------
  # Train the HARF model (single-threaded)
  # -------------------------------
  message("Training HARF model for ", instance$data_name, "...")
  # Test with the first 50 features and cell type for faster training
  # Split data into training (70%) and testing (30%) sets to assess intergeanrability
  # instance <- list(data = as.data.frame(fread(org_hember_dt_files["manno"]))) 
  # instance$evidence <- FALSE
  # train_idx <- sample(seq_len(nrow(instance$data)), size = floor(0.7 * nrow(instance$data)), replace = FALSE)
  # instance$data <- instance$data[, c(1:200, which(colnames(instance$data) == "tumor_stage"))]
  instance$data <- as.data.frame(instance$data)
  train_idx <- instance$train_idx
  start_time <- Sys.time()
  harf_model <- h_arf(
    omx_data = instance$data[train_idx, -which(colnames(instance$data) == "tumor_stage")],
    cli_lab_data = data.frame(tumor_stage = instance$data$tumor_stage[train_idx]),
    feature_ordering = colnames(instance$data),
    parallel = FALSE,   
    verbose = TRUE,
    target = "tumor_stage",
    ...
  )
  
  training_time_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  
  # -------------------------------
  # Prepare for data synthesis
  # -------------------------------
  cols_features <- setdiff(colnames(instance$data), c("tumor_stage"))
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
      n_synth = if (instance$evidence) 1 else nrow(instance$data[train_idx, ]),
      evidence = if (instance$evidence) data.frame(tumor_stage = instance$data$tumor_stage[train_idx]) else NULL,
      parallel = FALSE,   
      verbose = TRUE
    )
    
    iter_end <- Sys.time()
    real_train <- as.data.table(instance$data[train_idx, ])
    synth_classical_data <- as.data.table(synth_single_cell)
    # Compute performance measures
    print("I am here 1...")
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
    print("I am here 2...")
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
    print("I am here 3...")
    MMD_rbk <- tryCatch(
      MMD(real_train, synth_classical_data),
      error = function(e) {
        print(e)
        data.table(MMD.metric = NA_real_, 
                   MMD.pair = NA_real_, 
                   MMD.result = NA_real_)} 
    )
    
    
    # Retrieving training and testing data
    test_data <- instance$data[-train_idx, ]
    real_train <- as.data.frame(instance$data[train_idx, ])
    # Train RF on original training indices
    rf_org <- ranger(
      x = real_train[, -which(colnames(real_train) == "tumor_stage")],
      y = as.factor(real_train$tumor_stage),
      num.trees = 1000, 
      min.node.size = 5,
      probability = TRUE
    )
    # Train RF on synthetic data
    synth_single_cell <- as.data.frame(synth_single_cell)
    rf_synth <- ranger(
      x = synth_single_cell[ , -which(colnames(synth_single_cell) == "tumor_stage")],
      y = as.factor(synth_single_cell$tumor_stage),
      num.trees = 1000,
      min.node.size = 5,
      probability = TRUE
    )
    # Predict on test data using both models
    pred_rf_org <- predict(rf_org, data = test_data[, -which(colnames(test_data) == "tumor_stage")])$predictions[ , 1]
    pred_rf_synth <- predict(rf_synth, data = test_data[, -which(colnames(test_data) == "tumor_stage")])$predictions[ , 1]
    # Evaluate AUC for RF
    auc_rf_org <- roc(test_data$tumor_stage, pred_rf_org)$auc
    auc_rf_synth <- roc(test_data$tumor_stage, pred_rf_synth)$auc
    auc_rf_diff <- auc_rf_org - auc_rf_synth
    
    # Train Lasso on original training indices
    x_org_train <- model.matrix(tumor_stage ~ . - 1,
                                data = real_train)
    y_org_train <- real_train$tumor_stage  # Keep as factor
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
    y_test <- as.numeric(test_data$tumor_stage == "Late")  # 1 = Late, 0 = Early
    x_test <- model.matrix(tumor_stage ~ . - 1, data = test_data)
    
    # Predict probabilities
    pred_lasso_org <- predict(lasso_org_model, newx = x_test, s = "lambda.min", type = "response")
    
    # Compute AUC
    auc_lasso_org <- roc(y_test, as.vector(pred_lasso_org))$auc
    
    # Fit lasso on synthetic data
    x_synth_train <- model.matrix(tumor_stage ~ . - 1, data = synth_single_cell)
    y_synth_train <- as.numeric(synth_single_cell$tumor_stage == "Late")  # 1 = Late, 0 = Early
    lasso_synth_model <- cv.glmnet(
      x_synth_train,
      y_synth_train,
      alpha = 1,
      family = "binomial")
    # Predict probabilities on test set
    pred_lasso_synth <- predict(lasso_synth_model, newx = x_test, s = "lambda.min", type = "response")
    auc_lasso_synth <- roc(y_test, as.vector(pred_lasso_synth))$auc
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
arf_ds_pred <- function(data, job, instance, ...) {
  # --------------------------------------
  # Train the ARF model (single-threaded)
  # --------------------------------------
  train_idx <- instance$train_idx
  message("Training ARF model for ", instance$data_name, "...")
  start_time <- Sys.time()
  
  classical_arf <- adversarial_rf(
    x = instance$data[train_idx, ],
    prune = TRUE,
    delta = 0,
    verbose = TRUE,
    parallel = FALSE,   
    ...
  )
  
  classical_forde <- forde(
    classical_arf,
    instance$data[train_idx, ],
    parallel = FALSE
  )
  
  training_time_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  # -------------------------------
  # Prepare for data synthesis
  # -------------------------------
  cols_features <- setdiff(colnames(instance$data), "tumor_stage")
  n_iter <- 1  # change to 100 for full run
  results_list <- vector("list", n_iter)
  for (i in 1:n_iter) {
    message("Synthesizing data, iteration ", i, "...")
    iter_start <- Sys.time()
    
    synth_classical_data <- forge(
      classical_forde,
      n_synth = if (instance$evidence) 1 else nrow(instance$data[train_idx, ]),
      evidence = if (instance$evidence) data.frame(tumor_stage = instance$data$tumor_stage[train_idx]) else NULL,
      parallel = FALSE   
    )
    
    iter_end <- Sys.time()
    real_train <- as.data.table(instance$data[train_idx, ])
    synth_classical_data <- as.data.table(synth_classical_data)
    # Compute performance measures
    print("I am here 1...")
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
    print("I am here 2...")
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
    print("I am here 3...")
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
    test_data <- instance$data[-train_idx, ]
    
    # Train RF on original training indices
    rf_org <- ranger(
      x = instance$data[train_idx, -which(colnames(instance$data) == "tumor_stage")],
      y = as.factor(instance$data$tumor_stage[train_idx]),
      num.trees = 1000,
      min.node.size = 5,
      probability = TRUE
    )
    # Train RF on synthetic data
    rf_synth <- ranger(
      x = synth_classical_data[ , -which(colnames(synth_classical_data) == "tumor_stage")],
      y = as.factor(synth_classical_data$tumor_stage),
      num.trees = 1000,
      min.node.size = 5,
      probability = TRUE
    )
    # Predict on test data using both models
    pred_rf_org <- predict(rf_org, data = test_data[, -which(colnames(test_data) == "tumor_stage")])$predictions[ , 1]
    pred_rf_synth <- predict(rf_synth, data = test_data[, -which(colnames(test_data) == "tumor_stage")])$predictions[ , 1]
    # Evaluate AUC for RF
    auc_rf_org <- roc(test_data$tumor_stage, pred_rf_org)$auc
    auc_rf_synth <- roc(test_data$tumor_stage, pred_rf_synth)$auc
    auc_rf_diff <- auc_rf_org - auc_rf_synth
    # Train Lasso on original training indices
    x_org_train <- model.matrix(tumor_stage ~ . - 1, data = instance$data[train_idx, ])
    y_org_train <- instance$data$tumor_stage[train_idx]  # Keep as factor
    y_org_train_numeric <- as.numeric(y_org_train == "Late")  # 1 = Late, 0 = Early
    # Fit Lasso on original data
    lasso_org_model <- cv.glmnet(
      x_org_train,
      y_org_train_numeric,
      alpha = 1,
      family = "binomial"
    )
    # Prepare test set
    x_test <- model.matrix(tumor_stage ~ . - 1, data = test_data)
    y_test <- as.numeric(test_data$tumor_stage == "Late")  # 1 = Late, 0 = Early
    # Predict probabilities
    pred_lasso_org <- predict(lasso_org_model, newx = x_test, s = "lambda.min", type = "response")
    auc_lasso_org <- roc(y_test, as.vector(pred_lasso_org))$auc
    # Fit lasso on synthetic data
    x_synth_train <- model.matrix(tumor_stage ~ . - 1, data = synth_classical_data)
    y_synth_train <- as.numeric(synth_classical_data$tumor_stage == "Late")  # 1 = Late, 0 = Early
    lasso_synth_model <- cv.glmnet(
      x_synth_train,
      y_synth_train,
      alpha = 1,
      family = "binomial"
    )
    # Predict probabilities on test set
    pred_lasso_synth <- predict(lasso_synth_model, newx = x_test, s = "lambda.min", type = "response")
    auc_lasso_synth <- roc(y_test, as.vector(pred_lasso_synth))$auc
    
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
      time = as.numeric(difftime(iter_end, iter_start, units = "mins")) + training_time_minutes,
      Synthesizer = "ARF"
    )
    # Clean up memory
    rm(synth_classical_data)
    gc()
  }
  return(rbindlist(results_list))
}

