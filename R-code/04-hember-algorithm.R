# HARF synthesizer function
# ==========================
# Function to synthhesize Hemberger et al. datasets with harf

harf_synthesizer <- function(data, job, instance, ...) {
  # -------------------------------
  # Train the HARF model (single-threaded)
  # -------------------------------
  message("Training HARF model for ", instance$data_name, "...")
  # Test with the first 50 features and cell type for faster training
  # instance$data <- instance$data[, c(1:50, which(colnames(instance$data) == "cell_type"))]
  start_time <- Sys.time()
  
  harf_model <- h_arf(
    omx_data = instance$data[, -which(colnames(instance$data) == "cell_type")],
    cli_lab_data = data.frame(cell_type = instance$data$cell_type),
    feature_ordering = colnames(instance$data),
    parallel = FALSE,   
    verbose = TRUE,
    ...
  )
  
  training_time_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  
  # -------------------------------
  # Prepare for data synthesis
  # -------------------------------
  cols_features <- setdiff(colnames(instance$data), "cell_type")
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
      n_synth = if (instance$evidence) 1 else nrow(instance$data),
      evidence = if (instance$evidence) data.frame(cell_type = instance$data$cell_type) else NULL,
      parallel = FALSE,   
      verbose = TRUE
    )
    
    iter_end <- Sys.time()
    
    real_train <- as.data.table(instance$data)
    
    # Compute performance measures
    UVD <- univariate_distance(
      real_train = real_train[, ..cols_features],
      syn = synth_single_cell[, ..cols_features]
    )
    
    CD <- tryCatch(
      fastCor_dist_measure(real_train = real_train, syn = synth_single_cell),
      error = function(e) {
        print(e)
        data.table(CD.metric = NA_real_, 
                   CD.metric_info = NA_real_, 
                   CD.pair = NA_real_, 
                   CD.result = NA_real_)
      }
    )
    
    MMD_rbk <- tryCatch(
      MMD(real_train = real_train, syn = synth_single_cell),
      error = function(e) {
        print(e)
        data.table(MMD.metric = NA_real_, 
                   MMD.metric_info = NA_real_, 
                   MMD.pair = NA_real_, 
                   MMD.result = NA_real_)
      }
    )
    
    ari_nmi_org <- tryCatch(
      cluster_and_eval(sc_data = instance$data),
      error = function(e) {
        print(e)
        c(ARI = NA_real_, NMI = NA_real_)
      }
    )
    
    ari_nmi_harf <- tryCatch(
      cluster_and_eval(sc_data = synth_single_cell),
      error = function(e) {
        print(e)
        c(ARI = NA_real_, NMI = NA_real_)
      }
    )
    ari_nmi_diff <- abs(ari_nmi_org - ari_nmi_harf)
    
    results_list[[i]] <- data.table(
      Data = instance$data_name,
      iteration = i,
      UVD = UVD,
      CD = CD,
      MMD_rbk = MMD_rbk,
      ARI_ORG = ari_nmi_org["ARI"],
      NMI_ORG = ari_nmi_org["NMI"],
      ARI_HARF = ari_nmi_harf["ARI"],
      NMI_HARF = ari_nmi_harf["NMI"],
      ARI_DIFF = ari_nmi_diff["ARI"],
      NMI_DIFF = ari_nmi_diff["NMI"],
      time = as.numeric(difftime(iter_end, iter_start, units = "mins")) + training_time_minutes
    )
    # Clean up memory
    rm(synth_single_cell)
    gc()
  }
  
  # Combine results
  return(rbindlist(results_list))
}


# ARF synthesizer function
# ========================
# Function to synthesize all Hemberger et al. datasets with ARF
# See harf_synthesizer for argument details.
arf_synthesizer <- function(data, job, instance, ...) {
  # -------------------------------
  # Train the ARF model (single-threaded)
  # -------------------------------
  message("Training ARF model for ", instance$data_name, "...")
  start_time <- Sys.time()
  
  classical_arf <- adversarial_rf(
    x = instance$data,
    prune = TRUE,
    delta = 0,
    verbose = TRUE,
    parallel = FALSE,   
    ...
  )
  
  classical_forde <- forde(
    classical_arf,
    instance$data
  )
  
  training_time_minutes <- as.numeric(difftime(Sys.time(), Sys.time(), units = "mins"))
  
  # -------------------------------
  # Prepare for data synthesis
  # -------------------------------
  cols_features <- setdiff(colnames(instance$data), "cell_type")
  n_iter <- 1  # change to 100 for full run
  results_list <- vector("list", n_iter)
  
  # -------------------------------
  # Synthesize data & compute performance measures
  # -------------------------------
  for (i in 1:n_iter) {
    message("Synthesizing data, iteration ", i, "...")
    iter_start <- Sys.time()
    
    synth_classical_data <- forge(
      classical_forde,
      n_synth = if (instance$evidence) 1 else nrow(instance$data),
      evidence = if (instance$evidence) data.frame(cell_type = instance$data$cell_type) else NULL,
      parallel = FALSE   
    )
    
    iter_end <- Sys.time()
    real_train <- as.data.table(instance$data)
    synth_classical_data <- as.data.table(synth_classical_data)
    
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
      fastCor_dist_measure(real_train, synth_classical_data),
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
                   MMD.metric_info = NA_real_, 
                   MMD.pair = NA_real_, 
                   MMD.result = NA_real_)} 
    )
    
    ari_nmi_org <- tryCatch(
      cluster_and_eval(sc_data = instance$data),
      error = function(e) {
        print(e)
        c(ARI = NA_real_, NMI = NA_real_)
      }
    )
    
    ari_nmi_arf <- tryCatch(
      cluster_and_eval(sc_data = synth_classical_data),
      error = function(e) {
        print(e)
        c(ARI = NA_real_, NMI = NA_real_)
      }
    )
    ari_nmi_diff <- abs(ari_nmi_org - ari_nmi_arf)
    
    
    results_list[[i]] <- data.table(
      Data = instance$data_name,
      iteration = i,
      UVD = UVD,
      CD = CD,
      MMD_rbk = MMD_rbk,
      ARI_ORG = ari_nmi_org["ARI"],
      NMI_ORG = ari_nmi_org["NMI"],
      ARI_ARF = ari_nmi_arf["ARI"],
      NMI_ARF = ari_nmi_arf["NMI"],
      ARI_DIFF = ari_nmi_diff["ARI"],
      NMI_DIFF = ari_nmi_diff["NMI"],
      time = as.numeric(difftime(iter_end, iter_start, units = "mins")) + training_time_minutes
    )
    
    # Clean up memory
    rm(synth_classical_data)
    gc()
  }
  
  # Combine results
  return(rbindlist(results_list))
}


# ================================================
# HARF synthesizer function for interchangeability
# ================================================
harf_interchangeable <- function(data, job, instance, ...) {
  # -------------------------------
  # Train the HARF model (single-threaded)
  # -------------------------------
  message("Training HARF model for ", instance$data_name, "...")
  # Test with the first 50 features and cell type for faster training
  # Split data into training (70%) and testing (30%) sets to assess intergeanrability
  # instance <- list(data = as.data.frame(fread(org_hember_dt_files["manno"]))) 
  # instance$evidence <- FALSE
  # train_idx <- sample(seq_len(nrow(instance$data)), size = floor(0.7 * nrow(instance$data)), replace = FALSE)
  # instance$data <- instance$data[, c(1:200, which(colnames(instance$data) == "cell_type"))]
  instance$data <- as.data.frame(instance$data)
  train_idx <- instance$train_idx
  start_time <- Sys.time()
  harf_model <- h_arf(
    omx_data = instance$data[train_idx, -which(colnames(instance$data) == "cell_type")],
    cli_lab_data = data.frame(cell_type = instance$data$cell_type[train_idx]),
    feature_ordering = colnames(instance$data),
    parallel = FALSE,   
    verbose = TRUE,
    ...
  )
  
  training_time_minutes <- as.numeric(difftime(Sys.time(), start_time, units = "mins"))
  
  # -------------------------------
  # Prepare for data synthesis
  # -------------------------------
  cols_features <- setdiff(colnames(instance$data), "cell_type")
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
      n_synth = if (instance$evidence) 1 else nrow(instance$data),
      evidence = if (instance$evidence) data.frame(cell_type = instance$data$cell_type) else NULL,
      parallel = FALSE,   
      verbose = TRUE
    )
    
    iter_end <- Sys.time()
    synth_single_cell <- as.data.frame(synth_single_cell)
    
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
      fastCor_dist_measure(real_train, synth_classical_data),
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
                   MMD.metric_info = NA_real_, 
                   MMD.pair = NA_real_, 
                   MMD.result = NA_real_)} 
    )
    
    
    # Clustering test data with SC3 to evaluate interchangeability
    test_data <- instance$data[-train_idx, ]
    test_clusters <- cluster_it(test_data)
    
    # Train RF on original training indices
    rf_org <- ranger(
      x = instance$data[train_idx, -which(colnames(instance$data) == "cell_type")],
      y = as.factor(instance$data$cell_type[train_idx]),
      num.trees = 1000, 
      min.node.size = 5
    )
    # Train RF on synthetic data
    rf_synth <- ranger(
      x = synth_single_cell[ , -which(colnames(synth_single_cell) == "cell_type")],
      y = as.factor(synth_single_cell$cell_type),
      num.trees = 1000,
      min.node.size = 5
    )
    # Predict on test data using both models
    pred_org <- predict(rf_org, data = test_data[, -which(colnames(test_data) == "cell_type")])$predictions
    pred_synth <- predict(rf_synth, data = test_data[, -which(colnames(test_data) == "cell_type")])$predictions
    # Evaluate ARI and NMI
    ari_org <- adjustedRandIndex(test_clusters, pred_org)
    ari_synth <- adjustedRandIndex(test_clusters, pred_synth)
    ari_diff <- ari_org - ari_synth
    nmi_org <- NMI(test_clusters, pred_org)
    nmi_synth <- NMI(test_clusters, pred_synth)
    nmi_diff <- nmi_org - nmi_synth
    ari_nmi_org <- c(ARI = ari_org, NMI = nmi_org)
    ari_nmi_syn <- c(ARI = ari_synth, NMI = nmi_synth)
    ari_nmi_diff <- c(ARI = ari_diff, NMI = nmi_diff)
    
    
    results_list[[i]] <- data.table(
      Data = instance$data_name,
      iteration = i,
      UVD = UVD,
      CD = CD,
      MMD_rbk = MMD_rbk,
      ARI_ORG = ari_nmi_org["ARI"],
      NMI_ORG = ari_nmi_org["NMI"],
      ARI_HARF = ari_nmi_syn["ARI"],
      NMI_HARF = ari_nmi_syn["NMI"],
      ARI_DIFF = ari_nmi_diff["ARI"],
      NMI_DIFF = ari_nmi_diff["NMI"],
      time = as.numeric(difftime(iter_end, iter_start, units = "mins")) + training_time_minutes
    )
    # Clean up memory
    rm(synth_single_cell)
    gc()
  }
  
  # Combine results
  return(rbindlist(results_list))
}

# ================================================
# ARF synthesizer function for interchangeability
# ================================================
arf_interchangeable <- function(data, job, instance, ...) {
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
  
  training_time_minutes <- as.numeric(difftime(Sys.time(), Sys.time(), units = "mins"))
  # -------------------------------
  # Prepare for data synthesis
  # -------------------------------
  cols_features <- setdiff(colnames(instance$data), "cell_type")
  n_iter <- 1  # change to 100 for full run
  results_list <- vector("list", n_iter)
  for (i in 1:n_iter) {
    message("Synthesizing data, iteration ", i, "...")
    iter_start <- Sys.time()
    
    synth_classical_data <- forge(
      classical_forde,
      n_synth = if (instance$evidence) 1 else nrow(instance$data[train_idx, ]),
      evidence = if (instance$evidence) data.frame(cell_type = instance$data[train_idx, ]$cell_type) else NULL,
      parallel = FALSE   
    )
    
    iter_end <- Sys.time()
    synth_classical_data <- as.data.frame(synth_classical_data)
    
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
      fastCor_dist_measure(real_train, synth_classical_data),
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
                   MMD.metric_info = NA_real_, 
                   MMD.pair = NA_real_, 
                   MMD.result = NA_real_)} 
    )
    
    # Clustering test data with SC3 to evaluate interchangeability
    test_data <- instance$data[-train_idx, ]
    test_clusters <- cluster_it(test_data)
    
    # Train RF on original data
    rf_org <- ranger(
      x = instance$data[train_idx, -which(colnames(instance$data) == "cell_type")],
      y = as.factor(instance$data$cell_type[train_idx]),
      num.trees = 1000, 
      min.node.size = 5
    )
    # Train RF on synthetic data
    rf_synth <- ranger(
      x = synth_classical_data[, -which(colnames(synth_classical_data) == "cell_type")],
      y = as.factor(synth_classical_data$cell_type),
      num.trees = 1000,
      min.node.size = 5
    )
    # Predict on test data using both models
    pred_org <- predict(rf_org, data = test_data[, -which(colnames(test_data) == "cell_type")])$predictions
    pred_synth <- predict(rf_synth, data = test_data[, -which(colnames(test_data) == "cell_type")])$predictions
    # Evaluate ARI and NMI
    ari_org <- adjustedRandIndex(test_clusters, pred_org)
    ari_synth <- adjustedRandIndex(test_clusters, pred_synth)
    ari_diff <- ari_org - ari_synth
    nmi_org <- NMI(test_clusters, pred_org)
    nmi_synth <- NMI(test_clusters, pred_synth)
    nmi_diff <- nmi_org - nmi_synth
    ari_nmi_org <- c(ARI = ari_org, NMI = nmi_org)
    ari_nmi_syn <- c(ARI = ari_synth, NMI = nmi_synth)
    ari_nmi_diff <- c(ARI = ari_diff, NMI = nmi_diff)
    return(data.table(
      Data = instance$data_name,
      iteration = i,
      UVD = UVD,
      CD = CD,
      MMD_rbk = MMD_rbk,
      ARI_ORG = ari_nmi_org["ARI"],
      NMI_ORG = ari_nmi_org["NMI"],
      ARI_ARF = ari_nmi_syn["ARI"],
      NMI_ARF = ari_nmi_syn["NMI"],
      ARI_DIFF = ari_nmi_diff["ARI"],
      NMI_DIFF = ari_nmi_diff["NMI"],
      time = as.numeric(difftime(iter_end, iter_start, units = "mins")) + training_time_minutes
    ))
    rm(synth_classical_data)
    gc()
  }
  return(rbindlist(results_list))
}
