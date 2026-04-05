# HARF synthesizer function
# ==========================
# Function to synthhesize Hemberger et al. datasets with harf

harf_synthesizer <- function(data, job, instance, ...) {
  # -------------------------------
  # Train the HARF model (single-threaded)
  # -------------------------------
  message("Training HARF model for ", instance$data_name, "...")
  # Test with the first 50 features and cell type for faster training
  instance$data <- instance$data[, c(1:50, which(colnames(instance$data) == "cell_type"))]
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
      error = function(e) NA_real_
    )
    
    MMD_rbk <- tryCatch(
      MMD(real_train = real_train, syn = synth_single_cell),
      error = function(e) NA_real_
    )
    
    ari_nmi_org <- tryCatch(
      cluster_and_eval(sc_data = instance$data),
      error = function(e) c(ARI = NA_real_, NMI = NA_real_)
    )

    ari_nmi_harf <- tryCatch(
      cluster_and_eval(sc_data = synth_single_cell),
      error = function(e) c(ARI = NA_real_, NMI = NA_real_)
    )
    print(ari_nmi_org)
    print(ari_nmi_harf)
    ari_nmi_diff <- abs(ari_nmi_org - ari_nmi_harf)
    
    results_list[[i]] <- data.table(
      Data = instance$data_name,
      iteration = i,
      UVD = UVD,
      CD = CD,
      MMD_rbk = MMD_rbk,
      ARI = ari_nmi_diff["ARI"],
      NMI = ari_nmi_diff["NMI"],
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
      error = function(e) NA_real_
    )
    
    CD <- tryCatch(
      fastCor_dist_measure(real_train, synth_classical_data),
      error = function(e) NA_real_
    )
    MMD_rbk <- tryCatch(
      MMD(real_train, synth_classical_data),
      error = function(e) NA_real_
    )
    
    ari_nmi_org <- tryCatch(
      cluster_and_eval(sc_data = instance$data),
      error = function(e) c(ARI = NA_real_, NMI = NA_real_)
    )
    
    ari_nmi_arf <- tryCatch(
      cluster_and_eval(sc_data = synth_classical_data),
      error = function(e) c(ARI = NA_real_, NMI = NA_real_)
    )
    ari_nmi_diff <- abs(ari_nmi_org - ari_nmi_arf)
       
    
    results_list[[i]] <- data.table(
      Data = instance$data_name,
      iteration = i,
      UVD = UVD,
      CD = CD,
      MMD_rbk = MMD_rbk,
      ARI = ari_nmi_diff["ARI"],
      NMI = ari_nmi_diff["NMI"],
      time = as.numeric(difftime(iter_end, iter_start, units = "mins")) + training_time_minutes
    )
    
    # Clean up memory
    rm(synth_classical_data)
    gc()
  }
  
  # Combine results
  return(rbindlist(results_list))
}