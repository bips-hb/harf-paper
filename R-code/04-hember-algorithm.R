# HARF synthesizer function
# ==========================
# Function to synthhesize Hemberger et al. datasets with harf
harf_synthesizer <- function (
    data,
    job,
    instance,
    ...
) {
  # Parallelize
  ncores <- detectCores() - 1
  cl <- makeCluster(ncores)
  registerDoParallel(cores = ncores)
  # Train the HARF model
  message("Training HARF model for ", instance$data_name, "...")
  start_time <- Sys.time()
  harf_model <- h_arf(
    omx_data = instance$data[, -which(colnames(instance$data) == "cell_type")],
    cli_lab_data = data.frame(cell_type = instance$data$cell_type),
    feature_ordering = colnames(instance$data),
    parallel = TRUE,
    verbose = TRUE,
    ...
  )
  # Repeat synthesize data and evaluate performance measures for 100 iterations
  estimated_measures <- NULL
  for (i in 1:10) {
    # Only 10 iterations for testing, change to 100 for final submission
    message("synthesizing data for iteration ", i, "...\n")
    synth_single_cell <- h_forge(
      harf_obj = harf_model,
      n_synth = if (instance$evidence) {1} else {nrow(instance$data)},
      evidence = if (instance$evidence) {data.frame(cell_type = instance$data$cell_type)} else {NULL},
      parallel = TRUE,
      verbose = TRUE
    )
    end_time <- Sys.time()
    # Evaluate performance measures
    real_train <- as.data.table(instance$data)
    UVD <- univariate_distance(real_train = real_train[ , which(colnames(instance$data) != "cell_type"), with = FALSE],
                               syn = synth_single_cell[ , which(colnames(synth_single_cell) != "cell_type"), with = FALSE])
    CD <- fastCor_dist_measure(real_train = real_train, syn = synth_single_cell)
    MMD_rbk <- MMD(real_train = real_train, syn = synth_single_cell)
    estimated_measures <- rbind(estimated_measures,
                                c(UVD = UVD,
                                  CD = CD,
                                  MMD_rbk = MMD_rbk,
                                  time = as.numeric(difftime(end_time, start_time, units = "mins")),
                                  iteration = i))
    rm(synth_single_cell)
    gc()
  }
  stopCluster(cl)
  return(data.table(Data = instance$data_name, estimated_measures))
}


# ARF synthesizer function
# ========================
# Function to synthesize all Hemberger et al. datasets with ARF
arf_synthesizer <- function (
    data,
    job,
    instance,
    ...
) {
  # Parallelize
  ncores <- detectCores() - 1
  cl <- makeCluster(ncores)
  registerDoParallel(cores = ncores)
  # Train the HARF model
  message("Training ARF model for ", instance$data_name, "...")
  start_time <- Sys.time()
  classical_arf <- adversarial_rf(
    x = instance$data,
    prune = TRUE,
    delta = 0,
    verbose = TRUE,
    parallel = TRUE,
    ...
  )
  # forde
  message("synthesizing data for iteration ", i, "...\n")
  classical_forde <- forde(
    classical_arf,
    instance$data
  )
  # Repeat synthesize data and evaluate performance measures for 100 iterations
  estimated_measures <- NULL
  for (i in 1:10) {
    # Only 10 iterations for testing, change to 100 for final submission
    message("synthesizing data for iteration ", i, "...\n")
    synth_classical_data <- forge(
      classical_forde,
      n_synth = if (instance$evidence) {1} else {nrow(instance$data)},
      evidence = if (instance$evidence) {data.frame(cell_type = instance$data$cell_type)} else {NULL},
      parallel = TRUE
    )
     end_time <- Sys.time()
     # Evaluate performance measures
    real_train <- as.data.table(instance$data)
    UVD <- univariate_distance(real_train = real_train[ , which(colnames(instance$data) != "cell_type"), with = FALSE],
    syn = synth_classical_data[ , which(colnames(synth_classical_data) != "cell_type"), with = FALSE])
    CD <- fastCor_dist_measure(real_train = real_train, syn = synth_classical_data)
    MMD_rbk <- MMD(real_train = real_train, syn = synth_classical)
    estimated_measures <- rbind(estimated_measures,
                                c(UVD = UVD,
                                  CD = CD,
                                  MMD_rbk = MMD_rbk,
                                  time = as.numeric(difftime(end_time, start_time, units = "mins")),
                                  iteration = i))
    rm(synth_classical_data)
    gc()
  }
  stopCluster(cl)
  return(data.table(Data = instance$data_name, estimated_measures))
}
