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
  registerDoParallel(cores = ncores <- detectCores() - 1)
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
  for (i in 1:100) {
    message("synthesizing data for iteration ", i, "...\n")
    synth_single_cell <- h_forge(
      harf_obj = harf_model,
      n_synth = nrow(instance$data),
      evidence = if (instance$evidence) {
        data.frame(cell_type = instance$data$cell_type)
      } else {NULL},
      parallel = TRUE,
      verbose = TRUE
    )
    end_time <- Sys.time()
    # Evaluate performance measures
    estimated_measures <- NULL
    synth_single_cell <- as.data.frame(synth_single_cell)
    UVD <- univariate_distance(real_train = instance$data[ , which(colnames(instance$data) != "cell_type"), with = FALSE],
                               syn = synth_single_cell[ , which(colnames(synth_single_cell) != "cell_type"), with = FALSE])
    CD <- fastCor_dist_measure(real_train = instance$data, syn = synth_single_cell)
    MMD_rbk <- MMD(real_train = instance$data, syn = synth_single_cell)
    estimated_measures <- rbind(estimated_measures,
                                c(UVD = UVD,
                                  CD = CD,
                                  MMD_rbk = MMD_rbk,
                                  time = as.numeric(difftime(end_time, start_time, units = "mins")),
                                  iteration = i))
  }
  stopCluster(cl)
  return(data.table(Data = instance$data_name, estimated_measures))
}


# ARF synthesizer function
# ========================
# Function to synthesize all Hemberger et al. datasets with ARF
arf_synthesizer <- function (
    dt_name = "",
    org_file_paths,
    syn_file_prfxs,
    parallel = FALSE,
    verbose = FALSE,
    i = 1
) {
  org_file_path <- org_file_paths[dt_name]
  syn_file_prfx <- syn_file_prfxs[dt_name]
  # Read original data
  org_dt <- as.data.frame(fread(org_file_path, check.names = FALSE))
  colnames(org_dt) <- make.names(colnames(org_dt), unique = TRUE)
  n <- nrow(org_dt)
  # Train ARF models
  message("Training ARF model for ", basename(org_file_path), "...\n")
  start_time <- Sys.time()
  classical_arf <- adversarial_rf(
    x = org_dt,
    num_trees = 10,
    min_node_size = 5,
    prune = TRUE,
    delta = 0,
    verbose = verbose
  )
  # forde
  message("Forde for ", basename(org_file_path), "...\n")
  classical_forde <- forde(
    classical_arf,
    org_dt
  )
  # Unconditional synthesis for classical ARF
  message("Synthesizing data for ", basename(org_file_path), "...\n")
  synth_classical_data <- forge(
    classical_forde,
    n_synth = n
  )
  end_time <- Sys.time()
  runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
  fil_nm <- sprintf("%s_arf%02d.csv", syn_file_prfx, i)
  message("Writing synthesized data to ", fil_nm, "...\n")
  fwrite(
    synth_classical_data,
    file = fil_nm
  )
  return(c(dt_name = dt_name,
           runtime = runtime_minutes,
           file_name = fil_nm))
}

# Synpop synthesizer function
# ===========================
# TODO: Implement Synpop synthesizer function
