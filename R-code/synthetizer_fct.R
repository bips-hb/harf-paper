# HARF synthesizer function
# ==========================
# Function to synthhesize Hemberger et al. datasets with harf
harf_synthesizer <- function (
    dt_name = "",
    org_file_paths,
    syn_file_prfxs,
    chunck_size = 10,
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
  # Train the HARF model
  message("Training HARF model for ", basename(org_file_path), "...")
  start_time <- Sys.time()
  harf_model <- h_arf(
    omx_data = org_dt[, -which(colnames(org_dt) == "cell_type")],
    cli_lab_data = data.frame(cell_type = org_dt$cell_type),
    parallel = parallel,
    chunck_size = chunck_size,
    num_clusters = floor((ncol(org_dt) - 1) / chunck_size),
    verbose = verbose
  )
  rm(org_dt)
  # Synthesize data
  message("Synthesizing data for ", basename(org_file_path), "...")
  synth_single_cell <- h_forge(
    harf_obj = harf_model,
    n_synth = n,
    evidence = NULL,
    parallel = parallel,
    verbose = verbose
  )
  end_time <- Sys.time()
  runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
  fil_nm <- sprintf("%s_harf%02d.csv", syn_file_prfx, i)
  message("Writing synthesized data to ", fil_nm, "...\n")
  fwrite(
    synth_single_cell,
    file = fil_nm
  )
  return(c(dt_name = dt_name, runtime = runtime_minutes, file_name = fil_nm))
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
