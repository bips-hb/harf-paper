# HARF synthesizer function
# ==========================
# Function to synthhesize TCGA datasets with harf
harf_tcga_synthesizer <- function (
    dt_name = "",
    org_file_paths,
    syn_file_prfxs,
    chunck_size = 10,
    parallel = FALSE,
    verbose = FALSE,
    N_synth = 10
) {
  org_file_path <- org_file_paths[dt_name]
  syn_file_prfx <- syn_file_prfxs[dt_name]
  # Read original data
  org_dt <- as.data.frame(fread(org_file_path, check.names = FALSE))
  org_dt$patientID <- NULL
  colnames(org_dt) <- make.names(colnames(org_dt), unique = TRUE)
  n <- nrow(org_dt)
  print(dim(org_dt))
  # Train the HARF model
  message("Training HARF model for ", basename(org_file_path), "...")
  start_time <- Sys.time()
  harf_model <- h_arf(
    omx_data = org_dt[, -which(names(org_dt) %in%
                                 c("years_to_birth",
                                   "age_at_diagnosis",
                                   "gender",
                                   "tumor_stage"))],
    cli_lab_data = org_dt[, c("years_to_birth",
                              "age_at_diagnosis",
                              "gender",
                              "tumor_stage")],
    parallel = parallel,
    chunck_size = chunck_size,
    num_clusters = floor((ncol(org_dt) - 1) / chunck_size),
    verbose = verbose
  )
  # Synthesize data
  message("Synthesizing data for ", basename(org_file_path), "...")
  # Generate N_synth datasets
  run_time <- lapply(1:N_synth, function(i) {
    message("Synthesis iteration ", i, " of ", N_synth, "...\n")
    # Conditional resampling on tumor_stage
    sub_harf_synth_list <- lapply(unique(org_dt$tumor_stage), function(stage_level) {
      message("Synthesizing for tumor_stage: ", stage_level)
      sub_harf_synth <- h_forge(
        harf_obj = harf_model,
        n_synth = sum(org_dt$tumor_stage == stage_level),
        evidence = data.frame(tumor_stage = stage_level),
        verbose = TRUE,
        parallel = FALSE
      )
      return(sub_harf_synth)
    })
    synth_tumor <- do.call(rbind, sub_harf_synth_list)

    end_time <- Sys.time()
    runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
    fil_nm <- sprintf("%s_harf%02d.csv", syn_file_prfx, i)
    message("Writing synthesized data to ", fil_nm, "...\n")
    fwrite(
      synth_tumor,
      file = fil_nm
    )
    return(c(dt_name = dt_name, runtime = runtime_minutes, file_name = fil_nm))
  })
  return(do.call(rbind, run_time))
}


# ARF synthesizer function
# ========================
# Function to synthesize all Hemberger et al. datasets with ARF
arf_tcga_synthesizer <- function (
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
  org_dt$patientID <- NULL
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
  # Conditional resampling on tumor_stage
  sub_harf_synth_list <- lapply(unique(org_dt$tumor_stage), function(stage_level) {
    message("Synthesizing for tumor_stage: ", stage_level)
    sub_harf_synth <- forge(
      classical_forde,
      n_synth = 1,
      evidence = data.frame(tumor_stage = arg_dt$tumor_stage),
      verbose = TRUE,
      parallel = FALSE
    )
    return(sub_harf_synth)
  })
  synth_tumor <- do.call(rbind, sub_harf_synth_list)
  end_time <- Sys.time()
  runtime_minutes <- as.numeric(difftime(end_time, start_time, units = "mins"))
  fil_nm <- sprintf("%s_arf%02d.csv", syn_file_prfx, i)
  message("Writing synthesized data to ", fil_nm, "...\n")
  fwrite(
    synth_tumor,
    file = fil_nm
  )
  return(c(dt_name = dt_name,
           runtime = runtime_minutes,
           file_name = fil_nm))
}

# Synpop synthesizer function
# ===========================
# TODO: Implement Synpop synthesizer function
