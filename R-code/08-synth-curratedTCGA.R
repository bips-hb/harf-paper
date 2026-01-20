# This file contains code to synthesize the Hemberger et al. datasets.
# We synthesize 10 datasets based on the preprocessed Hemberger et al. data.
library(this.path)
source(file.path(this.dir(), "00-library-and-setup.R"))
source(file.path(this.dir(), "synthetizer_tcga_fct.R"))

# Synthesize TCGA datasets
# ====================================

# HARF synthesizer chunk sizes for TCGA datasets
# **********************************************

harf_chunck_tcga <- c(
  brca = 10,
  luad = 10,
  lusc = 10,
  kirc = 10,
  coad = 10
)
set.seed(1245)
syn_harf_tcga_list <- lapply(names(orig_tcga_data_files), function (nam) {
    message("Synthesizing TCGA dataset ", nam, " with HARF...\n")
    harf_tcga_synthesizer(
      dt_name = nam,
      org_file_paths = orig_tcga_data_files,
      syn_file_prfxs = syn_tcga_dt_prefixes,
      chunck_size = harf_chunck_tcga[nam],
      parallel = TRUE,
      verbose = TRUE,
      N_synth = 10
    )
})
# Save HARF synthesis results
saveRDS(
  syn_harf_tcga_list,
  file = file.path(
    syn_tcga_data_dir,
    "synth_tcga_harf_results.rds"
  )
)
