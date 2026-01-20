# This file contains code to synthesize the Hemberger et al. datasets.
# We synthesize 10 datasets based on the preprocessed Hemberger et al. data.
library(this.path)
source(file.path(this.dir(), "00-library-and-setup.R"))
source(file.path(this.dir(), "synthetizer_fct.R"))


# Synthesize Hemberger et al. datasets
# ====================================

# HARF synthesizer chunk sizes for Hemberger et al. datasets
# **********************************************************

harf_chunck_hember <- c(
  lake = 10,
  manno = 10,
  li = 10,
  patel = 10,
  deng = 10
)
set.seed(1245)
syn_harf_hember_list <- lapply(names(org_hember_dt_files), function (nam) {
  lapply(1:5, function (i) {
    message("Synthesizing Hemberger dataset ", nam,
            ", iteration ", i, " with HARF...\n")
    harf_synthesizer(
      dt_name = nam,
      org_file_paths = org_hember_dt_files,
      syn_file_prfxs = syn_hember_dt_prefixes,
      chunck_size = harf_chunck_hember[nam],
      parallel = TRUE,
      verbose = TRUE,
      i = i
    )
  })
})
# Save HARF synthesis results
saveRDS(
  syn_harf_hember_list,
  file = file.path(
    syn_hember_dt_dir,
    "synth_hember_harf_results.rds"
  )
)
# Synthesize Hemberger et al. datasets with ARF
# *********************************************
set.seed(5421)
syn_arf_hember_list <- lapply(names(org_hember_dt_files), function (nam) {
  lapply(1:5, function (i) {
    message("Synthesizing Hemberger dataset ", nam,
            ", iteration ", i, " with ARF...\n")
    arf_synthesizer(
      dt_name = nam,
      org_file_paths = org_hember_dt_files,
      syn_file_prfxs = syn_hember_dt_prefixes,
      parallel = TRUE,
      verbose = FALSE,
      i = i
    )
  })
})
