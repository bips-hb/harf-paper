# This file contains code to synthesize the SIMLR datasets.
# We synthesize 10 datasets based on the preprocessed SIMLR data.
library(this.path)
source(file.path(this.dir(), "00-library-and-setup.R"))
source(file.path(this.dir(), "synthetizer_fct.R"))


# Synthesize SIMLR datasets
# ====================================

# HARF synthesizer chunk sizes for SIMLR datasets
# **********************************************************

harf_chunck_simlr <- c(
  mecs = 10,
  kolod = 10,
  pollen = 10,
  usoskin = 10,
  zelsel = 10
)

syn_harf_simlr_list <- lapply(names(org_simlr_dt_proc_files), function (nam) {
  lapply(1:5, function (i) {
    message("Synthesizing Hemberger dataset ", nam,
            ", iteration ", i, " with HARF...\n")
    harf_synthesizer(
      org_file_path = org_hember_dt_files[nam],
      syn_file_prfx = syn_hember_dt_prefixes[nam],
      i = i,
      chunck_size = harf_chunck_hember[nam],
      parallel = FALSE,
      verbose = FALSE
    )
  })
})

# Synthesize SIMLR datasets with ARF
# *********************************************
syn_arf_simlr_list <- lapply(names(org_simlr_dt_proc_files), function (nam) {
  lapply(1:5, function (i) {
    message("Synthesizing SIMLR dataset ", nam,
            ", iteration ", i, " with ARF...\n")
    arf_synthesizer(
      org_file_path = org_hember_dt_files[nam],
      syn_file_prfx = syn_hember_dt_prefixes[nam],
      parallel = FALSE,
      verbose = FALSE
    )
  })
})
