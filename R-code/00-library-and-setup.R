# Libraries
# =================================
if (FALSE) {
  # Use this to install required packages
  devtools::install_github("bips-hp/harf")
  install.packages("devtools")
  install.packages("this.path")
  install.packages("R.utils")
  install.packages("data.table")
  install.packages("rsvd")
  install.packages("Rtsne")
  install.packages("cowplot")
  if (!require("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
  BiocManager::install("SingleCellExperiment")
  BiocManager::install("scater")
  install.packages("ggplot2")
  install.packages("corrplot")
  install.packages("doParallel")
}
library("harf")
library("this.path")
library("R.utils")
library("data.table")
library("SingleCellExperiment")
library("scater")
library("rsvd")
library("Rtsne")
library("cowplot")
library("ggplot2")
library("corrplot")
library("doParallel")

# Setups
# =================================
r_code_dir <- dirname(this.path())

# Original data directory
# ================================
org_data_dir <- file.path(dirname(r_code_dir), "data/original")
org_hember_dt_dir <- file.path(org_data_dir, "hember-lab")
# org_baron_dt_dir <- file.path(org_hember_dt_dir, "brain-baron")
# org_fan_dt_dir <- file.path(org_hember_dt_dir, "embryo-fan-mouse")
org_lake_dt_dir <- file.path(org_hember_dt_dir, "brain-lake")
org_manno_dt_dir <- file.path(org_hember_dt_dir, "brain-manno")
org_li_dt_dir <- file.path(org_hember_dt_dir, "tissue-li")
org_patel_dt_dir <- file.path(org_hember_dt_dir, "tissue-patel")
org_deng_dt_dir <- file.path(org_hember_dt_dir, "embryo-deng-mouse")
org_simlr_data_dir <- file.path(org_data_dir, "SIMLR")
org_tcga_tgex_data_dir <- file.path(org_data_dir, "TCGA+GTEX")

# Preprocessed data
# org_baron_dt_file <- file.path(org_baron_dt_dir, "processed_brain_baron_data.csv")
# org_fan_dt_file <- file.path(org_fan_dt_dir, "E-MTAB-3321.processed.1.zip")
org_lake_dt_file <- file.path(org_lake_dt_dir, "processed_brain_lake_data.csv")
org_manno_dt_file <- file.path(org_manno_dt_dir, "processed_brain_manno_data.csv")
org_li_dt_file <- file.path(org_li_dt_dir, "processed_tissue_li_data.csv")
org_patel_dt_file <- file.path(org_patel_dt_dir, "processed_tissue_patel_data.csv")
org_deng_dt_file <- file.path(org_deng_dt_dir, "embryo_deng_mouse_processed_data.csv")

org_hember_dt_files <- c(
  # baron = org_baron_dt_file,
  # fan = org_fan_dt_file,
  lake = org_lake_dt_file,
  manno = org_manno_dt_file,
  li = org_li_dt_file,
  patel = org_patel_dt_file,
  deng = org_deng_dt_file
)
org_simlr_dt_file <- dir(org_simlr_data_dir, pattern = "*.RData",
                         full.names = TRUE)
org_tcga_tgex_dt_file <- dir(org_tcga_tgex_data_dir,
                              pattern = "*.rds", full.names = TRUE)


# Synthetic data directory
# ================================
syn_data_dir <- file.path(dirname(r_code_dir), "data/synthetic")
syn_hember_dt_dir <- file.path(syn_data_dir, "hember-lab")
syn_simlr_data_dir <- file.path(syn_data_dir, "SIMLR")
syn_tcga_tgex_data_dir <- file.path(syn_data_dir, "TCGA+GTEX")
dir.create(syn_data_dir, showWarnings = FALSE)
dir.create(syn_hember_dt_dir, showWarnings = FALSE)
dir.create(syn_simlr_data_dir, showWarnings = FALSE)
dir.create(syn_tcga_tgex_data_dir, showWarnings = FALSE)
# syn_baron_dt_dir <- file.path(syn_hember_dt_dir, "brain-baron")
# syn_fan_dt_dir <- file.path(syn_hember_dt_dir, "embryo-fan-mouse")
syn_lake_dt_dir <- file.path(syn_hember_dt_dir, "brain-lake")
syn_manno_dt_dir <- file.path(syn_hember_dt_dir, "brain-manno")
syn_li_dt_dir <- file.path(syn_hember_dt_dir, "tissue-li")
syn_patel_dt_dir <- file.path(syn_hember_dt_dir, "tissue-patel")
syn_deng_dt_dir <- file.path(syn_hember_dt_dir, "embryo-deng-mouse")
# dir.create(syn_baron_dt_dir, showWarnings = FALSE)
# dir.create(syn_fan_dt_dir, showWarnings = FALSE)
dir.create(syn_lake_dt_dir, showWarnings = FALSE)
dir.create(syn_manno_dt_dir, showWarnings = FALSE)
dir.create(syn_li_dt_dir, showWarnings = FALSE)
dir.create(syn_patel_dt_dir, showWarnings = FALSE)
dir.create(syn_deng_dt_dir, showWarnings = FALSE)
# Prepare file prefixes for synthetic datasets
# syn_baron_dt_prefix <- file.path(syn_baron_dt_dir, "synth_baron")
# syn_fan_dt_prefix <- file.path(syn_fan_dt_dir, "synth_fan")
syn_lake_dt_prefix <- file.path(syn_lake_dt_dir, "synth_lake")
syn_manno_dt_prefix <- file.path(syn_manno_dt_dir, "synth_manno")
syn_li_dt_prefix <- file.path(syn_li_dt_dir, "synth_li")
syn_patel_dt_prefix <- file.path(syn_patel_dt_dir, "synth_patel")
syn_deng_dt_prefix <- file.path(syn_deng_dt_dir, "synth_deng")
syn_hember_dt_prefixes <- c(
  # baron = syn_baron_dt_prefix,
  # fan = syn_fan_dt_prefix,
  lake = syn_lake_dt_prefix,
  manno = syn_manno_dt_prefix,
  li = syn_li_dt_prefix,
  patel = syn_patel_dt_prefix,
  deng = syn_deng_dt_prefix
)
syn_simlr_dt_prefix <- file.path(syn_simlr_data_dir, "synth_simlr")
syn_tcga_tgex_dt_prefix <- file.path(syn_tcga_tgex_data_dir, "synth_tcga_tgex")

