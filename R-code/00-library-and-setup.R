# Libraries
# =================================
if (FALSE) {
  # Use this to install required packages
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
library("this.path")
library("R.utils")
library("data.table")
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
org_data_dir <- file.path(dirname(r_code_dir), "data/original")
org_hember_dt_dir <- file.path(org_data_dir, "hember-lab")
org_baron_dt_dir <- file.path(org_hember_dt_dir, "brain-baron")
org_lake_dt_dir <- file.path(org_hember_dt_dir, "brain-lake")
org_manno_dt_dir <- file.path(org_hember_dt_dir, "brain-manno")
org_li_dt_dir <- file.path(org_hember_dt_dir, "tissue-li")
org_patel_dt_dir <- file.path(org_hember_dt_dir, "tissue-patel")
org_simlr_data_dir <- file.path(org_data_dir, "SIMLR")
org_tcga_tgex_data_dir <- file.path(org_data_dir, "TCGA+GTEX")
# Synthetic data directory
syn_data_dir <- file.path(dirname(r_code_dir), "data/synthetic")
syn_hember_dt_dir <- file.path(syn_data_dir, "hember-lab")
syn_simlr_data_dir <- file.path(syn_data_dir, "SIMLR")
syn_tcga_tgex_data_dir <- file.path(syn_data_dir, "TCGA+GTEX")
dir.create(syn_data_dir, showWarnings = FALSE)
dir.create(syn_hember_dt_dir, showWarnings = FALSE)
dir.create(syn_simlr_data_dir, showWarnings = FALSE)
dir.create(syn_tcga_tgex_data_dir, showWarnings = FALSE)
# Preprocessed data
org_baron_dt_file <- file.path(org_baron_dt_dir, "processed_brain_baron_data.csv")
org_lake_dt_file <- file.path(org_lake_dt_dir, "processed_brain_lake_data.csv")
org_manno_dt_file <- file.path(org_manno_dt_dir, "processed_brain_manno_data.csv")
org_li_dt_file <- file.path(org_li_dt_dir, "processed_tissue_li_data.csv")
org_patel_dt_file <- file.path(org_patel_dt_dir, "processed_tissue_patel_data.csv")
org_hember_dt_files <- c(
  baron = org_baron_dt_file,
  lake = org_lake_dt_file,
  manno = org_manno_dt_file,
  li = org_li_dt_file,
  patel = org_patel_dt_file
)
org_simlr_dt_file <- dir(org_simlr_data_dir, pattern = "*.RData",
                         full.names = TRUE)
org_tcga_tgex_dt_file <- dir(org_tcga_tgex_data_dir,
                              pattern = "*.rds", full.names = TRUE)
