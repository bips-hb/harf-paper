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
  BiocManager::install(c("curatedTCGAData", "SummarizedExperiment"))
  install.packages("ggplot2")
  install.packages("corrplot")
  install.packages("doParallel")
}
library("harf")
library("arf")
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
library("curatedTCGAData")
library("SummarizedExperiment")
# Parallel backend
library(doParallel)
registerDoParallel(cores = 2)

# Register cores - Windows
# cl <- makeCluster(2)
# registerDoParallel(cl)

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

# Preprocessed Hember et al. data
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

# Raw SIMLR data file
org_simlr_mecs_file <- file.path(org_simlr_data_dir, "Test_1_mECS.RData")
org_simlr_kolod_file <- file.path(org_simlr_data_dir, "Test_2_Kolod.RData")
org_simlr_pollen_file <- file.path(org_simlr_data_dir, "Test_3_Pollen.RData")
org_simlr_usoskin_file <- file.path(org_simlr_data_dir, "Test_4_Usoskin.RData")
org_simlr_zeisel_file <- file.path(org_simlr_data_dir, "Test_5_Zeisel.RData")
org_simlr_dt_files <- c(
  mecs = org_simlr_mecs_file,
  kolod = org_simlr_kolod_file,
  pollen = org_simlr_pollen_file,
  usoskin = org_simlr_usoskin_file,
  zelsel = org_simlr_zeisel_file
)

# Processed SIMLR file
org_simlr_mecs_proc_file <- file.path(org_simlr_data_dir, "mecs_data.csv")
org_simlr_kolog_proc_file <- file.path(org_simlr_data_dir, "kolod_data.csv")
org_simlr_pollen_proc_file <- file.path(org_simlr_data_dir, "pollen_data.csv")
org_simlr_usoskin_proc_file <- file.path(org_simlr_data_dir, "usoskin_data.csv")
org_simlr_zeisel_proc_file <- file.path(org_simlr_data_dir, "zeisel_data.csv")

org_simlr_dt_proc_files <- c(
  mecs = org_simlr_mecs_proc_file,
  kolod = org_simlr_kolog_proc_file,
  pollen = org_simlr_pollen_proc_file,
  usoskin = org_simlr_usoskin_proc_file,
  zelsel = org_simlr_zeisel_proc_file
)

org_tcga_tgex_dt_file <- dir(org_tcga_tgex_data_dir,
                              pattern = "*.rds", full.names = TRUE)


# Synthetic Hember et al. data directory
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
dir.create(syn_lake_dt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(syn_manno_dt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(syn_li_dt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(syn_patel_dt_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(syn_deng_dt_dir, recursive = TRUE, showWarnings = FALSE)
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

# Synthetic SIMLR data directory
# ================================
syn_simlr_dt_prefix <- file.path(syn_simlr_data_dir, "synth_simlr")
syn_mecs_dt_prefix <- file.path(syn_simlr_data_dir, "synth_mecs")
syn_kolod_dt_prefix <- file.path(syn_simlr_data_dir, "synth_kolod")
syn_pollen_dt_prefix <- file.path(syn_simlr_data_dir, "synth_pollen")
syn_usoskin_dt_prefix <- file.path(syn_simlr_data_dir, "synth_usoskin")
syn_zeisel_dt_prefix <- file.path(syn_simlr_data_dir, "synth_zeisel")

syn_tcga_tgex_dt_prefix <- file.path(syn_tcga_tgex_data_dir, "synth_tcga_tgex")

# Prepare path to original currated TCGA data
orig_tcga_data_dir <- file.path(org_tcga_tgex_data_dir, "curatedTCGAData")
dir.create(orig_tcga_data_dir, showWarnings = FALSE)
orig_brca_data_dir <- file.path(orig_tcga_data_dir, "BRCA")
dir.create(orig_brca_data_dir, showWarnings = FALSE)
orig_luad_data_dir <- file.path(orig_tcga_data_dir, "LUAD")
dir.create(orig_luad_data_dir, showWarnings = FALSE)
orig_lusc_data_dir <- file.path(orig_tcga_data_dir, "LUSC")
dir.create(orig_lusc_data_dir, showWarnings = FALSE)
orig_kirc_data_dir <- file.path(orig_tcga_data_dir, "KIRC")
dir.create(orig_kirc_data_dir, showWarnings = FALSE)
orig_coad_data_dir <- file.path(orig_tcga_data_dir, "COAD")
dir.create(orig_coad_data_dir, showWarnings = FALSE)
orig_brca_data_file <- file.path(orig_brca_data_dir, "brca.txt")
orig_luad_data_file <- file.path(orig_luad_data_dir, "luad.txt")
orig_lusc_data_file <- file.path(orig_lusc_data_dir, "lusc.txt")
orig_kirc_data_file <- file.path(orig_kirc_data_dir, "kirc.txt")
orig_coad_data_file <- file.path(orig_coad_data_dir, "coad.txt")

orig_tcga_data_files <- c(
  brca = orig_brca_data_file,
  luad = orig_luad_data_file,
  lusc = orig_lusc_data_file,
  kirc = orig_kirc_data_file,
  coad = orig_coad_data_file
)

# Prepare path to synthetic currated TCGA data
syn_tcga_data_dir <- file.path(syn_tcga_tgex_data_dir, "curatedTCGAData")
dir.create(syn_tcga_data_dir, showWarnings = FALSE)
syn_brca_data_dir <- file.path(syn_tcga_data_dir, "BRCA")
dir.create(syn_brca_data_dir, showWarnings = FALSE)
syn_luad_data_dir <- file.path(syn_tcga_data_dir, "LUAD")
dir.create(syn_luad_data_dir, showWarnings = FALSE)
syn_lusc_data_dir <- file.path(syn_tcga_data_dir, "LUSC")
dir.create(syn_lusc_data_dir, showWarnings = FALSE)
syn_kirc_data_dir <- file.path(syn_tcga_data_dir, "KIRC")
dir.create(syn_kirc_data_dir, showWarnings = FALSE)
syn_coad_data_dir <- file.path(syn_tcga_data_dir, "COAD")
dir.create(syn_coad_data_dir, showWarnings = FALSE)
syn_brca_dt_prefix <- file.path(syn_brca_data_dir, "synth_brca")
syn_luad_dt_prefix <- file.path(syn_luad_data_dir, "synth_luad")
syn_lusc_dt_prefix <- file.path(syn_lusc_data_dir, "synth_lusc")
syn_kirc_dt_prefix <- file.path(syn_kirc_data_dir, "synth_kirc")
syn_coad_dt_prefix <- file.path(syn_coad_data_dir, "synth_coad")
syn_tcga_dt_prefixes <- c(
  brca = syn_brca_dt_prefix,
  luad = syn_luad_dt_prefix,
  lusc = syn_lusc_dt_prefix,
  kirc = syn_kirc_dt_prefix,
  coad = syn_coad_dt_prefix
)
