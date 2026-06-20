source("../00-library-and-setup.R")


# Set the path to original AD data
ad_dir <- "/huels_lab/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/data/original/AD"
metab_file <- file.path(ad_dir, "metabolomics_final.txt")
ad_results_dir <- file.path(res_dir, "ad")
dir.create(ad_results_dir, showWarnings = FALSE)
