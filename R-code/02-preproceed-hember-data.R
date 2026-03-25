# This script preprocesses the Hemberger et al. datasets such that they can be
# passed to h_arf as data.frames.
source(file.path(r_code_dir, "01-prepare-hember-data/create_sce_hember.R"))
# source(file.path(r_code_dir, "01-prepare-hember-data/01-baron.R"),
#        chdir = TRUE)
# source(file.path(r_code_dir, "01-prepare-hember-data/02-baron-processing.R"),
#        chdir = TRUE)
source(file.path(r_code_dir, "01-prepare-hember-data/01-lake.R"),
       chdir = TRUE)
source(file.path(r_code_dir, "01-prepare-hember-data/02-lake-processing.R"),
       chdir = TRUE)
source(file.path(r_code_dir, "01-prepare-hember-data/01-manno.R"),
       chdir = TRUE)
source(file.path(r_code_dir, "01-prepare-hember-data/02-manno-processing.R"),
       chdir = TRUE)
source(file.path(r_code_dir, "01-prepare-hember-data/01-li.R"),
       chdir = TRUE)
source(file.path(r_code_dir, "01-prepare-hember-data/02-li-processing.R"),
       chdir = TRUE)
source(file.path(r_code_dir, "01-prepare-hember-data/01-patel.R"),
       chdir = TRUE)
source(file.path(r_code_dir, "01-prepare-hember-data/02-patel-processing.R"),
       chdir = TRUE)
# source(file.path(r_code_dir, "01-prepare-hember-data/01-deng.R"),
#        chdir = TRUE)
# source(file.path(r_code_dir, "01-prepare-hember-data/02-deng-processing.R"),
#        chdir = TRUE)
