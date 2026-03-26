library(batchtools)

source(file.path(r_code_dir, "00-library-and-setup.R"))

# 1.  Prepare registry for Hemberger et al. datasets synthesis with HARF
makeClusterFunctionsSlurm(template = "~/batchtools/batchtools.slurm.tmpl")
template <- "~/batchtools/batchtools.slurm.tmpl"
partition <- "week-long-cpu"

reg <- makeExperimentRegistry(
  file.dir = "/huels_lab/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/registry/harf-paper",
  conf.file = "~/batchtools/batchtools.conf.R",
  packages = character(0L),
  source = c(
    file.path(r_code_dir, "00-library-and-setup.R"),
    file.path(perf_dir, "utils.R"),
    file.path(perf_dir, "evaluation_functions.R")
  ),
  seed = 123
)

# 2. Add problems and algorithms to registry
addProblem(name = "lake", 
           data = list(file_name = org_hember_dt_files["lake"], data_names = "lake"), 
           fun = create_single_cell_data,
           registry = reg)
addProblem(name = "manno", 
           data = list(file_name = org_hember_dt_files["manno"], data_names = "manno"), 
           fun = create_single_cell_data,
           registry = reg)
addProblem(name = "li",
           data = list(file_name = org_hember_dt_files["li"], data_names = "li"), 
           fun = create_single_cell_data,
           registry = reg)
addProblem(name = "patel",
           data = list(file_name = org_hember_dt_files["patel"], data_names = "patel"),
           fun = create_single_cell_data,
           registry = reg)

# 3. Add algorithms to registry
addAlgorithm(name = "harf_synthesizer", fun = harf_synthesizer, registry = reg)

# 4. Parameter design
hember_pdes <- data.frame(evidence = c(FALSE, TRUE))

# 5. Algorithm design
hember_ades <- expand.grid(
  num_trees = 10,
  chunck_size = c(5, 10, 15, 20, 25),
  num_btwn_pcs = c(2, 3, 4, 5)
)

# 6. Add experiments to registry
addExperiments(reg = reg,
               prob.designs = hember_pdes,
               algo.designs = hember_ades, repls = 10)
summarizeExperiments()

# 7. Test before submitting to cluster
id1 = head(findExperiments(prob.name = "lake", algo.name = "harf_synthesizer"), 1)
