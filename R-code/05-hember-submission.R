library(batchtools)

source(file.path(r_code_dir, "00-library-and-setup.R"))
source(file.path(r_code_dir, "03-hember-problem.R"))
source(file.path(r_code_dir, "04-hember-algorithm.R"))

# 1.  Prepare registry for Hemberger et al. datasets synthesis with HARF
makeClusterFunctionsSlurm(template = "~/batchtools/batchtools.slurm.tmpl")
template <- "~/batchtools/batchtools.slurm.tmpl"
partition <- "week-long-cpu"

unlink(file.path(reg_dir, "hember"), recursive = TRUE)
reg <- makeExperimentRegistry(
  file.dir = file.path(reg_dir, "hember"),
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
           reg = reg)
addProblem(name = "manno", 
           data = list(file_name = org_hember_dt_files["manno"], data_names = "manno"), 
           fun = create_single_cell_data,
           reg = reg)
addProblem(name = "li",
           data = list(file_name = org_hember_dt_files["li"], data_names = "li"), 
           fun = create_single_cell_data,
           reg = reg)
addProblem(name = "patel",
           data = list(file_name = org_hember_dt_files["patel"], data_names = "patel"),
           fun = create_single_cell_data,
           reg = reg)

# 3. Add algorithms to registry
addAlgorithm(name = "harf_synthesizer", fun = harf_synthesizer, reg = reg)
addAlgorithm(name = "arf_synthesizer", fun = arf_synthesizer, reg = reg)

# 4. Parameter design
pdes <- data.frame(evidence = c(FALSE, TRUE))
pdes <- list(lake = pdes, manno = pdes, li = pdes, patel = pdes)
# 5. Algorithm design
ades <- expand.grid(
  num_trees = 10,
  chunck_size = c(100),
  num_btwn_pcs = c(2, 3, 4, 5)
)
ades <- list(harf_synthesizer = ades, arf_synthesizer = ades)
# 6. Add experiments to registry
addExperiments(reg = reg,
               prob.designs = pdes,
               algo.designs = ades,
               repls = 10)
summarizeExperiments()

# 7. Test before submitting to cluster
id1 = head(findExperiments(prob.name = "lake", algo.name = "harf_synthesizer"), 1)
testJob(id = id1, reg = reg)

id2 = head(findExperiments(prob.name = "lake", algo.name = "arf_synthesizer"), 1)
# testJob(id = id2, reg = reg)

ids <- findExperiments()
submitJobs(reg = reg, 
           ids = ids[ , chunk := chunk(job.id, chunk.size = 800)],
           resources = list(walltime = "12:00:00",
                                       memory = "5g",
                                        ncpus = 5,
                                       partition = partition))
waitForJobs(reg = reg)
