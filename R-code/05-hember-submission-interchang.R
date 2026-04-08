library(batchtools)

source(file.path(r_code_dir, "00-library-and-setup.R"))
source(file.path(r_code_dir, "03-hember-problem.R"))
source(file.path(r_code_dir, "04-hember-algorithm.R"))

# 1.  Prepare registry for Hemberger et al. datasets synthesis with HARF
makeClusterFunctionsSlurm(template = "~/batchtools/batchtools.slurm.tmpl")
partition <- "day-long-cpu"

unlink(file.path(reg_dir, "hember-interchange"), recursive = TRUE)
reg <- makeExperimentRegistry(
  file.dir = file.path(reg_dir, "hember-interchange"),
  conf.file = "~/batchtools/batchtools.conf.R",
  packages = character(0L),
  work.dir = "/home/ckuetef/projects/harf-paper/R-code",
  source = c(
    file.path(r_code_dir, "00-library-and-setup.R"),
    file.path(perf_dir, "utils.R"),
    file.path(perf_dir, "evaluation_functions.R"),
    file.path(perf_dir, "cluster-and-eval.R")
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
addAlgorithm(name = "harf_synthesizer", fun = harf_interchangeable, reg = reg)
addAlgorithm(name = "arf_synthesizer", fun = arf_interchangeable, reg = reg)

# 4. Parameter design
pdes <- data.frame(evidence = c(FALSE, TRUE))
pdes <- list(lake = pdes, manno = pdes, li = pdes, patel = pdes)
# 5. Algorithm design
ades <- expand.grid(
  num_trees = 10,
  chunck_size = c(25, 50, 75, 100),
  num_btwn_pcs = c(2, 3)
)
ades <- list(harf_synthesizer = ades, arf_synthesizer = data.frame(num_trees = 10))
# 6. Add experiments to registry
addExperiments(reg = reg,
               prob.designs = pdes,
               algo.designs = ades,
               repls = 100)
summarizeExperiments()

# 7. Test before submitting to cluster
id1 = head(findExperiments(prob.name = "lake", algo.name = "harf_synthesizer"), 175)[1, ]
# testJob(id = id1, reg = reg)

# id2 = head(findExperiments(prob.name = "lake", algo.name = "arf_synthesizer"), 1)
# testJob(id = id2, reg = reg)

ids <- findExperiments()
ids <- ids[, chunk := chunk(job.id, chunk.size = 400)]
submitJobs(reg = reg, ids = ids,
           resources = list(walltime = "10:50:00",
                            memory = 1024 * 3,
                            ncpus = 1,
                            partition = partition,
                            ntasks = 1,
                            name = "harf_hember_inter"))
waitForJobs(reg = reg)
