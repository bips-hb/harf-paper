library(batchtools)

source(file.path(r_code_dir, "00-library-and-setup.R"))
source(file.path(r_code_dir, "09-tcga-problem.R"))
source(file.path(r_code_dir, "10-tcga-algorithm.R"))

# 1.  Prepare registry for Hemberger et al. datasets synthesis with HARF
template <- "~/batchtools/batchtools.slurm.tmpl"
makeClusterFunctionsSlurm(template = template)
partition <- "day-long-cpu"

# unlink(file.path(reg_dir, "hember-interchange2"), recursive = TRUE)
reg <- makeExperimentRegistry(
  file.dir = file.path(reg_dir, "tcga-downstreem-prediction"),
  conf.file = "~/batchtools/batchtools.conf.R",
  packages = c("glmnet", "pROC", "data.table"),
  work.dir = "/home/ckuetef/projects/harf-paper/R-code",
  source = c(
    file.path(r_code_dir, "00-library-and-setup.R"),
    file.path(perf_dir, "utils.R")
  ),
  seed = 123
)

# 2. Add problems and algorithms to registry
addProblem(name = "luad", 
           data = list(file_name = orig_tcga_data_files["luad"], data_names = "luad"), 
           fun = create_tcga_data,
           reg = reg)
addProblem(name = "lusc", 
           data = list(file_name = orig_tcga_data_files["lusc"], data_names = "lusc"), 
           fun = create_tcga_data,
           reg = reg)
addProblem(name = "kirc",
           data = list(file_name = orig_tcga_data_files["kirc"], data_names = "kirc"), 
           fun = create_tcga_data,
           reg = reg)
addProblem(name = "coad",
           data = list(file_name = orig_tcga_data_files["coad"], data_names = "coad"),
           fun = create_tcga_data,
           reg = reg)

# 3. Add algorithms to registry
addAlgorithm(name = "harf_ds_pred", fun = harf_ds_pred, reg = reg)
addAlgorithm(name = "arf_ds_pred", fun = arf_ds_pred, reg = reg)

# 4. Parameter design
pdes <- data.frame(evidence = c(FALSE, TRUE))
pdes <- list(luad = pdes, lusc = pdes, kirc = pdes, coad = pdes)
# 5. Algorithm design
ades <- expand.grid(
  num_trees = 10,
  chunck_size = seq(5, 50, by = 5),
  num_btwn_pcs = c(2)
)
ades <- list(harf_ds_pred = ades, arf_ds_pred = data.frame(num_trees = 10))
# 6. Add experiments to registry
addExperiments(reg = reg,
               prob.designs = pdes,
               algo.designs = ades,
               repls = 100)
summarizeExperiments()

# 7. Test before submitting to cluster
id1 = head(findExperiments(prob.name = "luad", algo.name = "harf_ds_pred"), 175)[1, ]
testJob(id = id1, reg = reg)

id2 = head(findExperiments(prob.name = "lusc", algo.name = "arf_ds_pred"), 1)
# testJob(id = id2, reg = reg)

ids <- findExperiments(reg = reg)
ids <- ids[, chunk := chunk(job.id, chunk.size = 1000)]
submitJobs(reg = reg, ids = ids,
           resources = list(walltime = "4:50:00",
                            memory = 1024 * 4,
                            ncpus = 1,
                            partition = "short-cpu",
                            ntasks = 1,
                            name = "tcga-downstreem-prediction"))
waitForJobs(reg = reg)
