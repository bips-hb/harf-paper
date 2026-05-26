library(batchtools)

source(file.path(r_code_dir, "00-library-and-setup.R"))
source(file.path(r_code_dir, "13-AD/01-ad-problem.R"))
source(file.path(r_code_dir, "13-AD/02-ad-algorithm.R"))

# 1.  Prepare registry for Hemberger et al. datasets synthesis with HARF
template <- "~/batchtools/batchtools.slurm.tmpl"
makeClusterFunctionsSlurm(template = template)
partition <- "day-long-cpu"

unlink(file.path(reg_dir, "ad"), recursive = TRUE)
reg <- makeExperimentRegistry(
  file.dir = file.path(reg_dir, "ad"),
  conf.file = "~/batchtools/batchtools.conf.R",
  packages = character(0L),
  work.dir = "/home/ckuetef/projects/harf-paper/R-code/13-AD",
  source = c(
    file.path(r_code_dir, "00-library-and-setup.R"),
    file.path(perf_dir, "utils.R"),
    file.path(perf_dir, "evaluation_functions.R"),
    file.path(perf_dir, "cluster-and-eval.R")
  ),
  seed = 123
)

# 2. Add problems and algorithms to registry
addProblem(name = "ad", 
           data = list(file_name = metab_file), 
           fun = create_ad_data,
           reg = reg)
# 3. Add algorithms to registry
addAlgorithm(name = "harf_ad_pred", fun = harf_ad_pred, reg = reg)
addAlgorithm(name = "arf_ad_pred", fun = arf_ad_pred, reg = reg)

# 4. Parameter design
pdes <- expand.grid(evidence = c(FALSE, TRUE), prop_synth = c(1.1, 1.2, 1.3, 1.4, 1.5))
pdes <- list(ad = pdes)
# 5. Algorithm design
ades <- expand.grid(
  num_trees = 10,
  chunck_size = seq(5, 50, by = 5),
  num_btwn_pcs = c(2)
)
ades <- list(harf_ad_pred = ades, arf_ad_pred = data.frame(num_trees = 10))
# 6. Add experiments to registry
addExperiments(reg = reg,
               prob.designs = pdes,
               algo.designs = ades,
               repls = 20)
summarizeExperiments()

# 7. Test h-ARF before submitting to cluster
id1 = head(findExperiments(prob.name = "ad", algo.name = "harf_ad_pred"), 1)[1, ]
testJob(id = id1, reg = reg)
# 8. Test ARF before submitting to cluster
id2 = head(findExperiments(prob.name = "ad", algo.name = "arf_ad_pred"), 1)[1, ]
testJob(id = id2, reg = reg)

ids <- findExperiments(reg = reg)
ids <- ids[, chunk := chunk(job.id, chunk.size = 500)]
submitJobs(reg = reg, ids = ids,
           resources = list(walltime = "4:50:00",
                            memory = 1024 * 4,
                            ncpus = 1,
                            partition = "short-cpu",
                            ntasks = 1,
                            name = "ad_prediction"))
waitForJobs(reg = reg)
