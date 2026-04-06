library(batchtools)

source(file.path(r_code_dir, "00-library-and-setup.R"))
source(file.path(r_code_dir, "03-hember-problem.R"))
source(file.path(r_code_dir, "04-hember-algorithm.R"))

# 1.  Prepare registry for Hemberger et al. datasets synthesis with HARF
makeClusterFunctionsSlurm(template = "~/batchtools/batchtools.slurm.tmpl")
partition <- "day-long-cpu"

unlink(file.path(reg_dir, "hember"), recursive = TRUE)
reg <- makeExperimentRegistry(
  file.dir = file.path(reg_dir, "hember"),
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
addAlgorithm(name = "harf_synthesizer", fun = harf_synthesizer, reg = reg)
addAlgorithm(name = "arf_synthesizer", fun = arf_synthesizer, reg = reg)

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
ids <- ids[, chunk := chunk(job.id, chunk.size = 500)]
submitJobs(reg = reg, ids = ids,
           resources = list(walltime = "10:50:00",
                            memory = 1024 * 15,
                            ncpus = 1,
                            partition = partition,
                            ntasks = 1,
                            name = "harf_hember"))
waitForJobs(reg = reg)

# ============================
# Resubmission for failed jobs
# ============================
reg <- loadRegistry(file.dir = file.path(reg_dir, "hember"),
                    writeable = TRUE,
                    conf.file = "~/batchtools/batchtools.conf.R")
ids_expired <- findExpired(reg = reg)
ids_error <- findErrors(reg = reg)
ids_failed <- rbind(ids_expired, ids_error)
if (length(ids_failed) > 0) {
  message("Resubmitting failed jobs...")
  submitJobs(reg = reg, ids = ids_failed, 
             resources = list(walltime = "10:59:00",
                              memory = 1024 * 10,
                              ncpus = 1,
                              partition = partition))
}

# 8. Collect results
ids <- findExperiments(algo.name = "harf_synthesizer", reg = reg)
ids_done <- findDone(ids, reg = reg)
# Load and flat job paramters
job_pars <- getJobPars(ids_done, reg = reg)
job_pars_algo <- rbindlist(job_pars$algo.pars, idcol = "job.id", fill = TRUE)
job_pars_prob <- rbindlist(job_pars$prob.pars, idcol = "job.id", fill = TRUE)
job_pars_DT <- merge(job_pars_algo, job_pars_prob, by = "job.id")
# Load and flat results for HARF
res_harf <- reduceResultsList(ids = ids_done, reg = reg)
res_harf_DT <- rbindlist(res_harf, idcol = "job.id", fill = TRUE)
res_harf_DT$algorithm <- "HARF"
res_harf_DT <- merge(res_harf_DT, job_pars_DT, by = "job.id")
saveRDS(res_harf_DT, file.path(res_hember_dt_dir, "harf_hember_results.rds"))

# Load and flat results for ARF
ids <- findExperiments(algo.name = "arf_synthesizer", reg = reg)
ids_done <- findDone(ids, reg = reg)
res_arf <- reduceResultsList(ids = ids_done, reg = reg)
res_arf_DT <- rbindlist(res_arf, idcol = "job.id")
res_arf_DT$algorithm <- "ARF"
res_arf_DT <- merge(res_arf_DT, job_pars_DT, by = "job.id")
res_arf_DT$chunck_size <- 0
saveRDS(res_arf_DT, file.path(res_hember_dt_dir, "arf_hember_results.rds"))

# Rbind HARF and ARF results
res_all_DT <- rbind(res_harf_DT, res_arf_DT)
saveRDS(res_all_DT, file.path(res_hember_dt_dir,
                              "harf_arf_hember_results.rds"))



