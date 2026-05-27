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
  chunck_size = seq(50, 50, by = 5),
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

# ============================
# Resubmission for failed jobs
# ============================
reg <- loadRegistry(file.dir = file.path(reg_dir, "ad"),
                    writeable = TRUE,
                    conf.file = "~/batchtools/batchtools.conf.R")
ids_expired <- findExpired(reg = reg)
ids_error <- findErrors(reg = reg)
ids_failed <- findNotDone(reg = reg)
if (length(ids_failed) > 0) {
  message("Resubmitting failed jobs...")
  ids_failed <- ids_failed[, chunk := chunk(job.id, chunk.size = nrow(ids_failed))]
  submitJobs(reg = reg, ids = ids_failed, 
             resources = list(walltime = "10:59:00",
                              memory = 1024 * 5,
                              ncpus = 1,
                              partition = partition))
}

# 8. Collect results
ids <- findExperiments(algo.name = "harf_ad_pred", reg = reg)
ids_done <- findDone(ids, reg = reg)
# Load and flat job paramters
job_pars <- getJobPars(ids_done, reg = reg)
job_pars_algo <- rbindlist(job_pars$algo.pars, idcol = "job.id", fill = TRUE)
job_pars_prob <- rbindlist(job_pars$prob.pars, idcol = "job.id", fill = TRUE)
job_pars_DT <- merge(job_pars_algo, job_pars_prob, by = "job.id")
# Load and flat results for HARF
res_harf <- reduceResultsList(ids = ids_done, reg = reg)
res_harf_DT <- rbindlist(lapply(res_harf, function (res) {
  #Rename columns with prefix  MMD_rbk.MMD* by MMD_rbk*
  setnames(res, gsub("MMD_rbk.MMD", "MMD_rbk", names(res)))
  return(res)
}), idcol = "job.id", fill = TRUE)
res_harf_DT$algorithm <- "HARF"
res_harf_DT <- merge(res_harf_DT, job_pars_DT, by = "job.id")
saveRDS(res_harf_DT, file.path(res_hember_dt_dir, "harf_ad.rds"))

# Load and flat results for ARF
ids <- findExperiments(algo.name = "arf_ad_pred", reg = reg)
ids_done <- findDone(ids, reg = reg)
job_pars <- getJobPars(ids_done, reg = reg)
job_pars_algo <- rbindlist(job_pars$algo.pars, idcol = "job.id", fill = TRUE)
job_pars_prob <- rbindlist(job_pars$prob.pars, idcol = "job.id", fill = TRUE)
job_pars_DT <- merge(job_pars_algo, job_pars_prob, by = "job.id")
res_arf <- reduceResultsList(ids = ids_done, reg = reg)
res_arf_DT <- rbindlist(res_arf, idcol = "job.id", fill = TRUE)
res_arf_DT <- rbindlist(lapply(res_arf, function (res) {
  #Rename columns with prefix  MMD_rbk.MMD* by MMD_rbk*
  setnames(res, gsub("MMD_rbk.MMD", "MMD_rbk", names(res)))
  return(res)
}), idcol = "job.id", fill = TRUE)
res_arf_DT$algorithm <- "ARF"
res_arf_DT <- merge(res_arf_DT, job_pars_DT, by = "job.id")
res_arf_DT$chunck_size <- 0
saveRDS(res_arf_DT, file.path(res_hember_dt_dir, "arf_ad_pred.rds"))

# Rbind HARF and ARF results
# Rename ARI_HARF and ARI_ARF by ARI and NMI_HARF and NMI_ARF by NMI
res_harf_DT[, `:=`(ARI_SYN = ARI_HARF, NMI_SYN = NMI_HARF)]
res_arf_DT[, `:=`(ARI_SYN = ARI_ARF, NMI_SYN = NMI_ARF)]
res_harf_DT$ARI_HARF <- NULL
res_harf_DT$NMI_HARF <- NULL
res_arf_DT$ARI_ARF <- NULL
res_arf_DT$NMI_ARF <- NULL
res_all_DT <- rbind(res_harf_DT, res_arf_DT, fill = TRUE)
saveRDS(res_all_DT, file.path(res_hember_dt_dir, "harf_arf_ad_pred.rds"))
