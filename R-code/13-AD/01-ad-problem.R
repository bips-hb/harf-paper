#' Generate TCGA data for predictive modeling
#'
#' @param data Path the currated TCGA data
#' @param job Current job
#' @param evidence To be passed to ARF and h-HARF algorithms
#'
#' @returns A data.frame containing the original data, which will be used as the problem for Batchexperiment.
#'
create_ad_data <- function(data, job, evidence = FALSE, prop_synth = 1) {
  
  ad_metab <- data.table::fread(data$file_name, check.names = FALSE)
  ad_metab$race <- NULL
  ad_metab$ApoE_bi <- ifelse(ad_metab$ApoE_bi == 2, 1L, 0L)
  metab_clin_feats <- c("age_at_death", "sex", "ApoE_bi", "Braak_bin3")
  metab_feats <- grep("*Meta*", colnames(ad_metab), value = TRUE)
  
  ad_metab <- ad_metab[, c(metab_clin_feats, metab_feats), with = FALSE]
  ad_metab <- ad_metab[complete.cases(ad_metab)]
  
  ad_metab$Braak_bin3 <- factor(ad_metab$Braak_bin3, levels = c(0, 1))
  
  train_idx <- caret::createDataPartition(ad_metab$Braak_bin3, p = 0.7, list = FALSE)[,1]
  
  train <- ad_metab[train_idx]
  test  <- ad_metab[-train_idx]
  
  # -----------------------
  # FEATURE SELECTION
  # -----------------------
  feat_var <- apply(train[, ..metab_feats], 2, var, na.rm = TRUE)
  variance_threshold <- quantile(feat_var, 0.5, na.rm = TRUE)
  # Choose only 100 for testing, but in practice we can use all above the threshold
  # selected_metab_feats <- names(sort(feat_var, decreasing = TRUE))[1:100]
  selected_metab_feats <- names(feat_var[feat_var >= variance_threshold])
  train_feats <- c(metab_clin_feats, selected_metab_feats)
  
  # NOW SUBSET CONSISTENTLY
  train <- train[, ..train_feats]
  test  <- test[, ..train_feats]
  
  # -----------------------
  # SAFE SCALING
  # -----------------------
  means <- sapply(train[, ..selected_metab_feats], mean)
  sds   <- sapply(train[, ..selected_metab_feats], sd)
  
  for (f in selected_metab_feats) {
    train[[f]] <- (train[[f]] - means[[f]]) / sds[[f]]
    test[[f]]  <- (test[[f]]  - means[[f]]) / sds[[f]]
  }
  
  age_mean <- mean(train$age_at_death)
  age_sd   <- sd(train$age_at_death)
  
  train[, age_at_death := (age_at_death - age_mean) / age_sd]
  test[, age_at_death := (age_at_death - age_mean) / age_sd]
  
  list(
    train_data = train,
    test_data  = test,
    metab_clin_feats = metab_clin_feats,
    metab_feats = selected_metab_feats,
    evidence = evidence,
    prop_synth = prop_synth
  )
}