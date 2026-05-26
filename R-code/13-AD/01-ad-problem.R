#' Generate TCGA data for predictive modeling
#'
#' @param data Path the currated TCGA data
#' @param job Current job
#' @param evidence To be passed to ARF and h-HARF algorithms
#'
#' @returns A data.frame containing the original data, which will be used as the problem for Batchexperiment.
#'
create_ad_data <- function (
    data,
    job,
    evidence = FALSE,
    prop_synth = 1
) {
  # Read data.table
  ad_metab <- fread(data$file_name, check.names = FALSE)
  # Exclude race == 1
  ad_metab$race <- NULL
  metab_clin_feats <- c("age_at_death",
                        "sex",
                        "ApoE_bi",
                        "Braak_bin3")
  metab_feats <- grep(pattern = "*Meta*", colnames(ad_metab), value = TRUE)
  # Also scale age
  ad_metab[, age_at_death := scale(age_at_death)]
  train_feats <- c(metab_clin_feats, metab_feats)
  ad_metab <- ad_metab[, ..train_feats]
  ad_metab <- ad_metab[complete.cases(ad_metab)]
  train_indices <- caret::createDataPartition(ad_metab$Braak_bin3, p = 0.7, list = FALSE)[ , "Resample1"]
  # Scale metabolomics features and choose the 10% top variable features based on training data
  ad_metab_train <- ad_metab[train_indices, ]
  ad_metab_train[, (metab_feats) := lapply(.SD, scale), .SDcols = metab_feats]
  var_threshold <- quantile(apply(ad_metab_train[, ..metab_feats], 2, var), probs = 0.9)
  selected_metab_feats <- metab_feats[apply(ad_metab_train[, ..metab_feats], 2, var) >= var_threshold]
  selected_metab_feats <- selected_metab_feats
  train_feats <- c(metab_clin_feats, selected_metab_feats)
  # Scale the selected metabolomics features in the test set using the training set parameters
  ad_metab_test <- ad_metab[-train_indices, ..train_feats]
  ad_metab_test[, (selected_metab_feats) := lapply(.SD, scale), .SDcols = selected_metab_feats]
  
  return(list(train_data = ad_metab_train[, ..train_feats],
              test_data = ad_metab_test[, ..train_feats],
              metab_clin_feats = metab_clin_feats, 
              metab_feats = selected_metab_feats,
              evidence = evidence,
              prop_synth = prop_synth,
              test_idx = setdiff(seq_len(nrow(ad_metab_train)), train_indices)
  ))
}
