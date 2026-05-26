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
    evidence = FALSE
) {
  # Read data.table
  ad_metab <- fread(data$file_name, check.names = FALSE)
  # Exclude race == 1
  ad_metab <- ad_metab[ad_metab$race != 1, ]
  ad_metab$race <- NULL
  metab_clin_feats <- c("age_at_death",
                        "sex",
                        "ApoE_bi",
                        "Braak_bin3")
  metab_feats <- grep(pattern = "*Meta*", colnames(ad_metab), value = TRUE)
  train_feats <- c(metab_clin_feats, metab_feats)
  train_indices <- caret::createDataPartition(ad_metab$Braak_bin3, p = 0.7, list = FALSE)
  return(list(data = ad_metab[, ..train_feats],
              metab_clin_feats = metab_clin_feats, 
              metab_feats = metab_feats,
              evidence = evidence,
              train_idx = train_indices,
              test_idx = setdiff(seq_len(nrow(ad_metab)), train_indices)
  ))
}
