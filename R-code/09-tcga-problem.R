#' Generate TCGA data for predictive modeling
#'
#' @param data Path the currated TCGA data
#' @param job Current job
#' @param evidence To be passed to ARF and h-HARF algorithms
#'
#' @returns A data.frame containing the original data, which will be used as the problem for Batchexperiment.
#'
create_tcga_data <- function (
    data,
    job,
    evidence = FALSE
) {
  # Read data.table
  org_dt <- fread(data$file_name, check.names = FALSE)
  suppressWarnings({
    org_dt$patientID <- NULL
  })
  colnames(org_dt) <- make.names(colnames(org_dt), unique = TRUE)
  char_cols <- names(which(sapply(org_dt, is.character)))
  org_dt[, (char_cols) := lapply(.SD, as.factor), .SDcols = char_cols]
  org_dt[ , ] <- lapply(org_dt, function(x)
    if (is.factor(x)) as.integer(x) else x
  )
  org_dt <- as.data.frame(org_dt)
  # Reduce dataset for 100 variables for testing
   if (ncol(org_dt) > 100) {
    set.seed(123)
    selected_cols <- sample(colnames(org_dt)[!colnames(org_dt) %in% c("years_to_birth",
                                                                    "tumor_stage")], 100)
    org_dt <- org_dt[, c(selected_cols, "years_to_birth", "tumor_stage")]
   }
  # Scale numerical variables
  num_cols <- names(which(sapply(org_dt, is.numeric)))
  org_dt[, num_cols] <- scale(org_dt[, num_cols])
  # Randomly select 10 genes to build an artificial tumor_stage balanced dataset
  selected_genes <- sample(colnames(org_dt)[!colnames(org_dt) %in% c("years_to_birth",
                                                                     "tumor_stage")], floor(0.005 * ncol(org_dt)))
  # Add age as effect variable to selected genes
  selected_genes <- c(selected_genes, "years_to_birth")
  beta <- runif(length(selected_genes), -2, 2)
  # Use logistic distribution to create a more realistic two-class problem
  logit <- as.matrix(org_dt[, selected_genes]) %*% beta
  prob <- 1 / (1 + exp(-logit))
  org_dt$tumor_stage <- ifelse(prob > 0.5, "Late", "Early")
  org_dt$tumor_stage <- factor(org_dt$tumor_stage,
                               levels = c("Early", "Late"))
  n <- nrow(org_dt)
  return(list(data = org_dt,
              file_name = data$file_name, 
              data_name = data$data_name,
              evidence = evidence,
              train_idx = sample(seq_len(nrow(org_dt)), size = floor(0.7 * nrow(org_dt)), replace = FALSE)
  ))
}
