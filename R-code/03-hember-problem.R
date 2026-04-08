#' Generate data as problem for Batchexperiment
#'
#' @param data Batchexperiment formal argument
#' @param job Batchexperiment formal argument
#' @param file_name Path to data file
#' @param data_name Name of the dataset
#'
#' @returns A data.frame containing the original data, which will be used as the problem for Batchexperiment.
#'
create_single_cell_data <- function (
    data,
    job,
    evidence = FALSE
) {
  # Read data.table
  org_dt <- as.data.frame(fread(data$file_name, check.names = FALSE))
  colnames(org_dt) <- make.names(colnames(org_dt), unique = TRUE)
  # Use 100 gene expressions and cell type for testing
  # org_dt <- org_dt[, c(1:100, which(colnames(org_dt) == "cell_type"))]
  n <- nrow(org_dt)
  return(list(data = org_dt,
              file_name = data$file_name, 
              data_name = data$data_name,
              evidence = evidence,
              train_idx = sample(seq_len(nrow(org_dt)), size = floor(0.7 * nrow(org_dt)), replace = FALSE)
              )
         )
}