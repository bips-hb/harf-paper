cluster_and_eval <- function(sc_data) {
  sc_data <- as.data.frame(sc_data) 
  
  # Create SingleCellExperiment object
  sce <- SingleCellExperiment::SingleCellExperiment(
    assays = list(counts = t(as.matrix(sc_data[ , - which(colnames(sc_data)  == "cell_type")])))
  )
  SingleCellExperiment::logcounts(sce) <- SingleCellExperiment::counts(sce) # Log-normalization
  sce$cell_type <- sc_data$cell_type
  # Perform clustering using SC3 algorithm
  rowData(sce)$feature_symbol <- paste0("Gene", seq_len(nrow(sce)))
  sce <- SC3::sc3(sce, ks = length(unique(sce$cell_type)), gene_filter = FALSE, n_cores = 1)
  # Evaluation metrics
  ari <- adjustedRandIndex(colData(sce)[, 1], colData(sce)[, 2])
  nmi <- NMI(colData(sce)[, 1], colData(sce)[, 2])
  
  #  Return as a list
  return(c(ARI = ari, NMI = nmi))
}
