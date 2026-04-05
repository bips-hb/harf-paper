

cluster_and_eval <- function(sc_data) {
  # ️Identify the cell type column
  cell_type_var <- sapply(sc_data, is.character) | sapply(sc_data, is.factor)
  cell_type_col <- which(cell_type_var)[1]          # Take the first factor/character column
  true_labels <- sc_data[[cell_type_col]]
  
  #  Create SingleCellExperiment object
  counts_matrix <- t(as.matrix(sc_data[, -cell_type_col]))
  sce <- SingleCellExperiment(assays = list(counts = counts_matrix))
  
  #  Log-normalization
  logcounts(sce) <- counts(sce)              # safer than raw counts
  
  #  Store cell types
  sce$cell_type <- true_labels
  
  #  PCA (replace rpca with prcomp if rpca is not available)
  pc_res <- prcomp(t(logcounts(sce)), scale. = TRUE)
  projected_data <- pc_res$x %*% pc_res$rotation   # rotated PCs
  
  # Clustering using hierarchical clustering
  dist_matrix <- dist(projected_data)
  hclust_res <- hclust(dist_matrix, method = "ward.D2")
  predicted_clusters <- cutree(hclust_res, k = length(unique(sce$cell_type)))
  
  # Evaluation metrics
  ari <- adjustedRandIndex(predicted_clusters, sce$cell_type)
  nmi <- NMI(predicted_clusters, sce$cell_type)
  
  #  Return as a list
  return(c(ARI = ari, NMI = nmi))
}
