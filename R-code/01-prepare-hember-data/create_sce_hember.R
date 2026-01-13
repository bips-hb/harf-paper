library(SingleCellExperiment)
library(scater)

create_sce_from_counts <- function(counts, colData, rowData = NULL) {
  if(is.null(rowData)) {
    sceset <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), 
                                   colData = colData)
  } else {
    sceset <- SingleCellExperiment(assays = list(counts = as.matrix(counts)), 
                                   colData = colData,
                                   rowData = rowData)
  }
  # this function writes to logcounts slot
  exprs(sceset) <- log1p(scater::calculateCPM(sceset, size.factors = NULL))
  # use gene names as feature symbols
  rowData(sceset)$feature_symbol <- rownames(sceset)
  # remove features with duplicated names
  if(is.null(rowData)) {
    sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]
  }
  # QC
  # isSpike(sceset, "ERCC") <- grepl("^ERCC-", rownames(sceset))
  if (sum(grepl("^ERCC-", rownames(sceset)))) {
    altExp(sceset, "ERCC") <- grepl("^ERCC-", rownames(sceset))
    sceset <- calculateQCMetrics(sceset, feature_controls = list("ERCC" = altExp(sceset, "ERCC")))
  } 
  return(sceset)
}

create_sce_from_normcounts <- function(normcounts, colData, rowData = NULL) {
  if(is.null(rowData)) {
    sceset <- SingleCellExperiment(assays = list(normcounts = as.matrix(normcounts)), 
                                   colData = colData)
  } else {
    sceset <- SingleCellExperiment(assays = list(normcounts = as.matrix(normcounts)), 
                                   colData = colData,
                                   rowData = rowData)
  }
  logcounts(sceset) <- log1p(normcounts(sceset))
  # use gene names as feature symbols
  rowData(sceset)$feature_symbol <- rownames(sceset)
  # remove features with duplicated names
  if(is.null(rowData)) {
    sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]
  }
  # QC
  if (sum(grepl("^ERCC-", rownames(sceset)))) {
    altExp(sceset, "ERCC") <- grepl("^ERCC-", rownames(sceset))
    sceset <- calculateQCMetrics(sceset, feature_controls = list("ERCC" = altExp(sceset, "ERCC")))
  } 
  return(sceset)
}

create_sce_from_logcounts <- function(logcounts, colData, rowData = NULL) {
  if(is.null(rowData)) {
    sceset <- SingleCellExperiment(assays = list(logcounts = as.matrix(logcounts)), 
                                   colData = colData)
  } else {
    sceset <- SingleCellExperiment(assays = list(logcounts = as.matrix(logcounts)), 
                                   colData = colData,
                                   rowData = rowData)
  }
  # use gene names as feature symbols
  rowData(sceset)$feature_symbol <- rownames(sceset)
  # remove features with duplicated names
  if(is.null(rowData)) {
    sceset <- sceset[!duplicated(rowData(sceset)$feature_symbol), ]
  }
  # QC
  if (sum(grepl("^ERCC-", rownames(sceset)))) {
    altExp(sceset, "ERCC") <- grepl("^ERCC-", rownames(sceset))
    sceset <- calculateQCMetrics(sceset, feature_controls = list("ERCC" = altExp(sceset, "ERCC")))
  } 
  return(sceset)
}
