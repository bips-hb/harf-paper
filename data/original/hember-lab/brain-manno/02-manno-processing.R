setwd(dirname(this.path::this.path()))
# Preprocess Baron et al. human pancreas data

cell_data <- fread("manno_not_filtered.txt")
cell_annotation <- data.table(cell_type = cell_data$cell_type)
cell_data$cell_type <- NULL
# Remove constant genes
cell_data <- cell_data[apply(cell_data, 1, function(x) var(as.numeric(x), na.rm = TRUE)) != 0, ]
# lop1 transformation
cell_data <- as.data.frame(cell_data)
# Create singleCellExperiment object
sce <- SingleCellExperiment(assays = list(counts = t(cell_data)))
sce$cell_type <- cell_annotation$cell_type
# Frequency filtering (FRQ): filter cells genes expressed in less than 6% of the cells
freq_threshold <- 0.06
gene_freq <- rowSums(assay(sce) > 0) / ncol(sce)
sce <- sce[gene_freq >= freq_threshold, ]
# Highly expressed genes (HiE): select genes with expression above the 90th percentile

n_genes <- nrow(counts(sce))
top_n <- ceiling(0.10 * n_genes)

hie_matrix <- apply(counts(sce), 2, function(x) {
  idx <- order(x, decreasing = TRUE)[seq_len(top_n)]
  hie <- logical(length(x))
  hie[idx] <- TRUE
  hie
})
rownames(hie_matrix) <- rownames(counts)
min_cells <- ceiling(0.10 * ncol(sce))
keep_genes <- rowSums(hie_matrix) >= min_cells
sce <- sce[keep_genes, ]

# Save sce as data.table
sce_dt <- as.data.table(t(assay(sce)))
sce_dt[, cell_type := sce$cell_type]
fwrite(sce_dt, file = "processed_brain_manno_data.csv")

