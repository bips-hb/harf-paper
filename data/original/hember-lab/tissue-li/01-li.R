setwd(dirname(this.path::this.path()))
source("../create_sce.R")
### DATA
d <- read.csv("data.csv")

### ANNOTATIONS
genes <- unlist(lapply(strsplit(as.character(d[,1]), "_"), "[[", 2))
d <- d[!duplicated(genes), ]
rownames(d) <- genes[!duplicated(genes)]
d <- d[,2:ncol(d)]
# metadata
ann <- data.frame(cell_type1 = unlist(lapply(strsplit(colnames(d), "__"), "[[", 2)))
rownames(ann) <- colnames(d)

### SINGLECELLEXPERIMENT
sceset <- create_sce_from_counts(d, ann)
# Built a data.table
sceset_dt <- as.data.table(t(logcounts(sceset)))
sceset_dt$cell_type <- sceset$cell_type1
fwrite(x = sceset_dt,
       file = "li_not_filtered.txt", sep = "\t")
