setwd(dirname(this.path::this.path()))
source("../create_sce.R")
### DATA
d <- read.table("data.txt")
# select 5 patients
d <- d[,grepl("MGH26_", colnames(d)) |
       grepl("MGH264_", colnames(d)) |
       grepl("MGH28_", colnames(d)) |
       grepl("MGH29_", colnames(d)) |
       grepl("MGH30_", colnames(d)) |
       grepl("MGH31_", colnames(d))]

### ANNOTATIONS
patients <- unlist(lapply(strsplit(colnames(d), "_"), "[[", 1))
patients[patients == "MGH264"] <- "MGH26"
ann <- data.frame(cell_type1 = patients)
rownames(ann) <- colnames(d)

### SINGLECELLEXPERIMENT
sceset <- create_sce_from_logcounts(d, ann)
# Built a data.table
sceset_dt <- as.data.table(t(logcounts(sceset)))
sceset_dt$cell_type <- sceset$cell_type1
fwrite(x = sceset_dt,
       file = "patel_not_filtered.txt", sep = "\t")
