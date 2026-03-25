# setwd(org_li_dt_dir)
## Download the data from https://www.nature.com/articles/s41467-019-13025-7 and save it as "data.csv" in the same directory as this script.
system("bash ./li.sh")
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
colnames(sceset_dt) <- gsub("^'|'$", "", colnames(sceset_dt))
colnames(sceset_dt) <- gsub("-", "_", colnames(sceset_dt))
fwrite(x = sceset_dt,
       file = "li_not_filtered.txt", sep = "\t")
# setwd(dirname(this.dir()))
