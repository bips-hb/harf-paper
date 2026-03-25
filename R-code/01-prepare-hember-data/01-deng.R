# setwd(org_deng_dt_dir)
### DATA
deng_rpkms <- read.table("deng-rpkms.txt", check.names = F, header = T)
genes <- deng_rpkms[ , 1]
deng_rpkms <- as.matrix(deng_rpkms[ , 2:ncol(deng_rpkms)])
rownames(deng_rpkms) <- genes
cell_ids <- colnames(deng_rpkms)

### ANNOTATIONS
labs <- unlist(lapply(strsplit(cell_ids, "\\."), "[[", 1))
ann <- data.frame(cell_type2 = labs)
labs[labs == "zy" | labs == "early2cell"] = "zygote"
labs[labs == "mid2cell" | labs == "late2cell"] = "2cell"
labs[labs == "earlyblast" | labs == "midblast" | labs == "lateblast"] = "blast"
ann$cell_type1 <- labs
rownames(ann) <- cell_ids

### SINGLECELLEXPERIMENT
deng_rpkms <- create_sce_from_normcounts(deng_rpkms, ann)
# Built a data.table
deng_dt <- as.data.table(t(logcounts(deng_rpkms)))
deng_dt$cell_type <- deng_rpkms$cell_type1
colnames(sceset_dt) <- gsub("^'|'$", "", colnames(sceset_dt))
colnames(sceset_dt) <- gsub("-", "_", colnames(sceset_dt))
fwrite(x = deng_dt,
       file = "embryo_deng_mouse_data.txt", sep = "\t")
# setwd(dirname(this.dir()))
