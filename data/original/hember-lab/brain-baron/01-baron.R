setwd(dirname(this.path::this.path()))
source("../create_sce.R")
### DATA
# human1
h1 <- read.csv("GSM2230757_human1_umifm_counts.csv", header = T)
rownames(h1) <- h1[,1]
labels_h1 <- as.character(h1$assigned_cluster)
h1 <- h1[,4:ncol(h1)]
h1 <- t(h1)
# human2
h2 <- read.csv("GSM2230758_human2_umifm_counts.csv", header = T)
rownames(h2) <- h2[,1]
labels_h2 <- as.character(h2$assigned_cluster)
h2 <- h2[,4:ncol(h2)]
h2 <- t(h2)
# human3
h3 <- read.csv("GSM2230759_human3_umifm_counts.csv", header = T)
rownames(h3) <- h3[,1]
labels_h3 <- as.character(h3$assigned_cluster)
h3 <- h3[,4:ncol(h3)]
h3 <- t(h3)
# human4
h4 <- read.csv("GSM2230760_human4_umifm_counts.csv", header = T)
rownames(h4) <- h4[,1]
labels_h4 <- as.character(h4$assigned_cluster)
h4 <- h4[,4:ncol(h4)]
h4 <- t(h4)

# merge data
h <- cbind(h1, h2, h3, h4)
rm(h1, h2, h3, h4)
### ANNOTATIONS
# human
h_ann <- data.frame(
    human = c(
        rep(1, length(labels_h1)),
        rep(2, length(labels_h2)),
        rep(3, length(labels_h3)),
        rep(4, length(labels_h4))
    ),
    cell_type1 = c(labels_h1, labels_h2, labels_h3, labels_h4))
rownames(h_ann) <- colnames(h)


### SINGLECELLEXPERIMENT
h_sceset <- create_sce_from_counts(h, h_ann)
# Built a data.table
sceset_dt <- as.data.table(t(counts(h_sceset)))
sceset_dt$cell_type <- h_sceset$cell_type1

fwrite(x = sceset_dt,
       file = "baron_not_filtered.txt", sep = "\t")
