##### Human brain Manno et al. data preprocessing

## Download data from https://www.nature.com/articles/s41467-018-04997-2#Sec15 and save the three files "GSE76381_ESMoleculeCounts.txt", "GSE76381_EmbryoMoleculeCounts.txt", and "GSE76381_iPSMoleculeCounts.txt" in the same directory as this script.
system("bash ./manno.sh")
### DATA
fix_name <- function(name){
        tmp <- strsplit(name, "_")
        return(tmp[[1]][1])
}
fix_type <- function(x) {as.numeric(as.character(x))}

x1 = read.delim("GSE76381_ESMoleculeCounts.cef.txt", sep="\t", header=F)
cellid = x1[2,]
celltype1=x1[3,]
timepoint1=x1[4,]
x1 <- x1[-c(1:5),]
rownames(x1)<- x1[,1]
x1<-x1[,-c(1,2)]
colnames(x1) <- as.character(unlist(cellid))[-c(1,2)]
tmp_rownames <- rownames(x1);
x1 <- apply(x1, 2, fix_type)
rownames(x1) <- tmp_rownames
x2 = read.delim("GSE76381_EmbryoMoleculeCounts.cef.txt", sep="\t", header=F)
cellid = x2[2,]
celltype2=x2[3,]
timepoint2=x2[4,]
x2 <- x2[-c(1:5),]
rownames(x2)<- x2[,1]
x2<-x2[,-c(1,2)]
colnames(x2) <- as.character(unlist(cellid))[-c(1,2)]
tmp_rownames <- rownames(x2);
x2 <- apply(x2, 2, fix_type)
rownames(x2) <- tmp_rownames
x3 = read.delim("GSE76381_iPSMoleculeCounts.cef.txt", sep="\t", header=F)
cellid = x3[2,]
celltype3=x3[3,]
timepoint3=x3[4,]
x3 <- x3[-c(1:5),]
rownames(x3)<- x3[,1]
x3<-x3[,-c(1,2)]
colnames(x3) <- as.character(unlist(cellid))[-c(1,2)]
tmp_rownames <- rownames(x3);
x3 <- apply(x3, 2, fix_type)
rownames(x3) <- tmp_rownames

### ANNOTATIONS
all_genes <- sort(unique(c(rownames(x1), rownames(x2), rownames(x3))))
all_symbol <- sapply(all_genes, fix_name)
x1_order <- match(all_genes, rownames(x1))
x1 <- x1[x1_order,]
x1[is.na(x1)] <- 0
x2_order <- match(all_genes, rownames(x2))
x2 <- x2[x2_order,]
x2[is.na(x2)] <- 0
x3_order <- match(all_genes, rownames(x3))
x3 <- x3[x3_order,]
x3[is.na(x3)] <- 0
DATA <- cbind(x1, x2, x3)
rownames(DATA) <- all_genes
TYPE <- c(as.character(unlist(celltype1))[-c(1,2)], as.character(unlist(celltype2))[-c(1,2)], as.character(unlist(celltype3))[-c(1,2)])
AGE <- c(as.character(unlist(timepoint1))[-c(1,2)], as.character(unlist(timepoint2))[-c(1,2)], as.character(unlist(timepoint3))[-c(1,2)])
stuff <- matrix(unlist(strsplit(colnames(DATA), "-|_")), ncol=3, byrow=T)
WELL <- stuff[,3]
SOURCE <- rep(c("ESCs", "ventral midbrain", "iPSCs"), times=c(length(x1[1,]), length(x2[1,]), length(x3[1,])))
ANN <- data.frame(Species=rep("Homo sapiens", times=length(DATA[1,])), cell_type1=TYPE, Source=SOURCE, age=AGE, WellID=WELL, batch=paste(stuff[,1], stuff[,2]))
rownames(ANN) = colnames(DATA)

### SINGLECELLEXPERIMENT
sceset <- create_sce_from_counts(DATA, ANN)

# Built a data.table
sceset_dt <- as.data.table(t(counts(sceset)))
sceset_dt$cell_type <- sceset$cell_type1
colnames(sceset_dt) <- gsub("^'|'$", "", colnames(sceset_dt))
colnames(sceset_dt) <- gsub("-", "_", colnames(sceset_dt))
fwrite(x = sceset_dt,
       file = "manno_not_filtered.txt", sep = "\t")
