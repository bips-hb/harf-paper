library(HiClimR)
source(file.path(r_code_dir, "00-library-and-setup.R"))

# Compute pairwise correlations for each hember dataset
lake <- as.data.frame(fread(org_hember_dt_files["lake"]))
li <- as.data.frame(fread(org_hember_dt_files["li"]))
manno <- as.data.frame(fread(org_hember_dt_files["manno"]))
patel <- as.data.frame(fread(org_hember_dt_files["patel"]))

lake_cor <- fastCor(as.matrix(lake[, -which(colnames(lake) == "cell_type")]),
                    , upperTri = TRUE, nSplit = 2)
li_cor <- fastCor(as.matrix(li[, -which(colnames(li) == "cell_type")]),
                  , upperTri = TRUE, nSplit = 2)
manno_cor <- fastCor(as.matrix(manno[, -which(colnames(manno) == "cell_type")]),
                     , upperTri = TRUE, nSplit = 2)
patel_cor <- fastCor(as.matrix(patel[, -which(colnames(patel) == "cell_type")]),
                      , upperTri = TRUE, nSplit = 2)

allcors <- c(lake = mean(lake_cor, na.rm = TRUE),
             li = mean(li_cor, na.rm = TRUE),
             manno = mean(manno_cor, na.rm = TRUE),
             patel = mean(patel_cor, na.rm = TRUE))
print(round(allcors, 3))

# Compute pairwise correlations for each TCGA dataset
luad <- as.data.frame(fread(orig_tcga_data_files["luad"]))
lusc <- as.data.frame(fread(orig_tcga_data_files["lusc"]))
kirc <- as.data.frame(fread(orig_tcga_data_files["kirc"]))
coad <- as.data.frame(fread(orig_tcga_data_files["coad"]))
luad_cor <- fastCor(as.matrix(luad[, -which(colnames(luad) %in% c("gender", "tumor_stage", "years_to_birth"))]),
                     , upperTri = TRUE, nSplit = 2)
lusc_cor <- fastCor(as.matrix(lusc[, -which(colnames(lusc) %in% c("gender", "tumor_stage", "years_to_birth"))]),
                     , upperTri = TRUE, nSplit = 2)
kirc_cor <- fastCor(as.matrix(kirc[, -which(colnames(kirc) %in% c("gender", "tumor_stage", "years_to_birth"))]),
                     , upperTri = TRUE, nSplit = 2)
coad_cor <- fastCor(as.matrix(coad[, -which(colnames(coad) %in% c("gender", "tumor_stage", "years_to_birth"))]),
                     , upperTri = TRUE, nSplit = 2)
allcors_tcga <- c(luad = mean(luad_cor, na.rm = TRUE),
                  lusc = mean(lusc_cor, na.rm = TRUE),
                  kirc = mean(kirc_cor, na.rm = TRUE),
                  coad = mean(coad_cor, na.rm = TRUE))
print(round(allcors_tcga, 3))
