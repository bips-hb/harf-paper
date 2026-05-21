library(HiClimR)
source(file.path(r_code_dir, "00-library-and-setup.R"))

lake <- as.data.frame(fread(org_hember_dt_files["lake"]))
li <- as.data.frame(fread(org_hember_dt_files["li"]))
manno <- as.data.frame(fread(org_hember_dt_files["manno"]))
patel <- as.data.frame(fread(org_hember_dt_files["patel"]))

# Compute pairwise correlations for each dataset
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
print(allcors)
