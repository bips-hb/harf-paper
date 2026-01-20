filter_data <- function(expr_mat, frq_threshold = 0.06) {
  # Filter genes that are expressed in less than frq_threshold fraction of cells
  expr_mat_filtered <- expr_mat[, colSums(expr_mat > 0) > floor(frq_threshold * nrow(expr_mat))]
  return(expr_mat_filtered)
}

org_file_path <- org_simlr_dt_files[1]
load(org_file_path)
mecs <- as.data.table(filter_data(t(Test_1_mECS$in_X)))
mecs$cell_type <- Test_1_mECS$true_labs$V1
fwrite(mecs, file = org_simlr_dt_proc_files[1])

org_file_path <- org_simlr_dt_files[2]
load(org_file_path)
kolod <- as.data.table(filter_data(t(Test_2_Kolod$in_X)))
kolod$cell_type <- Test_2_Kolod$true_labs$V1
fwrite(kolod, file = org_simlr_dt_proc_files[2])

org_file_path <- org_simlr_dt_files[3]
load(org_file_path)
pollen <- as.data.table(filter_data(t(Test_3_Pollen$in_X)))
pollen$cell_type <- Test_3_Pollen$true_labs$V1
fwrite(pollen, file = org_simlr_dt_proc_files[3])

org_file_path <- org_simlr_dt_files[4]
load(org_file_path)
usoskin <- as.data.table(filter_data(t(Test_4_Usoskin$in_X)))
usoskin$cell_type <- Test_4_Usoskin$true_labs$V1
fwrite(usoskin, file = org_simlr_dt_proc_files[4])

org_file_path <- org_simlr_dt_files[5]
load(org_file_path)
zelsel <- as.data.table(filter_data(t(Zelsel$in_X)))
zelsel$cell_type <- Zelsel$true_labs[ , "V1"]
fwrite(zeisel, file = org_simlr_dt_proc_files[5])

