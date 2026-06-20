# GitHub official repository:
# https://github.com/hemberg-lab/scRNA.seq.datasets/tree/master
source(file.path(dirname(this.dir()), "00-library-and-setup.R"))

# Download Hemberger/Lake et al. cell type labels
# ===============================================
# URLs
lake_url <- "http://genome-tech.ucsd.edu/public/Lake_Science_2016/Lake-2016_Gene_TPM.dat.gz"
lake_anno_url <- "http://genome-tech.ucsd.edu/public/Lake_Science_2016/Lake-2016_Gene_TPM_Sample-annotation.txt"
# Destination paths
lake_file <- file.path(org_lake_dt_dir, "Lake-2016_Gene_TPM.dat.gz")
lake_anno_file <- file.path(org_lake_dt_dir,
                            "Lake-2016_Gene_TPM_Sample-annotation.txt")
# Download the files
# ==================
download.file(lake_url, lake_file, mode = "wb")
download.file(lake_anno_url, lake_anno_file, mode = "wb")
# Unzip the TPM data
gunzip(lake_file, remove = TRUE, overwrite = TRUE)  # removes the .gz after extraction
# Clean up: rename unzipped file
file.rename(sub("\\.gz$", "", lake_file),
            file.path(org_lake_dt_dir, "Lake-2016_Gene_TPM.dat"))

# Download Hemberger/Manno et al. data
# ====================================
# http://genome-tech.ucsd.edu/ZhangLab/index.php/data/epigenomics-and-transcriptomics/sns/
# Lake et al. (2016) Neuronal subtypes and diversity revealed by single-nucleus RNA sequencing of the human brain. Science. 352 (6293): 1586-1590
# List of files to download
files <- c(
  "GSE76381_ESMoleculeCounts.cef.txt.gz",
  "GSE76381_EmbryoMoleculeCounts.cef.txt.gz",
  "GSE76381_MouseAdultDAMoleculeCounts.cef.txt.gz",
  "GSE76381_MouseEmbryoMoleculeCounts.cef.txt.gz",
  "GSE76381_iPSMoleculeCounts.cef.txt.gz"
)
base_url <- "ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE76nnn/GSE76381/suppl/"
for (f in files) {
  dest <- file.path(org_manno_dt_dir, f)
  # Download
  download.file(paste0(base_url, f), destfile = dest, mode = "wb")
  # Unzip
  gunzip(dest, remove = TRUE, overwrite = TRUE)  # remove = TRUE deletes the .gz file
}
# Clean up: rename unzipped files
for (f in files) {
  unzipped_file <- sub("\\.gz$", "", f)
  file.rename(file.path(org_manno_dt_dir, unzipped_file),
              file.path(org_manno_dt_dir, sub("\\.cef\\.txt$", ".txt", unzipped_file)))
}

# Download Hemberger/Li et al. data
# =================================
# File paths
li_gz_file <- file.path(org_li_dt_dir, "data.csv.gz")
li_csv_file <- file.path(org_li_dt_dir, "data.csv")
# Download the file
download.file(
  url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FCell%5FLine%5FCOUNT%2Ecsv%2Egz",
  destfile = li_gz_file,
  mode = "wb"
)
# Unzip
gunzip(li_gz_file, destname = li_csv_file, remove = TRUE, overwrite = TRUE)

# Download Hemberger/Patel et al. data
# ====================================
# File paths
patel_gz_file <- file.path(org_patel_dt_dir, "data.csv.gz")
patel_csv_file <- file.path(org_patel_dt_dir, "data.csv")
# Download the file
download.file(
  url = "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE81861&format=file&file=GSE81861%5FCell%5FLine%5FCOUNT%2Ecsv%2Egz",
  destfile = patel_gz_file,
  mode = "wb"
)
# Unzip
gunzip(patel_gz_file, destname = patel_csv_file, remove = TRUE, overwrite = TRUE)

