source(file.path(dirname(this.dir()), "00-library-and-setup.R"))

download_tcga_dt <- function(disease_code, file_nms) {
  message("Downloading and processing data for disease: ", disease_code)
  disease_code <- tolower(disease_code)
  # Downloading
  se <- curatedTCGAData(
    diseaseCode = disease_code,
    assays = "RNASeq2GeneNorm",
    version = "2.1.1",
    dry.run = FALSE
  )
  # Extract the assay data
  X <- assay(se[[1]])
  X <- X[complete.cases(X), ]
  X <- t(X)
  # Extract the clinical data
  if (disease_code %in% c("brca")) {
    clinical_var <- c("patientID",
                      "years_to_birth",
                      "Age.at.Initial.Pathologic.Diagnosis",
                      "gender",
                      "AJCC.Stage")
    clinical_data <- se@colData[ , clinical_var]
    colnames(clinical_data) <- c("patientID",
                                 "years_to_birth",
                                 "age_at_diagnosis",
                                 "gender",
                                 "tumor_stage")
    clinical_data <- clinical_data[!clinical_data$tumor_stage %in% c("[Not Available]", "Stage X"), ]
    stage <- clinical_data$tumor_stage
    clinical_data <- clinical_data[!clinical_data$tumor_stage %in% c("[Not Available]", "Stage X"), ]
    stage <- clinical_data$tumor_stage
    stage_bin <- ifelse(grepl("^Stage (I|II)", stage), "Early", "Late")
  }
  if (disease_code %in% c("luad", "lusc", "coad")) {
    gender <- "gender"
    if (disease_code %in% c("luad", "lusc")) {
      age_at_diagnosis <- "Age.at.diagnosis"
    }
    if (disease_code == "lusc") {
      age_at_diagnosis <- "Age.at.diagnosis"
    }
    if (disease_code == "coad") {
      age_at_diagnosis <- "age_at_initial_pathologic_diagnosis"
      gender <- "gender.x"
    }
    clinical_var <- c("patientID",
                      "years_to_birth",
                      age_at_diagnosis,
                      gender,
                      "pathologic_stage")
    clinical_data <- se@colData[ , clinical_var]
    colnames(clinical_data) <- c("patientID",
                                 "years_to_birth",
                                 "age_at_diagnosis",
                                 "gender",
                                 "tumor_stage")
    clinical_data <- clinical_data[!clinical_data$tumor_stage %in% c("Not Reported", "Stage X"), ]
    stage <- clinical_data$tumor_stage
    stage_bin <- ifelse(
      grepl("^stage (i|ii)(|a|b|c)?$", stage), "Early",
      ifelse(grepl("^stage (iii|iv)(a|b|c)?$", stage), "Late", NA)
    )
  }
  if (disease_code %in% c("kirc")) {
    clinical_var <- c("patientID",
                      "years_to_birth",
                      "patient.age_at_initial_pathologic_diagnosis",
                      "gender",
                      "pathologic_stage")
    clinical_data <- se@colData[ , clinical_var]
    colnames(clinical_data) <- c("patientID",
                                 "years_to_birth",
                                 "age_at_diagnosis",
                                 "gender",
                                 "tumor_stage")
    clinical_data <- clinical_data[!clinical_data$tumor_stage %in% c("Not Reported", "Stage X"), ]
    stage <- clinical_data$tumor_stage
    stage_bin <- ifelse(
      grepl("stage (i|ii)$", stage), "Early", "Late")
  }
  stage_bin <- factor(stage_bin, levels = c("Early", "Late"))
  clinical_data$tumor_stage <- stage_bin
  # Merge expression and clinical data
  sample_barcodes <- rownames(X)
  patient_ids <- substr(sample_barcodes, 1, 12)
  X <- X[!duplicated(patient_ids), ]
  vars <- apply(X, 2, var)
  thr  <- quantile(vars, 0.2)   # remove lowest 20%
  X <- X[, vars > thr]
  patient_ids <- patient_ids[!duplicated(patient_ids)]
  X <- X[order(patient_ids), ]
  X_DT <- as.data.table(X)
  X_DT$patientID <- patient_ids
  X_clinical <- merge(X_DT, clinical_data, by = "patientID")
  X_clinical$patientID <- NULL
  X_clinical <- X_clinical[complete.cases(X_clinical), ]
  X_clinical <- as.data.frame(X_clinical)
  X_clinical$tumor_stage <- as.factor(X_clinical$tumor_stage)
  # Save data
  fwrite(
    X_clinical,
    file = file_nms[disease_code],
    sep = "\t"
  )
  return(c(disease = disease_code, n = nrow(X_clinical), p = ncol(X_clinical)))
}

diseases <- c(
  "brca",  # Breast invasive carcinoma
  "luad",  # Lung adenocarcinoma
  "lusc",  # Lung squamous cell carcinoma
  "coad",  # Colon adenocarcinoma
  "kirc"   # Kidney renal clear cell carcinoma
)
lapply(
  diseases,
  download_tcga_dt,
  file_nms = orig_tcga_data_files
)
