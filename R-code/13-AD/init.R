source("../00-library-and-setup.R")
# library(data.table)
# library(caret)
# library(glmnet)
# library(ranger)
# library(pROC)
# library(ggplot2)

# ad_dir <- "R://ahuels/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/AD/"
ad_dir <- "/huels_lab/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/AD"
metab_file <- file.path(ad_dir, "metabolomics_final.txt")
ad_results_dir <- file.path(res_dir, "ad")
dir.create(ad_results_dir, showWarnings = FALSE)
# ad_metab <- fread(metab_file, header = TRUE, data.table = FALSE)

# Exclude race == 1
# ad_metab <- ad_metab[ad_metab$race != 1, ]
# ad_metab$race <- NULL
# 
# metab_clin_feats <- c("age_at_death",
#                       "sex",
#                       "ApoE_bi",
#                       "Braak_bin3")
# metab_feats <- grep(pattern = "*Meta*", colnames(ad_metab), value = TRUE)
# 
# clin_data <- ad_metab[, metab_clin_feats]
# metab_data <- ad_metab[, metab_feats]
# 
# # Split data into training and testing set, with stratification by Braak_bin3
# #set.seed(123)
# train_indices <- caret::createDataPartition(clin_data$Braak_bin3, p = 0.7, list = FALSE)
# train_clin <- clin_data[train_indices, ]
# train_metab <- metab_data[train_indices, ]
# test_clin <- clin_data[-train_indices, ]
# test_metab <- metab_data[-train_indices, ]
# 
# # Combine clinical and metabolomics data for training and testing using glmnet
# train_data <- cbind(train_clin, train_metab)
# test_data <- cbind(test_clin, test_metab)
# # train_data$ApoE_bi <- NULL
# # test_data$ApoE_bi <- NULL
# 
# #
# 
# # ========================================
# # Random forest
# # ========================================
# rf_model <- ranger::ranger(
#   x = train_data[, -which(colnames(train_data) == "Braak_bin3")],
#   y = as.factor(train_data$Braak_bin3),
#   num.trees = 5000,
#   mtry = floor(sqrt(ncol(train_data) - 1)),
#   min.node.size = 5,
#   probability = TRUE,
# )
# print(rf_model$prediction.error)
# 
# rf_pred <- predict(rf_model,
#                    data = test_data[, -which(colnames(test_data) == "Braak_bin3")],
#                    type = "response")$predictions
# rf_auc <- pROC::roc(test_data$Braak_bin3, rf_pred[, 1])$auc
# print(paste("Random Forest AUC:", rf_auc))
# 
# library(mltools)
# threshold_rf <- coords(rf_auc , "best", ret = "threshold")
# mcc_rf <- mltools::mcc(preds = ifelse(rf_pred[, 2] > threshold_rf$threshold,
#                                       1, 0), actuals = test_data$Braak_bin3)
# print(paste("Random Forest MCC:", mcc_rf))
# 
# # ========================================
# # Lasso regression
# # ========================================
# x_train <- as.matrix(train_data[, -which(colnames(train_data) == "Braak_bin3")])
# y_train <- as.factor(train_data$Braak_bin3)
# lasso_model <- glmnet::cv.glmnet(
#   x = x_train,
#   y = y_train,
#   alpha = 1,
#   family = "binomial",
#   nfolds = 5
# )
# # Print not null model coefficients at lambda.min
# print(coef(lasso_model, s = "lambda.min")[coef(lasso_model, s = "lambda.min")[, 1] != 0, ])
# print(lasso_model$lambda.min)
# lasso_pred <- predict(lasso_model, newx = as.matrix(test_data[, -which(colnames(test_data) == "Braak_bin3")]), s = "lambda.min", type = "response")
# lasso_auc <- pROC::roc(test_data$Braak_bin3, as.vector(lasso_pred))#$auc
# print(paste("Lasso AUC:", lasso_auc))
# 
# threshold_lasso <- coords(lasso_auc , "best", ret = "threshold")
# 
# mcc_lasso <- mltools::mcc(preds = ifelse(as.vector(lasso_pred) > threshold_lasso$threshold, 1, 0),
#                           actuals = test_data$Braak_bin3)
# print(paste("Lasso MCC:", mcc_lasso))
# 
# 
# # ========================================
# # Generate synthetic data with h-HARF
# # ========================================
# train_data <- as.data.frame(train_data)
# train_data$Braak_bin3 <- as.factor(train_data$Braak_bin3)
# train_data <- train_data[complete.cases(train_data), ]
# harf_model <- h_arf(
#   omx_data = train_data[ , metab_feats[1:100]],
#   cli_lab_data = train_data[ , metab_clin_feats],
#   feature_ordering = c(metab_feats[1:100], metab_clin_feats),
#   parallel = FALSE,   
#   verbose = TRUE,
#   target = "Braak_bin3",
#   chunck_size = 5
# )
# 
# pc_data <- harf_model$meta_features
# # ranger model using the meta-features
# rf_harf_model <- ranger::ranger(
#   x = pc_data[, -which(colnames(pc_data) == "Braak_bin3")],
#   y = as.factor(pc_data$Braak_bin3),
#   num.trees = 5000,
#   min.node.size = 5,
#   classification = TRUE
# )
# rf_harf_model$prediction.error
# 
# # Generate synthetic data with h_forge
# harf_syn <- h_forge(
#   harf_obj = harf_model,
#   n_synth = nrow(train_data),
#   evidence = NULL,
#   parallel = FALSE,   
#   verbose = TRUE
# )
# 
# 
# # ranger model using the original features
# rf_orig_model <- ranger::ranger(
#   x = train_data[ , -which(colnames(harf_syn) == "Braak_bin3")],#[, c(metab_clin_feats, metab_feats[1:100])],
#   y = as.factor(train_data$Braak_bin3),
#   num.trees = 5000,
#   min.node.size = 5,
#   probability = TRUE
# )
# rf_orig_model$prediction.error
# 
# # ranger model using the synthetic data
# harf_syn <- as.data.frame(harf_syn)
# rf_syn_model <- ranger::ranger(
#   x = harf_syn[, -which(colnames(harf_syn) == "Braak_bin3")],
#   y = as.factor(harf_syn$Braak_bin3),
#   num.trees = 5000,
#   min.node.size = 5,
#   probability = TRUE
# )
# rf_syn_model$prediction.error
# 
# # Predict on test set using the original model
# rf_orig_pred <- predict(
#   rf_orig_model,
#   data = test_data,#[, c(metab_clin_feats, metab_feats[1:100])],
#   type = "response")$predictions
# rf_orig_auc <- pROC::roc(test_data$Braak_bin3, rf_orig_pred[, 1])$auc
# print(paste("Random Forest Original AUC:", rf_orig_auc))
# threshold_rf_orig <- coords(rf_orig_auc , "best", ret = "threshold")
# mcc_rf_orig <- mltools::mcc(preds = ifelse(rf_orig_pred[, 2] > threshold_rf_orig$threshold, 1, 0), actuals = test_data$Braak_bin3)
# print(paste("Random Forest Original MCC:", mcc_rf_orig))
# 
# # Predict on test set using the synthetic model
# rf_syn_pred <- predict(
#   rf_syn_model,
#   data = test_data[, c(metab_clin_feats, metab_feats[1:100])],
#   type = "response")$predictions
# rf_syn_auc <- pROC::roc(test_data$Braak_bin3, rf_syn_pred
#                          [, 1])$auc
# print(paste("Random Forest Synthetic AUC:", rf_syn_auc))
# threshold_rf_syn <- coords(rf_syn_auc , "best", ret = "threshold")
# mcc_rf_syn <- mltools::mcc(preds = ifelse(rf_syn_pred[, 2] > threshold_rf_syn$threshold, 1, 0), actuals =
#                          test_data$Braak_bin3)
# print(paste("Random Forest Synthetic MCC:", mcc_rf_syn))
