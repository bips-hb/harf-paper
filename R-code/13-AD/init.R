source("../00-library-and-setup.R")
# library(data.table)
# library(caret)
# library(glmnet)
# library(ranger)
# library(pROC)
# library(ggplot2)

#ad_dir <- "R://ahuels/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/AD/"
ad_dir <- "/huels_lab/AIRCO/01_projects/019_adrc_bb_prediction_machine_learning/harf-paper/AD"
metab_file <- file.path(ad_dir, "metabolomics_final.txt")
# ad_metab <- fread(metab_file, header = TRUE, data.table = FALSE)
# 
# # Exclude race == 1
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
# set.seed(123)
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
