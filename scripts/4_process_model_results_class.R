## Processing results subgroups
rm(list=ls())
library(dplyr)
library(ggplot2)
source("functions.R")

path_true <- 'CKD/output_XGB_class_CKD_2020_12_18__22-07-15'
path_permuted <- 'CKD/output_XGB_class_CKD_2020_12_21__21-58-35_PERMUTED'
data_path <- 'CKD/input_data'
labels <- c("CKD", "Healthy")

compared_to_permuted_class(path_true, path_permuted)
plot_feature_importance_class(path_true, 30)
plot_features_tests_class(data_path, path_true, top_n=15, labels)
plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)
plot_feature_importance_class_microbiome(path_true, top_n=20)


path_true <- 'CKD/WithoutDM/output_XGB_class_CKDnoDM_2020_12_19__11-53-08'
path_permuted <- 'CKD/WithoutDM/output_XGB_class_CKDNoDM_2020_12_21__16-46-24_PERMUTED'
data_path <- 'CKD/WithoutDM/input_data'
labels <- c("CKD", "Healthy")

compared_to_permuted_class(path_true, path_permuted)
plot_feature_importance_class(path_true, 30)
plot_features_tests_class(data_path, path_true, top_n=15, labels)
plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)
plot_feature_importance_class_microbiome(path_true, top_n=20)


path_true <- 'CKDDM/output_XGB_class_CKDDM_2020_12_19__13-54-29'
path_permuted <- 'CKDDM/output_XGB_class_CKDDM_2020_12_22__09-38-25_PERMUTED'
data_path <- 'CKDDM/input_data'
labels <- c("CKD+T2DM", "CKD-T2DM")

compared_to_permuted_class(path_true, path_permuted)
plot_feature_importance_class(path_true, 30)
plot_features_tests_class(data_path, path_true, top_n=15, labels)
plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)
plot_feature_importance_class_microbiome(path_true, top_n=20)

path_true <- 'Urine/CKD/output_XGB_class_UrineCKD_2021_01_11__14-46-07'
#path_permuted <- 'Urine/CKD/output_XGB_class_CKDDM_2020_12_22__09-38-25_PERMUTED'
data_path <- 'Urine/CKD/input_data'
labels <- c("CKD", "Healthy")

#compared_to_permuted_class(path_true, path_permuted)
plot_feature_importance_class(path_true, 30)
plot_features_tests_class(data_path, path_true, top_n=15, labels)
plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)
plot_feature_importance_class_urine(path_true, top_n=20)

path_true <- 'Urine/CKD/WithoutDM/output_XGB_class_UrineCKDsens_2021_01_11__22-53-32'
#path_permuted <- 'Urine/CKD/output_XGB_class_CKDDM_2020_12_22__09-38-25_PERMUTED'
data_path <- 'Urine/CKD/WithoutDM/input_data'
labels <- c("CKD", "Healthy")

#compared_to_permuted_class(path_true, path_permuted)
plot_feature_importance_class(path_true, 30)
plot_features_tests_class(data_path, path_true, top_n=15, labels)
plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)
plot_feature_importance_class_urine(path_true, top_n=20)


path_true <- 'Urine/CKDDM/output_XGB_class_UrineCKDDM_2021_01_12__10-30-18'
#path_permuted <- 'Urine/CKD/output_XGB_class_UrineCKDDM_2021_01_12__10-30-18_PERMUTED'
data_path <- 'Urine/CKDDM/input_data'
labels <- c("CKD+T2DM", "CKD-T2DM")

#compared_to_permuted_class(path_true, path_permuted)
plot_feature_importance_class(path_true, 30)
plot_features_tests_class(data_path, path_true, top_n=15, labels)
plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)
plot_feature_importance_class_urine(path_true, top_n=20)
