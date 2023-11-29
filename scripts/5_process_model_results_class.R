## Processing machine learning results

library(dplyr)
library(tidyr)
library(ggplot2)
source("scripts/functions.R")

## Plasma metabolites
path_true <- 'CKD_plasma/output_XGB_class_nonD_CKD_plasma_2023_11_28__23-15-43'
path_permuted <- 'CKD_plasma/output_XGB_class_nonD_CKD_plasma_2023_11_29__00-17-53_PERMUTED'
data_path <- 'CKD_plasma/input_data'
labels <- c("CKD", "No CKD")

compared_to_permuted_class(path_true, path_permuted)
plot_feature_importance_class(path_true, 20)
plot_features_tests_class(data_path, path_true, top_n=20, labels)
plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)
plot_feature_importance_class_plasma(path_true, top_n=20)

## Urine metabolites
path_true <- 'CKD_urine_kreat/output_XGB_class_nonD_CKD_urinekreat_2023_11_28__23-15-36'
path_permuted <- 'CKD_urine_kreat/output_XGB_class_nonD_CKD_urinekreat_2023_11_29__00-22-34_PERMUTED'
data_path <- 'CKD_urine_kreat/input_data'
labels <- c("CKD", "No CKD")

compared_to_permuted_class(path_true, path_permuted)
plot_feature_importance_class(path_true, 20)
plot_features_tests_class(data_path, path_true, top_n=20, labels)
plot_features_tests_top(data_path, path_true, top_n=20, nrow=4, labels)
plot_feature_importance_class_urine(path_true, top_n=20)

