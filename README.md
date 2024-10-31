HELIUS: Metabolome profiles in chronic kidney disease
--------------------------------------------------------
Author: Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Summary


## File structure
The scripts in this repo include the following:  


## Descriptives
CKD subjects and controls were selected from the HELIUS cohort based on the following criteria:  
- ACR-KDIGO albuminuria stage A1 (controls) or A2 (CKD)  
- Availability of plasma, urine and fecal samples  

## Machine learning analyses
For the machine learning analyses, a XGBoost algorithm was used in a nested cross-validated design with 100 iterations. The script for these analyses can be found in the `XGBoost.py` file. For each model, the study population was divided in a training (80%) and test (20%) set. For each iteration, the model was trained on the training set on which the hyperparameters were optimized in a 5-fold CV. The resulting model was tested on the test set. The parameter grid can be found in `param_grid_medium.json` in the scripts/xgboost folder. Each model resulted in a ROC and feature importance of predictors, averaged over 200 iterations.

In the `process_model_results_class.R` file, the functions from the `functions.R` file are used to generate feature importance plots, and violin plots for best predictors to assess differences in concentrations between groups.

The following models were used:  
- Plasma metabolite profiles non-diabetic CKD versus healthy controls: AUC 0.59
- Urine metabolite profiles non-diabetic CKD versus healthy controls: AUC 0.65

### Plasma metabolites
![ROC CKD](results/ROC_CKD.png){width=400px}

![Violin plots best predictors CKD][ckd_model]

### Urine metabolites
![ROC CKD](results/ROC_CKDsens.png){width=400px}

![Violin plots best predictors sensitivity model][ckdsens_model]


## Regression analyses with highest ranked predictors



## Overlap plasma and urine metabolites


[ckd_model]: results/ckd_boxplots.svg "Best predictors CKD"
[ckdsens_model]: results/ckdsens_boxplots.svg "Best predictors CKD without T2DM"
[ckddm_model]: results/ckddm_boxplots.svg "Best predictors CKD+T2DM"
[corr_met]: results/correlations/ckd_cor.svg "Correlations CKD predictors"

