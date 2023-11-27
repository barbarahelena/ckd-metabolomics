---
title: "README"
author: "Barbara Verhaar"
date: "12/21/2020"
output:
  html_document:
    df_print: paged
always_allow_html: yes
---

HELIUS: Metabolome profiles in chronic kidney disease
--------------------------------------------------------
Author: Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Summary
These analyses aimed to assess differences in plasma metabolites in chronic kidney disease (CKD) versus healthy controls; and differences between subjects with CKD and diabetes type 2 (T2DM) and subjects with CKD without T2DM.  

In summary, machine learning classification models were used to find the best predicting metabolites for CKD. The AUC for the model to discriminate between CKD and healthy controls was 0.69. The AUC to discriminate CKD+T2DM and CKD-T2DM was 0.80. 

Regression models were used to assess associations between CKD or DM and these best predictors while adjusting for most important confounders. Linear models appeared not to be very useful, the best predictors were not associated with albuminuria in a linear model. Logistic regression for CKD versus controls (and for the second question, T2DM vs no T2DM) showed that many of the associations remained after adjusting for important confounders.

## File structure
The scripts in this repo include the following:  
- `XGBeast.py`: python script for XGBoost machine learning model  
- `barbara.json`: parameter grid  
- `functions.R`: functions for processing results  
- `process_model_results_class.R`: generating plots of model results  
- `corrplots.R`: correlation and network plots showing correlations between best metabolites  
- `glmmodels.R`: linear models (outcome: albuminuria) and logistic regression models (outcome: CKD or DM)  

## Descriptives
CKD subjects and controls were selected from the HELIUS cohort based on the following criteria:  
- ACR-KDIGO albuminuria stage A1 (controls) or A2 (CKD)  
- Availability of plasma, urine and fecal samples  

The descriptive characteristics of these subjects can be found in the table below.
```{r, echo=FALSE, warning=FALSE, message=FALSE}
library(dplyr)
library(tableone)
library(kableExtra)
helius <- readRDS("data/heliusdf.RDS")
table1 <- helius %>% 
    select(Age, Age_cat, Sex, Ethnicity, BMI, CurrSmoking, CVD, DM, HT, 
           AntiHT, SBP, DBP, MDRD, ACR_KDIGO, HbA1C, LDL, Trig, FramRisk, group) %>% 
    CreateTableOne(data=., strata = 'group', addOverall = TRUE, test = FALSE)
kableone(table1) %>% kable_paper()

```

## Machine learning analyses
For the machine learning analyses, a XGBoost algorithm was used in a nested cross-validated design with 100 iterations. The script for these analyses can be found in the `XGBoost.py` file. For each model, the study population was divided in a training (80%) and test (20%) set. For each iteration, the model was trained on the training set on which the hyperparameters were optimized in a 5-fold CV. The resulting model was tested on the test set. The parameter grid can be found in `barbara.json`. Each model resulted in a ROC and feature importance of predictors, averaged over 100 iterations.

In the `process_model_results_class.R` file, the functions from the `functions.R` file are used to generate feature importance plots, and violin plots for best predictors to assess differences in concentrations between groups.

The following models were used:  
- CKD versus healthy controls: AUC 0.69  
- Sensitivity analyses without 45 T2DM subjects: AUC 0.59  
- CKD+T2DM versus CKD-T2DM: AUC 0.80  

### CKD model
![ROC CKD](results/ROC_CKD.png){width=400px}

![Violin plots best predictors CKD][ckd_model]

### CKD versus controls - no T2DM subjects
![ROC CKD](results/ROC_CKDsens.png){width=400px}

![Violin plots best predictors sensitivity model][ckdsens_model]

### CKD+T2DM versus CKD-T2DM
![ROC CKD](results/ROC_CKDDM.png){width=400px}

![Violin plots best predictors CKD+T2DM versus CKD-T2DM][ckddm_model]


## Correlations between best-predicting metabolites
We could also look at the correlations between the best predicting metabolites, to see if there are clusters.

### Correlations between metabolites model CKD versus controls
There are three groups that stand out, looking at the plot below:  
1. Formylated/acetylated amino acids: N-formylmethionine, N-acetylvaline and N-acetylalanine  
2. Amino acids with a carboxyethyl-group  
3. Oxalate and glycerate  

![Correlation plot main model][corr_met]

### Correlations between metabolites model CKD+T2DM versus CKD-T2DM
There are two main clusters in the plot below:  
1. Glucose and carboxy-amino acid  
2. Glycerophospholipid  
3. Acyl cholines  

![Correlation plot T2DM vs no T2DM model][corr_metdm]

### Network plot from CKD-controls model
In the network plot below, the same correlations between best predicting metabolites are illustrated with the sub pathways that they are part of. There are many different sub-pathways and not very clear clusters. However, glucose has a central position in this plot.

![Network plot][network_dm]

## Logistic regression
With logistic regression, we can assess the association between CKD (binary outcome) and the metabolites, while adjusting for confounders. The odds ratios (OR) represent the odds of having CKD with 1 SD increase in the metabolite.

![Logistic regression][ckd_logreg]

## Linear models
There are also scripts to generate linear models with the best predictors and the outcome albuminuria (log-transformed)  (see `glmmodels.R`). These linear models are only performed in the CKD-groups, since albuminuria is almost absent in healthy controls. These linear models showed almost no significant associations.

![Linear regression][ckd_linreg]

## Conclusions
The AUC for the main model (predicting CKD versus no CKD) was reasonable, but got worse in the sensitivity analysis without T2DM subjects. The AUC for discriminating CKD+T2DM from CKD-T2DM was best: 0.80. 

In the models that included T2DM subjects, glucose, other saccharides and glycosylated products have a central role. However, there are other metabolites among the best predictors that do not cluster with glucose.


[ckd_model]: results/ckd_boxplots.svg "Best predictors CKD"
[ckdsens_model]: results/ckdsens_boxplots.svg "Best predictors CKD without T2DM"
[ckddm_model]: results/ckddm_boxplots.svg "Best predictors CKD+T2DM"
[corr_met]: results/correlations/ckd_cor.svg "Correlations CKD predictors"
[corr_metdm]: results/correlations/ckddm_cor.svg "Correlations T2DM in CKD predictors"
[network_dm]: results/correlations/ckd_network.svg "Network plot T2DM in CKD predictors"
[ckd_logreg]: results/logreg/ckd_logreg.svg "Logistic regression CKD"
[ckd_linreg]: results/linreg/ckd_linreg.svg "Linear regression CKD"
