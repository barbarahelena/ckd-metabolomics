#!/bin/bash
#SBATCH -c 24
#SBATCH --mem=16G
#SBATCH --time=12:00:00
#SBATCH -p 'genoa'
eval "$(conda shell.bash hook)"
conda activate xgb
python /projects/0/prjs0784/ckd-metabolomics/scripts/xgboost/XGBeast.py \
    -name nonD_CKD_plasma \
    -path /projects/0/prjs0784/ckd-metabolomics/CKD_plasma \
    -x class \
    -n 200 \
    -rand_seed 1234 \
    -t 23 \
    -param /projects/0/prjs0784/ckd-metabolomics/scripts/xgboost/param_medium.json 
python /projects/0/prjs0784/ckd-metabolomics/scripts/xgboost/XGBeast.py \
    -name nonD_CKD_plasma \
    -path /projects/0/prjs0784/ckd-metabolomics/CKD_plasma \
    -x class \
    -n 200 \
    -rand_seed 1234 \
    -t 23 \
    -permute \
    -param /projects/0/prjs0784/ckd-metabolomics/scripts/xgboost/param_medium.json 