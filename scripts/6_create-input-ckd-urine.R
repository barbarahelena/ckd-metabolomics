## Create XGB input files for (non-diabetic) CKD vs control model urine
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
rm(list=ls())
library(dplyr)

## Functions to make tab-delim data files for ML models
write_y <- function(x, name_y, data_path){
    if(missing(name_y)){
        cat('\n\nYou need to provide a name for the y data file!\n')
    }
    if(!name_y %in% c('y_binary.txt', 'y_reg.txt')){
        cat('\nThe file name is not compatible with XGBeast!\n' )
    }
    if(any(is.na(x))){
        cat('\nThere are missing values in the outcome data!\n')
    }
    data_path <- file.path(data_path, 'input_data')
    if(!exists(data_path)) dir.create(data_path)
    write.table(x, file = file.path(data_path, name_y), row.names = F, col.names = F, sep = '\t', quote = F)
}

write_data <- function(x, data_path){
    x <- as.matrix(x)
    if(any(is.na(x))){
        cat('There are missing values in the input data!\n')
    }
    data_path <- file.path(data_path, 'input_data')
    if(!exists(path)) dir.create(data_path)
    write.table(x, file.path(data_path, 'X_data.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
    write.table(colnames(x), file.path(data_path,'feat_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
    write.table(rownames(x), file.path(data_path,'subject_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
}


## Open datasets
umet <- readRDS("data/urine_metabolites_kreat.RDS") # urine metabolites adjusted for kreatinine
df <- readRDS('data/helius_ckd_xgb.RDS') # IDs and CKD yes/no
umet <- umet[which(rownames(umet) %in% df$ID),] # filter out DKD

df <- df[match(rownames(umet), df$ID), ] # put IDs of CKD dataframe in same order as metdf
rownames(umet)
df$ID
all(df$ID == rownames(umet)) # TRUE

## Make input data for classification model using functions
path <- 'CKD_urine_kreat'
if(!exists(path)) dir.create(path)
write_data(umet, path)
y <- as.data.frame(df$CKD)
y
write_y(y, name_y = 'y_binary.txt', path)
