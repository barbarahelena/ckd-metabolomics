# Create XGB input files for CKD vs healthy models

library(dplyr)
rm(list=ls())
setwd("~/Documents/VUmc/machine-learning")

# Functions to make tab-delim data files for ML models

# write y / predicted outcome
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
  dir.create(data_path)
  write.table(x, file = file.path(data_path, name_y), row.names = F, col.names = F, sep = '\t', quote = F)
}

# write x / determinants
write_data <- function(x, data_path){
  x <- as.matrix(x)
  if(any(is.na(x))){
    cat('There are missing values in the input data!\n')
  }
  data_path <- file.path(data_path, 'input_data')
  dir.create(data_path)
  write.table(x, file.path(data_path, 'X_data.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
  write.table(colnames(x), file.path(data_path,'feat_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
  write.table(rownames(x), file.path(data_path,'subject_ids.txt'), row.names = F, col.names = F, sep = '\t', quote = F)
}

# Create dataset metabolites without xenometabolites
## Open metabolite data
dd <- readxl::read_xlsx('data/HELIUS_EDTA_plasma_metabolites.xlsx')
met <- dd$Metabolite
dd$Metabolite <- NULL
dd2 <- as.data.frame(t(as.matrix(dd)))
dim(dd)
colnames(dd2)
colnames(dd2) <- met
rownames(dd2) <- paste0('S', rownames(dd2))

## Open infofile with categories of metabolites
infomet <- readxl::read_xlsx('data/Info_plasma_metabolites_b.xlsx')
colnames(infomet)
infomet <- infomet %>% select(metabolite = `BIOCHEMICAL`, sup = `SUPER PATHWAY`, sub = `SUB PATHWAY`)
infomet$sup <- as.factor(infomet$sup)
summary(infomet$sup)

## Select non-xenobiotics
noxeno <- infomet$metabolite[which(infomet$sup!='Xenobiotics')]
noxeno[1:5]
colnames(dd2)
ddd <- dd2[,which(colnames(dd2) %in% noxeno)]
dim(ddd)


# Create files for CKD main model
## Open HELIUS BP dataframe with BP and residuals data
setwd("~/Documents/VUmc/machine-learning")
d <- readRDS('data/helius_ckd.RDS')
head(d)
d$ID<-paste0('S', d$ID)
d <- d %>% select(ID, CKD)

## Put d and dd in same sequence of IDs
d <- d[match(rownames(ddd), d$ID), ]

## Check that outcome subject ids match metabolite subjects ids
all(d$ID == rownames(ddd)) # TRUE
d$ID
rownames(ddd)

## Make input data for classification model using functions
setwd('CKD-metabolites/helius-metabolomics-ckd/')
path <- 'CKD'
dir.create(path)
write_data(ddd, path)
y <- as.data.frame(d$CKD)
y
write_y(y, name_y = 'y_binary.txt', path)



# Make data files sensitivity model without DM
## Open HELIUS BP dataframe with BP and residuals data
setwd("~/Documents/VUmc/machine-learning")
d <- readRDS('data/helius_ckd.RDS')
head(d)
d$ID<-paste0('S', d$ID)
d <- d %>% filter(DM==0) %>% select(ID, CKD)
dim(d)

# Put d and dd in same sequence of IDs
rownames(ddd)
metdm <- ddd %>% filter(rownames(ddd) %in% d$ID)
d <- d[match(rownames(metdm), d$ID), ]

# check that outcome subject ids match metabolite subjects ids
all(d$ID == rownames(metdm)) # TRUE
d$ID
rownames(metdm)

# make input data for classification model 
setwd('CKD-metabolites/helius-metabolomics-ckd/')
path <- 'CKD/WithoutDM'
dir.create(path)
write_data(metdm, path)
y <- as.data.frame(d$CKD)
y
write_y(y, name_y = 'y_binary.txt', path)



# Make data files models CKD+T2DM and -T2DM
## Open HELIUS BP dataframe with BP and residuals data
setwd("~/Documents/VUmc/machine-learning")
d <- readRDS('data/helius_ckddm.RDS')
head(d)
table(d$group)
d$ID<-paste0('S', d$ID)
dim(d)

# Put d and dd in same sequence of IDs
rownames(ddd)
metckd <- ddd %>% filter(rownames(ddd) %in% d$ID)
d <- d[match(rownames(metckd), d$ID), ]

# check that outcome subject ids match metabolite subjects ids
all(d$ID == rownames(metckd)) # TRUE
d$ID
rownames(metckd)

# make input data for classification model 
setwd('CKD-metabolites/helius-metabolomics-ckd/')
path <- 'CKDDM'
dir.create(path)
write_data(metckd, path)
y <- as.data.frame(d$group)
y
table(y$`d$group`)
write_y(y, name_y = 'y_binary.txt', path)
