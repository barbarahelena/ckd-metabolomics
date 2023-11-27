### Table 1 
### Barbara Verhaar - b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(tableone)

## Open datasets
helius <- readRDS("data/helius_complete.RDS")
helius_ckd <- readRDS("data/helius_ckd.RDS")
helius_dkd <- readRDS("data/helius_dkd.RDS")

## Path for results
path <- "results"
if(!dir.exists(path)){dir.create(path)}
path <- "results/tables"
if(!dir.exists(path)){dir.create(path)}

## Table 1 complete
table1 <- helius %>% 
    select(Age, Age_cat, Sex, Ethnicity, BMI, CurrSmoking, CVD, DM, HT, 
           AntiHT, SBP, DBP, MDRD, ACR_KDIGO, HbA1C, LDL, Trig, FramRisk, group) %>% 
    CreateTableOne(data=., strata = 'group', addOverall = TRUE, test = TRUE) %>% 
    print(nonnormal=c("FramRisk", "Trig"))
openxlsx::write.xlsx(c(rownames(table1),table1), file.path(path, "table1_complete.xlsx"))

## Table 1 CKD
table1 <- helius_ckd %>% 
    select(Age, Age_cat, Sex, Ethnicity, BMI, CurrSmoking, CVD, DM, HT, 
           AntiHT, SBP, DBP, MDRD, ACR_KDIGO, HbA1C, LDL, Trig, FramRisk, CKD_group) %>% 
    CreateTableOne(data=., strata = 'CKD_group', addOverall = TRUE, test = TRUE) %>% 
    print(nonnormal=c("FramRisk", "Trig"))
openxlsx::write.xlsx(c(rownames(table1),table1), file.path(path, "table1_ckd.xlsx"))

## Table 1 DKD
table1 <- helius_dkd %>% 
    select(Age, Age_cat, Sex, Ethnicity, BMI, CurrSmoking, CVD, DM, HT, 
           AntiHT, SBP, DBP, MDRD, ACR_KDIGO, HbA1C, LDL, Trig, FramRisk, CKD_group) %>% 
    CreateTableOne(data=., strata = 'CKD_group', addOverall = TRUE, test = TRUE) %>% 
    print(nonnormal=c("FramRisk", "Trig"))
openxlsx::write.xlsx(c(rownames(table1),table1), file.path(path, "table1_dkd.xlsx"))
