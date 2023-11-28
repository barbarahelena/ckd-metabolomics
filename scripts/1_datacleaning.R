### Data cleaning CKD project
### Barbara Verhaar - b.j.verhaar@amsterdamumc.nl

## Libraries
library(haven)
library(rio)
library(tidyverse)

## Open datasets
heliusData <- read_sav('data/201104_HELIUS data Barbara Verhaar.sav')
metabolites <- rio::import('data/HELIUS_EDTA_plasma_metabolites.xlsx')
urinemet <- rio::import('data/HELIUS_Urine_metabolites.xlsx')

infomet <- rio::import('data/Info_plasma_metabolites_b.xlsx')
infoumet <- rio::import('data/Info_urine_metabolites.xlsx')

## Data cleaning
colnames(heliusData)
str(heliusData$H1_Lab_UitslagMIAL)
helius <- heliusData %>% 
    select(ID=Heliusnr, Age=H1_lft, Sex=H1_geslacht, Ethnicity=H1_EtnTotaal, 
           Smoking=H1_Roken, BMI=H1_LO_BMI, WHR=H1_LO_WHR, SBP=H1_LO_GemBPSysZit, 
           DBP=H1_LO_GemBPDiaZit, HT=H1_HT_SelfBPMed, DM=H1_Diabetes_SelfGlucMed,
           DMMed=H1_Diabetesmiddelen, Claudicatio=H1_CI_Rose, PossInf=H1_possINF_Rose, 
           AP=H1_AP_Rose, CVD=H1_CVD_Rose, AntiHT=H1_Antihypertensiva, 
           MDRD=H1_MDRD_eGFR, CKDEPI=H1_CKDEPI_eGFR, CKDStage=H1_CKDEPI_stage, 
           MetSyn=H1_MetSyn_MetabolicSyndrome, LDL=H1_Lab_uitslagRLDL, Trig=H1_Lab_UitslagTRIG,
           FramRisk=H1_Fram_CVD, SCORENL=H1_SCORE_CVDmort_NL,
           ACR_KDIGO=H1_ACR_KDIGO, Microalb=H1_Microalbuminurie, 
           HbA1C=H1_Lab_UitslagIH1C, Kreat=H1_Lab_UitslagKREA_HP, Alb=H1_Lab_UitslagMIAL) %>% 
    mutate(
        ID = paste0('S', ID),
        group = case_when(
            ACR_KDIGO==1 ~ 'Controls',
            DM==1 & ACR_KDIGO==2 ~ 'CKD+T2DM',
            DM==0 & ACR_KDIGO==2 ~ 'CKD-T2DM'
        ),
        group = factor(group, levels = c('Controls', 'CKD+T2DM', 'CKD-T2DM')),
        CKD_group = case_when(
            group =='Controls' ~ 'Controls',
            group == 'CKD+T2DM' ~ 'DKD',
            group == 'CKD-T2DM' ~ 'CKD'
        ),
        CKD_group = factor(CKD_group, levels = c('Controls', 'CKD', 'DKD')),
        Sex = as_factor(Sex, levels = c("labels"), ordered = FALSE),
        Sex = fct_recode(Sex, "Male" = "man", "Female" = "vrouw"),
        Ethnicity = as_factor(Ethnicity, levels = c("labels"), ordered = FALSE),
        Ethnicity = fct_recode(Ethnicity, "African Surinamese"="Creools", "Ghanaian"="Ghanees",
                               "South-Asian Surinamese"="Hind", "Dutch"="NL"),
        Smoking = as_factor(Smoking, levels = c("labels"), ordered = FALSE),
        CurrSmoking = fct_recode(Smoking, "Yes" = "Ja", "No" = "Nee, ik heb nooit gerookt", "No" = "Nee, maar vroeger wel"),
        CurrSmoking = fct_infreq(CurrSmoking),
        HT = as_factor(HT, levels = c("labels"), ordered = FALSE),
        DM = as_factor(DM, levels = c("labels"), ordered = FALSE),
        DM = fct_recode(DM, "Diabetes" = "Ja", "No diabetes" = "Nee"),
        DMMed = as_factor(DMMed, levels = c("labels"), ordered = FALSE),
        Claudicatio = as_factor(Claudicatio, levels = c("labels"), ordered = FALSE),
        PossInf = as_factor(PossInf, levels = c("labels"), ordered = FALSE),
        AP = as_factor(AP, levels = c("labels"), ordered = FALSE),
        CVD = as_factor(CVD, levels = c("labels"), ordered = FALSE),
        AntiHT = as_factor(AntiHT, levels = c("labels"), ordered = FALSE),
        CKDStage = as_factor(CKDStage, levels = c("labels"), ordered = FALSE),
        MetSyn=as_factor(MetSyn, levels = c("labels"), ordered = FALSE),
        ACR_KDIGO = as_factor(ACR_KDIGO, levels = c("labels"), ordered = FALSE),
        Microalb = as_factor(Microalb, levels = c("labels"), ordered = FALSE),
        EthnAfr = ifelse(Ethnicity == "Hind" | Ethnicity == "Dutch", "Non-African", "African"),
        Age_cat = ifelse(Age<=50, "Young", "Old"),
        BMI_cat = ifelse(BMI<25, "BMI<25", "BMI>=25")
    ) %>% 
    mutate_if(is.numeric, as.numeric) %>% # to make spss numerics, numerics again (strange thing to do)
    droplevels()

dim(helius)
colnames(helius)
helius %>% group_by(CKD_group) %>% 
    summarise(mean(Alb))

helius_ckd <- helius %>% filter(DM == "No diabetes") %>% droplevels(.)

helius_dkd <- helius %>% filter(!(DM == "No diabetes" & CKD_group == "CKD")) %>% 
    mutate(CKD_group = fct_recode(CKD_group, "DKD" = "CKD")) %>% 
    droplevels(.)
dim(helius_dkd)
table(helius_dkd$CKD_group)

## Save datasets
saveRDS(helius, "data/helius_complete.RDS")
saveRDS(helius_ckd, "data/helius_ckd.RDS")
saveRDS(helius_dkd, "data/helius_dkd.RDS")

## Save ML set with CKD-group and ID
table(helius_ckd$CKD_group) ## 200 controls, 124 CKD
ckd <- helius_ckd %>% select(ID, CKD_group) %>% 
    mutate(CKD = ifelse(CKD_group=="CKD", 1, 0)) 
table(ckd$CKD) ## 200 controls and 124 CKD
head(ckd)
saveRDS(ckd, 'data/helius_ckd_xgb.RDS')

## Create dataset metabolites without xenometabolites
infomet <- infomet %>% select(met=BIOCHEMICAL, sup=`SUPER PATHWAY`, sub=`SUB PATHWAY`)
noxeno <- infomet$met[which(infomet$sup!='Xenobiotics')]
saveRDS(noxeno, "data/infomet.RDS")

u_infomet <- infoumet %>% select(met=BIOCHEMICAL, sup=`SUPER PATHWAY`, sub=`SUB PATHWAY`)
u_noxeno <- u_infomet$met[which(u_infomet$sup!='Xenobiotics')]
saveRDS(u_noxeno, "data/infomet_urine.RDS")

mets <- metabolites[which(metabolites$Metabolite %in% noxeno),]
colnames(mets)[2:ncol(mets)] <- paste0('S', colnames(mets)[2:ncol(mets)]) # add letter to ID, safety reasons
rownames(mets) <- mets$Metabolite
mets$Metabolite <- NULL
mets <- t(mets)
saveRDS(mets, "data/plasma_metabolites.RDS")

u_mets <- urinemet[which(urinemet$Metabolite %in% u_noxeno),]
colnames(u_mets)[2:ncol(u_mets)] <- paste0('S', colnames(u_mets)[2:ncol(u_mets)]) # add letter to ID, safety reasons
rownames(u_mets) <- u_mets$Metabolite
u_mets$Metabolite <- NULL 
u_mets <- t(u_mets)
saveRDS(u_mets, "data/urine_metabolites.RDS")

summary(rownames(mets) %in% rownames(u_mets)) # overlap urine and plasma metabolites 422 of 722

u_mets_kreat <- u_mets / u_mets[,'creatinine']
all(u_mets[,1] / u_mets[,'creatinine'] == u_mets_kreat[,1])
saveRDS(u_mets, "data/urine_metabolites_kreat.RDS")

