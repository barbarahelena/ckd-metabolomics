### Data cleaning CKD project
### Barbara Verhaar - b.j.verhaar@amsterdamumc.nl

## Libraries
packages <- c("rio", "haven", "dplyr", "tableone", "tidyverse", "kableExtra", 
              "naniar", "ggpubr", "ggsci", "gridExtra", "clipr")
pacman::p_load(packages, character.only = TRUE)

## Open datasets
dataDir <- 'data' 
infomet <- rio::import(file.path(dataDir,'Info_plasma_metabolites_b.xlsx'))
metabolites <- rio::import(file.path(dataDir,'HELIUS_EDTA_plasma_metabolites.xlsx'))
infoumet <- rio::import(file.path(dataDir,'Info_urine_metabolites.xlsx'))
urinemet <- rio::import(file.path(dataDir,'HELIUS_Urine_metabolites.xlsx'))
heliusData <- read_sav(file.path(dataDir, "201104_HELIUS data Barbara Verhaar.sav"))

## Data cleaning
colnames(heliusData)
str(heliusData$H1_Lab_UitslagMIAL)
helius <- heliusData %>% 
    select(Heliusnr, Age=H1_lft, Sex=H1_geslacht, Ethnicity=H1_EtnTotaal, 
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
        group = case_when(
            ACR_KDIGO==1 ~ 'Controls',
            DM==1 & ACR_KDIGO==2 ~ 'CKD+T2DM',
            DM==0 & ACR_KDIGO==2 ~ 'CKD-T2DM'
        ),
        group = factor(group, levels = c('Controls', 'CKD+T2DM', 'CKD-T2DM')),
        CKD_group = case_when(
            group =='Controls' ~ 'Controls',
            group == 'CKD+T2DM' ~ 'CKD',
            group == 'CKD-T2DM' ~ 'CKD'
        ),
        CKD_group = factor(CKD_group, levels = c('Controls', 'CKD')),
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
    mutate_if(is.numeric, as.numeric) %>% # strange thing to do - reason is the spss source file
    droplevels()

dim(helius)
colnames(helius)
helius %>% group_by(CKD_group) %>% 
    summarise(mean(Alb))

saveRDS(helius, "data/heliusdf.RDS")

infomet <- infomet %>% select(met=BIOCHEMICAL, sup=`SUPER PATHWAY`, sub=`SUB PATHWAY`)
saveRDS(infomet, "data/infomet.RDS")

u_infomet <- infoumet %>% select(met=BIOCHEMICAL, sup=`SUPER PATHWAY`, sub=`SUB PATHWAY`)
saveRDS(u_infomet, "data/infomet_urine.RDS")

noxeno <- infomet$met[which(infomet$sup!='Xenobiotics')]
mets <- metabolites[which(metabolites$Metabolite %in% noxeno),]

u_noxeno <- u_infomet$met[which(u_infomet$sup!='Xenobiotics')]
u_mets <- umetabolites[which(umetabolites$Metabolite %in% u_noxeno),]

summary(mets$Metabolite %in% u_mets$Metabolite)

## Save RDS with CKD-group and ID
table(helius$CKD_group) ## 200 controls and 169 CKD
ckd <- helius %>% select(ID=Heliusnr, CKD_group, DM) %>% 
    mutate(CKD = ifelse(CKD_group=="CKD", 1, 0),
           DM = ifelse(DM=="Diabetes", 1, 0)) %>% 
    select(ID, CKD, DM)
table(ckd$CKD) ## 200 controls and 169 CKD
table(ckd$DM) ## 324 controls and 45 DM
head(ckd)
saveRDS(ckd, 'data/helius_ckd.RDS')

## Save RDS CKD- vs CKD-T2DM
table(helius$group) ## 200 - 45 - 124
levels(helius$group)
ckddm <- helius %>% select(ID=Heliusnr, group) %>% 
    filter(group!="Controls") %>% 
    mutate(group = ifelse(group=="CKD+T2DM", 1, 0))
table(ckddm$group) ## 124 vs 45
head(ckddm)
saveRDS(ckddm, 'data/helius_ckddm.RDS')
