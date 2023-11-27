### Data cleaning CKD project
### Barbara Verhaar - b.j.verhaar@amsterdamumc.nl

## Libraries
packages <- c("rio", "haven", "dplyr", "tableone", "tidyverse", "kableExtra", 
              "naniar", "ggpubr", "ggsci", "gridExtra", "clipr")
pacman::p_load(packages, character.only = TRUE)

## Functions
theme_Publication <- function(base_size=12, base_family="sans") {
    library(grid)
    library(ggthemes)
    (theme_foundation(base_size=base_size, base_family=base_family)
        + theme(plot.title = element_text(face = "bold", family = 'Helvetica',
                                          size = rel(1.0), hjust = 0.5),
                text = element_text(family = 'Helvetica'),
                panel.background = element_rect(colour = NA),
                plot.background = element_rect(colour = NA),
                panel.border = element_rect(colour = NA),
                axis.title = element_text(face = "bold",size = rel(1)),
                axis.title.y = element_text(angle=90,vjust =2),
                axis.title.x = element_text(vjust = -0.2),
                axis.text = element_text(), 
                axis.line.x = element_line(colour="black"),
                axis.ticks.y = element_blank(), axis.line.y = element_blank(), 
                axis.ticks = element_line(),
                panel.grid.major = element_line(colour="#f0f0f0"),
                panel.grid.minor = element_blank(),
                legend.key = element_rect(colour = NA),
                legend.position = "bottom",
                # legend.direction = "horizontal",
                legend.key.size= unit(0.2, "cm"),
                legend.spacing  = unit(0, "cm"),
                # legend.title = element_text(face="italic"),
                plot.margin=unit(c(10,5,5,5),"mm"),
                strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                strip.text = element_text(face="bold")
        ))
    
} 

## Open datasets
dataDir <- 'data' 
infomet <- rio::import(file.path(dataDir,'Info_plasma_metabolites_b.xlsx'))
metabolites <- rio::import(file.path(dataDir,'HELIUS_EDTA_plasma_metabolites.xlsx'))
infoumet <- rio::import(file.path(dataDir,'Info_urine_metabolites.xlsx'))
urinemet <- rio::import(file.path(dataDir,'HELIUS_Urine_metabolites.xlsx'))
heliusData <- read_sav(file.path(dataDir, "201104_HELIUS data Barbara Verhaar.sav"))

## Path for results
path <- "results"
if(!dir.exists(path)){dir.create(path)}
path <- "results/descriptives"
if(!dir.exists(path)){dir.create(path)}

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

## Table 1
table1 <- helius %>% 
    select(Age, Age_cat, Sex, Ethnicity, BMI, CurrSmoking, CVD, DM, HT, 
           AntiHT, SBP, DBP, MDRD, ACR_KDIGO, HbA1C, LDL, Trig, FramRisk, group) %>% 
    CreateTableOne(data=., strata = 'group', addOverall = TRUE, test = FALSE) %>% 
    print(nonnormal=c("FramRisk", "Trig"))
#kableone(table1) %>%  kable_classic()
#write_clip(kableone(table1))
openxlsx::write.xlsx(c(rownames(table1),table1), file.path(path, "table1.xlsx"))


## PCA plots of metabolites
dim(metabolites) # 890 - 370
mat <- metabolites[,-1]
rownames(mat) <- metabolites[,1]

infomet$sup <- as.factor(infomet$sup)
summary(infomet$sup)
noxeno <- infomet$met[which(infomet$sup!='Xenobiotics')]
mat2 <- mat[which(rownames(mat) %in% noxeno),]

pca <- prcomp(t(mat2), center = T, scale = T)
sum <- summary(pca)

df <- as.data.frame(pca$x[, 1:2])
head(df)
df$ID <- rownames(df)
df$CKD <- helius$CKD_group[match(df$ID, helius$Heliusnr)]
df$CKDT2DM <- helius$group[match(df$ID, helius$Heliusnr)]
df$Sex <- helius$Sex[match(df$ID, helius$Heliusnr)]
df$HT <- helius$HT[match(df$ID, helius$Heliusnr)]
df$DM <- helius$DM[match(df$ID, helius$Heliusnr)]
df$Age_cat <- helius$Age_cat[match(df$ID, helius$Heliusnr)]
df$Ethnicity <- helius$Ethnicity[match(df$ID, helius$Heliusnr)]
df$BMI_cat <- helius$BMI_cat[match(df$ID, helius$Heliusnr)]

plA <- ggplot(df, aes(x=PC1, y=PC2, color=Age_cat)) +
    geom_point()+
    ggtitle('Age') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x=str_c('PC1 ', round(sum$importance[2,1]*100, 1), '%'), 
         y=str_c('PC2 ', round(sum$importance[2,2]*100,1), '%')) +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank())

plB <- ggplot(df, aes(x=PC1, y=PC2, color=Sex)) +
    geom_point() +
    ggtitle('Sex') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x=str_c('PC1 ', round(sum$importance[2,1]*100, 1), '%'), 
         y=str_c('PC2 ', round(sum$importance[2,2]*100,1), '%')) +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank())

plC <- ggplot(df, aes(x=PC1, y=PC2, color=Ethnicity)) +
    geom_point() +
    ggtitle('Ethnicity') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x=str_c('PC1 ', round(sum$importance[2,1]*100, 1), '%'), 
         y=str_c('PC2 ', round(sum$importance[2,2]*100,1), '%')) +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank())

plD <- ggplot(df, aes(x=PC1, y=PC2, color=BMI_cat)) +
    geom_point()+
    ggtitle('BMI') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x=str_c('PC1 ', round(sum$importance[2,1]*100, 1), '%'), 
         y=str_c('PC2 ', round(sum$importance[2,2]*100,1), '%')) +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank())

pdf(file.path(path,"PCAs_subgroups.pdf"), width = 10, height = 10)
grid.arrange(plA, plB, plC, plD, nrow=2, top=textGrob("PCA metabolites for subgroups", gp=gpar(font=2)))
dev.off()
# svg(file.path(path,"PCAs_subgroups.svg"), width = 10, height = 10)
# grid.arrange(plA, plB, plC, plD, nrow=2, top=textGrob("PCA metabolites for subgroups", gp=gpar(font=2)))
# dev.off()

plE <- ggplot(df, aes(x=PC1, y=PC2, color=CKDT2DM)) +
    geom_point()+
    scale_color_lancet() +
    ggtitle('PCA metabolites - CKD groups') +
    theme_Publication() +
    stat_ellipse() +
    labs(x=str_c('PC1 ', round(sum$importance[2,1]*100, 1), '%'), 
         y=str_c('PC2 ', round(sum$importance[2,2]*100,1), '%')) +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank())

plF <- ggplot(df, aes(x=PC1, y=PC2, color=CKD)) +
    geom_point()+
    ggtitle('PCA metabolites - CKD') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x=str_c('PC1 ', round(sum$importance[2,1]*100, 1), '%'), 
         y=str_c('PC2 ', round(sum$importance[2,2]*100,1), '%')) +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank())

pdf(file.path(path,"PCAs_ckdgroups.pdf"), width = 12, height = 6)
grid.arrange(plE, plF, nrow=1, top=textGrob("PCA metabolites for CKD groups", gp=gpar(font=2)))
dev.off()
# svg(file.path(path,"PCAs_ckdgroups.svg"), width = 12, height = 6)
# grid.arrange(plE, plF, nrow=1, top=textGrob("PCA metabolites for CKD groups", gp=gpar(font=2)))
# dev.off()


## PCA urine metabolites
dim(urinemet) # 1075 - 368
mat <- urinemet[,-1]
rownames(mat) <- urinemet[,1]

u_infomet$sup <- as.factor(u_infomet$sup)
summary(u_infomet$sup)
noxeno <- u_infomet$met[which(u_infomet$sup!='Xenobiotics')]
mat2 <- mat[which(rownames(mat) %in% noxeno),]

pca <- prcomp(t(mat2), center = T, scale = T)
sum <- summary(pca)

df <- as.data.frame(pca$x[, 1:2])
head(df)
df$ID <- rownames(df)
df$CKD <- helius$CKD_group[match(df$ID, helius$Heliusnr)]
df$CKDT2DM <- helius$group[match(df$ID, helius$Heliusnr)]
df$Sex <- helius$Sex[match(df$ID, helius$Heliusnr)]
df$HT <- helius$HT[match(df$ID, helius$Heliusnr)]
df$DM <- helius$DM[match(df$ID, helius$Heliusnr)]
df$Age_cat <- helius$Age_cat[match(df$ID, helius$Heliusnr)]
df$Ethnicity <- helius$Ethnicity[match(df$ID, helius$Heliusnr)]
df$BMI_cat <- helius$BMI_cat[match(df$ID, helius$Heliusnr)]

plA <- ggplot(df, aes(x=PC1, y=PC2, color=Age_cat)) +
    geom_point()+
    ggtitle('Age') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x=str_c('PC1 ', round(sum$importance[2,1]*100, 1), '%'), 
         y=str_c('PC2 ', round(sum$importance[2,2]*100,1), '%')) +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank())

plB <- ggplot(df, aes(x=PC1, y=PC2, color=Sex)) +
    geom_point() +
    ggtitle('Sex') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x=str_c('PC1 ', round(sum$importance[2,1]*100, 1), '%'), 
         y=str_c('PC2 ', round(sum$importance[2,2]*100,1), '%')) +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank())

plC <- ggplot(df, aes(x=PC1, y=PC2, color=Ethnicity)) +
    geom_point() +
    ggtitle('Ethnicity') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x=str_c('PC1 ', round(sum$importance[2,1]*100, 1), '%'), 
         y=str_c('PC2 ', round(sum$importance[2,2]*100,1), '%')) +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank())

plD <- ggplot(df, aes(x=PC1, y=PC2, color=BMI_cat)) +
    geom_point()+
    ggtitle('BMI') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x=str_c('PC1 ', round(sum$importance[2,1]*100, 1), '%'), 
         y=str_c('PC2 ', round(sum$importance[2,2]*100,1), '%')) +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank(), legend.title = element_blank())

pdf(file.path(path,"Urine_PCAs_subgroups.pdf"), width = 10, height = 10)
grid.arrange(plA, plB, plC, plD, nrow=2, top=textGrob("PCA urine metabolites for subgroups", gp=gpar(font=2)))
dev.off()
# svg(file.path(path,"PCAs_subgroups.svg"), width = 10, height = 10)
# grid.arrange(plA, plB, plC, plD, nrow=2, top=textGrob("PCA metabolites for subgroups", gp=gpar(font=2)))
# dev.off()

plE <- ggplot(df, aes(x=PC1, y=PC2, color=CKDT2DM)) +
    geom_point()+
    scale_color_lancet() +
    ggtitle('PCA metabolites - CKD groups') +
    theme_Publication() +
    stat_ellipse() +
    labs(x=str_c('PC1 ', round(sum$importance[2,1]*100, 1), '%'), 
         y=str_c('PC2 ', round(sum$importance[2,2]*100,1), '%')) +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank())

plF <- ggplot(df, aes(x=PC1, y=PC2, color=CKD)) +
    geom_point()+
    ggtitle('PCA metabolites - CKD') +
    scale_color_lancet() +
    theme_minimal() +
    stat_ellipse() +
    labs(x=str_c('PC1 ', round(sum$importance[2,1]*100, 1), '%'), 
         y=str_c('PC2 ', round(sum$importance[2,2]*100,1), '%')) +
    theme_Publication() +
    theme(axis.line.x = element_blank(), axis.ticks.x=element_blank())

pdf(file.path(path,"Urine_PCAs_ckdgroups.pdf"), width = 12, height = 6)
grid.arrange(plE, plF, nrow=1, top=textGrob("PCA urine metabolites for CKD groups", gp=gpar(font=2)))
dev.off()
# svg(file.path(path,"PCAs_ckdgroups.svg"), width = 12, height = 6)
# grid.arrange(plE, plF, nrow=1, top=textGrob("PCA metabolites for CKD groups", gp=gpar(font=2)))
# dev.off()


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
