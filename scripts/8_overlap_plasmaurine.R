## Linear models
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
packages <- c("rio", "haven", "dplyr", "tidyverse", "corrr", "corrplot", "Hmisc", "igraph", 
              "ggraph", "network", "naniar", "ggpubr", "ggsci", "gridExtra",
              "RColorBrewer", "plotly")
pacman::p_load(packages, character.only = TRUE)

splitgroup <- function(df, infodf, cat){
    group <- infodf %>% filter(sup == cat) %>% select(metabolite)
    df <- df %>% select(any_of(group$metabolite))
}

## Opening data
helius <- readRDS("data/heliusdf.RDS")
summary <- rio::import('data/Info_plasma_metabolites_b.xlsx')
infomet <- summary %>% select(metabolite=BIOCHEMICAL, sup=`SUPER PATHWAY`, sub=`SUB PATHWAY`)
infomet$sup <- as.factor(infomet$sup)
plasmamet <- rio::import('data/HELIUS_EDTA_plasma_metabolites.xlsx')
urinemet <- rio::import('data/HELIUS_urine_metabolites.csv')
rownames(plasmamet) <- plasmamet$Metabolite
plasmamet$Metabolite <- NULL
plasmamet2 <- as.data.frame(t(as.matrix(plasmamet)))

rownames(urinemet) <- urinemet$Heliusnr
urinemet[,1:5] <- NULL

u_sel <- urinemet[,which(names(urinemet) %in% names(plasmamet2))]
u_sel <- u_sel %>% rownames_to_column() %>% 
    mutate_at(c(1:557),as.numeric) %>% 
    column_to_rownames()
u_sel <- u_sel[order(row.names(u_sel)),]

p_sel <- plasmamet2[,which(names(plasmamet2) %in% names(urinemet))]
p_sel <- p_sel %>% filter(rownames(.) %in% rownames(u_sel)) %>% rownames_to_column() %>% 
    mutate_all(as.numeric) %>% 
    column_to_rownames()
p_sel <- p_sel[order(row.names(p_sel)),]

p <- lapply(levels(infomet$sup), splitgroup, df=p_sel, infodf=infomet)
names(p) <- levels(infomet$sup)

u <- lapply(levels(infomet$sup), splitgroup, df=u_sel, infodf=infomet)
names(u) <- levels(infomet$sup)

## Correlation principal components
noxeno <- infomet$met[which(infomet$sup!='Xenobiotics')]
pmat <- as.matrix(p_sel[,which(colnames(p_sel) %in% noxeno)])
ppca <- prcomp(pmat, center = T)
psum <- summary(ppca)
df1 <- as.data.frame(psum$x)[,1:5]
colnames(df1) <- str_c("P_", colnames(df1))

umat <- as.matrix(u_sel[,which(colnames(u_sel) %in% noxeno)])
upca <- prcomp(umat, center = T)
usum <- summary(upca)
df2 <- as.data.frame(usum$x)[,1:5]
colnames(df2) <- str_c("U_", colnames(df2))

df <- merge(df1, df2, by="row.names")
df$Row.names <- NULL

## Corrplot
cor <- cor(df1, df2, method="spearman")
corrplot(cor,  type = "upper", tl.col = "black",
         order = "original", tl.srt = 45, insig="blank", sig.level = 0.05, 
         #p.mat=cor$P,
         method="color", mar=c(0,0,1,0), addCoef.col = "black", diag = T, 
         addgrid.col = "grey")





