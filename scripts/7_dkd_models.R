## Regression models
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
library(tidyverse)
library(ggpubr)
library(ggsci)
library(RColorBrewer)
library(broom)
library(aplot)

## Functions
glm_prep <- function(df, metabolites, helius){
    df <- df %>% select(Metabolite=FeatName, RelFeatImp) %>% 
        arrange(-RelFeatImp) %>% slice(1:20)
    metabolites <- as.data.frame(t(metabolites)) %>% mutate(Metabolite = rownames(.))
    df_met <- metabolites %>% inner_join(., df, by='Metabolite') %>% 
        arrange(-RelFeatImp)
    rownames(df_met) <- df_met$Metabolite
    df_met$Metabolite <- NULL
    df_met$RelFeatImp <- NULL
    df_met <- as.data.frame(t(as.matrix(df_met)))
    df_met$ID <- rownames(df_met)
    helius_sub <- left_join(helius, df_met, by='ID')
    return(helius_sub)
}

log_group_dkd <- function(df, dfname, writetable = FALSE, figure = FALSE, plottitle){
    theme_Publication <- function(base_size=12, base_family="sans") {
        library(grid)
        library(ggthemes)
        (theme_foundation(base_size=base_size, base_family=base_family)
            + theme(plot.title = element_text(face = "bold",
                                              size = rel(1.0), hjust = 0.5),
                    #family = 'Helvetica'
                    text = element_text(),
                    panel.background = element_rect(colour = NA),
                    plot.background = element_rect(colour = NA),
                    panel.border = element_rect(colour = NA),
                    axis.title = element_text(face = "bold",size = rel(0.8)),
                    axis.title.y = element_text(angle=90,vjust =2),
                    axis.title.x = element_text(vjust = -0.2),
                    axis.text = element_text(), 
                    axis.line.x = element_line(colour="black"),
                    axis.ticks.x = element_line(),
                    axis.ticks.y = element_blank(),
                    panel.grid.major = element_line(colour="#f0f0f0"),
                    panel.grid.minor = element_blank(),
                    legend.key = element_rect(colour = NA),
                    legend.position = "right",
                    # legend.direction = "horizontal",
                    legend.key.size= unit(0.2, "cm"),
                    legend.spacing  = unit(0, "cm"),
                    # legend.title = element_text(face="italic"),
                    plot.margin=unit(c(5,5,5,5),"mm"),
                    strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                    strip.text = element_text(face="bold"),
                    plot.caption = element_text(face = "italic", size=rel(0.6))
            ))
    } 
    
    if(writetable == FALSE & figure == FALSE){
        print("No output set!")
    }    else{
        ## Variable selection
        dfsub <- df %>% select(CKD_group, Age, Sex, BMI, HT, Ethnicity, FramRisk, tail(names(.), 20)) %>% 
            mutate(CKD_log=ifelse(CKD_group=="Controls", 0, 1))
        
        ## Models
        res_log <- c()
        for (i in c((ncol(dfsub)-19):ncol(dfsub)-1)){
            dfsub$met <- NULL    
            dfsub$met <- log10(dfsub[,i][[1]])
            m0 <- glm(CKD_log ~ scale(met), data = dfsub, family = "binomial")
            m1 <- glm(CKD_log ~ scale(met) + Age + Sex + BMI + HT, data=dfsub, family="binomial")
            
            metname <- colnames(dfsub)[i]
            m0 <- tidy(m0, conf.int=T)[2,]
            m1 <- tidy(m1, conf.int = T)[2,]
            
            resRow <- cbind(metname, exp(m0$estimate), exp(m0$conf.low), exp(m0$conf.high), m0$p.value,
                            exp(m1$estimate), exp(m1$conf.low), exp(m1$conf.high), m1$p.value)
            colnames(resRow) <- c("Metabolite", 
                                  "m0-est", "m0-l95", "m0-u95", "m0-p", 
                                  "m1-est", "m1-l95", "m1-u95", "m1-p")
            res_log <- rbind(res_log, resRow)
            dfsub$met <- NULL 
        }
        
        reslog <- as.data.frame(res_log)
        afronden2 <- function(x) return(as.numeric(format(round(x, 2),2)))
        afronden5 <- function(x) return(as.numeric(format(round(x, 5),5)))
        reslog2 <- reslog %>% 
            mutate_at(c(2:9), as.character) %>% 
            mutate_at(c(2:9), as.numeric) %>% 
            mutate_at(c(2:4, 6:8), afronden2) %>% 
            mutate_at(c(5,9), afronden5) %>% 
            mutate(
                `m0-q` = p.adjust(`m0-p`, 'fdr'),
                `m1-q` = p.adjust(`m1-p`, 'fdr')
            )
        
        ## Output
        if(writetable == TRUE){
            openxlsx::write.xlsx(reslog2, file.path("results/logreg/DKD", str_c(dfname,"_logreg.xlsx")))
        }
        if(figure == TRUE){
            labslist <- c("Unadjusted", "Adjusted")
            reslong <- reslog2 %>% 
                pivot_longer(c(2:9), names_to=c("model", "cat"), 
                             names_prefix="m", 
                             names_sep='-',
                             values_to="value") %>% 
                pivot_wider(names_from = cat, values_from = value) %>% 
                mutate(model = factor(model, levels = c("0", "1"), 
                                      labels = labslist),
                       Metabolite = factor(Metabolite, levels = colnames(dfsub)[8:27]),
                       Metabolite = fct_rev(Metabolite))
            
            ylab <- "OR for DKD per log10 increase"
            colors <- c(pal_jco()(4)[1], pal_jco()(4)[4])
            pl <- ggplot(reslong, aes(x=Metabolite,y=est, color=model)) +
                geom_hline(yintercept = 1, color = "grey40") +
                geom_point(position=position_dodge(-0.5)) +
                geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.5)) +
                expand_limits(y=0)+
                scale_y_log10(n.breaks = 6)+
                theme_Publication()+
                labs(title = plottitle, x = "", y = ylab) +
                scale_color_manual(values = colors) +
                coord_flip()
            return(pl)
        } 
        
        
    }
}

## Opening HELIUS file, metabolite data and best predictor files
helius_dkd <- readRDS("data/helius_dkd.RDS")
infomet_plasma <- readRDS("data/infomet.RDS") 
infomet_urine <- readRDS("data/infomet_urine.RDS")
plasma_metabolites <- readRDS("data/plasma_metabolites.RDS")
urine_metabolites <- readRDS("data/urine_metabolites.RDS")
best_ckd_plasma <- rio::import('CKD_plasma/output_XGB_class_nonD_CKD_plasma_2023_11_28__23-15-43/feature_importance.txt')
best_ckd_urine <- rio::import('CKD_urine_kreat/output_XGB_class_nonD_CKD_urinekreat_2023_11_28__23-15-36/feature_importance.txt')

## Preparation for models
pl_dkd <- glm_prep(best_ckd_plasma, plasma_metabolites, helius_dkd)
ur_dkd <- glm_prep(best_ckd_urine, urine_metabolites, helius_dkd)

## Logistic regression DKD
plot_p <- log_group_dkd(pl_dkd, "plasma_dkd", writetable = TRUE, figure = TRUE, 
                        "Plasma metabolites") 
plot_u <- log_group_dkd(ur_dkd, "urine_dkd", writetable = TRUE, figure = TRUE,
                        "Urine metabolites")

pl <- ggarrange(plot_p, plot_u, labels = c("A", "B"))

ggsave(pl, filename = "results/logreg/DKD/plasmaurine_logreg.pdf", device = "pdf", width = 16)
ggsave(pl, filename = "results/logreg/DKD/plasmaurine_logreg.svg", device = "svg", width = 16)
