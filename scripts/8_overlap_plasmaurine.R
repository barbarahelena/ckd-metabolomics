## Overlap best predictors plasma / urine

packages <- c("rio", "haven", "dplyr", "tidyverse", "corrr", "corrplot", "Hmisc", "igraph", 
              "ggraph", "network", "naniar", "ggpubr", "ggsci", "gridExtra",
              "RColorBrewer")
pacman::p_load(packages, character.only = TRUE)

overlap <- function(best_pl, best_u, modelname){
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
    pl1 <- best_pl %>% 
        select(Metabolite=FeatName, RelFeatImp) %>% 
        arrange(-RelFeatImp) %>% 
        slice(1:100)
    
    u1 <- best_u %>% 
        select(Metabolite=FeatName, RelFeatImp) %>% 
        arrange(-RelFeatImp) %>% 
        slice(1:100)
    
    pl2 <- best_pl %>% 
        select(Metabolite=FeatName, RelFeatImp) %>% 
        arrange(-RelFeatImp) %>% 
        slice(1:20)
    
    u2 <- best_u %>% 
        select(Metabolite=FeatName, RelFeatImp) %>% 
        arrange(-RelFeatImp) %>% 
        slice(1:20)
    
    plu <- inner_join(pl1, u1, by = "Metabolite", suffix = c("Plasma", "Urine")) %>% 
        arrange(-RelFeatImpPlasma)
    colnames(plu) <- c("Metabolite", "Plasma", "Urine")
    plu2 <- pivot_longer(plu, cols = 2:3)
    
    plot <- ggplot(plu2, aes(x = fct_reorder(Metabolite, value), y = value, fill = name)) + 
        theme_Publication() +
        geom_bar(data=subset(plu2, name == "Plasma"), aes(y = value*-1), stat = "identity", alpha = 0.8) + 
        geom_bar(data=subset(plu2, name == "Urine"), stat = "identity", alpha = 0.8) + 
        scale_y_continuous(limits = c(-100, 100), breaks=seq(-100,100,20),labels=abs(seq(-100,100,20))) + 
        coord_flip() + 
        scale_fill_manual(values = rev(pal_lancet()(2))) + 
        labs(title = "Overlapping best predictors (top 100)", x="",
             y="Relative importance (%)" ) +
        theme(legend.title = element_blank())
    ggsave(str_c('results/overlap100_', modelname ,'.pdf'), height = 7, width = 10)
    
    u2$unique <- ifelse(u2$Metabolite %in% pl1$Metabolite, "in both models", "only in this model")
    pl2$unique <- ifelse(pl2$Metabolite %in% u1$Metabolite, "in both models", "only in this model")
    pl2$unique <- as.factor(pl2$unique)
    u2$unique <- as.factor(u2$unique)
    
    plot1 <-  ggplot(data=u2, aes(y=RelFeatImp, x=fct_reorder(Metabolite, RelFeatImp), fill=unique)) +
        scale_fill_manual(values = c("in both models"=pal_lancet()(2)[1], "only in this model"="gray")) +
        theme_Publication()+
        geom_bar(stat="identity", alpha = 0.8)+
        coord_flip() +
        labs(title="Urine", y="Relative Importance (%)", x="", caption = "") +
        theme(legend.title = element_blank(),
              legend.position = "right")
    
    plot2 <-  ggplot(data=pl2, aes(y=RelFeatImp, x=fct_reorder(Metabolite, RelFeatImp), fill=unique)) +
        scale_fill_manual(values = c("in both models"=pal_lancet()(2)[2], "only in this model"="gray")) +
        theme_Publication()+
        geom_bar(stat="identity", alpha = 0.8)+
        coord_flip() +
        labs(title="Plasma", y="Relative Importance (%)", x="", caption = "") +
        theme(legend.title = element_blank(),
              legend.position = "right")
    
    p <- ggarrange(plot, ggarrange(plot1, plot2, nrow = 2, labels = c("B", "C")),
                   ncol = 2, labels = "A")
    annotate_figure(p, bottom = text_grob("Figure B and C show the top 20 predictors; these are marked as \"in both models\" 
    if this predictor appears in the top 100 predictors of the other category (plasma or urine)",
                                          hjust = 1, x=1, face = "italic", size = 10))
    ggsave(str_c('results/overlap_', modelname,'.pdf'), width=18, height=8.5)
}


## CKD vs HC 
best_pl <- rio::import('Plasma/CKD/output_XGB_class_CKD_2020_12_18__22-07-15/feature_importance.txt')
best_u <- rio::import('Urine/CKD/output_XGB_class_UrineCKD_2021_01_11__14-46-07/feature_importance.txt')
overlap(best_pl, best_u, "CKDvsHC")

## CKD-T2DM vs HC 
best_pl <- rio::import('Plasma/CKD/WithoutDM/output_XGB_class_CKDnoDM_2020_12_19__11-53-08/feature_importance.txt')
best_u <- rio::import('Urine/CKD/WithoutDM/output_XGB_class_UrineCKDsens_2021_01_11__22-53-32/feature_importance.txt')
overlap(best_pl, best_u, "CKD-T2DMvsHC")

## CKD+T2DM vs CKD-T2DM
best_pl <- rio::import('Plasma/CKDDM/output_XGB_class_CKDDM_2020_12_19__13-54-29/feature_importance.txt')
best_u <- rio::import('Urine/CKDDM/output_XGB_class_UrineCKDDM_2021_01_12__10-30-18/feature_importance.txt')
overlap(best_pl, best_u, "CKD+T2DMvsCKD-T2DM")


