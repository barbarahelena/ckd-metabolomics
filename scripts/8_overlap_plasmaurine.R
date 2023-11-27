## Linear models
## Barbara Verhaar, b.j.verhaar@amsterdamumc.nl

## Libraries
packages <- c("rio", "dplyr", "tidyverse", "Hmisc", "ggpubr", "ggsci", "RColorBrewer", "broom", "aplot")
pacman::p_load(packages, character.only = TRUE)

## Functions
distribution <- function(df, col){
    theme_Empty <- function() {
        theme_minimal() +
            theme(axis.title = element_blank(), axis.ticks.x = element_blank(), axis.line.x = element_blank(), 
                  axis.text.x = element_blank(), axis.text.y = element_blank(), panel.grid = element_blank(), 
                  axis.ticks.y = element_blank(),
                  plot.margin=unit(c(4,0,0,0), "mm")
            )
    }
    pl <- ggplot(df) + 
        geom_density(aes(df[,col]), fill = pal_nejm()(5)[5], color = pal_nejm()(5)[5]) +      
        theme_Empty()
}

lm_prep <- function(df, metabolites, helius){
    df <- df %>% select(Metabolite=FeatName, RelFeatImp) %>% 
        arrange(-RelFeatImp) %>% slice(1:20)
    
    df_met <- metabolites %>% inner_join(., df, by='Metabolite') %>% 
        arrange(-RelFeatImp)
    rownames(df_met) <- df_met$Metabolite
    df_met$Metabolite <- NULL
    df_met$RelFeatImp <- NULL
    df_met <- as.data.frame(t(as.matrix(df_met)))
    head(df_met)
    df_met$ID <- as.integer(rownames(df_met))
    helius$ID <- helius$Heliusnr
    helius_sub <- left_join(helius, df_met, by='ID')
    return(helius_sub)
}

log_group <- function(df, dfname, writetable = FALSE, figure = FALSE, contrast = "CKD"){
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
            dfsub <- df %>% select(CKD_group, Age, Sex, BMI, HT, DM, Ethnicity, FramRisk, tail(names(.), 20)) %>% 
                mutate(CKD_log=ifelse(CKD_group=="Controls", 0, 1),
                       DM_log=ifelse(DM=="No diabetes", 0, 1))
        
            ## Models
            res_log <- c()
            for (i in c(9:28)){
                dfsub$met <- NULL    
                dfsub$met <- dfsub[,i]
                if(contrast=="CKD"|contrast=="Sens"){
                    m0 <- glm(CKD_log ~ scale(met), data = dfsub, family = "binomial")
                    m1 <- glm(CKD_log ~ scale(met) + Age + Sex + BMI, data = dfsub, family = "binomial")
                    m2 <- glm(CKD_log ~ scale(met) + Age + Sex + BMI + HT + DM, data=dfsub, family="binomial")
                    m3 <- glm(CKD_log ~ scale(met) + Age + Sex + BMI + HT + DM + Ethnicity, data=dfsub, family="binomial")
                }else{
                    m0 <- glm(DM_log ~ scale(met), data = dfsub, family = "binomial")
                    m1 <- glm(DM_log ~ scale(met) + Age + Sex + BMI, data = dfsub, family = "binomial")
                    m2 <- glm(DM_log ~ scale(met) + Age + Sex + BMI + HT, data=dfsub, family = "binomial")
                    m3 <- glm(DM_log ~ scale(met) + Age + Sex + BMI + HT + Ethnicity, data=dfsub, family = "binomial")
                }
                
                metname <- colnames(dfsub)[i]
                m0 <- tidy(m0, conf.int=T)[2,]
                m1 <- tidy(m1, conf.int=T)[2,]
                m2 <- tidy(m2, conf.int = T)[2,]
                m3 <- tidy(m3, conf.int = T)[2,]
                
                resRow <- cbind(metname, exp(m0$estimate), exp(m0$conf.low), exp(m0$conf.high), m0$p.value,
                                exp(m1$estimate), exp(m1$conf.low), exp(m1$conf.high), m1$p.value,
                                exp(m2$estimate), exp(m2$conf.low), exp(m2$conf.high), m2$p.value,
                                exp(m3$estimate), exp(m3$conf.low), exp(m3$conf.high), m3$p.value)
                colnames(resRow) <- c("Metabolite", 
                                      "m0-est", "m0-l95", "m0-u95", "m0-p", 
                                      "m1-est", "m1-l95", "m1-u95", "m1-p",
                                      "m2-est", "m2-l95", "m2-u95", "m2-p",
                                      "m3-est", "m3-l95", "m3-u95", "m3-p")
                res_log <- rbind(res_log, resRow)
                dfsub$met <- NULL 
            }
            
            reslog <- as.data.frame(res_log)
            afronden2 <- function(x) return(as.numeric(format(round(x, 2),2)))
            afronden5 <- function(x) return(as.numeric(format(round(x, 5),5)))
            reslog2 <- reslog %>% 
                mutate_at(c(2:17), as.character) %>% 
                mutate_at(c(2:17), as.numeric) %>% 
                mutate_at(c(2:4, 6:8, 10:12, 14:17), afronden2) %>% 
                mutate_at(c(5,9,13,17), afronden5)
            
            ## Output
            if(figure == TRUE){
                if(contrast=="CKD"){
                    labslist <- c("Unadjusted", "Age, Sex, BMI", "+Diabetes, Hypertension", "+Ethnicity")
                }else{
                    labslist <- c("Unadjusted", "Age, Sex, BMI", "+Hypertension", "+Ethnicity")
                }
                reslong <- reslog2 %>% 
                    pivot_longer(c(2:17), names_to=c("model", "cat"), 
                                 names_prefix="m", 
                                 names_sep='-',
                                 values_to="value") %>% 
                    pivot_wider(names_from = cat, values_from = value) %>% 
                    mutate(model = factor(model, levels = c("0", "1", "2", "3"), 
                                          labels = labslist),
                           Metabolite = factor(Metabolite, levels = colnames(dfsub)[9:28]),
                           Metabolite = fct_rev(Metabolite))
                
                ylab <- ifelse(contrast=="DM", "OR for T2DM per SD increase", "OR for CKD per SD increase")
                colors <- c(pal_jco()(4)[1], pal_jco()(4)[2], pal_jco()(4)[4], pal_jco()(5)[5])
                pl <- ggplot(reslong, aes(x=Metabolite,y=est, color=model)) +
                    geom_hline(yintercept = 1, color = "grey40") +
                    geom_point(position=position_dodge(-0.8)) +
                    geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.8)) +
                    expand_limits(y=0)+
                    theme_Publication()+
                    labs(title = "Logistic regression: best predictors", x = "", y = ylab) +
                    scale_color_manual(values = colors) +
                    coord_flip()
                #ggsave(pl, filename = str_c("results/logreg/", dfname, "_logreg.pdf"), device = "pdf", width = 8)
                #ggsave(pl, filename = str_c("results/logreg/", dfname, "_logreg.svg"), device = "svg", width = 8)
                return(pl)
            } 
            
            if(writetable == TRUE){
                openxlsx::write.xlsx(reslog2, file.path("results/logreg", str_c(dfname,"_logreg.xlsx")))
            }
    }
}

lin_model <- function(df, dfname, writetable = FALSE, figure = FALSE, contrast = "CKD"){
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
        if(contrast=="CKD"|contrast=="Sens"){
            df <- df %>% filter(CKD_group=="CKD")
        }
        if(contrast=="DM"){
            df <- df %>% filter(DM=="Diabetes")
        }
        dfsub <- df %>% select(Alb, Age, Sex, BMI, HT, DM, Ethnicity, FramRisk, tail(names(.), 20)) %>% 
            mutate(logAlb = log(Alb+1))
        hist(dfsub$logAlb)
        
        ## Models
        res_lin <- c()
        for (i in c(9:28)){
            dfsub$met <- NULL    
            dfsub$met <- dfsub[,i]
            if(contrast=="CKD"|contrast=="Sens"){
                m0 <- glm(logAlb ~ scale(met), data = dfsub, family = "gaussian")
                m1 <- glm(logAlb ~ scale(met) + Age + Sex + BMI, data = dfsub, family = "gaussian")
                m2 <- glm(logAlb ~ scale(met) + Age + Sex + BMI + HT + DM, data = dfsub, family = "gaussian")
                m3 <- glm(logAlb ~ scale(met) + Age + Sex + BMI + HT + DM + Ethnicity, data=dfsub, family="gaussian")
            }else{
                m0 <- glm(logAlb ~ scale(met), data = dfsub, family = "gaussian")
                m1 <- glm(logAlb ~ scale(met) + Age + Sex + BMI, data = dfsub, family = "gaussian")
                m2 <- glm(logAlb ~ scale(met) + Age + Sex + BMI + HT, data=dfsub, family = "gaussian")
                m3 <- glm(logAlb ~ scale(met) + Age + Sex + BMI + HT + Ethnicity, data=dfsub, family = "gaussian")
            }
            
            metname <- colnames(dfsub)[i]
            m0 <- tidy(m0, conf.int=T)[2,]
            m1 <- tidy(m1, conf.int=T)[2,]
            m2 <- tidy(m2, conf.int = T)[2,]
            m3 <- tidy(m3, conf.int = T)[2,]
            
            resRow <- cbind(metname, m0$estimate, m0$conf.low, m0$conf.high, m0$p.value,
                            m1$estimate, m1$conf.low, m1$conf.high, m1$p.value,
                            m2$estimate, m2$conf.low, m2$conf.high, m2$p.value,
                            m3$estimate, m3$conf.low, m3$conf.high, m3$p.value)
            colnames(resRow) <- c("Metabolite", 
                                  "m0-est", "m0-l95", "m0-u95", "m0-p", 
                                  "m1-est", "m1-l95", "m1-u95", "m1-p",
                                  "m2-est", "m2-l95", "m2-u95", "m2-p",
                                  "m3-est", "m3-l95", "m3-u95", "m3-p")
            res_lin <- rbind(res_lin, resRow)
            dfsub$met <- NULL 
        }
        
        reslin <- as.data.frame(res_lin)
        afronden2 <- function(x) return(as.numeric(format(round(x, 2),2)))
        afronden5 <- function(x) return(as.numeric(format(round(x, 5),5)))
        reslin2 <- reslin %>% 
            mutate_at(c(2:17), as.character) %>% 
            mutate_at(c(2:17), as.numeric) %>% 
            mutate_at(c(2:4, 6:8, 10:12, 14:17), afronden2) %>% 
            mutate_at(c(5,9,13,17), afronden5)
        
        ## Output
        if(figure == TRUE){
            if(contrast=="CKD"){
                labslist <- c("Unadjusted", "Age, Sex, BMI", "+Diabetes, Hypertension", "+Ethnicity")
            }else{
                labslist <- c("Unadjusted", "Age, Sex, BMI", "+Hypertension", "+Ethnicity")
            }
            reslong <- reslin2 %>% 
                pivot_longer(c(2:17), names_to=c("model", "cat"), 
                             names_prefix="m", 
                             names_sep='-',
                             values_to="value") %>% 
                pivot_wider(names_from = cat, values_from = value) %>% 
                mutate(model = factor(model, levels = c("0", "1", "2", "3"), 
                                      labels = labslist),
                       Metabolite = factor(Metabolite, levels = colnames(dfsub)[9:28]),
                       Metabolite = fct_rev(Metabolite))
            
            colors <- c(pal_jco()(4)[1], pal_jco()(4)[2], pal_jco()(4)[4], pal_jco()(5)[5])
            pl <- ggplot(reslong, aes(x=Metabolite,y=est, color=model)) +
                geom_hline(yintercept = 0, color = "grey40") +
                geom_point(position=position_dodge(-0.8)) +
                geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.8)) +
                theme_Publication()+
                labs(title = "Linear regression: best predictors - albuminuria", x = "", 
                     y = "Difference in albuminuria per SD increase metabolite") +
                scale_color_manual(values = colors) +
                coord_flip()
            ggsave(pl, filename = str_c("results/linreg/", dfname, "_linreg.pdf"), device = "pdf", width = 8)
            ggsave(pl, filename = str_c("results/linreg/", dfname, "_linreg.svg"), device = "svg", width = 8)
        } 
        
        if(writetable == TRUE){
            openxlsx::write.xlsx(reslin2, file.path("results/linreg", str_c(dfname,"_linreg.xlsx")))
        }
    }       
}

log_strata <- function(df, dfname, eth, writetable = FALSE, figure = TRUE, contrast = "CKD"){
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
        dfsub <- df %>% filter(Ethnicity==eth) %>% 
            select(CKD_group, Age, Sex, BMI, HT, DM, Ethnicity, FramRisk, tail(names(.), 20)) %>% 
            mutate(CKD_log=ifelse(CKD_group=="Controls", 0, 1),
                   DM_log=ifelse(DM=="No diabetes", 0, 1)) 
        
        ## Models
        res_log <- c()
        #for (a in levels(Ethnicity)){
            #dfsub <- dfsub %>% filter(Ethnicity==a)
                for (i in c(9:28)){
                    dfsub$met <- NULL    
                    dfsub$met <- dfsub[,i]
                    if(contrast=="CKD"|contrast=="Sens"){
                        m0 <- glm(CKD_log ~ scale(met), data = dfsub, family = "binomial")
                        m1 <- glm(CKD_log ~ scale(met) + Age + Sex + BMI, data = dfsub, family = "binomial")
                        m2 <- glm(CKD_log ~ scale(met) + Age + Sex + BMI + HT, data=dfsub, family="binomial") 
                            #sommige groepen hebben geen DM subjects
                    }else{
                        m0 <- glm(DM_log ~ scale(met), data = dfsub, family = "binomial")
                        m1 <- glm(DM_log ~ scale(met) + Age + Sex + BMI, data = dfsub, family = "binomial")
                        m2 <- glm(DM_log ~ scale(met) + Age + Sex + BMI + HT, data=dfsub, family = "binomial")
                    }
                    
                    metname <- colnames(dfsub)[i]
                    m0 <- tidy(m0, conf.int=T)[2,]
                    m1 <- tidy(m1, conf.int=T)[2,]
                    m2 <- tidy(m2, conf.int = T)[2,]
                    
                    resRow <- cbind(metname, exp(m0$estimate), exp(m0$conf.low), exp(m0$conf.high), m0$p.value,
                                    exp(m1$estimate), exp(m1$conf.low), exp(m1$conf.high), m1$p.value,
                                    exp(m2$estimate), exp(m2$conf.low), exp(m2$conf.high), m2$p.value)
                    colnames(resRow) <- c("Metabolite", 
                                          "m0-est", "m0-l95", "m0-u95", "m0-p", 
                                          "m1-est", "m1-l95", "m1-u95", "m1-p",
                                          "m2-est", "m2-l95", "m2-u95", "m2-p")
                    res_log <- rbind(res_log, resRow)
                    dfsub$met <- NULL 
                }
                
                reslog <- as.data.frame(res_log)
                afronden2 <- function(x) return(as.numeric(format(round(x, 2),2)))
                afronden5 <- function(x) return(as.numeric(format(round(x, 5),5)))
                reslog2 <- reslog %>% 
                    mutate_at(c(2:13), as.character) %>% 
                    mutate_at(c(2:13), as.numeric) %>% 
                    mutate_at(c(2:4, 6:8, 10:12), afronden2) %>% 
                    mutate_at(c(5,9,13), afronden5)
                
                ## Output
                if(figure == TRUE){
                    if(contrast=="CKD"){
                        labslist <- c("Unadjusted", "Age, Sex, BMI", "+Diabetes, Hypertension")
                    }else{
                        labslist <- c("Unadjusted", "Age, Sex, BMI", "+Hypertension")
                    }
                    reslong <- reslog2 %>% 
                        pivot_longer(c(2:13), names_to=c("model", "cat"), 
                                     names_prefix="m", 
                                     names_sep='-',
                                     values_to="value") %>% 
                        pivot_wider(names_from = cat, values_from = value) %>% 
                        mutate(model = factor(model, levels = c("0", "1", "2"), 
                                              labels = labslist),
                               Metabolite = factor(Metabolite, levels = colnames(dfsub)[9:28]),
                               Metabolite = fct_rev(Metabolite))
                    
                    ylab <- ifelse(contrast=="DM", "OR for T2DM per SD increase", "OR for CKD per SD increase")
                    colors <- c(pal_jco()(4)[1], pal_jco()(4)[2], pal_jco()(4)[4])
                    pl <- ggplot(reslong, aes(x=Metabolite,y=est, color=model)) +
                        geom_hline(yintercept = 1, color = "grey40") +
                        geom_point(position=position_dodge(-0.8)) +
                        geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.8)) +
                        expand_limits(y=0)+
                        theme_Publication()+
                        labs(title = str_c("Logistic regression: ", str_to_lower(eth)), x = "", y = ylab) +
                        scale_color_manual(values = colors) +
                        coord_flip()
                    #ggsave(pl, filename = str_c("results/logreg/", str_to_lower(eth), dfname, "_logreg.pdf"), device = "pdf", width = 8)
                    #ggsave(pl, filename = str_c("results/logreg/", a, dfname, "_logreg.svg"), device = "svg", width = 8)
                } 
                
                if(writetable == TRUE){
                    openxlsx::write.xlsx(reslog2, file.path("results/logreg", str_c(dfname,str_to_lower(eth),"_logreg.xlsx")))
                }
                return(pl)
        #}
    }
}

log_strata2 <- function(df, dfname, writetable = FALSE, figure = TRUE, contrast = "CKD"){
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
                    axis.text.y = element_text(size=rel(0.7)), 
                    axis.line.x = element_line(colour="black"),
                    axis.ticks.x = element_line(),
                    axis.ticks.y = element_blank(),
                    panel.grid.major = element_line(colour="#f0f0f0"),
                    panel.grid.minor = element_blank(),
                    legend.key = element_rect(colour = NA),
                    legend.position = "right",
                    #legend.direction = "horizontal",
                    legend.text = element_text(size=rel(0.8), hjust = 0),
                    legend.key.size= unit(0.2, "cm"),
                    legend.spacing  = unit(0, "cm"),
                    legend.title = element_blank(),
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
        dfsub <- df %>% select(CKD_group, Age, Sex, BMI, HT, DM, Ethnicity, FramRisk, tail(names(.), 20)) %>% 
            mutate(CKD_log=ifelse(CKD_group=="Controls", 0, 1),
                   DM_log=ifelse(DM=="No diabetes", 0, 1)) 
       
        ## Models
        res_log <- c()
        if(dfname=="ckddm"){
            dfsub <- dfsub %>% filter(Ethnicity=="African Surinamese"|Ethnicity=="South-Asian Surinamese")
            dfsub$Ethnicity <- droplevels(dfsub$Ethnicity)
        }
        for (a in levels(dfsub$Ethnicity)){
            dftemp <- dfsub %>% filter(Ethnicity==a)
            for (i in c(9:28)){
                dftemp$met <- NULL    
                dftemp$met <- dftemp[,i]
                m1 <- glm(CKD_log ~ scale(met) + Age + Sex + BMI, data = dftemp, family = "binomial")
                metname <- colnames(dftemp)[i]
                m1 <- tidy(m1, conf.int=T)[2,]
                resRow <- cbind(a, metname, exp(m1$estimate), exp(m1$conf.low), exp(m1$conf.high), m1$p.value)
                colnames(resRow) <- c("Ethnicity","Metabolite", "est", "l95", "u95", "p")
                res_log <- rbind(res_log, resRow)
            }
        }
        reslog <- as.data.frame(res_log)
        
        afronden2 <- function(x) return(as.numeric(format(round(x, 2),2)))
        afronden5 <- function(x) return(as.numeric(format(round(x, 5),5)))
        reslog2 <- reslog %>% 
             mutate_at(c(3:6), as.character) %>% 
             mutate_at(c(3:6), as.numeric) %>% 
             mutate_at(c(3:5), afronden2) %>% 
             mutate_at(6, afronden5)
        
        ## Output
        if(figure == TRUE){
            labslist <- levels(dfsub$Ethnicity)
            reslog2 <- reslog2 %>%
                mutate(Metabolite = factor(Metabolite, levels = colnames(dfsub)[9:28]),
                       Metabolite = fct_rev(fct_inorder(Metabolite)))

            ylab <- ifelse(contrast=="DM", "OR for T2DM per SD increase", "OR for CKD per SD increase")
            colors <- pal_jco()(4)
            pl1 <- reslog2 %>% filter(Metabolite %in% levels(Metabolite)[11:20]) %>% 
                            ggplot(., aes(x=Metabolite,y=est, color=Ethnicity)) +
                                geom_hline(yintercept = 1, color = "grey40") +
                                geom_point(position=position_dodge(-0.8)) +
                                geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.8)) +
                                expand_limits(y=0)+
                                theme_Publication()+
                                labs(x = "", y = ylab) +
                                scale_color_manual(values = colors) +
                                coord_flip(ylim = c(0,10))+
                                theme(legend.position = "none")
            pl2 <- reslog2 %>% filter(Metabolite %in% levels(Metabolite)[1:10]) %>% 
                            ggplot(., aes(x=Metabolite,y=est, color=Ethnicity)) +
                                geom_hline(yintercept = 1, color = "grey40") +
                                geom_point(position=position_dodge(-0.8)) +
                                geom_errorbar(aes(ymin=l95,ymax=u95,width=.3), position=position_dodge(-0.8)) +
                                expand_limits(y=0)+
                                theme_Publication()+
                                labs(x = "", y = ylab) +
                                scale_color_manual(values = colors) +
                                coord_flip(ylim=c(0,10)) 
            pl <- ggarrange(pl1, pl2, widths = c(5,7))
            annotate_figure(pl, top = text_grob("Logistic regression: strata ethnicity", face = "bold"))
            ggsave(str_c("results/logreg/strataeth_", dfname, ".pdf"), width = 12, height = 6)
        }

        # if(writetable == TRUE){
        #     openxlsx::write.xlsx(reslog2, file.path("results/logreg", str_c(dfname,"ethstrata_logreg.xlsx")))
        # }
        # return(pl)
        # }
    }
}

scatter_strata <- function(df, dfname){
    scatter <- function(df, dfname, metabolite){
        theme_Publication <- function(base_size=11, base_family="sans") {
            library(grid)
            library(ggthemes)
            (theme_foundation(base_size=base_size, base_family=base_family)
                + theme(plot.title = element_text(face = "bold",
                                                  size = rel(0.8), hjust = 0.5),
                        text = element_text(),
                        panel.background = element_rect(colour = NA),
                        plot.background = element_rect(colour = NA),
                        panel.border = element_rect(colour = NA),
                        axis.title = element_text(face = "bold",size = rel(0.7)),
                        axis.title.y = element_text(angle=90,vjust = 2),
                        axis.title.x = element_text(vjust = -0.2),
                        axis.text = element_text(size = rel(0.6)), 
                        axis.line = element_line(colour="black"),
                        axis.ticks = element_line(),
                        panel.grid.major = element_line(colour="#f0f0f0"),
                        panel.grid.minor = element_blank(),
                        legend.key = element_rect(colour = NA),
                        legend.position = "none",
                        plot.margin=unit(c(5,2,2,2),"mm"),
                        strip.background=element_rect(colour="#f0f0f0",fill="#f0f0f0"),
                        strip.text = element_text(face="bold", size = rel(0.5))
                ))
            
        } 
        df <- df %>% select(CKD_group, metabolite, Ethnicity)
        names(df)[2] <- 'Metabolite'
        if(dfname=="ckddm"){
            labels <- c("CKD+T2DM", "CKD-T2DM")
        }else{
            labels <- c("CKD", "Controls")
        }
        comps <- list(c(labels[1],labels[2]))
        cl <-  c(pal_jco()(2)[1], pal_jco()(2)[2])
        
        pl <- ggplot(df, aes(x=CKD_group, y=Metabolite, fill=CKD_group)) + 
            geom_violin(alpha=0.7, trim = TRUE) +
            scale_fill_manual(values = cl)+
            geom_boxplot(width=0.1, fill="white")+
            theme_Publication()+
            theme(legend.position = 'none')+
            ggpubr::stat_compare_means(comparisons = comps, paired = F, size=rel(0.5))+
            xlab('Group')+
            ylab('Metabolite concentration')+
            ggtitle(metabolite)+
            facet_wrap(~Ethnicity)
    }
    
    plots <- lapply(tail(colnames(df), 20), scatter, df = df, dfname=dfname)
    ggarrange(plotlist=plots, ncol = 5, nrow=4)
    ggsave(str_c("results/logreg/scatter_",dfname,".pdf"), device = "pdf", width=20, height=20)
}

## Opening HELIUS file, metabolite data and best predictor files
helius <- readRDS("../data/heliusdf.RDS")
summary <- rio::import('../data/Info_plasma_metabolites_b.xlsx')
infomet <- summary %>% select(metabolite=BIOCHEMICAL, sup=`SUPER PATHWAY`, sub=`SUB PATHWAY`)
metabolites <- rio::import('../data/HELIUS_EDTA_plasma_metabolites.xlsx')
best_ckd <- rio::import('Plasma/CKD/output_XGB_class_CKD_2020_12_18__22-07-15/feature_importance.txt')
best_ckdsens <- rio::import('Plasma/CKD/WithoutDM/output_XGB_class_CKDnoDM_2020_12_19__11-53-08/feature_importance.txt')
best_ckddm <- rio::import('Plasma/CKDDM/output_XGB_class_CKDDM_2020_12_19__13-54-29/feature_importance.txt')

## Preparation for models
ckd <- lm_prep(best_ckd, metabolites, helius)
ckdsens <- lm_prep(best_ckdsens, metabolites, helius)
ckddm <- lm_prep(best_ckddm, metabolites, helius)

## Logistic regression
log_group(ckd, "ckd", writetable = TRUE, figure = TRUE, contrast = "CKD") # For table AND figure
log_group(ckdsens, "ckdsens", writetable = TRUE, figure = TRUE, contrast = "Sens")
log_group(ckddm, "ckddm", writetable = TRUE, figure = TRUE, contrast = "DM")

## Linear models
lin_model(ckd, "ckd", writetable = TRUE, figure = TRUE, contrast = "CKD")
lin_model(ckdsens, "ckdsens", writetable = TRUE, figure = TRUE, contrast = "Sens")
lin_model(ckddm, "ckddm", writetable = TRUE, figure = TRUE, contrast = "DM")

## Strata
plots <- lapply(levels(ckdsens$Ethnicity), log_strata, df = ckdsens, dfname="ckdsens", contrast = "Sens")
ggarrange(plots[[1]], plots[[2]], plots[[3]], plots[[4]], ncol = 2, nrow=2, labels = c("A", "B", "C", "D"))
ggsave("results/logreg/logreg_eth.pdf", device = "pdf", width=14, height=10)

## Strata in 1 fig
log_strata2(ckd, dfname="ckd")
log_strata2(ckdsens, dfname="ckdsens")
log_strata2(ckddm, dfname="ckddm")

## Scatter plots
scatter_strata(ckd, dfname="ckd")
scatter_strata(ckdsens, dfname="ckdsens")
scatter_strata(ckddm, dfname="ckddm")


## Files urine
umetabolites <- rio::import('data/HELIUS_urine_metabolites.xlsx')
ubest_ckd <- rio::import('Urine/CKD/output_XGB_class_UrineCKD_2021_01_11__14-46-07/feature_importance.txt')
ubest_ckdsens <- rio::import('Urine/CKD/WithoutDM/output_XGB_class_UrineCKDsens_2021_01_11__22-53-32/feature_importance.txt')
ubest_ckddm <- rio::import('Urine/CKDDM/output_XGB_class_UrineCKDDM_2021_01_12__10-30-18/feature_importance.txt')

## Preparation for models
uckd <- lm_prep(ubest_ckd, umetabolites, helius)
uckdsens <- lm_prep(ubest_ckdsens, umetabolites, helius)
uckddm <- lm_prep(ubest_ckddm, umetabolites, helius)

## Logistic regression
log_group(uckd, "uckd", writetable = TRUE, figure = TRUE, contrast = "CKD") # For table AND figure
log_group(uckdsens, "uckdsens", writetable = TRUE, figure = TRUE, contrast = "Sens")
log_group(uckddm, "uckddm", writetable = TRUE, figure = TRUE, contrast = "DM")

## Linear models
lin_model(uckd, "uckd", writetable = TRUE, figure = TRUE, contrast = "CKD")
lin_model(uckdsens, "uckdsens", writetable = TRUE, figure = TRUE, contrast = "Sens")
lin_model(uckddm, "uckddm", writetable = TRUE, figure = TRUE, contrast = "DM")

## Strata in 1 fig
log_strata2(uckd, dfname="uckd")
log_strata2(ckdsens, dfname="uckdsens")
log_strata2(ckddm, dfname="uckddm")

## Scatter plots
scatter_strata(ckd, dfname="ckd")
scatter_strata(ckdsens, dfname="ckdsens")
scatter_strata(ckddm, dfname="ckddm")



