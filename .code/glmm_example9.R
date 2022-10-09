## ----setup, include=FALSE, warnings=FALSE, message=FALSE----------------------
knitr::opts_chunk$set(echo = TRUE)


## ----libraries, results='markdown', eval=TRUE, message=FALSE, warning=FALSE----
library(car)       #for regression diagnostics
library(broom)     #for tidy output
library(broom.mixed) #for tidy output
library(ggfortify) #for model diagnostics
library(sjPlot)    #for outputs
library(knitr)     #for kable
library(effects)   #for partial effects plots
library(ggeffects) #for effects plots in ggplotjk
library(emmeans)   #for estimating marginal means
library(MASS)      #for glm.nb
library(MuMIn)     #for AICc
library(tidyverse) #for data wrangling
library(DHARMa)    #for assessing dispersion etc
library(glmmTMB)    #for glmmTMB
library(performance) #for diagnostic plots
library(see)         #for diagnostic plots
library(ordinal)    #for ordinal models


## ----readDataP, results='markdown', eval=TRUE---------------------------------
hughes = read_csv('../data/hughes_full.csv', trim_ws=TRUE)
glimpse(hughes)


## ----processDataP, results='markdown', eval=TRUE------------------------------
hughes = hughes %>%
    mutate(fYear=factor(Year),
           Score=ifelse(Score==5,4,Score),
           oScore = factor(Score, ordered=TRUE),
           nScore = as.numeric(factor(Score, ordered=TRUE)),
           SectorThree=factor(SectorThree, levels=c('North','Central','South')),
           fReef=factor(ReefID),
           nReef=as.numeric(fReef))
# now make a version that is just 2016
hughes <- hughes %>% filter(fYear==2016) %>%
    dplyr::select(REEF=ReefID, HABITAT=Habitat, SECTOR=SectorThree, SCORE=Score)
write_csv(hughes, file='../data/hughes.csv')
## hughes.colors = c('#FFFFFF', rev(heat.colors(length(levels(hughes$oScore))))[-1])


## ----readData, results='markdown', eval=TRUE----------------------------------
hughes = read_csv('../data/hughes.csv', trim_ws=TRUE)
glimpse(hughes)


## ----processData, results='markdown', eval=TRUE-------------------------------
hughes = hughes %>%
    mutate(oSCORE = factor(SCORE, ordered=TRUE),
           HABITAT = factor(HABITAT),
           SECTOR=factor(SECTOR, levels=c('North','Central','South')),
           REEF=factor(REEF))
## hughes.colors = c('#FFFFFF', rev(heat.colors(length(levels(hughes$oSCORE))))[-1])


## ----EDA2, results='markdown', eval=TRUE, fig.width=10, fig.height=10---------
# Scatterplot
hughes %>%
    ggplot(aes(y = oSCORE, x = HABITAT)) +
    geom_point(position = position_jitter()) +
    facet_wrap(~SECTOR)
# Line plot - too messy to make much of this...
hughes %>%
    group_by(SECTOR, REEF, HABITAT) %>%
    summarise(SCORE = mean(SCORE)) %>%
    ungroup() %>%
    ggplot(aes(y = SCORE, x = as.numeric(HABITAT), group = REEF)) +
    geom_blank(aes(x = HABITAT)) +
    geom_line() +
    facet_grid(~SECTOR)


hughes %>%
    group_by(SECTOR, HABITAT, oSCORE) %>%
    summarise(n = n()) %>%
    ungroup() %>%
    group_by(SECTOR, HABITAT) %>%
    mutate(prop = n/sum(n)) %>%
    mutate(oSCORE = factor(oSCORE, levels = rev(levels(oSCORE)))) ->
    hughes.sum
    
## hughes.sum <- hughes %>%
##     count(SECTOR,HABITAT,oSCORE) %>%
##     group_by(SECTOR, HABITAT) %>%
##     mutate(prop=prop.table(n),
##            oSCORE=factor(oSCORE, levels=rev(levels(oSCORE))))

hughes.sum %>% head

ggplot(data=hughes.sum, aes(y=prop, x=HABITAT)) +
    geom_bar(stat='Identity', aes(fill=oSCORE), color='black') +
    facet_grid(~SECTOR) +
    ## scale_fill_manual('Bleaching score', values=rev(hughes.colors) ) +
    scale_fill_manual('Bleaching score', values=c(heat.colors(5)[-5], '#FFFFFF') ) +
    scale_y_continuous('Proportion of Reef', expand=c(0,0))+
    theme_bw() + 
    theme(panel.spacing.y=unit(10,'pt'))


## ----mod1Fit, results='markdown', eval=FALSE, cache=TRUE, results='hide'------
## library(ordinal)
## hughes.clmm=ordinal::clmm(oSCORE ~ HABITAT*SECTOR+(1|REEF), data=hughes)
## hughes.clmm=clmm(oSCORE ~ HABITAT*SECTOR+(1|REEF), data=hughes)
## hughes.clmm1=ordinal::clmm(oSCORE ~ HABITAT*SECTOR+(HABITAT|REEF), data=hughes)
## AIC(hughes.clmm, hughes.clmm1)
## allEffects(hughes.clmm) %>% plot
## allEffects(hughes.clmm1) %>% plot
## ## plot_model(hughes.clmm, type='eff', terms=c('HABITAT','SECTOR'))
## ## plot_model(hughes.clmm, type='eff', terms=c('SECTOR','HABITAT'))
## ## ggpredict(hughes.clmm) %>% plot
## hughes.clmm1 %>% ggpredict(c('HABITAT', 'SECTOR')) %>% plot
## hughes.clmm1 %>% ggemmeans(~HABITAT|SECTOR) %>% plot
## 
## ## predict(hughes.clmm)
## ## model=hughes.clmm
## ## coefs <- c(model$beta, unlist(model$ST))
## #autoplot(hughes.clmm1, what='qq')
## ## simulateResiduals(hughes.clmm, plot=TRUE)
## ## hughes.clmm=ordinal::clmm(oScore ~ Habitat*fYear*SectorThree+(1|ReefID), data=hughes)
## summary(hughes.clmm)
## summary(hughes.clmm1)
## ## Habitat F is associated with a lower prob of higher bleaching - but not significant
## ## Habitat L is associated with a lower prob of higher bleaching (In the North)
## exp(-1.655)
## ## Habitat U is associated with higher prob of higher bleaching (in north)
## exp(0.99968)
## ## Central is associated with lower prob of higher bleaching (in habitat c)
## exp(-2.29)
## ## South is associated with a lower prob of higher bleaching
## exp(-7.51)
## # Significant interactions due to F:Central (Lower) and F:Southern (lower)
## 
## exp(-9.1435)
## exp(-0.7224)
## exp(-0.3714)
## 
## emmeans(hughes.clmm, ~oSCORE|HABITAT+SECTOR, mode='prob')
## emmeans(hughes.clmm1, ~oSCORE|HABITAT+SECTOR, mode='prob')
## emmeans(hughes.clmm, ~HABITAT|SECTOR, mode='mean.class')
## emmeans(hughes.clmm1, ~HABITAT|SECTOR, mode='mean.class')
## emmeans(hughes.clmm, ~HABITAT|SECTOR, mode='mean.class') %>% pairs()
## emmeans(hughes.clmm1, ~HABITAT|SECTOR, mode='mean.class') %>% pairs()
## 
## ## newdata = emmeans(hughes.clmm, ~ HABITAT|SECTOR, mode='mean.class') %>% as.data.frame %>%
## newdata = emmeans(hughes.clmm1, ~ HABITAT|SECTOR, mode='mean.class') %>% as.data.frame %>%
##     mutate(across(c(mean.class, asymp.LCL, asymp.UCL), function(x) x-1))
## newdata
## ScoreBoundaries = data.frame(Score=factor(0:4), ymin=c(0:4), ymax=c(1:5))
## ggplot(newdata) +
##     geom_blank(aes(y=mean.class, x=HABITAT)) +
##     geom_hline(yintercept=1, linetype='dashed', size=0.1) +
##     geom_hline(yintercept=2, linetype='dashed', size=0.1) +
##     geom_hline(yintercept=3, linetype='dashed', size=0.1) +
##     #geom_rect(data=ScoreBoundaries, aes(ymin=ymin, ymax=ymax, xmin=-Inf, xmax=Inf, fill=Score), alpha=0.2) +
##     geom_pointrange(aes(y=mean.class, x=HABITAT, ymin=asymp.LCL, ymax=asymp.UCL)) +
##     facet_grid(~SECTOR) +
##     scale_y_continuous('Bleaching score', breaks=(0:4), labels=0:4, limits=c(0,4),expand=c(0,0)) +
##     theme_bw() +
##     theme(panel.spacing.y=unit(10,'pt'))
## 
## 
## ## Pairwise for habitat
## ## newdata = emmeans(hughes.clmm, ~HABITAT|SECTOR, mode='mean.class') %>%
## newdata = emmeans(hughes.clmm1, ~HABITAT|SECTOR, mode='mean.class') %>%
##     pairs() %>% confint() %>% as.data.frame()
## 
## ggplot(newdata) +
##     geom_hline(yintercept=0) +
##     geom_pointrange(aes(y=estimate, x=contrast, ymin=asymp.LCL, ymax=asymp.UCL, color=SECTOR),
##                     position=position_dodge(width=0.5)) +
##     facet_grid(~SECTOR) +
##     coord_flip() +
##     scale_y_continuous('Effect size')+
##     theme_bw()
## 
## ## Pairwise for year
## 
## ## newdata = emmeans(hughes.clmm, pairwise~ fYear|Habitat+SectorThree, mode='mean.class')$contrasts %>%
## ##                                                                                   confint %>% as.data.frame
## ## ggplot(newdata) +
## ##     geom_hline(yintercept=0) +
## ##     geom_pointrange(aes(y=estimate, x=Habitat, ymin=asymp.LCL, ymax=asymp.UCL, color=SectorThree),
## ##                     position=position_dodge(width=0.5)) +
## ##     facet_grid(~SectorThree) +
## ##     coord_flip() +
## ##     scale_y_continuous('Effect size')+
## ##     theme_bw()
## 


## ----name, results='markdown', eval=FALSE-------------------------------------
## hughes %>% dplyr::filter(fYear=='2016', SectorThree=='North') %>%
##     group_by(Habitat) %>%
##     count(oScore) %>%
##     group_by(Habitat) %>%
##     mutate(Prop = n/sum(n))
##     ggplot() + geom_point(aes(oScore
## 
## hughes.clmm=ordinal::clmm(oScore ~ Habitat+(1|ReefID),
##                           data=hughes %>% dplyr::filter(fYear=='2016', SectorThree=='North'))
## summary(hughes.clmm)
## emmeans(hughes.clmm, ~oScore|Habitat, mode='prob')
## emmeans(hughes.clmm, ~Habitat, mode='mean.class')
## emmeans(hughes.clmm, pairwise~Habitat, mode='mean.class')
## 

