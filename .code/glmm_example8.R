## ----setup, include=FALSE, warnings=FALSE, message=FALSE----------------------
knitr::opts_chunk$set(echo = TRUE)


## ----libraries, results='markdown', eval=TRUE, message=FALSE, warning=FALSE----
library(car)       #for regression diagnostics
library(broom)     #for tidy output
library(ggfortify) #for model diagnostics
library(sjPlot)    #for outputs
library(knitr)     #for kable
library(effects)   #for partial effects plots
library(emmeans)   #for estimating marginal means
library(MASS)      #for glm.nb
library(MuMIn)     #for AICc
library(tidyverse) #for data wrangling
library(DHARMa)    #for assessing dispersion etc
library(glmmTMB)    #for glmmTMB
library(performance) #for diagnostic plots
library(see)         #for diagnostic plots


## ----readData, results='markdown', eval=TRUE----------------------------------
elston <- read_csv('../data/elston.csv', trim_ws = TRUE)
elston %>% glimpse()


## ----EDA, results='markdown', eval=FALSE, hidden=FALSE------------------------
## ##Response = TICKS
## ## LOCATION (R)
## ## BROOD (R)
## ## ALTITUDE
## ## YEAR
## elston = elston %>%
##     mutate(fYear=factor(YEAR),
##            LOCATION=factor(LOCATION),
##            BROOD = factor(paste0(LOCATION,BROOD)))
## 
## ggplot(elston, aes(y=TICKS, x=ALTITUDE, color=fYear)) +
##   geom_smooth() +
##   geom_point() +
##   scale_y_log10()
## 
## elston.glmmTMB1a <- glmmTMB(TICKS ~ fYear*scale(ALTITUDE) + (1|LOCATION/BROOD),
##                           data=elston,
##                           family=poisson(link='log'),
##                           REML=TRUE)
## elston.glmmTMB1b <- glmmTMB(TICKS ~ fYear*scale(ALTITUDE) + (fYear|LOCATION/BROOD),
##                           data=elston,
##                           family=poisson(link='log'),
##                           REML=TRUE,
##                           control=glmmTMBControl(optimizer = optim,
##                                                  optArgs = list(method = 'BFGS')))
## elston.glmmTMB1c <- glmmTMB(TICKS ~ fYear*scale(ALTITUDE) + (scale(ALTITUDE)|LOCATION/BROOD),
##                           data=elston,
##                           family=poisson(link='log'),
##                           REML=TRUE,
##                           control=glmmTMBControl(optimizer = optim,
##                                                  optArgs = list(method = 'BFGS')))
## elston.glmmTMB1d <- glmmTMB(TICKS ~ fYear*scale(ALTITUDE) + (fYear*scale(ALTITUDE)|LOCATION/BROOD),
##                           data=elston,
##                           family=poisson(link='log'),
##                           REML=TRUE,
##                           control=glmmTMBControl(optimizer = optim,
##                                                  optArgs = list(method = 'BFGS')))
## AICc(elston.glmmTMB1a,  elston.glmmTMB1b, elston.glmmTMB1c, elston.glmmTMB1d)
## 
## 
## 
## plot_model(elston.glmmTMB1a, type='diag') %>% plot_grid
## performance::check_model(elston.glmmTMB1b)
## elston.resid <- elston.glmmTMB1b %>% simulateResiduals(plot=TRUE)
## 
## plot(allEffects(elston.glmmTMB1a),  multiline=TRUE,  ci.style='bands')
## ##plot_model(elston.glmmTMB, type='eff', terms=c('ALTITUDE', 'fYear'))
## 
## summary(elston.glmmTMB1a)
## tidy(elston.glmmTMB1a,  conf.int=TRUE,  exponentiate=TRUE)
## 
## emmeans(elston.glmmTMB1b,  pairwise~fYear|ALTITUDE,  type='response',
##         at=list(ALTITUDE= quantile(elston$ALTITUDE)))
## emmeans(elston.glmmTMB1b,  pairwise~fYear|ALTITUDE,  type='response',
##         at=list(ALTITUDE= quantile(elston$ALTITUDE)))$contrasts %>%
##   confint() %>%
##   as.data.frame %>%
##   ggplot(aes(y=ratio,  x=ALTITUDE,  color=contrast)) +
##   geom_hline(yintercept=1,  linetype='dashed') +
##   geom_pointrange(aes(ymin=lower.CL,  ymax=upper.CL),  position=position_dodge(width=0.2)) +
##   scale_y_log10() +
##   coord_flip()
## 
## elston.grid = with(elston,  list(ALTITUDE=modelr::seq_range(ALTITUDE,  n=100)))
## newdata = emmeans(elston.glmmTMB1b,  ~ALTITUDE|fYear, type='response',
##                   at = elston.grid) %>%
##   as.data.frame
## head(newdata)
## ggplot(newdata) +
##   geom_ribbon(aes(x=ALTITUDE, fill=fYear, ymin=lower.CL, ymax=upper.CL),  alpha=0.3) +
##   geom_line(aes(y=rate, x=ALTITUDE, color=fYear)) +
##   scale_y_log10()

