## ----setup, include=FALSE, warnings=FALSE, message=FALSE----------------------
knitr::opts_chunk$set(echo = TRUE)


## ----libraries, results='markdown', eval=TRUE, message=FALSE, warning=FALSE----
library(car)       #for regression diagnostics
library(broom)     #for tidy output
library(ggfortify) #for model diagnostics
library(sjPlot)    #for outputs
library(knitr)     #for kable
library(effects)   #for partial effects plots
library(DHARMa)    #for residual diagnostics
library(emmeans)   #for estimating marginal means
library(ggeffects) #for diagnostic plots in ggplotjk pR
library(MASS)      #for glm.nb
library(MuMIn)     #for AICc
library(pscl)      #for zero inflated models
library(glmmTMB)
library(nlme)
library(performance) #for residuals diagnostics
library(see)         #for plotting residuals
library(tidyverse) #for data wrangling
library(modelr)    #for auxillary modelling functions


## ----readData, results='markdown', eval=TRUE----------------------------------
fish = read_csv('../data/fish.csv', trim_ws=TRUE)
glimpse(fish)


## ----eda, results='markdown', eval=TRUE, hidden=TRUE--------------------------
fish <- fish %>%
  mutate(CAMPER = factor(CAMPER),
         CHILD=factor(CHILD))

ggplot(fish, aes(y=COUNT, x= CAMPER)) +
  geom_boxplot() +
  geom_point()
ggplot(fish, aes(y=COUNT, x= CAMPER)) +
  geom_boxplot() +
  geom_point()+
  scale_y_continuous(trans=scales::pseudo_log_trans())

ggplot(fish, aes(y=COUNT, x= CHILD, color=CAMPER)) +
  geom_point(position='jitter') +
  scale_y_continuous(trans=scales::pseudo_log_trans())

 

fish <- fish %>% mutate(CHILD=factor(CHILD))

ggplot(fish, aes(y=COUNT, x= CHILD)) + geom_boxplot() +
  scale_y_continuous(trans=scales::pseudo_log_trans())

ggplot(fish, aes(y=COUNT, x= CHILD, color=CAMPER)) + geom_boxplot() +
  scale_y_continuous(trans=scales::pseudo_log_trans())

ggplot(fish,  aes(x=COUNT)) + geom_bar() +
  scale_x_continuous(trans=scales::pseudo_log_trans()) +
  facet_wrap(~PERSONS,  scales='free')

## There is a issue in that there is no variation in response (e.g. all with 3 children caught no fish)


## ----fitModel, hidden=TRUE----------------------------------------------------
fish.glm <- glm(COUNT ~ CHILD * CAMPER,  data=fish,
                family=poisson(link='log'))


## ----ValidateModel, results='markdown', eval=TRUE, fig.width=7, fig.height=7, hidden=TRUE----
autoplot(fish.glm,  which=1:6)
performance::check_model(fish.glm)
performance::check_overdispersion(fish.glm)
performance::check_zeroinflation(fish.glm)
fish.resids <- simulateResiduals(fish.glm, plot=TRUE, integerResponse = TRUE)
testOverdispersion(fish.resids)
testZeroInflation(fish.resids)


## ----ValidateModel2, results='markdown', eval=TRUE, fig.width=7, fig.height=7, hidden=TRUE----

fish.glm1 <- glmmTMB(COUNT ~ CHILD * CAMPER,  ziformula = ~1,
                      data=fish,  family=poisson(link='log'))
fish.glm2 <- glmmTMB(COUNT ~ CHILD * CAMPER, ziformula= ~PERSONS, 
                      data=fish,  family=poisson(link='log'))

#fish.resids <- simulateResiduals(fish.glm2, integerResponse = TRUE)
fish.glm3 <- glmmTMB(COUNT ~ CHILD * CAMPER,  ziformula = ~PERSONS,
                      data=fish,  family=nbinom2(link='log'))

fish.glm4 <- glmmTMB(COUNT ~ CHILD + CAMPER,  ziformula = ~PERSONS,
                      data=fish,  family=nbinom2(link='log'))
AICc(fish.glm, fish.glm1, fish.glm2, fish.glm3, fish.glm4)
performance::check_model(fish.glm4)
fish.resids <- simulateResiduals(fish.glm4,  plot=TRUE, integerResponse = TRUE)
fish.resids <- simulateResiduals(fish.glm3,  plot=TRUE, integerResponse = TRUE)


## library(pscl)

## fish.glm1 <- zeroinfl(COUNT ~ CHILD * CAMPER | 1,
##                       data=fish,  dist='poisson')
## fish.glm2 <- zeroinfl(COUNT ~ CHILD * CAMPER | PERSONS,
##                       data=fish,  dist='poisson')
## #fish.resids <- simulateResiduals(fish.glm2, integerResponse = TRUE)
## fish.glm3 <- zeroinfl(COUNT ~ CHILD * CAMPER | PERSONS,
##                       data=fish,  dist='negbin')

## fish.glm4 <- zeroinfl(COUNT ~ CHILD + CAMPER | PERSONS,
##                       data=fish,  dist='negbin')
## AICc(fish.glm, fish.glm1, fish.glm2, fish.glm3, fish.glm4)




## autoplot(fish.glm4,  which=1:6)
## performance::check_model(fish.glm4)
## simulate(fish.glm4)


## ----partialplots, results='markdown', eval=TRUE, hidden=TRUE-----------------
ggemmeans(fish.glm4,  ~CHILD|CAMPER) %>% plot
plot_model(fish.glm4,  transform='I')

summary(fish.glm4)
                                        #tidy(fish.glm3,  exponentiate=TRUE)
exp(1.22)
exp(1.22+0.96)
exp(1.22-3.43)

exp(1.47)
exp(-1.52)
## when there is only one person in the group,
## - zeros are 4.3 time more likely to be false than true.
## - zeros are 'false' because the person did not go fishing (or lied)
## - hence it is 4.3 times more likely that a 0 is because they did not go fishing than they did go fishing.
## The more people in the group,  the less likely that a 0 is due to not going fishing
## - the odds decline by 78% (factor of 0.22) for each additional person in the group.


## ----mainEffects, results='markdown', eval=TRUE, hidden=TRUE------------------
emmeans(fish.glm4,  pairwise~CHILD,  type='response')

newdata = emmeans(fish.glm4,  ~CHILD|CAMPER,  type='response') %>%
  as.data.frame
newdata %>% head

ggplot(newdata,  aes(y=response,  x=CHILD,  color=CAMPER)) +
  geom_pointrange(aes(ymin=lower.CL,  ymax=upper.CL)) +
  theme_classic()

fish.grid <- with(fish, list(PERSONS = 1:4))
#emmeans(fish.glm4,  ~PERSONS,  at = fish.grid,  mode='zero',  lin.pred=FALSE,  type='response')
#emmeans(fish.glm4,  ~PERSONS,  at = fish.grid,  mode='prob0')
#emmeans(fish.glm4,  ~PERSONS,  mode='zero',  lin.pred=TRUE)  #link scale
#emmeans(fish.glm4,  ~PERSONS,  mode='prob0')
#https://cran.r-project.org/web/packages/emmeans/vignettes/models.html

