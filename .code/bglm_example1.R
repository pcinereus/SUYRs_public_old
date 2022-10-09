## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE,cache.lazy = FALSE, tidy='styler')


## ----libraries, results='markdown', eval=TRUE---------------------------------
library(tidyverse)  #for data wrangling etc
library(rstanarm)   #for fitting models in STAN
library(cmdstanr)   #for cmdstan
library(brms)       #for fitting models in STAN
library(standist)   #for exploring distributions
library(coda)       #for diagnostics
library(bayesplot)  #for diagnostics
library(ggmcmc)     #for MCMC diagnostics
library(DHARMa)     #for residual diagnostics
library(rstan)      #for interfacing with STAN
library(emmeans)    #for marginal means etc
library(broom)      #for tidying outputs
library(tidybayes)  #for more tidying outputs
library(HDInterval) #for HPD intervals
library(ggeffects)  #for partial plots
library(broom.mixed)#for summarising models
library(posterior)  #for posterior draws
library(ggeffects)  #for partial effects plots
library(patchwork)  #for multi-panel figures
library(bayestestR) #for ROPE
library(see)        #for some plots
theme_set(theme_grey()) #put the default ggplot theme back
source('helperFunctions.R')


## ----readData, results='markdown', eval=TRUE----------------------------------
fert <- read_csv('../data/fertilizer.csv', trim_ws = TRUE)
glimpse(fert)


## ----lm, results='markdown', eval=TRUE, hidden=TRUE---------------------------
summary(fert.lm <- lm(YIELD ~ FERTILIZER, data = fert))
sigma(fert.lm) 
summary(fert.lm <- lm(YIELD ~ scale(FERTILIZER, scale=FALSE), data = fert))
sigma(fert.lm)


## ----fitModel1a, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
fert.rstanarm <- stan_glm(YIELD ~ FERTILIZER,
                          data=fert,
                          iter = 5000,
                          warmup = 1000,
                          chains = 3,
                          thin = 5,
                          refresh = 0)


## ----fitModel1b, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
prior_summary(fert.rstanarm)


## ----fitModel1c, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
mean(fert$YIELD)
sd(fert$YIELD)
2.5*sd(fert$YIELD)


## ----fitModel1d, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
2.5*sd(fert$YIELD)/sd(fert$FERTILIZER)


## ----fitModel1e, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
1/sd(fert$YIELD)


## ----fitModel1f, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
fert.rstanarm1 <- update(fert.rstanarm,  prior_PD=TRUE)


## ----fitModel1g, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
ggpredict(fert.rstanarm1) %>% plot(add.data=TRUE)


## ----fitModel1h1, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, fig.width = 5, fig.height = 3----
standist::visualize("normal(164, 65)", xlim = c(0, 300)) +
standist::visualize("normal(0, 2)", xlim = c(-10, 10))
standist::visualize("student_t(3, 0, 65)",
                    "gamma(2,1)",
                    "exponential(0.016)",
                    xlim = c(-10, 50))

standist::visualize("student_t(3, 0, 65)",
                    "exponential(0.016)",
                    xlim = c(-10, 300))


## ----fitModel1h, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
fert.rstanarm2= stan_glm(YIELD ~ FERTILIZER, data=fert,
                         prior_intercept = normal(164, 65, autoscale=FALSE),
                         prior = normal(0, 2, autoscale=FALSE),
                         prior_aux = student_t(3, 0, 65, autoscale=FALSE),
                         prior_PD=TRUE, 
                         iter = 5000, warmup = 1000,
                         chains = 3, cores = 3,
                         thin = 5, refresh = 0
                         )


## ----fitModel1i, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
ggpredict(fert.rstanarm2) %>%
  plot(add.data=TRUE)


## ----fitModel1j, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
fert.rstanarm3= update(fert.rstanarm2,  prior_PD=FALSE)  


## ----modelFit1k, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
posterior_vs_prior(fert.rstanarm3, color_by='vs', group_by=TRUE,
                   facet_args=list(scales='free_y'))


## ----modelFit1l, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggpredict(fert.rstanarm3) %>% plot(add.data=TRUE)


## ----fitModel2a, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
fert.brm <- brm(bf(YIELD ~ FERTILIZER),
               data = fert,
               iter = 5000,
               warmup = 1000,
               chains = 3,
               thin = 5,
               refresh = 0,
               backend = 'cmdstan')


## ----fitModel2a2, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
fert.brm %>% names


## ----fitModel2b, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80), echo=2----
options(width=100)
prior_summary(fert.brm)
options(width=80)


## ----fitModel2c, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
median(fert$YIELD)
mad(fert$YIELD)


## ----student_t_dist, results='markdown', eval=TRUE----------------------------
standist::visualize("student_t(3,161.5,90.4)", xlim=c(-100,1000))


## ----student_t_dist1, results='markdown', eval=TRUE---------------------------
standist::visualize("student_t(3,0,90.4)", xlim=c(-10,500))


## ----normal_prior, results='markdown', eval=TRUE------------------------------
standist::visualize("normal(0, 2.5)", "student_t(3, 0, 2.5)", xlim = c(-50,50))


## ----fitModel2d, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
fert.brm1 <- brm(bf(YIELD ~ FERTILIZER),
                 data = fert,
                 prior = prior(student_t(3, 0, 2.5), class = 'b'), 
                 sample_prior = 'only', 
                 iter = 5000,
                 warmup = 1000,
                 chains = 3,
                 thin = 5,
                 backend = 'cmdstan',
                 refresh = 0)


## ----fitModel2e1, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
fert.brm1 %>% ggpredict() %>% plot(add.data=TRUE)


## ----fitModel2e2, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
fert.brm1 %>% ggemmeans(~FERTILIZER) %>% plot(add.data=TRUE)


## ----fitModel2e3, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
fert.brm1 %>% conditional_effects() %>%  plot(points=TRUE)


## ----priors, results='markdown', eval=TRUE------------------------------------
standist::visualize("normal(164, 65)", xlim = c(-100, 300))

standist::visualize("normal(0, 2)", xlim = c(-10, 10))

standist::visualize("student_t(3, 0, 65)",
                    "normal(0, 65)",
                    "gamma(2,0.05)",
                    xlim = c(0, 100))
    


## ----fitModel2h, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
priors <- prior(normal(164, 65),  class='Intercept') +
    prior(normal(0, 2.5), class='b') +
    prior(student_t(3, 0, 65),  class='sigma')

fert.brm2 = brm(bf(YIELD ~ FERTILIZER),
                data=fert,
                prior=priors,
                sample_prior = 'only', 
                iter = 5000,
                warmup = 1000,
                chains = 3, cores = 3,
                thin = 5,
                backend = 'cmdstan',
                refresh = 0)


## ----fitModel2i1, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
fert.brm2 %>%
    ggpredict() %>%
    plot(add.data=TRUE)


## ----fitModel2i2, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
fert.brm2 %>%
    ggemmeans(~FERTILIZER) %>%
    plot(add.data = TRUE)


## ----fitModel2i3, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
fert.brm2 %>%
    conditional_effects() %>%
    plot(points = TRUE)


## ----fitModel2j, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
fert.brm3 <- update(fert.brm2, sample_prior = 'yes', refresh = 0)
#OR 
fert.brm3 = brm(bf(YIELD ~ FERTILIZER),
                data=fert,
                prior=priors,
                sample_prior = 'yes', 
                iter = 5000,
                warmup = 1000,
                chains = 3, cores = 3,
                thin = 5,
                backend = 'cmdstan',
                refresh = 0)


## ----fitModel2j2, results='markdown', eval=TRUE, echo = FALSE, hidden=TRUE----
save(fert.brm3, file = '../ws/testing/fert.brm3.RData')


## ----fitModel2k1, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
fert.brm3 %>%
    ggpredict() %>%
    plot(add.data=TRUE)


## ----fitModel2k2, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
fert.brm3 %>%
    ggemmeans(~FERTILIZER) %>%
    plot(add.data = TRUE)


## ----fitModel2k3, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
fert.brm3 %>%
    conditional_effects() %>%
    plot(points = TRUE)


## ----posterior2, results='markdown', eval=TRUE--------------------------------
fert.brm3 %>% get_variables()
## fert.brm3 %>% hypothesis('Intercept=0', class='b') %>% plot
fert.brm3 %>% hypothesis('FERTILIZER = 0') %>% plot
## fert.brm3 %>% hypothesis('scaleFERTILIZERscaleEQFALSE=0') %>% plot
fert.brm3 %>% hypothesis('sigma = 0', class = '') %>% plot


## ----fitModel2k, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, fig.width = 8, fig.height = 4----
fert.brm3 %>% get_variables()
fert.brm3 %>% SUYR_prior_and_posterior()


## ----fitModel2l, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
fert.brm3 %>% standata()
fert.brm3 %>% stancode()


## ----INLApackages, results='markdown', eval=TRUE------------------------------
library(INLA)


## ----fitModel3a, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
fert.inla <- inla(YIELD ~ FERTILIZER,
                  data = fert,
                  family = 'gaussian',
                  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE)
                  )


## ----fitModel3a2, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
fert.inla %>% names


## ----fitModel3b0, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80), echo=TRUE----
inla.priors.used(fert.inla)


## ----fitModel3b3, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80)----
standist::visualize("gamma(1, 0.00005)", xlim=c(-100,100000))


## ----fitModel3b1, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80)----
standist::visualize("normal(0, 31)", xlim=c(-100,100))


## ----fitModel3b5, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80)----
standist::visualize("gamma(0.5, 0.032)", xlim=c(-5,100))


## ----fitModel3b9, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80)----
min(fert$YIELD)
median(fert$YIELD)
mad(fert$YIELD)


## ----fitModel3b8, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80)----
standist::visualize("normal(80, 271.32)", xlim=c(-700,1000))


## ----fitModel3b4, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, paged.print=FALSE,tidy.opts = list(width.cutoff = 80)----
standist::visualize("normal(0, 5.102)", xlim=c(-20,20))


## ----fitModel3b6, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
fert.inla1 <- inla(YIELD ~ FERTILIZER,
                  data = fert,
                  family = 'gaussian',
                  control.fixed = list(
                      mean.intercept = 80,
                      prec.intercept = 0.00001,
                      mean = 0,
                      prec = 0.0384),
                  control.family = list(hyper = list(prec = list(prior = "loggamma",
                                                                param = c(0.5, 0.032)))),
                  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE)
                  ) 


## ----fitModel3b7, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
## Fixed effects
rbind(Model.1 = fert.inla$summary.fixed, Model.2 = fert.inla1$summary.fixed)

## Family hyperpriors (variance)
rbind(Model.1 = fert.inla$summary.hyperpar, Model.2 = fert.inla1$summary.hyperpar)


## ----fitModel3d0, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
plot(fert.inla1,
     plot.fixed.effects = TRUE,
     plot.lincomb = FALSE,
     plot.random.effects = FALSE,
     plot.hyperparameters = TRUE,
     plot.predictor = FALSE,
     plot.q = FALSE,
     plot.cpo = FALSE,
     plot.prior = TRUE,
     single = FALSE)


## ----fitModel3d1, results='markdown', eval=TRUE, hidden=TRUE, fig.width = 6, fig.height = 3, cache=FALSE----
library(INLAutils)
plot_fixed_marginals(fert.inla1, priors = TRUE)

## ----fitModel3d1a, results='markdown', eval=TRUE, hidden=TRUE, fig.width = 3, fig.height = 3, cache=FALSE----
library(INLAutils)
plot_hyper_marginals(fert.inla1, CI = TRUE)


## ----fitModel3d1c, results='markdown', eval=TRUE, hidden=TRUE, fig.width = 6, fig.height = 6, cache=FALSE----
library(INLAutils)
autoplot(fert.inla1)


## ----fitModel3d2, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, fig.width=8----
marg.fix <- fert.inla1$marginals.fixed
(nm <- names(marg.fix))
p <- vector('list', length(nm))
for (wch in 1:length(nm)) {
    xy <- INLA:::inla.get.prior.xy(section = 'fixed', hyperid = nm[1] ,
                                   all.hyper = fert.inla1$all.hyper,
                                   range = 3*range(marg.fix[[wch]][,'x']))
    p[[wch]] <- ggplot() +
        geom_density(data=as.data.frame(marg.fix[[wch]]),
                     aes(x=x, y=y, fill='Posterior'), alpha=0.2,
                     stat='Identity') +
        geom_density(data=as.data.frame(xy), aes(x=x, y=y, fill='Prior'), alpha = 0.2,
                     stat='Identity')
}
patchwork::wrap_plots(p)


## ----modelValidation1a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation1b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
plot(fert.rstanarm3, plotfun='mcmc_trace')


## ----modelValidation1c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
plot(fert.rstanarm3, 'acf_bar')


## ----modelValidation1d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
plot(fert.rstanarm3, 'rhat_hist')


## ----modelValidation1e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
plot(fert.rstanarm3, 'neff_hist')


## ----Validation1f, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
plot(fert.rstanarm3, 'combo')
plot(fert.rstanarm3, 'violin')


## ----modelValidation1g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
stan_trace(fert.rstanarm3)


## ----modelValidation1h, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
stan_ac(fert.rstanarm3) 


## ----modelValidation1i, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
stan_rhat(fert.rstanarm3) 


## ----modelValidation1j, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
stan_ess(fert.rstanarm3)


## ----modelValidation1k, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
stan_dens(fert.rstanarm3, separate_chains = TRUE)


## ----modelValidation1l, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=9----
fert.ggs <- ggs(fert.rstanarm3)
ggs_traceplot(fert.ggs)


## ----modelValidation1m, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggs_autocorrelation(fert.ggs)


## ----modelValidation1n, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggs_Rhat(fert.ggs)


## ----modelValidation1o, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggs_effective(fert.ggs)


## ----modelValidation1p, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggs_crosscorrelation(fert.ggs)


## ----modelValidation1q, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggs_grb(fert.ggs)


## ----modelValidation2a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation2b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.brm3 %>% mcmc_plot(type = 'combo')
fert.brm3 %>% mcmc_plot(type = 'trace')
fert.brm3 %>% mcmc_plot(type = 'dens_overlay')


## ----modelValidation2c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.brm3 %>% mcmc_plot(type = 'acf_bar')


## ----modelValidation2d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.brm3 %>% mcmc_plot(type = 'rhat_hist')


## ----modelValidation2e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.brm3 %>% mcmc_plot(type = 'neff_hist')


## ----modelValidation2f, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.brm3 %>% mcmc_plot(type = 'combo')
fert.brm3 %>% mcmc_plot(type = 'violin')


## ----modelValidation2g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
stan_trace(fert.brm3$fit)
fert.brm3$fit %>% stan_trace()
fert.brm3$fit %>% stan_trace(inc_warmup = TRUE)


## ----modelValidation2h, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.brm3$fit %>% stan_ac() 


## ----modelValidation2i, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.brm3$fit %>% stan_rhat() 


## ----modelValidation2j, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.brm3$fit %>% stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.brm3$fit %>% stan_dens(separate_chains = TRUE)


## ----modelValidation2l, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=7----
fert.ggs <- fert.brm3$fit %>% ggs(inc_warmup = TRUE, burnin = TRUE)  # does not seem to ignore burnin?
fert.ggs %>% ggs_traceplot()
fert.ggs <- fert.brm3$fit %>% ggs(inc_warmup = FALSE, burnin = TRUE)  # does not seem to ignore burnin?
fert.ggs %>% ggs_traceplot()


## ----modelValidation2m, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=7----
fert.ggs %>% ggs_autocorrelation()


## ----modelValidation2n, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.ggs %>% ggs_Rhat()


## ----modelValidation2o, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.ggs %>% ggs_effective()


## ----modelValidation2p, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.ggs %>% ggs_crosscorrelation()


## ----modelValidation2q, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.ggs %>% ggs_grb()


## ----fitModel3d1b, results='markdown', eval=TRUE, hidden=TRUE, fig.width = 6, fig.height = 6, cache=FALSE----
fert.inla1.cpo <- data.frame(cpo = fert.inla1$cpo$cpo,
           pit = fert.inla1$cpo$pit,
           failure = fert.inla1$cpo$failure) %>%
    filter(failure == 0) %>%
    mutate(row = 1:n())

g1 <- fert.inla1.cpo %>%
    ggplot(aes(y = cpo, x = row)) +
    geom_point()
g2 <- fert.inla1.cpo %>%
    ggplot(aes(x = cpo)) +
    geom_histogram()
g3 <- fert.inla1.cpo %>%
    ggplot(aes(y = pit, x = row)) +
    geom_point()
g4 <- fert.inla1.cpo %>%
    ggplot(aes(x = pit)) +
    geom_histogram()

(g1 + g2)/(g3 + g4)


## ----modelValidation3a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation3b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
pp_check(fert.rstanarm3,  plotfun='dens_overlay')


## ----modelValidation3c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
pp_check(fert.rstanarm3, plotfun='error_scatter_avg')


## ----modelValidation3d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
pp_check(fert.rstanarm3, x=fert$FERTILIZER, plotfun='error_scatter_avg_vs_x')


## ----modelValidation3e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
pp_check(fert.rstanarm3, x=fert$FERTILIZER, plotfun='intervals')


## ----modelValidation3f, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
pp_check(fert.rstanarm3, x=fert$FERTILIZER,plotfun='ribbon')


## ----modelValidation3g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(fert.rstanarm3)


## ----modelValidation4a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
preds <- posterior_predict(fert.rstanarm3,  ndraws=250,  summary=FALSE)
fert.resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = fert$YIELD,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = 'gaussian')
plot(fert.resids)


## ----modelValidation5a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.brm3 %>% pp_check( type='dens_overlay', ndraws=100)


## ----modelValidation5c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.brm3 %>% pp_check(type='error_scatter_avg')


## ----modelValidation5d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.brm3 %>% pp_check(x='FERTILIZER', type='error_scatter_avg_vs_x')


## ----modelValidation5e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
fert.brm3 %>% pp_check(x='FERTILIZER', type='intervals')


## ----modelValidation5f, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
pp_check(fert.brm3, x='FERTILIZER',type='ribbon')


## ----modelValidation5g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(fert.brm3)


## ----modelValidation6a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
preds <- fert.brm3 %>% posterior_predict(ndraws=250,  summary=FALSE)
fert.resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = fert$YIELD,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = FALSE)
fert.resids %>% plot()


## ----fitModel3e0, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
fert.inla1$cpo$cpo
sum(fert.inla1$cpo$cop)
-mean(log(fert.inla1$cpo$cpo))


## ----fitModel3e1, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
fert.inla1$cpo$failure


## ----fitModel3e3, results='markdown', eval=TRUE, hidden=FALSE, cache=FALSE----
plot(fert.inla1,
     plot.fixed.effects = FALSE,
     plot.lincomb = FALSE,
     plot.random.effects = FALSE,
     plot.hyperparameters = FALSE,
     plot.predictor = FALSE,
     plot.q = FALSE,
     plot.cpo = TRUE,
     plot.prior = FALSE) 


## ----fitModel3d1b, results='markdown', eval=TRUE, hidden=TRUE, fig.width = 6, fig.height = 6, cache=FALSE----
fert.inla1.cpo <- data.frame(cpo = fert.inla1$cpo$cpo,
           pit = fert.inla1$cpo$pit,
           failure = fert.inla1$cpo$failure) %>%
    filter(failure == 0) %>%
    mutate(row = 1:n())

g1 <- fert.inla1.cpo %>%
    ggplot(aes(y = cpo, x = row)) +
    geom_point()
g2 <- fert.inla1.cpo %>%
    ggplot(aes(x = cpo)) +
    geom_histogram()
g3 <- fert.inla1.cpo %>%
    ggplot(aes(y = pit, x = row)) +
    geom_point()
g4 <- fert.inla1.cpo %>%
    ggplot(aes(x = pit)) +
    geom_histogram()

(g1 + g2)/(g3 + g4)


## ----fitModel3e4, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
ggplot_inla_residuals(fert.inla1, observed = fert$YIELD)


## ----fitModel3e5, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE-----
ggplot_inla_residuals2(fert.inla1, observed = fert$YIELD)


## ----fitModel3f, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, fig.width = 8----
## draw 250 samples from the posterior
draws <- inla.posterior.sample(n=250, result = fert.inla)
## extract the latent predictions for the observed data
preds = t(sapply(draws, function(x) x$latent[1:nrow(fert)]))
## extract the first (family precision) hyperprior
preds.hyper <- sapply(draws, function(x) 1/sqrt(x$hyperpar[[1]]))
## use this to generate gaussian noise 
preds.hyper <- rnorm(length(preds.hyper), mean=0, sd=preds.hyper)
## add the noise to each prediction
preds <- sweep(preds, 1, preds.hyper, FUN = "+")
## generate the DHARMa residuals
fert.resids <- createDHARMa(simulatedResponse = t(preds),
                           observedResponse = fert$YIELD,
                           fittedPredictedResponse = apply(preds, 2, median),
                           integerResponse = FALSE)
fert.resids %>% plot()


## ----partialPlot1a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% ggpredict() %>% plot(add.data=TRUE)


## ----partialPlot1b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% ggemmeans(~FERTILIZER) %>% plot(add.data=TRUE)


## ----partialPlot1c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% fitted_draws(newdata=fert) %>%
 median_hdci() %>%
 ggplot(aes(x=FERTILIZER, y=.value)) +
 geom_ribbon(aes(ymin=.lower, ymax=.upper), fill='blue', alpha=0.3) + 
 geom_line()


## ----partialPlot2d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% conditional_effects() 
fert.brm3 %>% conditional_effects() %>% plot(points = TRUE)
fert.brm3 %>%
    conditional_effects(spaghetti = TRUE, ndraws = 200) %>%
    plot(points = TRUE) 


## ----partialPlot2a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% ggpredict() %>% plot(add.data=TRUE)


## ----partialPlot2b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% ggemmeans(~FERTILIZER) %>% plot(add.data=TRUE)


## ----partialPlot2c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>%
    fitted_draws(newdata = fert) %>%
    median_hdci() %>%
    ggplot(aes(x = FERTILIZER, y=.value)) +
    geom_ribbon(aes(ymin = .lower, ymax = .upper), fill = 'blue', alpha = 0.3) + 
    geom_line()


## ----partialPlot3a, results='markdown', eval=TRUE-----------------------------
plot(fert.inla1,
     plot.fixed.effects = FALSE,
     plot.lincomb = FALSE,
     plot.random.effects = FALSE,
     plot.hyperparameters = FALSE,
     plot.predictor =TRUE,
     plot.q = FALSE,
     plot.cpo = FALSE,
     plot.prior = FALSE)


## ----partialPlot3a1, results='markdown', eval=TRUE----------------------------
newdata <- fert %>% tidyr::expand(FERTILIZER = modelr::seq_range(FERTILIZER, n=100))
Xmat <- model.matrix(~FERTILIZER, newdata)

nms <- colnames(fert.inla1$model.matrix)
n <- sapply(nms, function(x) 0, simplify=FALSE)
draws <- inla.posterior.sample(n=250, result = fert.inla1, selection = n)
coefs <- t(sapply(draws, function(x) x$latent))
Fit = coefs %*% t(Xmat) %>%
    as.data.frame() %>%
    mutate(Rep=1:n()) %>%
    pivot_longer(cols=-Rep) %>%
    group_by(name) %>%
    median_hdci(value) %>%
    ungroup() %>%
    mutate(name = as.integer(as.character(name))) %>%
    arrange(name) 
newdata <- newdata %>%
    bind_cols(Fit)

newdata %>%
    ggplot(aes(x=FERTILIZER)) +
    geom_ribbon(aes(ymin=.lower, ymax=.upper), fill='orange', color=NA, alpha=0.3) +
    geom_line(aes(y=value), color='orange') +
    geom_point(data=fert, aes(y=YIELD)) +
    theme_classic()


## ----summariseModel1a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% summary()


## ----summariseModel1a1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
fert.sum <- summary(fert.rstanarm3)


## ----summariseModel1b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tidyMCMC(fert.rstanarm3$stanfit, estimate.method='median',  conf.int=TRUE,
         conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)

## ----summariseModel1b1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
fert.tidy <- tidyMCMC(fert.rstanarm3$stanfit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)


## ----summariseModel1dd, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3$stanfit %>%
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel1d2, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3$stanfit %>%
    summarise_draws(median,
                    ~HDInterval::hdi(.x, credMass = 0.9),
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel1m, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3$stanfit %>% as_draws_df()

fert.rstanarm3$stanfit %>%
    as_draws_df() %>%
    summarise_draws(median,
                    ~ HDInterval::hdi(.x),
                    rhat,
                    ess_bulk, ess_tail)


## ----summariseModel1d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% tidy_draws()


## ----summariseModel1c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.draw <- fert.rstanarm3 %>% gather_draws(`(Intercept)`, FERTILIZER, sigma)
## OR via regex
fert.draw <- fert.rstanarm3 %>% gather_draws(`.Intercept.*|FERT.*|sigma`,  regex=TRUE)
fert.draw


## ----summariseModel1c1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.draw %>% median_hdci


## ----summariseModel1c3, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
fert.gather <- fert.rstanarm3 %>% gather_draws(`(Intercept)`,FERTILIZER,sigma) %>%
  median_hdci


## ----summariseModel1c4, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
fert.rstanarm3 %>% 
  gather_draws(`(Intercept)`, FERTILIZER, sigma) %>% 
  ggplot() + 
  stat_halfeye(aes(x=.value,  y=.variable)) +
  facet_wrap(~.variable, scales='free')


## ----summariseModel1j, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% plot(plotfun='mcmc_intervals') 


## ----summariseModel1e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% spread_draws(`(Intercept)`, FERTILIZER, sigma)
# OR via regex
fert.rstanarm3 %>% spread_draws(`.Intercept.*|FERT.*|sigma`,  regex=TRUE)


## ----summariseModel1f, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% posterior_samples() %>% as_tibble()


## ----summariseModel1g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% bayes_R2() %>% median_hdci


## ----summariseModel1h, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
mcmcpvalue <- function(samp)
{
    ## elementary version that creates an empirical p-value for the
    ## hypothesis that the columns of samp have mean zero versus a
    ## general multivariate distribution with elliptical contours.
    
    ## differences from the mean standardized by the observed
    ## variance-covariance factor
    
    ## Note, I put in the bit for single terms
    if (length(dim(samp))==0) {
        std <- backsolve(chol(var(samp)),cbind(0, t(samp)) - mean(samp),transpose = TRUE)
        sqdist <- colSums(std * std)
        sum(sqdist[-1] > sqdist[1])/length(samp)
    }
    else {
        std <- backsolve(chol(var(samp)),cbind(0, t(samp)) - colMeans(samp),transpose = TRUE)
        sqdist <- colSums(std * std)
        sum(sqdist[-1] > sqdist[1])/nrow(samp)
    }
    
}

mcmcpvalue(as.matrix(fert.rstanarm3)[, c("FERTILIZER")])


## ----summariseModel2a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
fert.sum <- summary(fert.brm3)


## ----summariseModel2b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3$fit %>%
    tidyMCMC(estimate.method = 'median',
             conf.int = TRUE,
             conf.method = 'HPDinterval',
             rhat = TRUE, ess = TRUE)

## ----summariseModel2b1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
fert.tidy <- tidyMCMC(fert.brm3$fit, estimate.method='median',  conf.int=TRUE,  conf.method='HPDinterval',  rhat=TRUE, ess=TRUE)


## ----summariseModel2dd, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>%
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel2d2, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>%
    summarise_draws(median,
                    ~HDInterval::hdi(.x, credMass = 0.9),
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel2m, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% as_draws_df()


## ----summariseModel2d1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% tidy_draws()

fert.brm3 %>% get_variables()
fert.brm3$fit %>%
    tidy_draws() %>%
    dplyr::select(matches('^b_.*|^sigma'))  %>%
    summarise_draws(median,
                    HDInterval::hdi,
                    rhat, length, ess_bulk, ess_tail)


## ----summariseModel2d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>%
    tidy_draws() %>%
    gather_variables() %>%
    median_hdci()


## ----summariseModel2c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% get_variables()
fert.brm3 %>%
    gather_draws(b_Intercept, b_FERTILIZER, sigma) %>%
    median_hdci
## OR via regex
fert.brm3 %>%
    gather_draws(`b_.*|sigma`,  regex=TRUE) %>%
    median_hdci

## ----summariseModel2c1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
fert.gather <- fert.brm3 %>%
    gather_draws(b_Intercept, b_FERTILIZER, sigma) %>%
    median_hdci


## ----summariseModel2c2, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
fert.brm3 %>% 
  gather_draws(b_Intercept, b_FERTILIZER, sigma) %>% 
  ggplot() + 
  stat_halfeye(aes(x = .value,  y = .variable)) +
  facet_wrap(~.variable, scales = 'free')

fert.brm3 %>% 
  gather_draws(b_Intercept, b_FERTILIZER, sigma) %>% 
  ggplot() + 
  stat_halfeye(aes(x = .value,  y = .variable)) 


## ----summariseModel2e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% spread_draws(b_Intercept, b_FERTILIZER, sigma)
# OR via regex
fert.brm3 %>% spread_draws(`b_.*|sigma`,  regex=TRUE)


## ----summariseModel2j, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% mcmc_plot(type='intervals') 


## ----summariseModel2f, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% posterior_samples() %>% as_tibble()


## ----summariseModel2g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% bayes_R2(summary=FALSE) %>% median_hdci


## ----summariseModel2h, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
mcmcpvalue <- function(samp)
{
    ## elementary version that creates an empirical p-value for the
    ## hypothesis that the columns of samp have mean zero versus a
    ## general multivariate distribution with elliptical contours.
    
    ## differences from the mean standardized by the observed
    ## variance-covariance factor
    
    ## Note, I put in the bit for single terms
    if (length(dim(samp))==0) {
        std <- backsolve(chol(var(samp)),cbind(0, t(samp)) - mean(samp),transpose = TRUE)
        sqdist <- colSums(std * std)
        sum(sqdist[-1] > sqdist[1])/length(samp)
    }
    else {
        std <- backsolve(chol(var(samp)),cbind(0, t(samp)) - colMeans(samp),transpose = TRUE)
        sqdist <- colSums(std * std)
        sum(sqdist[-1] > sqdist[1])/nrow(samp)
    }
    
}

mcmcpvalue(as.matrix(fert.brm3)[, c("b_FERTILIZER")])


## ----summariseModel3a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.inla1 %>% summary()


## ----summariseModel3a1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
fert.inla1.sum <- summary(fert.inla1)


## ----posteriors3a1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=3----
brinla::bri.fixed.plot(fert.inla1)


## ----posteriors3a2, results='markdown', eval=TRUE, hidden=TRUE, fig.width=4, fig.height=3----
brinla::bri.hyperpar.plot(fert.inla1)


## ----posteriors3b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=4, fig.height=3----
## Get all the draws as a list
n <- 1000
draws <- inla.posterior.sample(n = n, result = fert.inla1, seed = 1)
## Redimension the list for the latent (fixed and fitted) values 
(latent.names <- draws[[1]]$latent %>% rownames)


## ----posteriors3c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=4, fig.height=3----
draws <- sapply(draws, function(x) x[['latent']]) %>%
    as.data.frame() %>%
    set_names(paste0('Draw',1:n)) %>%
    mutate(Parameter = latent.names) %>%
    dplyr::select(Parameter, everything())
draws


## ----posteriors3d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=4, fig.height=3----
fert.draws <- draws %>% pivot_longer(cols =- Parameter, names_to='Draws')
fert.draws


## ----posteriors3e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=4, fig.height=3----
fert.draws %>%
    filter(grepl('Intercept|FERTILIZER', Parameter)) %>%
    group_by(Parameter) %>%
    median_hdci(value)

fert.draws %>%
    filter(grepl('Intercept|FERTILIZER', Parameter)) %>%
    ggplot(aes(x=value, y= Parameter)) +
    geom_vline(xintercept = 0, linetype = 'dashed') + 
    stat_pointinterval(point_interval = median_hdci)
## or separate into panels
fert.draws %>%
    filter(grepl('Intercept|FERTILIZER', Parameter)) %>%
    ggplot(aes(x=value)) +
    geom_vline(xintercept = 0) +
    stat_pointinterval(point_interval = median_hdci) +
    facet_grid(~Parameter, scales = 'free_x')


## ----posteriors3f, results='markdown', eval=TRUE, hidden=TRUE, fig.width=4, fig.height=3----
fert.draws %>%
    filter(grepl('Intercept|FERTILIZER', Parameter)) %>%
    group_by(Parameter) %>%
    nest() %>%
    mutate(hpd = map(data, ~median_hdci(.x$value)),
           Density = map(data, function(.x) {d <- density(.x$value); data.frame(x=d$x, y=d$y)}),
           Quant = map2(Density, hpd, ~factor(findInterval(.x$x, .y[,c('ymin','ymax')])))
           ) %>%
    unnest(c(Density, Quant)) %>%
    ggplot(aes(x = x, y=y, fill = Quant)) +
    geom_vline(xintercept = 0) +
    geom_ribbon(aes(ymin=0, ymax = y)) +
    facet_wrap(~Parameter, scales = "free") 


## ----posteriors3g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=4, fig.height=3----
fert.draws %>%
    filter(grepl('Intercept|FERTILIZER', Parameter)) %>%
    ggplot(aes(x=value, y = Parameter, fill = stat(quantile))) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    ggridges::stat_density_ridges(geom = "density_ridges_gradient",
                                  calc_ecdf = TRUE,
                                  quantile_fun = function(x, probs) quantile(x,probs),
                                  quantiles = c(0.025, 0.975)) +
    scale_fill_viridis_d()
fert.draws %>%
    filter(grepl('Intercept|FERTILIZER', Parameter)) %>%
    ggplot(aes(x=value, y = Parameter, fill = stat(quantile))) +
    geom_vline(xintercept = 0, linetype = 'dashed') +
    ggridges::stat_density_ridges(geom = "density_ridges_gradient",
                                  calc_ecdf = TRUE,
                                  quantile_fun = function(x, probs) quantile(x,probs),
                                  quantiles = c(0.025, 0.975)) +
    scale_fill_viridis_d() +
    facet_grid(~Parameter, scales = 'free_x')


## ----posteriors3h, results='markdown', eval=TRUE, hidden=TRUE, fig.width=4, fig.height=3----
cor(fert.inla1$summary.fitted.values[,'mean'],
    fert$YIELD)^2


## ----predictions0a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
## establish a data set that defines the new data to predict against
newdata = data.frame(FERTILIZER = 110)


## ----predictions1b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% emmeans(~FERTILIZER,  at=newdata)


## ----predictions1a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% tidybayes::predicted_draws(newdata) %>% median_hdci()
## or
newdata %>% tidybayes::add_predicted_draws(fert.rstanarm3) %>% median_hdci
## or
fert.rstanarm3 %>% posterior_predict(newdata=newdata) %>%
  tidyMCMC(estimate.method='median', conf.int=TRUE, conf.method='HPDinterval')


## ----predictions1c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% epred_draws(newdata) %>% median_hdci()
## or
newdata %>% add_epred_draws(fert.rstanarm3) %>% median_hdci
## or
fert.rstanarm3 %>% posterior_epred(newdata=newdata) %>%
  tidyMCMC(estimate.method='median', conf.int=TRUE, conf.method='HPDinterval')


## ----predictions1d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% linpred_draws(newdata) %>% median_hdci()
## or
newdata %>% add_linpred_draws(fert.rstanarm3) %>% median_hdci
## or
fert.rstanarm3 %>% posterior_linpred(newdata=newdata) %>%
  tidyMCMC(estimate.method='median', conf.int=TRUE, conf.method='HPDinterval')


## ----predictions1h, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
## Establish an appropriate model matrix
Xmat = model.matrix(~FERTILIZER, data=newdata)
### get the posterior draws for the linear predictor
coefs <- fert.rstanarm3$stanfit %>%
    as_draws_df() %>%
    dplyr::select(c('(Intercept)','FERTILIZER')) %>%
    as.matrix()
fit <- coefs %*% t(Xmat)
fit %>% median_hdci
## or
coefs <- fert.rstanarm3 %>%
    tidy_draws() %>% 
    dplyr::select(c('(Intercept)','FERTILIZER')) %>%
    as.matrix()
fit <- coefs %*% t(Xmat)
fit %>% median_hdci
## or
coefs <- fert.rstanarm3$stanfit %>%
    as_draws_matrix() %>% 
    subset_draws(c("(Intercept)", "FERTILIZER"))
fit <- coefs %*% t(Xmat)
fit %>% median_hdci


## ----predictions2b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>%
    emmeans(~FERTILIZER, at = newdata)
## OR
fert.brm3 %>%
    emmeans(~FERTILIZER, at = newdata) %>%
    tidy_draws() %>%
    median_hdci()


## ----predictions2c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% tidybayes::predicted_draws(newdata) %>% median_hdci()
## or
newdata %>% tidybayes::add_predicted_draws(fert.brm3) %>% median_hdci
## or
fert.brm3 %>% posterior_predict(newdata=newdata) %>%
    as.mcmc() %>%
  tidyMCMC(estimate.method='median', conf.int=TRUE, conf.method='HPDinterval')


## ----predictions2d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% epred_draws(newdata) %>% median_hdci()
## or
newdata %>% add_epred_draws(fert.brm3) %>% median_hdci
## or
fert.brm3 %>% posterior_epred(newdata=newdata) %>%
    as.mcmc() %>%
  tidyMCMC(estimate.method='median', conf.int=TRUE, conf.method='HPDinterval')


## ----predictions2e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% linpred_draws(newdata) %>% median_hdci()
## or
newdata %>% add_linpred_draws(fert.brm3) %>% median_hdci
## or
fert.brm3 %>% posterior_linpred(newdata=newdata) %>%
    as.mcmc() %>% 
  tidyMCMC(estimate.method='median', conf.int=TRUE, conf.method='HPDinterval')


## ----predictions2h, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
## Establish an appropriate model matrix
Xmat = model.matrix(~FERTILIZER, data=newdata)
### get the posterior draws for the linear predictor
coefs <- fert.brm3 %>%
    as_draws_df() %>%
    dplyr::select(c('b_Intercept','b_FERTILIZER')) %>%
    as.matrix()
fit <- coefs %*% t(Xmat)
fit %>% median_hdci %>% bind_cols(newdata)
## or
coefs <- fert.brm3 %>%
    tidy_draws() %>% 
    dplyr::select(c('b_Intercept','b_FERTILIZER')) %>%
    as.matrix()
fit <- coefs %*% t(Xmat)
fit %>% median_hdci %>% bind_cols(newdata)
## or
coefs <- fert.brm3 %>%
    as_draws_matrix() %>% 
    subset_draws(c("b_Intercept", "b_FERTILIZER"))
fit <- coefs %*% t(Xmat)
fit %>% median_hdci %>% bind_cols(newdata)
## or
## Establish an appropriate model matrix
Xmat = model.matrix(~FERTILIZER, data=newdata)
### get the posterior draws for the linear predictor
coefs <- fert.brm3 %>% posterior_samples(pars=c('b_Intercept','b_FERTILIZER')) %>% as.matrix()
fit <- coefs %*% t(Xmat)
fit %>% median_hdci %>% bind_cols(newdata)


## ----predictions3a, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
## Expected values
Xmat <- model.matrix(~FERTILIZER, newdata)
nms <- colnames(fert.inla1$model.matrix)
sel <- sapply(nms, function(x) 0, simplify = FALSE)
n <- 1000
draws <- inla.posterior.sample(n = n, result = fert.inla1, selection = sel, seed = 1)
coefs <- t(sapply(draws, function(x) x$latent))
Fit <- coefs %*% t(Xmat) %>%
  as.data.frame() %>%
  mutate(Rep = 1:n()) %>%
  pivot_longer(cols = -Rep) %>%
  group_by(name) %>%
  median_hdci(value) %>%
  ungroup() %>%
  mutate(name = as.integer(as.character(name))) %>%
  arrange(name)
newdata.inla <- newdata %>%
  bind_cols(Fit)
newdata.inla

## Predictions
Fit <- coefs %*% t(Xmat) 
draws <- inla.posterior.sample(n = n, result = fert.inla1, seed = 1)
sigma <- sqrt(1/(sapply(draws, function(x) x$hyperpar)))
sigma <- sapply(sigma, function(x) rnorm(1,0,sigma))
Fit <- sweep(Fit, MARGIN = 1, sigma, FUN = "+") %>%
    as.data.frame() %>%
    mutate(Rep = 1:n()) %>%
    pivot_longer(cols = -Rep) %>%
    group_by(name) %>%
    median_hdci(value) %>%
    bind_cols(newdata) %>% 
    ungroup() %>%
    dplyr::select(FERTILIZER, everything(), -name)
Fit

## or
fun <- function(coefs = NA) {
    ## theta[1] is the precision
    return (Intercept + FERTILIZER * coefs[,'FERTILIZER'] +
            rnorm(nrow(coefs), sd = sqrt(1/theta[1])))
}
Fit <- inla.posterior.sample.eval(fun, draws, coefs = newdata) %>%
    as.data.frame() %>%
    bind_cols(newdata) %>%
    pivot_longer(cols = -FERTILIZER) %>%
    group_by(FERTILIZER) %>%
    median_hdci(value)
    
Fit



## ----predictions3b, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
fert.pred <- fert %>%
    bind_rows(newdata)
i.newdata <- (nrow(fert) +1):nrow(fert.pred)
fert.inla2 <- inla(YIELD ~ FERTILIZER,
  data = fert.pred,
  family = "gaussian",
  control.fixed = list(
    mean.intercept = 80,
    prec.intercept = 0.00001,
    mean = 0,
    prec = 0.0384
  ),
  control.family = list(hyper = list(prec = list(
    prior = "loggamma",
    param = c(0.5, 0.31)
  ))),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE)
)

fert.inla2$summary.fitted.values[i.newdata,]
## or on the link scale...
fert.inla2$summary.linear.predictor[i.newdata,]


## ----predictions3c, results='markdown', eval=TRUE, fig.width=8, fig.height=5----
Xmat <- model.matrix(~FERTILIZER, data=newdata)
lincomb <- inla.make.lincombs(as.data.frame(Xmat))

fert.inla3 <- inla(YIELD ~ FERTILIZER,
  data = fert,
  family = "gaussian",
  lincomb = lincomb,
  control.fixed = list(
    mean.intercept = 80,
    prec.intercept = 0.00001,
    mean = 0,
    prec = 0.0384
  ),
  control.family = list(hyper = list(prec = list(
    prior = "loggamma",
    param = c(0.5, 0.31)
  ))),
  control.compute = list(config = TRUE, dic = TRUE, waic = TRUE, cpo = TRUE)
)

fert.inla3$summary.lincomb.derived


## ----Probability1a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% hypothesis('FERTILIZER>0')


## ----Probability1b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,paged.print=FALSE----
fert.rstanarm3 %>% tidy_draws() %>% summarise(P=sum(FERTILIZER>0)/n())


## ----Probability1bb, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,paged.print=FALSE----
0.1 * sd(fert$YIELD)/sd(fert$FERTILIZER)
## Cannot pipe to rope if want to pipe rope to plot()
bayestestR::rope(fert.rstanarm3, parameters = 'FERTILIZER', range = c(-0.08, 0.08), ci_method = "HDI") 
bayestestR::rope(fert.rstanarm3, parameters = 'FERTILIZER', range = c(-0.08, 0.08), ci_method = "HDI") %>% plot()

bayestestR::equivalence_test(fert.rstanarm3, parameters = 'FERTILIZER', range = c(-0.08, 0.08), ci_method = "HDI")  
bayestestR::equivalence_test(fert.rstanarm3, parameters = 'FERTILIZER', range = c(-0.08, 0.08), ci_method = "HDI") %>% plot()



## ----Probability1c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5, echo=1:2----
newdata <- list(FERTILIZER=c(200, 100)) 
fert.rstanarm3 %>% emmeans(~FERTILIZER,  at=newdata) %>% pairs()
fert.mcmc <- fert.rstanarm3 %>% emmeans(~FERTILIZER,  at=newdata) %>% pairs() %>%
  as.data.frame()


## ----Probability1c2, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5, echo=1:2----
newdata <- list(FERTILIZER=c(200, 100)) 
fert.rstanarm3 %>% emmeans(~FERTILIZER,  at=newdata) %>% regrid(transform = "log") %>% pairs() %>% regrid()


## ----Probability1d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>% emmeans(~FERTILIZER,  at=newdata) %>%
    regrid(transform = "log") %>%
    pairs() %>%
    regrid() %>%
    tidy_draws() %>%
    summarise_draws(median,
                    HDInterval::hdi,
                    P= ~ sum(.x>1.5)/length(.x))

fert.mcmc <- fert.rstanarm3 %>% emmeans(~FERTILIZER,  at=newdata) %>% 
  tidy_draws() %>%
  rename_with(~str_replace(., 'FERTILIZER ', 'p')) %>%
  mutate(Eff=p200 - p100,
         PEff=100*Eff/p100)
fert.mcmc %>% head


## ----Probability1e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.mcmc %>% tidyMCMC(estimate.method='median',
                       conf.int=TRUE, conf.method='HPDinterval')


## ----Probability1f, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.mcmc %>% median_hdci(PEff)
fert.mcmc %>% median_hdci(PEff, Eff)


## ----Probability1g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5, paged.print=FALSE----
fert.mcmc %>% summarise(P=sum(PEff>50)/n())


## ----Probability1h, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.mcmc %>% hypothesis('PEff>50')


## ----Probability1i, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
 newdata = list(FERTILIZER=c(200, 100)) 
 fert.rstanarm3 %>%
     emmeans(~FERTILIZER,  at=newdata) %>%
     pairs() %>%
     tidy_draws() %>%
     summarise(across(contains('contrast'),
                      list(P = ~ sum(.>50)/n(),
                           HDCI = ~ median_hdci(.)),
                      .names = c('{.fn}')
                      )) %>%
     tidyr::unpack(HDCI)


## ----Probability1j, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
 newdata = list(FERTILIZER=c(200, 100)) 
 ## Simple
 fert.rstanarm3 %>%
     emmeans(~FERTILIZER,  at=newdata) %>%
     regrid(transform = 'log') %>%
     pairs() %>%
     regrid() 

 ## More advanced (both P and percent change)
 fert.mcmc <- fert.rstanarm3 %>%
     emmeans(~FERTILIZER,  at=newdata) %>%
     regrid(transform = 'log') %>%
     pairs() %>%
     regrid() %>%
     tidy_draws() %>%
     mutate(across(contains('contrast'), ~ 100*(. - 1)))

 fert.mcmc %>%
     summarise(across(contains('contrast'),
                      list(P = ~ sum(.>50)/n(),
                           HDCI = ~ median_hdci(.)),
                      .names = c('{.fn}')
                      )) %>%
     tidyr::unpack(HDCI)

 ## OR
 fert.mcmc %>%
     summarise_draws(median,
                     HDInterval::hdi,
                     P = ~ sum(.x>50)/length(.x))


## ----Probability1k, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.rstanarm3 %>%
     linpred_draws(newdata = as.data.frame(newdata)) %>%
     ungroup() %>%
     group_by(.draw) %>%
     summarise(Eff = -1 * diff(.linpred),
               PEff = 100*Eff/.linpred[2]) %>%
     ungroup() %>%
     mutate(P = sum(PEff>50)/n()) %>%
     pivot_longer(cols = -.draw) %>%
     group_by(name) %>%
     median_hdci()

 ##OR
 fert.rstanarm3 %>%
     epred_draws(newdata = as.data.frame(newdata)) %>%
     ungroup() %>%
     group_by(.draw) %>%
     summarise(Eff = -1 * diff(.epred),
               PEff = 100*Eff/.epred[2]) %>%
     ungroup() %>%
     mutate(P = sum(PEff>50)/n()) %>%
     pivot_longer(cols = -.draw) %>%
     group_by(name) %>%
     median_hdci()

 ##OR for prediction of individual values
 fert.rstanarm3 %>%
     predicted_draws(newdata = as.data.frame(newdata)) %>%
     ungroup() %>%
     group_by(.draw) %>%
     summarise(Eff = -1 * diff(.prediction),
               PEff = 100*Eff/.prediction[2]) %>%
     ungroup() %>%
     mutate(P = sum(PEff>50)/n()) %>%
     pivot_longer(cols = -.draw) %>%
     group_by(name) %>%
     median_hdci()


## ----Probability2a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% hypothesis('FERTILIZER > 0')
fert.brm3 %>% hypothesis('FERTILIZER > 0') %>% plot (ignore_prior= TRUE)
fert.brm3 %>% hypothesis('FERTILIZER > 0') %>% `[[`('samples') %>% median_hdci()


## ----Probability2b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,paged.print=FALSE----
fert.brm3 %>%
    gather_draws(b_FERTILIZER) %>%
    summarise(P=sum(.value>0)/n())
#OR
fert.brm3 %>%
    tidy_draws() %>%
    summarise(P=sum(b_FERTILIZER>0)/n())


## ----Probability2c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,paged.print=FALSE----
fert.brm3 %>% hypothesis('FERTILIZER>0.9')
fert.brm3 %>%
    hypothesis('FERTILIZER>0.9') %>%
    plot(ignore_prior = TRUE)
# This returns a list
fert.brm3 %>%
    hypothesis('FERTILIZER>0.9') %>%
    plot(ignore_prior = TRUE) %>%
    `[[`(1) +
    geom_vline(xintercept=0, linetype='dashed')

#OR
fert.brm3 %>%
    tidy_draws() %>%
    summarise(P=sum(b_FERTILIZER>0.7)/n())
fert.brm3 %>%
    tidy_draws() %>%
    ggplot(aes(x=b_FERTILIZER)) +
    geom_density(fill='orange') +
    geom_vline(xintercept=0.7, linetype='dashed')


## ----Probability2bb, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,paged.print=FALSE----
0.1 * sd(fert$YIELD)/sd(fert$FERTILIZER)
## Cannot pipe to rope if want to pipe rope to plot()
bayestestR::rope(fert.brm3, parameters = 'FERTILIZER', range = c(-0.08, 0.08), ci_method = "HDI") 
bayestestR::rope(fert.brm3, parameters = 'FERTILIZER', range = c(-0.08, 0.08), ci_method = "HDI") %>% plot()

bayestestR::equivalence_test(fert.brm3, parameters = 'FERTILIZER', range = c(-0.08, 0.08), ci_method = "HDI")  
bayestestR::equivalence_test(fert.brm3, parameters = 'FERTILIZER', range = c(-0.08, 0.08), ci_method = "HDI") %>% plot()



## ----Probability2c1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5, echo=1:2----
newdata <- list(FERTILIZER = c(200, 100)) 
fert.brm3 %>%
    emmeans(~FERTILIZER,  at = newdata) %>%
    pairs()
fert.mcmc <- fert.brm3 %>%
    emmeans(~FERTILIZER,  at = newdata) %>%
    pairs() %>%
    as.data.frame()


## ----Probability2c2, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5, echo=1:2----
newdata <- list(FERTILIZER=c(200, 100)) 
fert.brm3 %>%
    emmeans(~FERTILIZER,  at=newdata) %>%
    regrid(transform = "log") %>%
    pairs() %>%
    regrid()


## ----Probability2d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>% emmeans(~FERTILIZER,  at=newdata) %>%
    regrid(transform = "log") %>%
    pairs() %>%
    regrid() %>%
    tidy_draws() %>%
    summarise_draws(median,
                    HDInterval::hdi,
                    P= ~ sum(.x>1.5)/length(.x))

fert.mcmc <- fert.brm3 %>%
    emmeans(~FERTILIZER,  at = newdata) %>% 
    tidy_draws() %>%
    rename_with(~str_replace(., 'FERTILIZER ', 'p')) %>%
    mutate(Eff=p200 - p100,
           PEff=100*Eff/p100)
fert.mcmc %>% head


## ----Probability2e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.mcmc %>%
    tidyMCMC(estimate.method='median',
             conf.int=TRUE, conf.method='HPDinterval')


## ----Probability2f, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.mcmc %>% median_hdci(PEff)
fert.mcmc %>% median_hdci(PEff, Eff)


## ----Probability2ff, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.mcmc %>%
    ggplot() +
    geom_density(aes(x=PEff))


## ----Probability2g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5, paged.print=FALSE----
fert.mcmc %>% summarise(P=sum(PEff>50)/n())


## ----Probability2h, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.mcmc %>% hypothesis('PEff>50')


## ----Probability2i, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
newdata = list(FERTILIZER=c(200, 100)) 
fert.brm3 %>%
    emmeans(~FERTILIZER,  at=newdata) %>%
    pairs() %>%
    tidy_draws() %>%
    summarise(across(contains('contrast'),
                     list(P = ~ sum(.>80)/n(),
                          HDCI = ~ median_hdci(.)),
                     .names = c('{.fn}')
                     )) %>%
    tidyr::unpack(HDCI)


## ----Probability2j, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
newdata = list(FERTILIZER=c(200, 100)) 
## Simple
fert.brm3 %>%
    emmeans(~FERTILIZER,  at=newdata) %>%
    regrid(transform = 'log') %>%
    pairs() %>%
    regrid() 

## More advanced (both P and percent change)
fert.mcmc <- fert.brm3 %>%
    emmeans(~FERTILIZER,  at=newdata) %>%
    regrid(transform = 'log') %>%
    pairs() %>%
    regrid() %>%
    tidy_draws() %>%
    mutate(across(contains('contrast'), ~ 100*(. - 1)))

fert.mcmc %>%
    summarise(across(contains('contrast'),
                     list(P = ~ sum(.>50)/n(),
                          HDCI = ~ median_hdci(.)),
                     .names = c('{.fn}')
                     )) %>%
    tidyr::unpack(HDCI)



## ----Probability2k, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
fert.brm3 %>%
    linpred_draws(newdata = as.data.frame(newdata)) %>%
    ungroup() %>%
    group_by(.draw) %>%
    summarise(Eff = -1 * diff(.linpred),
              PEff = 100*Eff/.linpred[2]) %>%
    ungroup() %>%
    mutate(P = sum(PEff>50)/n()) %>%
    pivot_longer(cols = -.draw) %>%
    group_by(name) %>%
    median_hdci()

##OR
fert.brm3 %>%
    epred_draws(newdata = as.data.frame(newdata)) %>%
    ungroup() %>%
    group_by(.draw) %>%
    summarise(Eff = -1 * diff(.epred),
              PEff = 100*Eff/.epred[2]) %>%
    ungroup() %>%
    mutate(P = sum(PEff>50)/n()) %>%
    pivot_longer(cols = -.draw) %>%
    group_by(name) %>%
    median_hdci()

##OR for prediction of individual values
fert.brm3 %>%
    predicted_draws(newdata = as.data.frame(newdata)) %>%
    ungroup() %>%
    group_by(.draw) %>%
    summarise(Eff = -1 * diff(.prediction),
              PEff = 100*Eff/.prediction[2]) %>%
    ungroup() %>%
    mutate(P = sum(PEff>50)/n()) %>%
    pivot_longer(cols = -.draw) %>%
    group_by(name) %>%
    median_hdci()


## ----summaryFig1a, results='markdown', eval=TRUE, hidden=TRUE-----------------
fert.list = with(fert, list(FERTILIZER = seq(min(FERTILIZER), max(FERTILIZER), len=100)))
newdata = emmeans(fert.rstanarm3, ~FERTILIZER, at=fert.list) %>% as.data.frame
head(newdata)

ggplot(newdata, aes(y=emmean, x=FERTILIZER)) + 
geom_point(data=fert, aes(y=YIELD)) +
geom_line() + 
geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), fill='blue', alpha=0.3) +
scale_y_continuous('YIELD') +
scale_x_continuous('FERTILIZER') +
theme_classic()

## spaghetti plot
newdata = emmeans(fert.rstanarm3, ~FERTILIZER, at=fert.list) %>%
  gather_emmeans_draws()
newdata %>% head
ggplot(newdata,  aes(y=.value,  x=FERTILIZER)) +
  geom_line(aes(group=.draw),  alpha=0.01) +
  geom_point(data=fert,  aes(y=YIELD))


## ----summaryFig1b, results='markdown', eval=TRUE, hidden=TRUE-----------------
fert.list = with(fert, list(FERTILIZER = seq(min(FERTILIZER), max(FERTILIZER), len=100)))
newdata = emmeans(fert.rstanarm3, ~FERTILIZER, at=fert.list) %>% as.data.frame
head(newdata)

ggplot(newdata, aes(y=emmean, x=FERTILIZER)) + 
    geom_point(data=fert, aes(y=YIELD)) +
    geom_line() + 
    geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), fill='blue', alpha=0.3) +
    scale_y_continuous(expression(Grass~yield~(g.m^-2)), breaks = seq(50,300, by = 50)) +
    scale_x_continuous(expression(Fertilizer~concentration~(g.m^-2))) +
    theme_classic()

## spaghetti plot
newdata = emmeans(fert.rstanarm3, ~FERTILIZER, at=fert.list) %>%
    gather_emmeans_draws()
newdata %>% head
ggplot(newdata,  aes(y=.value,  x=FERTILIZER)) +
    geom_line(aes(group=.draw),  colour = 'orange', alpha=0.01) +
    geom_point(data=fert,  aes(y=YIELD)) +
    scale_y_continuous(expression(Grass~yield~(g.m^-2)), breaks = seq(50,300, by = 50)) +
    scale_x_continuous(expression(Fertilizer~concentration~(g.m^-2))) +
    theme_classic()


## ----summaryFig2a, results='markdown', eval=TRUE, hidden=TRUE-----------------
fert.grid <- with(fert, list(FERTILIZER = modelr::seq_range(FERTILIZER, n = 100)))
newdata <- fert.brm3 %>%
    emmeans(~FERTILIZER, at=fert.grid) %>%
    as.data.frame
head(newdata)

ggplot(newdata, aes(y=emmean, x=FERTILIZER)) + 
    geom_point(data=fert, aes(y=YIELD)) +
    geom_line() + 
    geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), fill='blue', alpha=0.3) +
    scale_y_continuous('YIELD') +
    scale_x_continuous('FERTILIZER') +
    theme_classic()

## spaghetti plot
newdata = emmeans(fert.brm3, ~FERTILIZER, at = fert.grid) %>%
    gather_emmeans_draws()
newdata %>% head
ggplot(newdata,  aes(y = .value,  x = FERTILIZER)) +
    geom_line(aes(group = .draw),  color = 'blue', alpha = 0.01) +
    geom_point(data = fert,  aes(y = YIELD)) +
    theme_classic()
  


## ----summaryFig2b, results='markdown', eval=TRUE, hidden=TRUE-----------------
fert.list = with(fert, list(FERTILIZER = seq(min(FERTILIZER), max(FERTILIZER), len=100)))
newdata = emmeans(fert.brm3, ~FERTILIZER, at=fert.list) %>% as.data.frame
head(newdata)

ggplot(newdata, aes(y=emmean, x=FERTILIZER)) + 
    geom_point(data=fert, aes(y=YIELD)) +
    geom_line() + 
    geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD), fill='blue', alpha=0.3) +
    scale_y_continuous(expression(Grass~yield~(g.m^-2)), breaks = seq(50,300, by = 50)) +
    scale_x_continuous(expression(Fertilizer~concentration~(g.m^-2))) +
    theme_classic()

## spaghetti plot
newdata = emmeans(fert.brm3, ~FERTILIZER, at=fert.list) %>%
    gather_emmeans_draws()
newdata %>% head
ggplot(newdata,  aes(y=.value,  x=FERTILIZER)) +
    geom_line(aes(group=.draw),  colour = 'orange', alpha=0.01) +
    geom_point(data=fert,  aes(y=YIELD)) +
    scale_y_continuous(expression(Grass~yield~(g.m^-2)), breaks = seq(50,300, by = 50)) +
    scale_x_continuous(expression(Fertilizer~concentration~(g.m^-2))) +
    theme_classic()

