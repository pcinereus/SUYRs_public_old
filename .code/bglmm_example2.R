## ----setup, include=FALSE, warnings=FALSE, message=FALSE----------------------
knitr::opts_chunk$set(echo = TRUE)


## ----libraries, results='markdown', eval=TRUE, message=FALSE, warning=FALSE----
library(rstanarm)   #for fitting models in STAN
library(brms)       #for fitting models in STAN
library(coda)       #for diagnostics
library(bayesplot)  #for diagnostics
library(ggmcmc)     #for MCMC diagnostics
library(rstan)      #for interfacing with STAN
library(emmeans)    #for marginal means etc
library(broom)      #for tidying outputs
library(DHARMa)     #for residual diagnostics
library(tidybayes)  #for more tidying outputs
library(ggeffects)  #for partial plots
library(broom.mixed)#for tidying MCMC outputs
library(tidyverse)  #for data wrangling etc
library(patchwork)  #for multiple plots
library(ggridges)   #for ridge plots 
source('helperFunctions.R')


## ----readData, results='markdown', eval=TRUE----------------------------------
norin <- read_csv('../data/norin.csv', trim_ws=TRUE)
glimpse(norin)


## ----eda1, results='markdown', eval=TRUE, hidden=TRUE-------------------------
norin <- norin %>% mutate(FISHID=factor(FISHID),
                         TRIAL=factor(TRIAL))


## ----eda2, results='markdown', eval=TRUE, hidden=TRUE-------------------------
ggplot(norin, aes(y=CHANGE, x=TRIAL)) + geom_boxplot()


## ----eda3, results='markdown', eval=TRUE, hidden=TRUE, fig.width=7, fig.height=4----
ggplot(norin, aes(y=CHANGE, x=SMR_contr, shape=TRIAL, color=TRIAL)) +
    geom_smooth(method='lm') + geom_point()
ggplot(norin, aes(y=CHANGE, x=SMR_contr, shape=TRIAL, color=TRIAL)) +
  geom_smooth() + geom_point()


## ----eda4, results='markdown', eval=TRUE, hidden=TRUE, fig.width=7, fig.height=4----
ggplot(norin, aes(y=CHANGE, x=as.numeric(FISHID), color=TRIAL)) +
    geom_point() + geom_line()


## ----eda5, results='markdown', eval=TRUE, hidden=TRUE, fig.width=7, fig.height=4----
ggplot(norin, aes(y=CHANGE, x=MASS, color=TRIAL)) +
  geom_point() +
  geom_smooth(method='lm')


## ----fitModel1a, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
norin.rstanarm <- stan_glmer(CHANGE ~  (1|FISHID)+TRIAL*SMR_contr+MASS,
                               data = norin,
                               family = gaussian(), 
                               iter = 5000,
                               warmup = 2000,
                               chains = 3,
                               thin = 5,
                               refresh = 0)


## ----fitModel1b, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
norin.rstanarm %>% prior_summary()


## ----fitModel1c, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
2.5*sd(norin$CHANGE)


## ----fitModel1d, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
2.5*sd(norin$CHANGE)/apply(model.matrix(~TRIAL*SMR_contr+MASS, norin)[, -1], 2, sd)


## ----fitModel1e, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
1/sd(norin$CHANGE)


## ----fitModel1f, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
norin.rstanarm1 <- update(norin.rstanarm,  prior_PD=TRUE)


## ----fitModel1g, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
norin.rstanarm1 %>%
  ggpredict(~SMR_contr*TRIAL) %>%
  plot(add.data=TRUE)
#OR
norin.rstanarm1 %>%
  ggemmeans(~SMR_contr*TRIAL) %>%
  plot(add.data=TRUE)


## ----fitModel1h, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
norin.rstanarm2 <- stan_glmer(CHANGE ~  (1|FISHID)+TRIAL*scale(SMR_contr)+scale(MASS),
                                data = norin,
                                family = gaussian(), 
                                prior_intercept = normal(17, 35, autoscale = FALSE),
                                prior = normal(0, 70, autoscale = FALSE),
                                prior_aux=rstanarm::exponential(0.03, autoscale = FALSE),
                                prior_covariance = decov(1, 1, 1, 1), 
                                prior_PD = TRUE, 
                                iter = 5000,
                                warmup = 1000,
                                chains = 3,
                                thin = 5,
                                refresh = 0
                                )


## ----fitModel1i, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
norin.rstanarm2 %>%
    ggpredict(~SMR_contr * TRIAL) %>%
    plot(add.data = TRUE)


## ----fitModel1j, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE, dependson='fitModel1h'----
norin.rstanarm3 <- update(norin.rstanarm2,  prior_PD=FALSE)


## ----modelFit1k, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=8----
get_variables(norin.rstanarm3)
posterior_vs_prior(norin.rstanarm3, color_by='vs', group_by=TRUE,
                   facet_args=list(scales='free_y'),
                   regex_pars = "^.Intercept|TRIAL|SMR_contr|MASS|sigma|Sigma")


## ----modelFit1l, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
norin.rstanarm3 %>%
    ggpredict(~SMR_contr * TRIAL) %>%
    plot(add.data = TRUE)
##OR
norin.rstanarm3 %>%
    ggemmeans(~SMR_contr * TRIAL) %>%
    plot(add.data = TRUE)


## ----fitModel2a, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE, paged.print=FALSE, tidy.opts = list(width.cutoff = 80), echo=c(-4,-6)----
norin.form <- bf(CHANGE ~  (1|FISHID)+TRIAL*SMR_contr+MASS,
                   family = gaussian() 
                   )
options(width=100)
norin.form %>% get_prior(data=norin)
options(width=80)
norin.brm <- brm(norin.form,
                  data=norin,
                  iter = 5000,
                  warmup = 1000,
                  chains = 3,
                  thin = 5,
                  refresh = 0)


## ----fitModel2h, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, fig.width = 10, fig.height = 7----
norin %>% 
    group_by(TRIAL) %>%
    summarise(median(CHANGE),
              mad(CHANGE))
norin %>% 
    group_by(TRIAL, FISHID) %>%
    summarise(median = median(CHANGE),
              MAD = mad(CHANGE)) %>%
    ungroup(FISHID) %>%
    summarise(sd(median))

sd(norin$CHANGE)/apply(model.matrix(~TRIAL*SMR_contr+MASS, norin)[, -1], 2, sd)

standist::visualize("normal(16,35)", xlim=c(-10,100))
standist::visualize("normal(0, 70)", "normal(0, 54)", xlim=c(-200,200))
standist::visualize("gamma(2, 1)", "gamma(35, 1)",
                    xlim=c(-10,100))


## ----fitModel2h1, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE------
priors <- prior(normal(16, 35), class = 'Intercept') +
    prior(normal(0, 70), class = 'b', coef = 'TRIALHypoxia') +
    prior(normal(0, 70), class = 'b', coef = 'TRIALLowSalinity') +
    prior(normal(0, 54), class = 'b', coef = 'SMR_contr') +
    prior(normal(0, 4), class = 'b', coef = 'MASS') +
    prior(normal(0, 13), class = 'b') +
    prior(gamma(35, 1), class = 'sigma') +
    prior(cauchy(0, 2.5), class = 'sd') 
norin.form <- bf(CHANGE ~ (1|FISHID) + TRIAL*SMR_contr + MASS,
                     family = gaussian()
                   )
norin.brm2 <- brm(norin.form, 
                  data = norin,
                  prior = priors,
                  sample_prior = 'only',
                  iter = 5000,
                  warmup = 1000,
                  chains = 3,
                  thin = 5,
                  refresh = 0
                  )


## ----partialPlot2h1a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.brm2 %>%
    ggpredict(~SMR_contr*TRIAL) %>%
    plot(add.data = TRUE)


## ----fitModel2h1b, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-----
norin.brm3 <- update(norin.brm2,  
                       sample_prior = 'yes',
                       control = list(adapt_delta = 0.99),
                       refresh = 0)
save(norin.brm3, file = '../ws/testing/norin.brm3')


## ----partialPlot2h1b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.brm3 %>%
    ggpredict(~SMR_contr*TRIAL) %>%
    plot(add.data = TRUE)


## ----posterior2h2, results='markdown', eval=TRUE------------------------------
norin.brm3 %>% get_variables()
norin.brm3 %>% hypothesis('TRIALHypoxia=0') %>% plot
norin.brm3 %>% hypothesis('SMR_contr=0') %>% plot


## ----posterior2h2a, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
norin.brm3 %>% SUYR_prior_and_posterior()


## ----fitModel2h3, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE------
priors <- prior(normal(16, 35), class = 'Intercept') +
    prior(normal(0, 70), class = 'b', coef = 'TRIALHypoxia') +
    prior(normal(0, 70), class = 'b', coef = 'TRIALLowSalinity') +
    prior(normal(0, 54), class = 'b', coef = 'SMR_contr') +
    prior(normal(0, 4), class = 'b', coef = 'MASS') +
    prior(normal(0, 13), class = 'b') +
    prior(gamma(35, 1), class = 'sigma') +
    prior(cauchy(0, 2.5), class = 'sd') +
    prior(lkj_corr_cholesky(1), class = 'L')
norin.form <- bf(CHANGE ~ (TRIAL|FISHID) + TRIAL*SMR_contr + MASS,
                     family = gaussian()
                   )
norin.brm4 <-  brm(norin.form, 
                  data = norin,
                  prior = priors,
                  sample_prior = 'yes',
                  iter = 5000,
                  warmup = 1000,
                  chains = 3,
                  thin = 5,
                  refresh = 0,
                  control = list(adapt_delta=0.99)
                  )
save(norin.brm4, file = '../ws/testing/norin.brm4')


## ----posterior2k, results='markdown', eval=TRUE-------------------------------
norin.brm4 %>% get_variables()
norin.brm4 %>% hypothesis('TRIALHypoxia=0') %>% plot
norin.brm4 %>% hypothesis('SMR_contr=0') %>% plot


## ----posterior2k1, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
norin.brm4 %>% SUYR_prior_and_posterior()


## ----posterior2k2, results='markdown', eval=TRUE, fig.width=10, fig.height=4----
norin.brm4 %>%
  posterior_samples %>%
  dplyr::select(-`lp__`) %>%
  pivot_longer(everything(), names_to = 'key') %>% 
  filter(!str_detect(key, '^r')) %>%
  mutate(Type = ifelse(str_detect(key, 'prior'), 'Prior', 'Posterior'),
         Class = case_when(
             str_detect(key, '(^b|^prior).*Intercept$') ~ 'Intercept',
             str_detect(key, 'b_TRIAL.*|prior_b_TRIAL.*') & str_detect(key, '.*\\:.*') ~ 'TRIAL',
             str_detect(key, 'b_SMR_contr|prior_b_SMR_contr') ~ 'SMR',
             str_detect(key, 'b_MASS|prior_b_MASS') ~ 'MASS',
             str_detect(key, '.*\\:.*|prior_b_.*\\:.*') ~ 'Interaction',
             str_detect(key, 'sd') ~ 'sd',
             str_detect(key, '^cor|prior_cor') ~ 'cor',
             str_detect(key, 'sigma') ~ 'sigma'),
         Par = str_replace(key, 'b_', '')) %>%
  ggplot(aes(x = Type,  y = value, color = Par)) +
  stat_pointinterval(position = position_dodge())+
  facet_wrap(~Class,  scales = 'free')



## ----fitModel2h3a, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-----
(l.1 <- norin.brm3 %>% loo())
(l.2 <- norin.brm4 %>% loo())
loo_compare(l.1, l.2)


## ----modelValidation2a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation2b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
norin.brm4 %>% mcmc_plot(type='trace')


## ----modelValidation2c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
norin.brm4 %>% mcmc_plot(type='acf_bar')


## ----modelValidation2d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
norin.brm4 %>% mcmc_plot(type='rhat_hist')


## ----modelValidation2e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
norin.brm4 %>% mcmc_plot(type='neff_hist')


## ----modelValidation2f, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
norin.brm4 %>% mcmc_plot(type='combo')
norin.brm4 %>% mcmc_plot(type='violin')


## ----modelValidation2g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
norin.brm4 %>% get_variables()
pars <- norin.brm4 %>% get_variables()
pars <- str_extract(pars, '^b_.*|^sigma$|^sd.*') %>% na.omit()

norin.brm4$fit %>%
    stan_trace(pars = pars)


## ----modelValidation2h, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
norin.brm4$fit %>%
    stan_ac(pars = pars)


## ----modelValidation2i, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
norin.brm4$fit %>% stan_rhat() 


## ----modelValidation2j, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
norin.brm4$fit %>% stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
norin.brm4$fit %>%
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation2l, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=7----
norin.ggs <- norin.brm4 %>% ggs(burnin = FALSE, inc_warmup = FALSE)
norin.ggs %>% ggs_traceplot()


## ----modelValidation2m, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=7----
ggs_autocorrelation(norin.ggs)


## ----modelValidation2n, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggs_Rhat(norin.ggs, scaling = 1.01)


## ----modelValidation2o, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggs_effective(norin.ggs)


## ----modelValidation2p, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggs_crosscorrelation(norin.ggs)


## ----modelValidation2q, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggs_grb(norin.ggs)


## ----modelValidation5a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
norin.brm4 %>% pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation5c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
norin.brm4 %>% pp_check(type = 'error_scatter_avg')


## ----modelValidation5e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
norin.brm4 %>% pp_check(group = 'TRIAL', type = 'intervals')
norin.brm3 %>% pp_check(group = 'TRIAL', type = 'intervals_grouped')
norin.brm3 %>% pp_check(group = 'TRIAL', type = 'violin_grouped')


## ----modelValidation5g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(norin.brm2)


## ----modelValidation6a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
preds <- norin.brm4 %>% posterior_predict(ndraws = 250,  summary = FALSE)
norin.resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = norin$CHANGE,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = FALSE)
plot(norin.resids, quantreg = FALSE)


## ----partialPlot2d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.brm4 %>%
    conditional_effects() %>%
    plot(ask = FALSE, points = TRUE, plot = FALSE) %>%
    wrap_plots()


## ----partialPlot2a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.brm4 %>%
    ggpredict(terms = c("SMR_contr", "TRIAL")) %>%
    plot(add.data = TRUE)


## ----partialPlot2b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.brm4 %>%
    ggemmeans(~SMR_contr*TRIAL) %>%
    plot(add.data = TRUE)


## ----summariseModel2a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.brm4 %>% summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
norin.sum <- summary(norin.brm4)


## ----summariseModel2i, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.brm4 %>% as_draws_df()
norin.brm4 %>%
  as_draws_df() %>%
  summarise_draws(
    "median",
    ~ HDInterval::hdi(.x),
    "rhat",
    "ess_bulk"
  )
## or if you want to exclude some parameters
norin.brm4 %>%
  as_draws_df() %>%
  summarise_draws(
    "median",
    ~ HDInterval::hdi(.x),
    "rhat",
    "ess_bulk"
  ) %>%
  filter(str_detect(variable, 'prior|^r_|^lp__', negate = TRUE)) 


## ----summariseModel2b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.brm4$fit %>%
    tidyMCMC(estimate.method = 'median',
             conf.int = TRUE,  conf.method = 'HPDinterval',
             rhat = TRUE, ess = TRUE)

## ----summariseModel2b1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
norin.tidy <- tidyMCMC(norin.brm4$fit, estimate.method='median',
                         conf.int=TRUE,  conf.method='HPDinterval',
                         rhat=TRUE, ess=TRUE)


## ----summariseModel2c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.brm4 %>% get_variables()
norin.draw <- norin.brm4 %>%
    gather_draws(`^b.Intercept$|b_.*|sd_.*|sigma`,  regex=TRUE)
norin.draw


## ----summariseModel2c1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.draw %>% median_hdci


## ----summariseModel2c3, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
norin.gather <- norin.brm4 %>%
    gather_draws(`b_Intercept|b_TREAT.*|sd_.*|sigma`,  regex=TRUE) %>%
  median_hdci


## ----summariseModel2c4, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
norin.brm4 %>%
    gather_draws(`b_Intercept|b_.*`, regex=TRUE) %>% 
    ggplot() +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                           .width = c(0.5, 0.8, 0.95), 
                           labels = scales::percent_format())
                           )), color='black') + 
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE) 

norin.brm4 %>% 
    gather_draws(`b_Intercept|b_.*`, regex=TRUE) %>% 
    ggplot() + 
    stat_halfeye(aes(x=.value,  y=.variable)) +
    facet_wrap(~.variable, scales='free')

norin.brm4 %>% 
    gather_draws(`b_Intercept|b_.*`, regex=TRUE) %>% 
    ggplot() + 
    stat_halfeye(aes(x=.value,  y=.variable)) +
    theme_classic()


## ----summariseModel2j, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.brm4$fit %>% plot(type='intervals') 


## ----summariseModel2ka, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
norin.brm4 %>% 
    gather_draws(`^b_.*`, regex=TRUE) %>% 
    filter(.variable != 'b_Intercept') %>%
    ggplot() + 
    stat_halfeye(aes(x=.value,  y=.variable)) +
    facet_wrap(~.variable, scales='free')

norin.brm4 %>% 
    gather_draws(`^b_.*`, regex=TRUE) %>% 
    filter(.variable != 'b_Intercept') %>%
    ggplot() + 
    stat_halfeye(aes(x=.value,  y=.variable)) +
    geom_vline(xintercept = 0, linetype = 'dashed')


## ----summariseModel2c7, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
norin.brm4 %>% 
    gather_draws(`^b_.*`, regex=TRUE) %>% 
    filter(.variable != 'b_Intercept') %>%
    ggplot() +  
    geom_density_ridges(aes(x=.value, y = .variable), alpha=0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed')
##Or in colour
norin.brm4 %>% 
    gather_draws(`^b_.*`, regex=TRUE) %>% 
    filter(.variable != 'b_Intercept') %>%
    ggplot() + 
    geom_density_ridges_gradient(aes(x=exp(.value),
                                     y = .variable,
                                     fill = stat(x)),
                                 alpha=0.4, colour = 'white',
                                 quantile_lines = TRUE,
                                 quantiles = c(0.025, 0.975)) +
    geom_vline(xintercept = 1, linetype = 'dashed') +
    scale_x_continuous(trans = scales::log2_trans()) +
    scale_fill_viridis_c(option = "C") 


## ----summariseModel2d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.brm4 %>% tidy_draws()


## ----summariseModel2e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.brm4 %>% spread_draws(`.*Intercept.*|^b_.*`,  regex=TRUE)


## ----summariseModel2f, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.brm4 %>% posterior_samples() %>% as_tibble()


## ----summariseModel2g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
norin.brm4 %>%
    bayes_R2(re.form = NA, summary=FALSE) %>%
    median_hdci
norin.brm4 %>%
    bayes_R2(re.form = ~(1|FISHID), summary=FALSE) %>%
    median_hdci
norin.brm4 %>%
    bayes_R2(re.form = ~(TRIAL|FISHID), summary=FALSE) %>%
    median_hdci


## ----fitModels, results='markdown', eval=FALSE, hidden=TRUE-------------------
## norin = norin %>% mutate(FISHID=factor(FISHID),
##                          TRIAL=factor(TRIAL))
## 
## ggplot(norin, aes(y=CHANGE, x=TRIAL)) + geom_boxplot()
## ggplot(norin, aes(y=CHANGE, x=SMR_contr, shape=TRIAL, color=TRIAL)) +
##     geom_smooth(method='lm') + geom_point()
## #ggplot(norin, aes(y=CHANGE, x=MASS, shape=TRIAL, color=TRIAL)) +
## #geom_smooth(method='lm') + geom_point()
## ggplot(norin, aes(y=CHANGE, x=as.numeric(FISHID), color=TRIAL)) +
##     geom_point() + geom_line()
## 
## #ggplot(norin, aes(y=MASS, x=TRIAL)) + geom_boxplot()
## ggplot(norin, aes(y=CHANGE, x=MASS, color=TRIAL)) + geom_point() + geom_smooth(method='lm')
## 
## norin.rstanarm = stan_glmer(CHANGE ~ (1|FISHID)+TRIAL*SMR_contr+MASS, data=norin,
##                             prior_PD=TRUE,
##                          iter=5000, warmup=2000, chains=3, thin=5, refresh=0)
## prior_summary(norin.rstanarm)
## 
## posterior_vs_prior(norin.rstanarm, color_by='vs', group_by=TRUE,
##                    facet_args=list(scales='free_y'), pars=c('(Intercept)'))
## ggpredict(norin.rstanarm, ~TRIAL*SMR_contr) %>% plot(add.data=TRUE)
## 
## norin.rstanarm %>% get_variables()
## plot(norin.rstanarm,  'mcmc_trace', regex_pars='^.Intercept|TRIAL|SMR|MASS|[sS]igma')
## plot(norin.rstanarm,  'mcmc_acf_bar', regex_pars='^.Intercept|TRIAL|SMR|MASS|[sS]igma')
## plot(norin.rstanarm,  'mcmc_rhat_hist', regex_pars='^.Intercept|TRIAL|SMR|MASS|[sS]igma')
## plot(norin.rstanarm,  'mcmc_neff_hist', regex_pars='^.Intercept|TRIAL|SMR|MASS|[sS]igma')
## 
## #norin.rstan1 = stan_glmer(CHANGE ~ (TRIAL|FISHID)+TRIAL*SMR_contr+MASS, data=norin,
## #                          iter=5000, warmup=2000, chains=3, thin=5, refresh=0, cores=3)
## norin.rstanarm1 = stan_glmer(CHANGE ~ (SMR_contr|FISHID) + TRIAL*SMR_contr+MASS, data=norin,
##                              prior_PD=FALSE,
##                           iter=5000, warmup=2000, chains=3, thin=5, refresh=0, cores=3)
## norin.rstanarm1 = update(norin.rstanarm1,  prior_PD=FALSE)
## 
## 
## 
## norin.rstanarm2 = stan_glmer(CHANGE ~ (TRIAL*SMR_contr|FISHID) + TRIAL*SMR_contr+MASS, data=norin,
##                              prior_PD=FALSE,
##                           iter=5000, warmup=2000, chains=3, thin=5, refresh=0, cores=3)
## 
## posterior_vs_prior(norin.rstanarm1, color_by='vs', group_by=TRUE,
##                    facet_args=list(scales='free_y'), pars=c('(Intercept)'))
## 
## ggpredict(norin.rstanarm1, ~TRIAL*SMR_contr) %>% plot(add.data=TRUE)
## 
## norin.rstanarm1 %>% get_variables()
## plot(norin.rstanarm1,  'mcmc_trace', regex_pars='^.Intercept|TRIAL|^SMR|MASS|[sS]igma')
## plot(norin.rstanarm1,  'mcmc_acf_bar', regex_pars='^.Intercept|TRIAL|^SMR|MASS|[sS]igma')
## plot(norin.rstanarm1,  'mcmc_rhat_hist', regex_pars='^.Intercept|TRIAL|^SMR|MASS|[sS]igma')
## plot(norin.rstanarm1,  'mcmc_neff_hist', regex_pars='^.Intercept|TRIAL|^SMR|MASS|[sS]igma')
## 
## 
## (l.1 <- loo(norin.rstanarm))
## (l.2 <- loo(norin.rstanarm1))
## loo_compare(l.1,  l.2)
## 
## 
## preds <- posterior_predict(norin.rstanarm,  nsamples=250,  summary=FALSE)
## norin.resids <- createDHARMa(simulatedResponse = t(preds),
##                             observedResponse = norin$CHANGE,
##                             fittedPredictedResponse = apply(preds, 2, median))
## plot(norin.resids)
## 
## 
## g=ggpredict(norin.rstanarm) %>% plot
## do.call('grid.arrange', g)
## 
## #ggemmeans(norin.rstan, ~TRIAL)
## 
## summary(norin.rstanarm)
## nms <- norin.rstanarm1 %>% get_variables()
## wch <- grep('^.Intercept|TRIAL|^SMR|[sS]igma', nms)
## tidyMCMC(norin.rstanarm$stanfit,conf.int=TRUE, conf.method='HPDinterval',
##          rhat=TRUE, ess=TRUE, pars=nms[wch], estimate.method='median')
## 
## tidyMCMC(norin.rstanarm1$stanfit,conf.int=TRUE, conf.method='HPDinterval',
##          rhat=TRUE, ess=TRUE, pars=nms[wch], estimate.method='median')
## 
## 
## norin.grid = with(norin, list(SMR_contr=seq(min(SMR_contr),max(SMR_contr), len=100)))
## newdata = emmeans(norin.rstanarm, ~TRIAL|SMR_contr, at=norin.grid) %>% as.data.frame
## head(newdata)
## ggplot(newdata, aes(y=emmean, x=SMR_contr, color=TRIAL)) +
##     geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD, fill=TRIAL), alpha=0.3,color=NA) +
##     geom_line()
## 
## norin.grid = with(norin, list(SMR_contr=c(min(SMR_contr),mean(SMR_contr),max(SMR_contr))))
## 
## emmeans(norin.rstan, pairwise~TRIAL|SMR_contr, at=norin.grid)
## 
## norin.em = emmeans(norin.rstan, pairwise~TRIAL|SMR_contr, at=norin.grid)$contrast %>%
##               gather_emmeans_draws() %>%
##               mutate(Fit=.value)
## norin.em
## norin.em %>% group_by(contrast) %>% median_hdci(Fit)
## norin.em %>% group_by(contrast, SMR_contr) %>% median_hdci(Fit)
## ## norin.em %>%
## ##     group_by(contrast) %>%
## ##     summarize(P=sum(Fit>0)/n())
## norin.em %>%
##     group_by(contrast, SMR_contr) %>%
##     summarize(P=sum(Fit>0)/n())
## 
## 
## bayes_R2(norin.rstanarm, re.form=NA) %>% median_hdi
## bayes_R2(norin.rstanarm, re.form=~(1|FISHID)) %>% median_hdi
## #bayes_R2(norin.rstan1, re.form=~(SMR_contr|FISHID)) %>% median_hdi
## 


## ----fitModels.brms, results='markdown', eval=FALSE, hidden=TRUE--------------
## norin = norin %>% mutate(FISHID=factor(FISHID),
##                          TRIAL=factor(TRIAL))
## 
## ggplot(norin, aes(y=CHANGE, x=TRIAL)) + geom_boxplot()
## ggplot(norin, aes(y=CHANGE, x=SMR_contr, shape=TRIAL, color=TRIAL)) +
##     geom_smooth(method='lm') + geom_point()
## ggplot(norin, aes(y=CHANGE, x=MASS, shape=TRIAL, color=TRIAL)) +
##     geom_smooth(method='lm') + geom_point()
## ggplot(norin, aes(y=CHANGE, x=as.numeric(FISHID), color=TRIAL)) +
##     geom_point() + geom_line()
## 
## ##ggplot(norin, aes(y=MASS, x=TRIAL)) + geom_boxplot()
## ##ggplot(norin, aes(y=CHANGE, x=MASS, color=TRIAL)) + geom_point() + geom_smooth(method='lm')
## 
## norin %>% group_by(TRIAL) %>%
##     summarise(median(CHANGE),
##               mad(CHANGE))
## priors <- prior(normal(50, 20), class='Intercept') +
##     prior(normal(0, 10), class='b') +
##     prior(gamma(2,1), class='sigma') +
##     prior(gamma(2,1), class='sd')
## 
## norin.form <- bf(CHANGE ~ (1|FISHID)+TRIAL*SMR_contr+MASS,
##                  family=gaussian)
## 
## norin.brm1 = brm(norin.form,
##                  data=norin,
##                  prior = priors,
##                  sample_prior = 'yes',
##                  iter=5000, warmup=2000,
##                  chains=3, thin=5, refresh=0)
## 
## norin.brm1 %>% get_variables()
## pars <- norin.brm1 %>% get_variables()
## wch <- grepl('^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd', pars, perl=TRUE)
## 
## stan_trace(norin.brm1$fit, pars=pars[wch])
## stan_ac(norin.brm1$fit, pars=pars[wch])
## stan_rhat(norin.brm1$fit, pars=pars[wch])
## stan_ess(norin.brm1$fit, pars=pars[wch])
## 
## ##mcmc_plot(norin.brms,  type='trace',
## ##          regex_pars='^b.*|sigma|^sd')
## ##mcmc_plot(norin.brms,  type='trace',
## ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
## ##mcmc_plot(norin.brms,  type='acf_bar',
## ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
## ##mcmc_plot(norin.brms,  type='rhat_hist',
## ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
## ##mcmc_plot(norin.brms,  type='neff_hist',
## ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
## 
## preds <- posterior_predict(norin.brm1,  nsamples=250,  summary=FALSE)
## norin.resids <- createDHARMa(simulatedResponse = t(preds),
##                             observedResponse = norin$CHANGE,
##                             fittedPredictedResponse = apply(preds, 2, median),
##                             integerResponse =FALSE)
## plot(norin.resids)
## #norin.rstan1 = stan_glmer(CHANGE ~ (TRIAL|FISHID)+TRIAL*SMR_contr+MASS, data=norin,
## #                          iter=5000, warmup=2000, chains=3, thin=5, refresh=0, cores=3)
## norin.form <- bf(CHANGE ~ (TRIAL|FISHID) + TRIAL*SMR_contr+MASS,
##                  family=gaussian)
## norin.brm2 = brm(norin.form, data=norin,
##                  prior = priors,
##                  sample_prior = 'yes',
##                  iter=5000, warmup=2000,
##                  chains=3, thin=5, refresh=0, cores=3,
##                  control=list(adapt_delta=0.99))
## 
## norin.brm2 %>% get_variables()
## 
## pars <- norin.brm2 %>% get_variables()
## ## wch <- grepl('^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd', pars, perl=TRUE)
## wch <- grepl('^b_.*|[sS]igma|^sd_.*', pars, perl=TRUE)
## 
## stan_trace(norin.brm2$fit, pars=pars[wch])
## stan_ac(norin.brm2$fit, pars=pars[wch])
## stan_rhat(norin.brm2$fit)#, pars=pars[wch])
## stan_ess(norin.brm2$fit)#, pars=pars[wch])
## ##mcmc_plot(norin.brm2,  type='trace',
## ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
## ##mcmc_plot(norin.brm2,  type='trace',
## ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
## ##mcmc_plot(norin.brm2,  type='acf_bar',
## ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
## ##mcmc_plot(norin.brm2,  type='rhat_hist',
## ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
## ##mcmc_plot(norin.brm2,  type='neff_hist',
## ##          regex_pars='^b.Intercept|TRIAL|SMR|MASS|[sS]igma|^sd')
## 
## (l.1 <- loo(norin.brm1))
## (l.2 <- loo(norin.brm2))
## loo_compare(l.1,  l.2)
## 
## 
## preds <- posterior_predict(norin.brm2,  nsamples=250,  summary=FALSE)
## norin.resids <- createDHARMa(simulatedResponse = t(preds),
##                             observedResponse = norin$CHANGE,
##                             fittedPredictedResponse = apply(preds, 2, median))
## plot(norin.resids)
## 
## g <- norin.brm2 %>%
##     conditional_effects() %>%
##     plot(points=TRUE, ask=FALSE)
## library(patchwork)
## g[[1]] + g[[2]] + g[[3]] + g[[4]]
## 
## 
## ##g=ggpredict(norin.brms1) %>% plot
## ##library(patchwork)
## ##g[[1]] + g[[2]] + g[[3]]
## 
## ##do.call('grid.arrange', g)
## 
## ggemmeans(norin.brm2, ~TRIAL) %>% plot
## 
## summary(norin.brm2)
## 
## tidyMCMC(norin.brm2$fit,conf.int=TRUE, conf.method='HPDinterval',
##          rhat=TRUE, ess=TRUE, estimate.method='median') %>%
##   slice(1:11)
## 
## pars <- norin.brm2 %>% get_variables()
## wch <- grep('^b.Intercept|TRIAL|^b.*SMR|[sS]igma|^sd', pars)
## tidyMCMC(norin.brms1$fit,conf.int=TRUE, conf.method='HPDinterval',
##          rhat=TRUE, ess=TRUE, pars=pars[wch], estimate.method='median')
## 
## bayes_R2(norin.brm2, re.form=NA,  summary=FALSE) %>%
##     median_hdci
## bayes_R2(norin.brm2, re.form=~(1|FISHID), summary=FALSE) %>%
##     median_hdci
## bayes_R2(norin.brm2, re.form=~(TRIAL|FISHID), summary=FALSE) %>%
##     median_hdci
## 
## emmeans(norin.brm2, pairwise~TRIAL)
## 
## 
## norin.em <- norin.brm2 %>%
##     emmeans(~TRIAL) %>%
##     pairs() %>%
##     gather_emmeans_draws() %>%
##     mutate(Fit=.value)
## 
## norin.em %>%
##   group_by(contrast) %>%
##   median_hdi()
## 
## norin.em %>%
##     ggplot() +
##     geom_vline(xintercept=0, linetype='dashed') +
##     stat_slab(aes(x=.value, y=contrast,
##                   fill = stat(ggdist::cut_cdf_qi(cdf,
##                             .width = c(0.5, 0.8, 0.95),
##                             labels = scales::percent_format())
##                             )), color='black') +
##     scale_fill_brewer('Interval', direction = -1, na.translate = FALSE) +
##     theme_bw()
## 
## norin.em %>%
##     group_by(contrast) %>%
##   summarize(P=sum(Fit>0)/n())
## 
## 
## norin.grid <- with(norin, list(SMR_contr=c(min(SMR_contr),
##                                            mean(SMR_contr),
##                                            max(SMR_contr))))
## 
## norin.em <- norin.brm2 %>%
##     emmeans(~TRIAL|SMR_contr, at=norin.grid) %>%
##     pairs() %>%
##     gather_emmeans_draws()
## 
## norin.em %>% head
## norin.em %>%
##     group_by(contrast, SMR_contr) %>%
##     median_hdi()
## 
## norin.em %>%
##     group_by(contrast, SMR_contr) %>%
##     summarize(P=sum(.value>0)/n())
## 
## norin.grid <- with(norin, list(SMR_contr=modelr::seq_range(SMR_contr, n=100)))
## newdata <- norin.brm2 %>%
##     emmeans(~SMR_contr|TRIAL, at=norin.grid) %>%
##     as.data.frame
## head(newdata)
## partial.obs <- norin %>%
##     mutate(Pred = predict(norin.brm2, re.form = NA, summary=TRUE)[,'Estimate'],
##            Resid = resid(norin.brm2)[,'Estimate'],
##            Obs = Pred + Resid)
## ggplot(newdata, aes(y=emmean, x=SMR_contr, color=TRIAL)) +
##     geom_point(data=partial.obs, aes(y=Obs)) +
##     ##geom_point(data=partial.obs, aes(y=CHANGE), shape=2) +
##     geom_ribbon(aes(ymin=lower.HPD, ymax=upper.HPD, fill=TRIAL), alpha=0.3,color=NA) +
##     geom_line()

