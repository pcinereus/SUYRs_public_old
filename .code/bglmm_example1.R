## ----setup, include=FALSE, warnings=FALSE, message=FALSE----------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE,cache.lazy = FALSE, tidy='styler')


## ----libraries, results='markdown', eval=TRUE, message=FALSE, warning=FALSE----
library(tidyverse)  #for data wrangling etc
library(rstanarm)   #for fitting models in STAN
library(cmdstanr)   #for cmdstan
library(brms)       #for fitting models in STAN
library(standist)   #for exploring distributions
library(HDInterval) #for HPD intervals
library(posterior)  #for posterior draws
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
library(patchwork)  #for multiple plots
library(ggridges)   #for ridge plots 
library(bayestestR) #for ROPE
library(see)        #for some plots
source('helperFunctions.R')


## ----readData, results='markdown', eval=TRUE----------------------------------
tobacco <- read_csv('../data/tobacco.csv', trim_ws = TRUE)
tobacco %>% glimpse()


## ----processData, results='markdown', eval=TRUE-------------------------------
tobacco <- tobacco %>% mutate(LEAF = factor(LEAF),
                             TREATMENT = factor(TREATMENT))
tobacco %>% head()


## ----tobaccoEDA2, results='markdown', eval=TRUE, hidden=TRUE------------------
ggplot(tobacco,  aes(y = NUMBER,  x = TREATMENT)) +
  geom_boxplot()


## ----tobaccoEDA3, results='markdown', eval=TRUE, hidden=TRUE------------------
ggplot(tobacco,  aes(y = NUMBER,  x = as.numeric(LEAF))) +
  geom_line(aes(linetype = TREATMENT))

## If we want to retain the original LEAF labels
ggplot(tobacco,  aes(y = NUMBER,  x = as.numeric(LEAF))) +
  geom_blank(aes(x = LEAF)) +
  geom_line(aes(linetype = TREATMENT))


## ----tobaccoEDA4, results='markdown', eval=TRUE, hidden=TRUE------------------
ggplot(tobacco,  aes(y = NUMBER,  x = TREATMENT,  group = LEAF)) +
  geom_point() +
  geom_line(aes(x = as.numeric(TREATMENT))) 


## ----fitModel1a, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
tobacco.rstanarm <- stan_glmer(NUMBER ~ (1|LEAF) + TREATMENT,
                               data = tobacco,
                               family = gaussian(), 
                               iter = 5000,
                               warmup = 2000,
                               chains = 3,
                               thin = 5,
                               refresh = 0)


## ----fitModel1b, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
tobacco.rstanarm %>% prior_summary()


## ----fitModel1c, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
2.5*sd(tobacco$NUMBER)


## ----fitModel1d, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
2.5*sd(tobacco$NUMBER)/sd(model.matrix(~TREATMENT, tobacco)[, 2])


## ----fitModel1e, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
1/sd(tobacco$NUMBER)


## ----fitModel1f, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
tobacco.rstanarm1 <- update(tobacco.rstanarm,  prior_PD=TRUE)


## ----fitModel1g, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
ggpredict(tobacco.rstanarm1) %>% plot(add.data=TRUE)


## ----fitModel1h, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-------
tobacco.rstanarm2 <- stan_glmer(NUMBER ~ (1|LEAF) + TREATMENT,
                                data = tobacco,
                                family = gaussian(), 
                                prior_intercept = normal(35, 7, autoscale = FALSE),
                                prior = normal(0, 13, autoscale = FALSE),
                                prior_aux=rstanarm::exponential(0.15, autoscale = FALSE),
                                prior_covariance = decov(1, 1, 1, 1), 
                                prior_PD = TRUE, 
                                iter = 5000,
                                warmup = 1000,
                                chains = 3,
                                thin = 5,
                                refresh = 0
                                )


## ----fitModel1i, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE------
tobacco.rstanarm2 %>%
    ggpredict() %>%
    plot(add.data = TRUE)


## ----fitModel1j, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE, dependson='fitModel1h'----
tobacco.rstanarm3 <- update(tobacco.rstanarm2,  prior_PD=FALSE)


## ----modelFit1k, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=8----
posterior_vs_prior(tobacco.rstanarm3, color_by='vs', group_by=TRUE,
                   facet_args=list(scales='free_y'))


## ----modelFit1l, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggemmeans(tobacco.rstanarm3,  ~TREATMENT) %>% plot(add.data=TRUE)
ggpredict(tobacco.rstanarm3,  ~TREATMENT) %>% plot(add.data=TRUE)


## ----fitModel2a, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE, paged.print=FALSE, tidy.opts = list(width.cutoff = 80), echo=c(-4,-6)----
tobacco.form <- bf(NUMBER ~ (1|LEAF) + TREATMENT,
                   family = gaussian() 
                   )
options(width=100)
tobacco.form %>% get_prior(data=tobacco)
options(width=80)
## tobacco.brm <- brm(tobacco.form,
##                   data=tobacco,
##                   iter = 5000,
##                   warmup = 1000,
##                   chains = 3,
##                   thin = 5,
##                   refresh = 0)


## ----fitModel2h, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, fig.width = 10, fig.height = 7----
tobacco %>% 
    group_by(TREATMENT) %>%
    summarise(median(NUMBER),
              mad(NUMBER))
sd(tobacco$NUMBER)   # 6.5
standist::visualize("normal(31,7)", xlim=c(-10,100))
standist::visualize("normal(0, 13)", xlim=c(-20,20))
standist::visualize(
              "student_t(3,0,6.5)",
              "gamma(2,1)",
                    "gamma(2,0.5)",
                    "gamma(5,0.1)",
                    "exponential(1)",
                    "exponential(0.15)",
                    "cauchy(0,6.5)",
                    "cauchy(0, 2.5)",  # since sqrt(6.5) = 2.5
                    "cauchy(0,1)",
                    xlim=c(-10,25))


## ----fitModel2h1, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE------
priors <- prior(normal(31,7), class = 'Intercept') +
    prior(normal(0, 13), class = 'b') +
    prior(student_t(3,0,6.5), class = 'sigma') +
    prior(cauchy(0, 2.5), class = 'sd') 
tobacco.form <- bf(NUMBER ~ (1|LEAF) + TREATMENT,
                     family = gaussian()
                   )
tobacco.brm2 <- brm(tobacco.form, 
                  data = tobacco,
                  prior = priors,
                  sample_prior = 'only',
                  iter = 5000,
                  warmup = 2500,
                  chains = 3, cores = 3,
                  thin = 5,
                  refresh = 0,
                  backend = "cmdstanr"
                  )


## ----partialPlot2h1a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm2 %>%
    ggpredict() %>%
    plot(add.data = TRUE)


## ----fitModel2h1b, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-----
tobacco.brm3 <- update(tobacco.brm2,  
                       sample_prior = 'yes',
                       control = list(adapt_delta = 0.99),
                       refresh = 0) 
save(tobacco.brm3, file = '../ws/testing/tobacco.brm3')


## ----partialPlot2h1b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm3 %>%
    ggpredict() %>%
    plot(add.data = TRUE)


## ----posterior2h2, results='markdown', eval=TRUE------------------------------
tobacco.brm3 %>% get_variables()
tobacco.brm3 %>% hypothesis('TREATMENTWeak=0') %>% plot


## ----posterior2h2a, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
tobacco.brm3 %>% SUYR_prior_and_posterior()


## ----fitModel2h3, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE------
priors <- prior(normal(31,7), class = 'Intercept') +
    prior(normal(0, 13), class = 'b') +
    prior(student_t(3,0,6.5), class = 'sigma') +
    prior(cauchy(0, 2.5), class = 'sd') +
    prior(lkj_corr_cholesky(1), class = 'cor')
tobacco.form <- bf(NUMBER ~ (TREATMENT|LEAF) + TREATMENT,
                     family = gaussian()
                   )
 
tobacco.brm4 <-  brm(tobacco.form, 
                  data = tobacco,
                  prior = priors,
                  sample_prior = 'yes',
                  iter = 5000,
                  warmup = 1000,
                  chains = 3,
                  thin = 5,
                  refresh = 0,
                  control = list(adapt_delta=0.99),
                  backend = "cmdstanr"
                  )
save(tobacco.brm4, file = '../ws/testing/tobacco.brm4')


## ----posterior2k, results='markdown', eval=TRUE-------------------------------
tobacco.brm4 %>% get_variables()
tobacco.brm4 %>% hypothesis('TREATMENTWeak=0') %>% plot


## ----posterior2k1, results='markdown', eval=TRUE, fig.width = 7, fig.height = 5----
tobacco.brm4 %>% SUYR_prior_and_posterior()


## ----posterior2k2, results='markdown', eval=TRUE, fig.width=10, fig.height=4----
tobacco.brm4 %>%
  posterior_samples %>%
  dplyr::select(-`lp__`) %>%
  pivot_longer(everything(), names_to = 'key') %>% 
  filter(!str_detect(key, '^r')) %>%
  mutate(Type = ifelse(str_detect(key, 'prior'), 'Prior', 'Posterior'),
         ## Class = ifelse(str_detect(key, 'Intercept'),  'Intercept',
         ##         ifelse(str_detect(key, 'b'),  'b', 'sigma')),
         Class = case_when(
             str_detect(key, '(^b|^prior).*Intercept$') ~ 'Intercept',
             str_detect(key, 'b_TREATMENT|prior_b') ~ 'TREATMENT',
             str_detect(key, 'sd') ~ 'sd',
             str_detect(key, '^cor|prior_cor') ~ 'cor',
             str_detect(key, 'sigma') ~ 'sigma'),
         Par = str_replace(key, 'b_', '')) %>%
  ggplot(aes(x = Type,  y = value, color = Par)) +
  stat_pointinterval(position = position_dodge())+
  facet_wrap(~Class,  scales = 'free')



## ----fitModel2h3a, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE-----
(l.1 <- tobacco.brm3 %>% loo()) 
(l.2 <- tobacco.brm4 %>% loo())
loo_compare(l.1, l.2)


## ----modelValidation2a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
available_mcmc()


## ----modelValidation2b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
tobacco.brm3 %>% mcmc_plot(type='trace')


## ----modelValidation2c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
tobacco.brm3 %>% mcmc_plot(type='acf_bar')


## ----modelValidation2d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
tobacco.brm3 %>% mcmc_plot(type='rhat_hist')


## ----modelValidation2e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
tobacco.brm2 %>% mcmc_plot(type='neff_hist')


## ----modelValidation2f, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
tobacco.brm3 %>% mcmc_plot(type='combo')
tobacco.brm3 %>% mcmc_plot(type='violin')


## ----modelValidation2g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
tobacco.brm3 %>% get_variables()
pars <- tobacco.brm3 %>% get_variables()
pars <- str_extract(pars, '^b_.*|^sigma$|^sd.*') %>% na.omit()

tobacco.brm3$fit %>%
    stan_trace(pars = pars)


## ----modelValidation2h, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
tobacco.brm3$fit %>%
    stan_ac(pars = pars)


## ----modelValidation2i, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
tobacco.brm3$fit %>% stan_rhat() 


## ----modelValidation2j, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
tobacco.brm2$fit %>% stan_ess()


## ----modelValidation2k, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
tobacco.brm3$fit %>%
    stan_dens(separate_chains = TRUE, pars = pars)


## ----modelValidation2l, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=7----
tobacco.ggs <- tobacco.brm3 %>% ggs(burnin = FALSE, inc_warmup = FALSE)
tobacco.ggs %>% ggs_traceplot()


## ----modelValidation2m, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=7----
ggs_autocorrelation(tobacco.ggs)


## ----modelValidation2n, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggs_Rhat(tobacco.ggs)


## ----modelValidation2o, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggs_effective(tobacco.ggs)


## ----modelValidation2p, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggs_crosscorrelation(tobacco.ggs)


## ----modelValidation2q, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
ggs_grb(tobacco.ggs)


## ----modelValidation5a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
available_ppc()


## ----modelValidation5b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
tobacco.brm3 %>% pp_check(type = 'dens_overlay', ndraws = 100)


## ----modelValidation5c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
tobacco.brm3 %>% pp_check(type = 'error_scatter_avg')


## ----modelValidation5e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
tobacco.brm3 %>% pp_check(group = 'TREATMENT', type = 'intervals')
tobacco.brm3 %>% pp_check(group = 'TREATMENT', type = 'intervals_grouped')
tobacco.brm3 %>% pp_check(group = 'TREATMENT', type = 'violin_grouped')


## ----modelValidation5g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=4----
#library(shinystan)
#launch_shinystan(tobacco.brm2)


## ----modelValidation6a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
preds <- tobacco.brm4 %>% posterior_predict(ndraws = 250,  summary = FALSE)
tobacco.resids <- createDHARMa(simulatedResponse = t(preds),
                            observedResponse = tobacco$NUMBER,
                            fittedPredictedResponse = apply(preds, 2, median),
                            integerResponse = FALSE)
plot(tobacco.resids, quantreg = FALSE)


## ----partialPlot2d, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm3 %>%
    conditional_effects() %>%
    plot(points = TRUE)


## ----partialPlot2a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm3 %>%
    ggpredict() %>%
    plot(add.data = TRUE)


## ----partialPlot2b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm3 %>%
    ggemmeans(~TREATMENT) %>%
    plot(add.data = TRUE) +
    geom_point(data = tobacco, aes(y = NUMBER, x = as.numeric(TREATMENT)))


## ----partialPlot2c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
Partial.obs <- tobacco.brm3$data %>%
    mutate(Pred = predict(tobacco.brm3)[,'Estimate'],
           Resid = resid(tobacco.brm3)[,'Estimate'],
           Obs = Pred + Resid)

tobacco.brm3 %>%
    fitted_draws(newdata = tobacco) %>%
    median_hdci() %>%
    ggplot(aes(x = TREATMENT, y = .value)) +
    geom_pointrange(aes(ymin = .lower, ymax = .upper)) + 
    geom_line() +
    geom_point(data = Partial.obs,  aes(y = Obs,  x = TREATMENT), color = 'red',
               position = position_nudge(x = 0.1)) +
    geom_point(data = tobacco,  aes(y = NUMBER,  x = TREATMENT),
               position = position_nudge(x = 0.05))

tobacco.brm3 %>%
    epred_draws(newdata = tobacco) %>%
    ggplot() +
    geom_violin(data = tobacco, aes(y = NUMBER, x = TREATMENT), fill = 'blue', alpha = 0.2) +
    geom_point(data = tobacco, aes(y = NUMBER, x = TREATMENT),
               position = position_jitter(width = 0.1, height = 0)) +
    geom_violin(aes(y = .epred, x = TREATMENT), fill = 'orange', alpha = 0.2) +
    theme_bw()


## ----summariseModel2a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm3 %>% summary()


## ----summariseModel2a1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5, echo=FALSE----
tobacco.sum <- summary(tobacco.brm3)


## ----summariseModel2i2, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm3 %>%
  summarise_draws(
    median,
    ~ HDInterval::hdi(.x),
    rhat,
    ess_bulk,
    ess_tail
  )

## or if you want to exclude some parameters
tobacco.brm3 %>%
  summarise_draws(
    median,
    ~ HDInterval::hdi(.x),
    rhat,
    ess_bulk,
    ess_tail
  ) %>%
  filter(str_detect(variable, 'prior|^r_|^lp__', negate = TRUE)) 


## ----summariseModel2i, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm3 %>% as_draws_df()
tobacco.brm3 %>%
  as_draws_df() %>%
  summarise_draws(
    median,
    ~ HDInterval::hdi(.x),
    rhat,
    ess_bulk,
    ess_tail
  )
## or if you want to exclude some parameters
tobacco.brm3 %>%
  as_draws_df() %>%
  summarise_draws(
    median,
    ~ HDInterval::hdi(.x),
    rhat,
    ess_bulk,
    ess_tail
  ) %>%
  filter(str_detect(variable, 'prior|^r_|^lp__', negate = TRUE)) 


## ----summariseModel2b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm3$fit %>%
    tidyMCMC(estimate.method = 'median',
             conf.int = TRUE,  conf.method = 'HPDinterval',
             rhat = TRUE, ess = TRUE)

## ----summariseModel2b1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
tobacco.tidy <- tidyMCMC(tobacco.brm3$fit, estimate.method='median',
                         conf.int=TRUE,  conf.method='HPDinterval',
                         rhat=TRUE, ess=TRUE)


## ----summariseModel2c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm3 %>% get_variables()
tobacco.draw <- tobacco.brm3 %>%
    gather_draws(`b.Intercept.*|b_TREAT.*|sd_.*|sigma`,  regex=TRUE)
tobacco.draw


## ----summariseModel2c1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.draw %>% median_hdci


## ----summariseModel2c3, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=FALSE----
tobacco.gather <- tobacco.brm3 %>%
    gather_draws(`b_Intercept.*|b_TREAT.*|sd_.*|sigma`,  regex=TRUE) %>%
  median_hdci


## ----summariseModel2c4, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
tobacco.brm3 %>%
    gather_draws(`b_Intercept.*|b_TREAT.*`, regex=TRUE) %>% 
    ggplot() +
    geom_vline(xintercept=0, linetype='dashed') +
    stat_slab(aes(x = .value, y = .variable,
                  fill = stat(ggdist::cut_cdf_qi(cdf,
                           .width = c(0.5, 0.8, 0.95), 
                           labels = scales::percent_format())
                           )), color='black') + 
    scale_fill_brewer('Interval', direction = -1, na.translate = FALSE) 

tobacco.brm3 %>% 
  gather_draws(`.Intercept.*|.*TREAT.*`, regex=TRUE) %>% 
  ggplot() + 
  stat_halfeye(aes(x=.value,  y=.variable)) +
  facet_wrap(~.variable, scales='free')

tobacco.brm3 %>% 
  gather_draws(`.Intercept.*|.*TREAT.*`, regex=TRUE) %>% 
  ggplot() + 
    stat_halfeye(aes(x=.value,  y=.variable)) +
    theme_classic()


## ----summariseModel2j, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm3$fit %>% plot(type='intervals') 


## ----summariseModel2ka, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
tobacco.brm3 %>% 
    gather_draws(`^b_.*`, regex=TRUE) %>% 
    filter(.variable != 'b_Intercept') %>%
    ggplot() + 
    stat_halfeye(aes(x=.value,  y=.variable)) +
    facet_wrap(~.variable, scales='free')

tobacco.brm3 %>% 
    gather_draws(`^b_.*`, regex=TRUE) %>% 
    filter(.variable != 'b_Intercept') %>%
    ggplot() + 
    stat_halfeye(aes(x=.value,  y=.variable)) +
    geom_vline(xintercept = 0, linetype = 'dashed')


## ----summariseModel2c7, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5,echo=TRUE----
tobacco.brm3 %>% 
    gather_draws(`^b_.*`, regex=TRUE) %>% 
    filter(.variable != 'b_Intercept') %>%
    ggplot() +  
    geom_density_ridges(aes(x=.value, y = .variable), alpha=0.4) +
    geom_vline(xintercept = 0, linetype = 'dashed')
##Or in colour 
tobacco.brm3 %>% 
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
tobacco.brm3 %>% tidy_draws()


## ----summariseModel2e, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm3 %>% spread_draws(`.*Intercept.*|.*TREAT.*`,  regex=TRUE)


## ----summariseModel2f, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm3 %>% posterior_samples() %>% as_tibble()


## ----summariseModel2g, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm3 %>%
    bayes_R2(re.form = NA, summary=FALSE) %>%
    median_hdci
tobacco.brm3 %>%
    bayes_R2(re.form = ~(1|LEAF), summary=FALSE) %>%
    median_hdci
## tobacco.brm3 %>%
##     bayes_R2(re.form = ~(TREATMENT|LEAF), summary=FALSE) %>%
##     median_hdci


## ----summariseModel2k, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
0.1 * sd(tobacco$NUMBER)
tobacco.brm3 %>% rope(range = c(-0.65, 0.65))
rope(tobacco.brm3, range = c(-0.65, 0.65)) %>% plot()

## Or based on fractional scale
tobacco.brm3 %>% emmeans(~TREATMENT) %>%
    gather_emmeans_draws() %>%
    group_by(.draw) %>%
    arrange(desc(TREATMENT)) %>%
    summarise(Diff = 100*(exp(diff(log(.value))) -1)) %>%
    rope(range = c(-10,10))



## ----predictions2a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
tobacco.brm3 %>% emmeans(~TREATMENT) %>%
    gather_emmeans_draws() %>%
    group_by(.draw) %>%
    arrange(desc(TREATMENT)) %>%
    summarise(Diff = 100*(exp(diff(log(.value))) -1)) %>%
    summarise_draws(
        median,
        ~ HDInterval::hdi(.x),
        rhat,
        ess_bulk,
        ess_tail
        )

## Or via gather and pivot
newdata <- tobacco.brm3 %>%
    emmeans(~TREATMENT) %>%
    gather_emmeans_draws() %>%
    pivot_wider(names_from=TREATMENT,values_from=.value) %>%
    mutate(Eff = Strong - Weak,
           PEff = 100*Eff/Weak)
newdata %>% median_hdci(PEff)
newdata %>% summarise(P = sum(PEff>0)/n())
newdata %>% summarise(P = sum(PEff>20)/n())
newdata %>%
    dplyr::select(-.chain, -.iteration) %>%
    hypothesis('PEff>20')

newdata <- tobacco.brm3 %>% emmeans(~TREATMENT) %>% as.data.frame
head(newdata)
ggplot(newdata, aes(y=emmean, x=TREATMENT)) +
    geom_pointrange(aes(ymin=lower.HPD, ymax=upper.HPD)) +
    theme_bw()

tobacco.brm3 %>%
    emmeans(~TREATMENT) %>%
    gather_emmeans_draws() %>%
    ggplot() +
    geom_density_ridges(aes(x = .value, y = TREATMENT), alpha = 0.5, fill = 'orange') +
    scale_x_continuous("Average number of lesions") + 
    theme_bw()



## ----predictions2b, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
newdat <- tobacco %>% tidyr::expand(TREATMENT)
newdata <- tobacco.brm3 %>%
    brms::posterior_epred(newdat, re_formula = NA) %>%
    as.data.frame() %>%
    rename_with(~as.character(newdat$TREATMENT)) %>% 
    mutate(Eff = Strong - Weak,
           PEff = 100*Eff/Weak)
head(newdata)
newdata %>% median_hdci(PEff)
newdata %>% summarise(P = sum(PEff>0)/n())
newdata %>% summarise(P = sum(PEff>20)/n())


## ----predictions2c, results='markdown', eval=TRUE, hidden=TRUE, fig.width=8, fig.height=5----
newdata <- tobacco.brm3 %>%
    emmeans(~TREATMENT) %>%
    pairs() %>%
    gather_emmeans_draws()
newdata %>% median_hdci()

## OR on percentage scale
newdata <- tobacco.brm3 %>%
    emmeans(~TREATMENT) %>%
    regrid(trans = 'log') %>%
    pairs() %>%
    regrid() %>% 
    gather_emmeans_draws() %>%
    mutate(.value = (.value - 1) * 100)
newdata %>% median_hdci()
newdata %>% summarise(P = sum(.value>0)/n())
newdata %>% summarise(P = sum(.value>20)/n())

