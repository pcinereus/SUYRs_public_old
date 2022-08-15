## SUYR_prior_and_posterior_old <- function(mod) {
##     dat <- mod$data
##     terms <- attr(dat, 'terms')
##     response <- all.vars(update(terms, .~ 1))
##     predictors <- all.vars(terms)[-1]

##     Xmat <- dat %>%
##         dplyr::select(any_of(predictors)) %>%
##         colMeans() %>%
##         as.matrix()
##     b <- mod %>%
##         as_draws_df() %>%
##         dplyr::select(starts_with('b_'),
##                       -contains('Intercept')) %>%
##         as.matrix()

##     scal <- as.vector(Xmat %*% t(b))

##     priors <- mod %>% get_variables() %>% str_subset("^prior_.*") 
##     pars <- mod %>% get_variables() %>%
##         str_subset(paste0("^b_(Intercept|",paste0(predictors,collapse="|"),")"))
##     aux <- priors %>% str_subset("prior_(Intercept|b)", negate = TRUE) %>%
##         str_remove("^prior_")
##     pars <- c(pars, aux)
##     ## ## nms <- mod %>% get_variables() %>% str_detect(pars)

##     mod.pp <- mod %>%
##         as_draws_df() %>%
##         select(any_of(c(pars, priors))) %>%
##         mutate(b_Intercept = b_Intercept + scal) %>%
##         pivot_longer(cols=everything(), names_to='key', values_to='value') %>% 
##         mutate(Type = ifelse(str_detect(key, 'prior'), 'Prior', 'Posterior'),
##                Parameter = ifelse(Type == 'Prior',
##                                   str_remove(key, "^prior_"),
##                                   str_remove(key, "^b_")
##                                   ),
##                Class = ifelse(Parameter %in% brms_names, 'b', Parameter))

##     return(
##         ggplot(data = NULL, aes(x=Type,  y=value)) +
##         stat_pointinterval(data = mod.pp %>% filter(Class != 'b' | Parameter == 'b')) +
##         stat_pointinterval(data = mod.pp %>% filter(Class == 'b' & Parameter != 'b'),
##                            aes(colour = Parameter))+
##         facet_wrap(~Class,  scales='free')
##         )
## }


## SUYR_prior_and_posterior <- function(mod) {
##     dat <- mod$data
##     terms <- attr(dat, 'terms')
##     response <- all.vars(update(terms, .~ 1))
##     predictors <- all.vars(terms)[-1]
##     rhs <- mod$formula %>% as.formula %>% brms:::str_rhs() #%>%
##         ## str_split('\\+', simplify = TRUE) %>%
##         ## str_trim()

##     Xmat <- model.matrix(as.formula(paste0('~',rhs)), dat)[,-1] %>%
##         colMeans() %>%
##         as.matrix()
    
##     ## Xmat <- dat %>%
##     ##     select(any_of(rhs)) %>%
##     ##     summarise(across(is.numeric, mean)) %>%
##     ##     as.matrix() 

##     b <- mod %>%
##         as_draws_df() %>%
##         dplyr::select(starts_with('b_'),
##                       -contains('Intercept')) %>%
##         as.matrix()
    
##     scal <- as.vector(t(Xmat) %*% t(b))
    
##     brms_names <- brms:::change_effects(brmsterms(mod$formula),
##                                         data = dat,
##                                         pars = variables(mod))
##     brms_names <- sapply(brms_names, function(x) str_remove(x$fnames, "b_"))
##     priors <- mod %>% get_variables() %>% str_subset("^prior_.*") 
##     pars <- mod %>% get_variables() %>%
##         str_subset(paste0("^b_(Intercept|",paste0(brms_names,collapse="|"),")"))
##     aux <- priors %>% str_subset("prior_(Intercept|b)", negate = TRUE) %>%
##         str_remove("^prior_")
##     pars <- c(pars, aux)

##     mod.pp <- mod %>%
##         as_draws_df() %>%
##         select(any_of(c(pars, priors))) %>%
##         mutate(b_Intercept = b_Intercept + scal) %>%
##         pivot_longer(cols=everything(), names_to='key', values_to='value') %>% 
##         mutate(Type = ifelse(str_detect(key, 'prior'), 'Prior', 'Posterior'),
##                Parameter = ifelse(Type == 'Prior',
##                                   str_remove(key, "^prior_"),
##                                   str_remove(key, "^b_")
##                                   ),
##                Class = ifelse(Parameter %in% brms_names, 'b', Parameter))

##     return(
##         ggplot(data = NULL, aes(x=Type,  y=value)) +
##         stat_pointinterval(data = mod.pp %>% filter(Class != 'b' | Parameter == 'b')) +
##         stat_pointinterval(data = mod.pp %>% filter(Class == 'b' & Parameter != 'b'),
##                            aes(colour = Parameter))+
##         facet_wrap(~Class,  scales='free')
##         )
## }

SUYR_prior_and_posterior <- function(mod) {
    dat <- mod$data
    terms <- attr(dat, 'terms')
    response <- all.vars(update(terms, .~ 1))
    predictors <- all.vars(terms)[-1]
    ## rhs <- mod$formula %>% as.formula %>% brms:::str_rhs()

    f <- mod$formula %>% as.formula %>% update(NULL ~.)
    rhs  <-
        deparse1(f) %>%
        str_remove("~") %>%
        ## paste(f[2],f[3],sep='~') %>%
        str_split("\\+") %>%
        unlist() %>%
        str_trim()
    rhs
    ## ## exclude any terms with a "|"
    ## rhs <-
    ##     rhs[-grep("\\|",rhs)]
    wch.rnd <- rhs[grep("\\|", rhs)]
    if (length(wch.rnd)>0) f <- update(f, paste("~ . -",wch.rnd))
    
    Xmat <- model.matrix(f, dat)[,-1] %>%
        as.matrix() %>% 
        colMeans() 
    ## if (length(Xmat)==1) Xmat <- Xmat %>% as.matrix()
    ## Xmat <- dat %>%
    ##     select(any_of(rhs)) %>%
    ##     summarise(across(everything(), mean)) %>%
    ##     as.matrix()

    b <- mod %>%
        as_draws_df() %>%
        dplyr::select(starts_with('b_'),
                      -contains('Intercept')) %>%
        as.matrix()
    
    scal <- as.vector(Xmat %*% t(b))
    
    ## fixed effects
    brms_names <- brms:::change_effects(brmsterms(mod$formula),
                                        data = dat,
                                        pars = variables(mod))
    brms_names <- sapply(brms_names, function(x) str_remove(x$fnames, "b_"))
    priors <- mod %>% get_variables() %>% str_subset("^prior_.*") 
    pars <- mod %>% get_variables() %>%
        str_subset(paste0("^b_(Intercept|",paste0(brms_names,collapse="|"),")"))
    
    ## auxillary
    aux <- priors %>% str_subset("prior_(Intercept|b)", negate = TRUE) %>%
        str_remove("^prior_")
    pars <- c(pars, aux)
    ## random effects
    if (length(wch.rnd)>0) {
        ran.pars <- brms:::change_re(mod$ranef, pars = variables(mod))[[1]]$fnames
        pars <- c(pars, ran.pars)
    }
    variables(mod)

    vars <- variables(mod)
    priors <- vars %>% str_subset("prior")
    all.pars <- priors %>% str_remove("prior_")
    fixed.pars <- vars %>% str_subset("^b_")
    other.pars <- all.pars %>% str_subset("^Intercept$|^b$", negate = TRUE)
    other.pars <- vars %>% str_subset(paste0("^", other.pars, collapse = '|'))
    pars <- c(fixed.pars, other.pars)
    
    ## coefs <- prior_summary(mod)$class %>% unique()
    ## coefs.regex <- paste0("^b_", coefs, collapse = "|")
    
    mod.pp <- mod %>%
        as_draws_df() %>%
        select(any_of(c(pars, priors))) %>%
        mutate(b_Intercept = b_Intercept + scal) %>%
        pivot_longer(cols=everything(), names_to='key', values_to='value') %>% 
        mutate(Type = ifelse(str_detect(key, 'prior'), 'Prior', 'Posterior'),
               Parameter = ifelse(Type == 'Prior',
                                  str_remove(key, "^prior_"),
                                  str_remove(key, "^b_")
                                  ),
               ## Parameter = ifelse(Type == 'Posterior',
               ##                     str_remove(Parameter, "__.*"),
               ##                    Parameter),
               Class = ifelse(Parameter %in% brms_names, 'b', Parameter),
               Class = ifelse(Type == 'Posterior', str_remove(Class, "__.*"), Class))

    return(
        ggplot(data = NULL, aes(x=Type,  y=value)) +
        stat_pointinterval(data = mod.pp %>% filter(Type == 'Prior')) +
        stat_pointinterval(data = mod.pp %>% filter(Type != 'Prior' & (Class != 'b' | Parameter == 'b')),
                           aes(colour = Parameter), position = position_dodge()) +
        stat_pointinterval(data = mod.pp %>% filter(Type != 'Prior' & (Class == 'b' & Parameter != 'b')),
                           aes(colour = Parameter), position = position_dodge())+
        facet_wrap(~Class,  scales='free')
        )
}
