
nBoot <- 10
fish.pred <- with(env,
                     expand.grid(DAS = modelr::seq_range(DAS, n = 100),
                                ALT = NA,
                                PEN = NA,
                                DEB = NA,
                                PH = NA,
                                DUR = NA,
                                PHO = NA,
                                NIT = NA,
                                AMM = NA,
                                OXY = NA,
                                DBO = NA)
                     )
fish.list <- vector('list', nBoot) 
fish.list
fish.sum <- vector('list', nBoot) 
for (i in 1:nBoot) {
    print(paste0('Boot number: ', i))
    ## Create random set
    fish.rnd <- env %>%
        sample_n(size = n(), replace=TRUE)
    ## Fit the trees
    fish.gbm = gbm(fish ~ DAS + ALT + PEN + DEB + PH + DUR + PHO + NIT + AMM +
              OXY + DBO,
              data=fish.rnd,
              distribution='poisson',
              var.monotone = c(1,-1,-1,1,0,1,0,0,0,0,0),
              n.minobsinnode = 2,
              n.trees=10000,
              interaction.depth=5,
              bag.fraction=0.5,
              shrinkage=0.001,
              train.fraction=1,
              cv.folds=3)
    ## Determine the best number of trees
    (best.iter = gbm.perf(fish.gbm,method='cv'))
    ## predict based on shell weight
    fit <- predict(fish.gbm, newdata = fish.pred, n.trees = best.iter) %>% exp()
    fish.list[[i]] <- data.frame(fish.pred, Boot = i, Fit = fit)
    ## relative influence
    fish.sum[[i]] <- summary(fish.gbm, n.trees = best.iter)
}
fish.fit <- do.call('rbind', fish.list)
fish.fit <- fish.fit %>%
    group_by(DAS, ALT) %>%
    ## summarise(Median = median(Fit),
    ##           Lower = quantile(Fit, p=0.025),
    ##           Upper = quantile(Fit, p=0.975))
    ggdist::median_hdci(Fit)       
g1 <- fish.fit %>% ggplot(aes(y=Fit, x=DAS, fill=ALT, color=ALT)) +
    geom_ribbon(aes(ymin=.lower, ymax=.upper), alpha=0.3, color=NA) +
    geom_line() +
    scale_fill_viridis_d() +
    scale_colour_viridis_d() +
    theme_classic()

fish.inf <- do.call('rbind', fish.sum)
fish.inf <- fish.inf %>%
    group_by(var) %>%
    ggdist::median_hdci(rel.inf)       

g2 <- fish.inf %>% ggplot(aes(y=var, x=rel.inf)) +
    geom_vline(xintercept=12.5, linetype='dashed') +
    geom_pointrange(aes(xmin=.lower, xmax=.upper)) +
    theme_classic()

g2 + patchwork::inset_element(g1, left=0.5, bottom=0.01, right=1, top=0.7)
g1 + patchwork::inset_element(g2, left=0.5, bottom=0.01, right=1, top=0.5)
