## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE)


## ----libraries, results='markdown', eval=TRUE---------------------------------
library(gbm)         #for gradient boosted models
library(car)
library(dismo)
library(pdp)
library(ggfortify)
library(randomForest)
library(tidyverse)
library(gridExtra)
library(patchwork)


## ----readData, results='markdown', eval=TRUE----------------------------------
fish = read_csv('../data/verneaux.fish.csv', trim_ws=TRUE)
glimpse(fish)

env = read_csv('../data/verneaux.env.csv', trim_ws=TRUE)
glimpse(env)


## ----processData, results='markdown', eval=TRUE-------------------------------
fish <- fish %>% rowSums


## ----EDA, results='markdown', eval=TRUE, hidden=TRUE, fig.width=15, fig.height=15----
car::scatterplotMatrix(cbind(fish, env))


## ----fitModel2, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE--------
fish.gbm = gbm(fish ~ DAS + ALT + PEN + DEB + PH + DUR + PHO + NIT + AMM +
              OXY + DBO,
              data=env,
              var.monotone = c(1,-1,-1,1,0,1,0,0,0,0,0),
              distribution='poisson',
              n.trees=10000,
              n.minobsinnode = 2,
              interaction.depth=5,
              bag.fraction=0.5,
              shrinkage=0.01,
              train.fraction=1,
              cv.folds=3)


## ----fitModel3, results='markdown', eval=TRUE, hidden=TRUE--------------------
(best.iter = gbm.perf(fish.gbm,method='OOB'))
(best.iter = gbm.perf(fish.gbm,method='cv'))


## ----fitModel4, results='markdown', eval=TRUE, hidden=TRUE, cache=TRUE--------
fish.gbm = gbm(fish ~ DAS + ALT + PEN + DEB + PH + DUR + PHO + NIT + AMM +
              OXY + DBO,
              data=env,
              var.monotone = c(1,-1,-1,1,0,1,0,0,0,0,0),
              distribution='poisson',
              n.trees=10000,
              n.minobsinnode = 2,
              interaction.depth=5,
              bag.fraction=0.5,
              shrinkage=0.001,
              train.fraction=1,
              cv.folds=3)


## ----fitModel5, results='markdown', eval=TRUE, hidden=TRUE--------------------
(best.iter = gbm.perf(fish.gbm,method='OOB'))
(best.iter = gbm.perf(fish.gbm,method='cv'))


## ----relativeInfluence1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=7, fig.height=10----
summary(fish.gbm, n.trees=best.iter)


## ----partialEffects1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6----
attr(fish.gbm$Terms,"term.labels")
plot(fish.gbm, 1, n.tree=best.iter)
plot(fish.gbm, 2, n.tree=best.iter)
fish.gbm %>%
    pdp::partial(pred.var='DAS',
                 n.trees=best.iter,
                 recursive=FALSE,
                 inv.link=exp) %>%
    autoplot()


## ----partialEffects2, results='markdown', eval=TRUE, hidden=TRUE, fig.width=10, fig.height=10----
nms <- attr(fish.gbm$Terms,"term.labels")
p <- vector('list', length(nms))
names(p) <- nms
for (nm in nms) {
  print(nm)
  p[[nm]] <- fish.gbm %>% pdp::partial(pred.var=nm,
                                 n.trees=best.iter,
                                 inv.link=exp,
                                 recursive=FALSE,
                                 type='regression') %>%
      autoplot() +
      ylim(18, 60)
}
patchwork::wrap_plots(p) 
#do.call('grid.arrange', p)


## ----partialEffects3, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE----
fish.gbm %>%
    pdp::partial(pred.var=c('DAS'),
                 n.trees=best.iter, recursive=FALSE, inv.link=exp) %>%
    autoplot()

fish.gbm %>%
    pdp::partial(pred.var=c('DAS','ALT'),
                 n.trees=best.iter, recursive=TRUE) %>%
    autoplot()

g1 = fish.gbm %>% pdp::partial(pred.var='DAS', n.trees=best.iter,
                            recursive=FALSE,inv.link=exp) %>%
    autoplot
g2 = fish.gbm %>% pdp::partial(pred.var='ALT', n.trees=best.iter,
                            recursive=FALSE,inv.link=exp) %>%
    autoplot


g1 + g2


## ----Accuracy1, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, fig.width=7, fig.height=7----
fish.acc <- env %>%
    bind_cols(Pred = predict(fish.gbm,
                             newdata=env,
                             n.tree=best.iter,
                             type='response'))

with(fish.acc,  cor(fish, Pred))


fish.acc %>%
  ggplot() +
  geom_point(aes(y=Pred,  x=fish))



## ----Interactions1, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, fig.width=7, fig.height=7----
attr(fish.gbm$Terms,"term.labels")
 
interact.gbm(fish.gbm, env,c(1,2), n.tree=best.iter)
interact.gbm(fish.gbm, env,c(1,10), n.tree=best.iter)
interact.gbm(fish.gbm, env,c(2,10), n.tree=best.iter)
interact.gbm(fish.gbm, env,c(1,2,10), n.tree=best.iter)

fish.gbm %>% pdp::partial(pred.var=c('DAS', 'OXY'),  n.trees=best.iter, recursive=FALSE) %>%
    autoplot
fish.gbm %>% pdp::partial(pred.var=c('DAS', 'OXY'),
                          n.trees=best.iter,
                          recursive=FALSE,
                          inv.link = exp) %>%
    autoplot
fish.gbm %>% pdp::partial(pred.var=c(1, 10),  n.trees=best.iter, recursive=FALSE) %>% autoplot
fish.gbm %>% pdp::partial(pred.var=c(2, 10),  n.trees=best.iter, recursive=FALSE) %>% autoplot



## ----Interactions2, results='markdown', eval=TRUE, hidden=TRUE, cache=FALSE, fig.width=7, fig.height=7----

fish.grid = plot(fish.gbm, c(1,10), n.tree=best.iter, return.grid=TRUE)
head(fish.grid)

ggplot(fish.grid, aes(y=DAS, x=OXY)) +
    geom_tile(aes(fill=y)) +
    geom_contour(aes(z=y)) +
    scale_fill_gradientn(colors=heat.colors(10))



## ----Interactions3, eval=TRUE, echo=FALSE-------------------------------------
terms <- attr(fish.gbm$Terms,"term.labels")
fish.int <- NULL
for (i in 1:(length(terms)-1)) {
    for (j in (i+1):length(terms)) {
        print(paste('i=',i, ' Name = ', terms[i]))
        print(paste('j=',j, ' Name = ', terms[j]))
        fish.int <- rbind(fish.int,
                             data.frame(Var1=terms[i], Var2=terms[j],
                                        "H.stat"=interact.gbm(fish.gbm, env,c(i,j),
                                                              n.tree=best.iter)
                                        ))
    }
}
fish.int %>% arrange(-H.stat)


## ----gbmstep1, eval=FALSE, hidden=TRUE, cache=FALSE---------------------------
## fish.gbm1 <- dismo::gbm.step(data=cbind(fish, env) %>% as.data.frame, gbm.x=2:11, gbm.y=1,
##                         tree.complexity=5,
##                         learning.rate=0.001,
##                       n.minobsinnode = 2,
##                         bag.fraction=0.5,
##                       n.train = 1,
##                         n.trees=10000,
##                         family='poisson')


## ----gbmstep2, eval=FALSE, hidden=TRUE, cache=FALSE---------------------------
## summary(abalone.gbm1)


## ----randomForest, results='markdown', eval=TRUE, hidden=TRUE-----------------
library(randomForest)
fish.rf = randomForest(fish ~ DAS + ALT + PEN + DEB + PH + DUR + PHO + NIT + AMM +
                           OXY + DBO,
                       data=env, importance=TRUE,
                       ntree=1000)
fish.imp = randomForest::importance(fish.rf)
## Rank by either:
## *MSE (mean decrease in accuracy)
## For each tree, calculate OOB prediction error.
## This also done after permuting predictors.
## Then average diff of prediction errors for each tree
## *NodePurity (mean decrease in node impurity)
## Measure of the total decline of impurity due to each
## predictor averaged over trees
100*fish.imp/sum(fish.imp)
varImpPlot(fish.rf)
## use brute force
fish.rf %>% pdp::partial('DAS') %>% autoplot


fish.rf.acc <- env %>%
    bind_cols(Pred = predict(fish.rf,
                             newdata=env))

with(fish.rf.acc,  cor(fish, Pred))


fish.rf.acc %>%
  ggplot() +
  geom_point(aes(y=Pred,  x=fish)) +
  geom_point(data = fish.acc, aes(y=Pred,  x=fish), colour = 'red')

