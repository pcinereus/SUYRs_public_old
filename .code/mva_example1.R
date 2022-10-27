## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE,cache.lazy = FALSE, tidy='styler')


## ----libraries, results='markdown', eval=TRUE---------------------------------
library(tidyverse)  #for data wrangling etc
library(vegan)
library(ggvegan)
library(ggrepel)
library(GGally)
library(corrplot)
library(EcolUtils)
library(car)
library(scales)
library(patchwork)


## ----readData, results='markdown', eval=TRUE----------------------------------
spider.abund <- read_csv(file = "../data/spider.abund.csv", trim_ws = TRUE)
spider.abund %>% glimpse()


## ----EDA, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.abund %>%
    cor %>%
    corrplot(type = 'upper',
             diag = FALSE)
## And now with axes arrange according to first princomp
spider.abund %>%
    cor %>%
    corrplot(type = 'upper',
             order = 'FPC',
             diag = FALSE)


## ----EDA1, results='markdown', cache = TRUE, eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.abund %>%
    ggpairs(lower = list(continuous = "smooth"),
            diag = list(continuous = "density"),
            axisLabels = "show")
## - clearly not normal

spider.abund^0.25 %>%
    ggpairs(lower = list(continuous = "smooth"),
            diag = list(continuous = "density"),
            axisLabels = "show")
## - still not normal


## ----DCA, results='markdown', cache = TRUE, eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.dca <- decorana(spider.abund^0.25)
spider.dca


## ----scaling, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.abund %>% head

## Standardize columns to maximums
spider.abund %>%
    head() %>% 
    decostand(method = "max", MARGIN = 2) 

## spider.abund %>%
##     head() %>% 
##     decostand(method = "range", MARGIN = 2) 

## spider.abund %>%
##     head() %>% 
##     decostand(method = "normalize", MARGIN = 2) 

## Standardize rows to totals
spider.abund %>%
    head() %>% 
    decostand(method = "total", MARGIN = 1) 

## Double standardization
spider.abund %>%
    head() %>% 
    wisconsin() 

## Transformations
spider.abund %>%
    sqrt() %>%
    head() %>% 
    wisconsin() 



## ----PCA1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.rda <- rda(spider.abund, scale=TRUE)
summary(spider.rda, display=NULL)


## ----screePlot1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
screeplot(spider.rda)
abline(a=1,b=0)


## ----biplot, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
## Quick and nasty ordination plots
biplot(spider.rda, scaling='species')
biplot(spider.rda, scaling='sites')


## ----biplot1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
## Quick and nasty ordination plots
pl<-vegan::ordiplot(spider.rda)
points(pl, "sites", pch=21, col="red", bg="yellow")
text(pl, "sites", col="red", cex=0.9)
text(pl, "species", col="blue", cex=0.9)
## line(pl, "sites")


## ----biplot3, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
## ggvegan provides ways to use ggplot
## library(devtools)
## devtools::install_github("gavinsimpson/ggvegan")
## library(ggvegan)
autoplot(spider.rda)
autoplot(spider.rda) + theme_bw()
autoplot(spider.rda,geom='text') + theme_bw()


## ----biplot4, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.rda.scores <- spider.rda %>%
    fortify()
spider.rda.scores

g <-
    ggplot(data = NULL, aes(y=PC2, x=PC1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=spider.rda.scores %>% filter(Score=='sites')) +
    geom_text(data=spider.rda.scores %>% filter(Score=='sites'),
              aes(label=Label), hjust=-0.2) +
    geom_segment(data=spider.rda.scores %>% filter(Score=='species'),
                 aes(y=0, x=0, yend=PC2, xend=PC1),
                 arrow=arrow(length=unit(0.3,'lines')), color='red') +
    ## geom_text(data=spider.rda.scores %>% filter(Score=='species'),
    ##           aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red') +
    geom_text_repel(data=spider.rda.scores %>% filter(Score=='species'),
                    aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red') +
    theme_bw()
g

                                        # Nice axes titles
eig <- eigenvals(spider.rda)

paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100 * eig[2]/sum(eig)))
g <- g +
    scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)',
                                                    100 * eig[2]/sum(eig))))+
    scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)',
                                                    100 * eig[1]/sum(eig))))

                                        #put a circle
circle.prob <- 0.68
## circle.prob <- 0.95
## circle.prob <- 0.95
r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(spider.rda$CA$u[,1:2]^2))^(1/4)
theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
circle <- data.frame(PC1 = r * cos(theta), PC2 = r * sin(theta))
g <- g + geom_path(data = circle, aes(y=PC2,x=PC1), color = muted('white'), size = 1/2, alpha = 1/3)
g
g1 <- g


## ----PCA2, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.std <- spider.abund %>%
    mutate(across(everything(), function(x) x^0.25)) %>%
    wisconsin()
spider.std
spider.rda <- rda(spider.std, scale=FALSE)
summary(spider.rda, display=NULL)


## ----biplot5, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.rda.scores <- spider.rda %>%
    fortify()
spider.rda.scores

g <-
    ggplot(data = NULL, aes(y=PC2, x=PC1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=spider.rda.scores %>% filter(Score=='sites')) +
    geom_text(data=spider.rda.scores %>% filter(Score=='sites'),
              aes(label=Label), hjust=-0.2) +
    geom_segment(data=spider.rda.scores %>% filter(Score=='species'),
                 aes(y=0, x=0, yend=PC2, xend=PC1),
                 arrow=arrow(length=unit(0.3,'lines')), color='red') +
    ## geom_text(data=spider.rda.scores %>% filter(Score=='species'),
    ##           aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red') +
    geom_text_repel(data=spider.rda.scores %>% filter(Score=='species'),
                    aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red') +
    theme_bw()
g

                                        # Nice axes titles
eig <- eigenvals(spider.rda)

paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100 * eig[2]/sum(eig)))
g <- g +
    scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)',
                                                    100 * eig[2]/sum(eig))))+
    scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)',
                                                    100 * eig[1]/sum(eig))))

                                        #put a circle
circle.prob <- 0.68
## circle.prob <- 0.95
## circle.prob <- 0.95
r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(spider.rda$CA$u[,1:2]^2))^(1/4)
theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
circle <- data.frame(PC1 = r * cos(theta), PC2 = r * sin(theta))
g <- g + geom_path(data = circle, aes(y=PC2,x=PC1), color = muted('white'), size = 1/2, alpha = 1/3)
g

g2 <- g


## ----PCA3, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.hel <- spider.abund %>%
     decostand(method = "hellinger")
spider.rda <- rda(spider.hel, scale=FALSE)
summary(spider.rda, display=NULL)


## ----biplot6, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.rda.scores <- spider.rda %>%
    fortify()
spider.rda.scores

g <-
    ggplot(data = NULL, aes(y=PC2, x=PC1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=spider.rda.scores %>% filter(Score=='sites')) +
    geom_text(data=spider.rda.scores %>% filter(Score=='sites'),
              aes(label=Label), hjust=-0.2) +
    geom_segment(data=spider.rda.scores %>% filter(Score=='species'),
                 aes(y=0, x=0, yend=PC2, xend=PC1),
                 arrow=arrow(length=unit(0.3,'lines')), color='red') +
    ## geom_text(data=spider.rda.scores %>% filter(Score=='species'),
    ##           aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red') +
    geom_text_repel(data=spider.rda.scores %>% filter(Score=='species'),
                    aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red') +
    theme_bw()
g

                                        # Nice axes titles
eig <- eigenvals(spider.rda)

paste(names(eig[2]), sprintf('(%0.1f%% explained var.)', 100 * eig[2]/sum(eig)))
g <- g +
    scale_y_continuous(paste(names(eig[2]), sprintf('(%0.1f%% explained var.)',
                                                    100 * eig[2]/sum(eig))))+
    scale_x_continuous(paste(names(eig[1]), sprintf('(%0.1f%% explained var.)',
                                                    100 * eig[1]/sum(eig))))

                                        #put a circle
circle.prob <- 0.68
## circle.prob <- 0.95
## circle.prob <- 0.95
r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(spider.rda$CA$u[,1:2]^2))^(1/4)
theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
circle <- data.frame(PC1 = r * cos(theta), PC2 = r * sin(theta))
g <- g + geom_path(data = circle, aes(y=PC2,x=PC1), color = muted('white'), size = 1/2, alpha = 1/3)
g
g3 <- g



## ----comparison, results='markdown', eval=TRUE, hidden=TRUE, fig.width=21, fig.height=6, hidden=TRUE----
g1 + g2 + g3


## ----readData1, results='markdown', eval=TRUE---------------------------------
spider.env <- read_csv(file = "../data/spider.env.csv", trim_ws = TRUE)
spider.env %>% glimpse()


## ----EDA4, results='markdown', eval=TRUE--------------------------------------
spider.env %>%
    as.data.frame() %>%
    ggpairs(lower = list(continuous = "smooth"),
            diag = list(continuous = "density"),
            axisLabels = "show")


## ----envfit, results='markdown', eval=TRUE------------------------------------
spider.envfit <- envfit(spider.rda, env = spider.env)
spider.envfit


## ----envfitPlot, results='markdown', eval=TRUE--------------------------------
spider.env.scores <- spider.envfit %>% fortify()
g <- g + 
    geom_segment(data=spider.env.scores,
                 aes(y=0, x=0, yend=PC2, xend=PC1),
                 arrow=arrow(length=unit(0.3,'lines')), color='blue') +
    geom_text(data=spider.env.scores,
              aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='blue')
g


## ----lmExtract, results='markdown', eval=TRUE---------------------------------
pc1 <- spider.rda.scores %>% filter(Score=='sites') %>% pull(PC1)
pc2 <- spider.rda.scores %>% filter(Score=='sites') %>% pull(PC2)


## ----lmVIF, results='markdown', eval=TRUE-------------------------------------
lm(1:nrow(spider.env) ~ soil.dry + bare.sand + fallen.leaves +
       moss + herb.layer + reflection, data=spider.env) %>%
    vif()
lm(1:nrow(spider.env) ~ herb.layer + fallen.leaves +
       bare.sand + moss, data=spider.env) %>%
    vif()


## ----lmfit, results='markdown', eval=TRUE-------------------------------------
lm(pc1 ~ herb.layer + fallen.leaves + bare.sand + moss, data=spider.env) %>%
    summary()
lm(pc2 ~ herb.layer + fallen.leaves + bare.sand + moss, data=spider.env) %>%
    summary()



## ----RDA1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.std <- spider.abund %>%
    mutate(across(everything(), function(x) x^0.25)) %>%
    wisconsin()
spider.std


## ----RDA, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.rda <- rda(spider.std ~ 
                      scale(herb.layer)+
                      scale(fallen.leaves) +
                      scale(bare.sand) +
                      scale(moss),
                  data=spider.env, scale=FALSE)


## ----RDA1a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
vif.cca(spider.rda)
summary(spider.rda, display=NULL)


## ----RDA2, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
goodness(spider.rda, display = "species")
goodness(spider.rda, display = "sites")


## ----RDA3, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
inertcomp(spider.rda)
## or proportionally
inertcomp(spider.rda, proportional = TRUE)


## ----RDA4, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
##overall test
anova(spider.rda)
#By axes
anova(spider.rda, by='axis')
#By covariate
anova(spider.rda, by='margin')


## ----RDA4a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
autoplot(spider.rda, geom='text') +
    theme_bw()


## ----RDA5, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
RsquareAdj(spider.rda)

