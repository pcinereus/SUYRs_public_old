## ----setup, include=FALSE-----------------------------------------------------
knitr::opts_chunk$set(echo = TRUE, message=FALSE, warning=FALSE,cache.lazy = FALSE, tidy='styler')


## ----libraries, results='markdown', eval=TRUE---------------------------------
library(tidyverse)  #for data wrangling etc
library(vegan)
library(ggvegan)
library(ggrepel)
library(GGally)
library(corrplot)
library(mvabund)
library(gllvm)
library(EcolUtils)
library(car)
library(glmmTMB)
library(scales)


## ----readData, results='markdown', eval=TRUE----------------------------------
spider.abund <- read_csv(file = "../data/spider.abund.csv", trim_ws = TRUE)
spider.abund %>% glimpse()


## ----CA1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.std <- spider.abund %>%
    decostand(method="total",MARGIN=2)
spider.std


## ----CA1a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.ca <- cca(spider.std, scale=FALSE)
summary(spider.ca, display=NULL)


## ----screePlot1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
screeplot(spider.ca)
sum(eigenvals(spider.ca))/length(eigenvals(spider.ca))
eigenvals(spider.ca)/sum(eigenvals(spider.ca))


## ----biplot3, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
## ggvegan provides ways to use ggplot
## library(devtools)
## devtools::install_github("gavinsimpson/ggvegan")
## library(ggvegan)
autoplot(spider.ca)
autoplot(spider.ca) + theme_bw()
autoplot(spider.ca,geom='text') + theme_bw()


## ----biplot4, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.ca.scores <- spider.ca %>%
    fortify()
spider.ca.scores

g <-
    ggplot(data = NULL, aes(y=CA2, x=CA1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=spider.ca.scores %>% filter(Score=='sites')) +
    geom_text(data=spider.ca.scores %>% filter(Score=='sites'),
              aes(label=Label), hjust=-0.2) +
    geom_segment(data=spider.ca.scores %>% filter(Score=='species'),
                 aes(y=0, x=0, yend=CA2, xend=CA1),
                 arrow=arrow(length=unit(0.3,'lines')), color='red') +
    ## geom_text(data=spider.ca.scores %>% filter(Score=='species'),
    ##           aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='red') +
    geom_text_repel(data=spider.ca.scores %>% filter(Score=='species'),
                    aes(y=CA2*1.1, x=CA1*1.1, label=Label), color='red') +
    theme_bw()
g

                                        # Nice axes titles
eig <- eigenvals(spider.ca)

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
r <- sqrt(qchisq(circle.prob, df = 2)) * prod(colMeans(spider.ca$CA$u[,1:2]^2))^(1/4)
theta <- c(seq(-pi, pi, length = 50), seq(pi, -pi, length = 50))
circle <- data.frame(CA1 = r * cos(theta), CA2 = r * sin(theta))
g <- g + geom_path(data = circle, aes(y=CA2,x=CA1), color = muted('white'), size = 1/2, alpha = 1/3)
g



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
spider.envfit <- envfit(spider.ca, env = spider.env)
spider.envfit


## ----envfitPlot, results='markdown', eval=TRUE--------------------------------
spider.env.scores <- spider.envfit %>% fortify()
g <- g + 
    geom_segment(data=spider.env.scores,
                 aes(y=0, x=0, yend=CA2, xend=CA1),
                 arrow=arrow(length=unit(0.3,'lines')), color='blue') +
    geom_text(data=spider.env.scores,
              aes(y=CA2*1.1, x=CA1*1.1, label=Label), color='blue')
g


## ----lmExtract, results='markdown', eval=TRUE---------------------------------
pc1 <- spider.ca.scores %>% filter(Score=='sites') %>% pull(CA1)
pc2 <- spider.ca.scores %>% filter(Score=='sites') %>% pull(CA2)


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



## ----CCA1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.std <- spider.abund %>%
    decostand(method="total",MARGIN=2)
spider.std


## ----CCA, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
spider.cca <- cca(spider.std ~ 
                      scale(herb.layer)+
                      scale(fallen.leaves) +
                      scale(bare.sand) +
                      scale(moss),
                  data=spider.env, scale=FALSE)


## ----CCA1a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
vif.cca(spider.cca)
summary(spider.cca, display=NULL)


## ----CCA2, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
goodness(spider.cca, display = "species")
goodness(spider.cca, display = "sites")


## ----CCA3, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
inertcomp(spider.cca)
## or proportionally
inertcomp(spider.cca, proportional = TRUE)


## ----CCA4, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
##overall test
anova(spider.cca)
#By axes
anova(spider.cca, by='axis')
#By covariate
anova(spider.cca, by='margin')


## ----CCA4a, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
autoplot(spider.cca, geom='text') +
    theme_bw()


## ----CCA5, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
RsquareAdj(spider.cca)

