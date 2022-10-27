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
macnally <- read.csv('../data/macnally_full.csv',strip.white=TRUE)
head(macnally)
macnally[1:5,1:5]


## ----PCoA1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
macnally.std <-
    macnally %>% dplyr::select(-HABITAT) %>%
    decostand(method="total",MARGIN=2)
macnally.std

macnally.dist = vegdist(macnally.std, method='bray')
macnally.capscale = capscale(macnally.dist~1, data=macnally$HABITAT)
macnally.capscale = capscale(macnally[,-1]~1, data=macnally$HABITAT)

summary(macnally.capscale, display=NULL)
plot(macnally.capscale)
autoplot(macnally.capscale, geom='text') + theme_bw()



## ----MDS1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
macnally.mds <- metaMDS(macnally[,-1], k=2,  plot=TRUE)
macnally.mds

#OR

macnally.std <- wisconsin(macnally[,c(-1)]^0.25)
macnally.dist <- vegdist(macnally.std,"bray")
macnally.mds1 <- metaMDS(macnally.std, k=2, plot=TRUE)
macnally.mds2 <- metaMDS(macnally.dist, k=2, plot=TRUE)



## ----Stressplot, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
stressplot(macnally.mds)


## ----biplot1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
autoplot(macnally.mds) + theme_bw()


## ----biplot2, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
macnally.mds.scores <- macnally.mds %>%
    fortify() %>%
    full_join(macnally %>%
              rownames_to_column(var='Label'))

g <-
    ggplot(data = NULL, aes(y=NMDS2, x=NMDS1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=macnally.mds.scores %>% filter(Score=='sites'),
               aes(color=HABITAT)) +
    geom_text(data=macnally.mds.scores %>% filter(Score=='sites'),
              aes(label=Label, color=HABITAT), hjust=-0.2) 

g


g + ggforce::geom_mark_ellipse(data=macnally.mds.scores %>% filter(Score=='sites'),
                      aes(y=NMDS2, x=NMDS1, fill=HABITAT), expand=0) 
g <- g + ggforce::geom_mark_hull(data=macnally.mds.scores %>% filter(Score=='sites'),
                      aes(y=NMDS2, x=NMDS1, fill=HABITAT), concavity = 10) 

g + theme_bw()


## ----envfit, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
Xmat <- model.matrix(~-1+HABITAT, data=macnally)
colnames(Xmat) <-gsub("HABITAT","",colnames(Xmat))
envfit <- envfit(macnally.mds, env=Xmat)
envfit

data.env.scores <- envfit %>% fortify()
g <- g + 
    geom_segment(data=data.env.scores,
                 aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
                 arrow=arrow(length=unit(0.3,'lines')), color='blue') +
    geom_text(data=data.env.scores,
              aes(y=NMDS2*1.1, x=NMDS1*1.1, label=Label), color='blue')
g + theme_bw()



## ----simper, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
simper(macnally.std, macnally$HABITAT)


## ----adonis, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
macnally <- macnally %>% mutate(HABITAT = factor(HABITAT))
macnally.dist <- vegdist(macnally[,-1], 'bray')
macnally.adonis <- adonis2(macnally.dist ~ HABITAT, data=macnally)
macnally.adonis

EcolUtils::adonis.pair(macnally.dist, macnally$HABITAT)



## ----betadisper, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
macnally.disp <- betadisper(macnally.dist, macnally$HABITAT)
boxplot(macnally.disp)
plot(macnally.disp)
anova(macnally.disp)
permutest(macnally.disp, pairwise = TRUE)
TukeyHSD(macnally.disp)


## ----betadisper1, results='markdown', eval=TRUE, hidden=TRUE, fig.width=6, fig.height=6, hidden=TRUE----
macnally.disp <- betadisper(macnally.dist, macnally$HABITAT, type="median",bias.adjust = TRUE)
boxplot(macnally.disp)
plot(macnally.disp)
anova(macnally.disp)
permutest(macnally.disp, pairwise = TRUE)
TukeyHSD(macnally.disp)

