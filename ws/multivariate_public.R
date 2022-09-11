#############################################################
##             RDA (Correlations, Euclidean distance)                
##            /  \                              
##Unconstrained   Constrained ---> ordination                    
##   (PCA)            (RDA)   ---> anova
##
##              CA (Chisq distance)
##             /  \
##Unconstrained   Constrained ---> ordination
## (CA)             (CCA)     ---> anova
##
##             PCoA (any distance)
##             /  \
##Unconstrained   Constrained ---> ordination
##                            ---> anova
##
##             dbRDA (any distance)
##             /  \
##Unconstrained   Constrained ---> ordination
##                            ---> anova
##
##Unconstrained  ---> ordination
##               ---> envfit (overlay enviromental data) (permutation test)
##               ---> lm/glm etc (response or predictor)
#############################################################
##     Dissimilarity
##            --> MDS      ---> ordination
##            --> bioenv   ---> variable importance       (perm test)
##            --> adonis2* ---> anova                     (perm test)
##            --> simper   ---> similarity percentages
##            --- betadisp ---> homogeneity of dispersion (perm test)
#############################################################
##     Model based ordination
##            ---> glmmTMB (via reduced rank / latent variable)
##            ---> gllvm (generalised latent variable models)
#############################################################
##     Model based
##            ---> manyglm ---> anova
##            ---> gllvm (generalized latent variable models)
#############################################################


library(tidyverse)
library(vegan)
library(GGally)
library(corrplot)
library(car)
library(mvabund)
library(scales)
library(ggvegan)
library(ggrepel)
library(glmmTMB)
library(gllvm)
library(EcoUtils)

## ---- PCA - spiders
data(spider)
glimpse(spider)

spider.abund <- spider$abund %>%
    as.data.frame()
write_csv(spider.abund,
          file = "../public/data/spider.abund.csv")
          
spider.abund <- read_csv(file = "../public/data/spider.abund.csv", trim_ws = TRUE)

glimpse(spider.abund)

{
    ## ---- EDA 1
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
    ## ----end
    ## ---- EDA 2
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
    ## ----end

    ## Conclusions:
    
    ## There is evidence of non-normality and non-linearity although
    ##we could attempt to normalize via sqrt transform, this is
    ##unlikely to fix all.vars linearity also not good. It is likely
    ##that these issues are the result of the sampling occuring over a
    ##larger scale than the natural range of the taxa.  This will
    ##result in some species being left skewed and others being right
    ##skewed. It will also result in non-linearity and the horseshoe
    ##effect.

    ## Technically, these techniques also assume homogeneity of
    ## variance.  Hwever, it is not possible to explore residuals from
    ## these techniques.
    
    ## Evidence of linearity and non-normality
    ## Solutions:
    ## Metric ordination
    ## 1. PCA - ignore at own peril
    ## 2. hellinger trasformation (tb-PCA)
    ## 3. CA (CCA)
    ## 4. Distance-based RDA (capscale)
    ## Non-metric
    ## 5. nMDS

    ## ---- PCA
    spider.std <- spider.abund %>%
        mutate(across(everything(), function(x) x^0.25)) %>%
        wisconsin()
    data.std
    spider.hel <- spider.abund %>%
        decostand(method = "hellinger")

    ## Run PCA - unconstrained axis rotations
    ## normally we would not scale, but for illustrative purposes...
    spider.rda <- rda(spider.abund, scale=TRUE)
    spider.rda <- rda(spider.std, scale=TRUE)
    spider.rda <- rda(spider.std, scale=FALSE)
    spider.rda <- rda(spider.hel, scale=FALSE)
    ## spider.rda <- rda(spider.std, scale=FALSE)
    summary(spider.rda, display=NULL)
    
    screeplot(spider.rda)
    abline(a=1,b=0)
    
    scores(spider.rda, choices=1:3, display='sites')
    scores(spider.rda, choices=1:3, display='species')
    
    ## Quick and nasty ordination plots
    biplot(spider.rda, scaling='species')
    biplot(spider.rda, scaling='sites')
    
    ## Quick and nasty ordination plots
    pl<-ordiplot(spider.rda)
    points(pl, "sites", pch=21, col="red", bg="yellow")
    text(pl, "sites", col="red", cex=0.9)
    text(pl, "species", col="blue", cex=0.9)
    ## line(pl, "sites")

    
    ## ggvegan provides ways to use ggplot
    ## library(devtools)
    ## devtools::install_github("gavinsimpson/ggvegan")
    ## library(ggvegan)
    autoplot(spider.rda)
    autoplot(spider.rda) + theme_bw()
    autoplot(spider.rda,geom='text') + theme_bw()

    ## Alternatively, we can extract the scores from the PCA/RDA etc
    ## into a dataframe.
    ## - fortify() is a general method used to extract/convert data into
    ##   a form that is suitable for ggplot.
    ## - the ggvegan package defines some fortify methods specifically
    ##   for analyses created with the vegan package

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
    
    ## ----end
    ## ---- Envfit
    X <- spider$x %>% as.data.frame()
    X %>%
        cor %>%
        corrplot(type = 'upper',
                 order = 'FPC',
                 diag = FALSE)
    X %>%
        as.data.frame() %>%
        ggpairs(lower = list(continuous = "smooth"),
                diag = list(continuous = "density"),
                axisLabels = "show")
    spider.env <- envfit(spider.rda, env = X)
    spider.env

    
    spider.env.scores <- spider.env %>% fortify()
    g <- g + 
        geom_segment(data=spider.env.scores,
                     aes(y=0, x=0, yend=PC2, xend=PC1),
                     arrow=arrow(length=unit(0.3,'lines')), color='blue') +
        geom_text(data=spider.env.scores,
                  aes(y=PC2*1.1, x=PC1*1.1, label=Label), color='blue')
    g

    ## ----end
    ## ---- lm
    pc1 <- spider.rda.scores %>% filter(Score=='sites') %>% pull(PC1)
    pc2 <- spider.rda.scores %>% filter(Score=='sites') %>% pull(PC2)
    
    lm(1:nrow(X) ~ soil.dry + bare.sand + fallen.leaves +
           moss + herb.layer + reflection, data=X) %>%
        vif()
    lm(1:nrow(X) ~ herb.layer + fallen.leaves + bare.sand + moss, data=X) %>%
        vif()
    lm(pc1 ~ herb.layer + fallen.leaves + bare.sand + moss, data=X) %>%
        summary()
    lm(pc2 ~ herb.layer + fallen.leaves + bare.sand + moss, data=X) %>%
        summary()

    ## ----end

}
## ----end

## ---- RDA spiders
{
    ## ---- RDA


    spider.rda <- rda(spider.std ~ 
                          scale(herb.layer)+
                          scale(fallen.leaves) +
                          scale(bare.sand) +
                          scale(moss),
                      data=X, scale=FALSE)
    spider.rda <- rda(spider.hel ~ 
                          scale(herb.layer)+
                          scale(fallen.leaves) +
                          scale(bare.sand) +
                          scale(moss),
                      data=X, scale=FALSE)
    summary(spider.rda, display=NULL)

    vif.cca(spider.rda)

    ## ----end
    ## ---- goodness of fit
    ## Goodness of fit test
    ## proportion of inertia accounted for by species
    ## up to chosen axes
    ## proportion can be assessed by either species or sites
    ## (depending on argument display)
    ## THis is the cumulative proportion of variance explained by each axis
    ## Often used to remove species from ordination (those not well fit)
    goodness(spider.rda)
    goodness(spider.rda, display = "sites")
    ## Proportion of inertia explained (decomposed) by constrained and unconstrained
    inertcomp(spider.rda)
    ## or proportionally
    inertcomp(spider.rda, proportional = TRUE)

    ## ----end
    ## ---- Anova
    ##overall test
    anova(spider.rda)
    anova(spider.rda, by='axis')
    anova(spider.rda, by='margin')
                                        #anova(spider.rda, by='margin', scope="Altitude")

    ## see the regression coefficients
    coef(spider.rda)

    RsquareAdj(spider.rda)

    screeplot(spider.rda)

    autoplot(spider.rda, geom='text')
    ## ----end
}

## ----end

## ---- CA spiders

spider.std <- spider$abund %>%
    decostand(method="total",MARGIN=2)
spider.std

spider.std %>% cor %>%
    corrplot(diag=FALSE, order='FPC')

spider.ca <- cca(spider.std, scale=FALSE)
summary(spider.ca, display=NULL)
#anova(spider.ca)

screeplot(spider.ca)
sum(eigenvals(spider.ca))/length(eigenvals(spider.ca))
eigenvals(spider.ca)/sum(eigenvals(spider.ca))
g <- autoplot(spider.ca, geom='text')
spider.env <- envfit(spider.ca, env=X)
spider.env
g + autoplot(spider.env)



## spider.dca <- decorana(spider.abund)
## summary(spider.dca, display="none")
## ordiplot(spider.dca)
## autoplot(spider.dca,  geom="text")
## ----end

## ---- bioenv spider

spider.dist <- vegdist(wisconsin(spider.abund^0.25),"bray")
spider.bioenv <- bioenv(spider.dist,
                      decostand(spider$x,"standardize"))
spider.bioenv
## ----end

## ---- adonis spider
## When covariates are continuous, adonis uses
## Distance-based redundancy analysis (dbRDA)
spider.dbrda <- dbrda(spider.dist~ soil.dry + bare.sand + fallen.leaves + moss,
                        data=spider$x)
#ordiplot(spider.dbrda)
autoplot(spider.dbrda)
autoplot(spider.dbrda,geom='text') + theme_bw()
spider.dbrda
eigenvals(spider.dbrda)
spider.eig <- eigenvals(spider.dbrda)
spider.eig/sum(spider.eig)

anova(spider.dbrda, by = "axis")
anova(spider.dbrda, by = "term")
anova(spider.dbrda, by = "margin")
## The above, really just call permutest
permutest(spider.dbrda)
permutest(spider.dbrda, by = "onedf")
permutest(spider.dbrda, by = "terms")


spider.adonis <- adonis(spider.dist ~ soil.dry + bare.sand + fallen.leaves + moss,
                        data=spider$x)
spider.adonis
## The above is sequentially
## Note, if we change the order (mainly important for interactions)
spider.adonis <- adonis(spider.dist ~ fallen.leaves + soil.dry + bare.sand + moss,
                        data=spider$x)
spider.adonis
## Adonis2 can do either sequential (as next)
spider.adonis <- adonis2(spider.dist ~ soil.dry + bare.sand + fallen.leaves + moss,
                        data=spider$x)
spider.adonis
## or marginal (next)
spider.adonis <- adonis2(spider.dist ~ soil.dry + bare.sand + fallen.leaves + moss,
                        data=spider$x, by = 'margin')
spider.adonis
## or overall
spider.adonis <- adonis2(spider.dist ~ soil.dry + bare.sand + fallen.leaves + moss,
                        data=spider$x, by = NULL)
spider.adonis

## spider.disp <- betadisper(spider.dist, group = spider$x$fallen.leaves)
## spider.disp

## boxplot(spider.disp)
## plot(spider.disp)
## anova(spider.disp)
## TukeyHSD(spider.disp)
## permutest(spider.disp, pairwise = TRUE)
## ----end

## ---- MDS spiders
## MUST READ IN THIS WAY..

spider.mds <- metaMDS(spider.abund, k=2,  plot=TRUE)
spider.mds

spider.std <- wisconsin(spider.abund^0.25)

#apply(macnally.std[,c(-1,-2)],2,max)
#apply(macnally.std[,c(-1,-2)],2,var, na.rm=TRUE)

#vegdist(macnally[,-1], method='bray')

spider.dist <- vegdist(spider.std,"bray")

spider.mds <- metaMDS(spider.std, k=2, plot=TRUE)
spider.mds <- metaMDS(spider.dist, k=2, plot=TRUE)
spider.mds <- metaMDS(spider.abund, k=2)


wascores(spider.mds$points, spider.abund)


spider.mds$stress

stressplot(spider.mds)

plot(macnally.mds)

## autoplot(macnally.mds)
spider.mds.scores <- spider.mds %>%
        fortify()

spider.env <- spider$x
g <-
    ggplot(data = NULL, aes(y=NMDS2, x=NMDS1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=spider.mds.scores %>%
                   filter(Score=='sites'),
               aes(color=spider.env$soil.dry)) +
    geom_text(data=spider.mds.scores %>%
                  filter(Score=='sites'),
              aes(label=Label,
                  color=spider.env$soil.dry), hjust=-0.2) +
    scale_colour_binned(type = "viridis")
g

spider.env <- spider.env %>%
    mutate(fSoil = cut(soil.dry, breaks = c(0, 1, 2, 4)))
g <-
    ggplot(data = NULL, aes(y=NMDS2, x=NMDS1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=spider.mds.scores %>%
                   filter(Score=='sites'),
               aes(color=spider.env$fSoil)) +
    geom_text(data=spider.mds.scores %>%
                  filter(Score=='sites'),
              aes(label=Label,
                  color=spider.env$fSoil), hjust=-0.2) 
    scale_colour_viridis_c()
g


spider.env <- spider.env %>%
    mutate(fFallen = cut(fallen.leaves, breaks = c(-1, 1, 3, 5)))
g <-
    ggplot(data = NULL, aes(y=NMDS2, x=NMDS1)) +
    geom_hline(yintercept=0, linetype='dotted') +
    geom_vline(xintercept=0, linetype='dotted') +
    geom_point(data=spider.mds.scores %>%
                   filter(Score=='sites'),
               aes(color=spider.env$fFallen)) +
    geom_text(data=spider.mds.scores %>%
                  filter(Score=='sites'),
              aes(label=Label,
                  color=spider.env$fFallen), hjust=-0.2) 
    scale_colour_viridis_c()
g


g + ggforce::geom_mark_ellipse(data=spider.mds.scores %>%
                                   filter(Score=='sites'),
                               aes(y=NMDS2, x=NMDS1, fill=spider.env$fFallen),
                               expand=0) 
g + ggforce::geom_mark_hull(data=macnally.mds.scores %>% filter(Score=='sites'),
                      aes(y=NMDS2, x=NMDS1, fill=HABITAT), expand=0) 
g + ggforce::geom_mark_hull(data=spider.mds.scores %>%
                                filter(Score=='sites'),
                            aes(y=NMDS2, x=NMDS1, fill=spider.env$fFallen),
                            expand=0, concavity = 10) 

Xmat <- model.matrix(~soil.dry + bare.sand + fallen.leaves + moss +
                         herb.layer + reflection, data=spider.env)
envfit <- envfit(spider.mds, env=Xmat)
envfit


spider.env.scores <- envfit %>% fortify()
g <- g + 
    geom_segment(data=spider.env.scores,
                 aes(y=0, x=0, yend=NMDS2, xend=NMDS1),
                 arrow=arrow(length=unit(0.3,'lines')), color='blue') +
    geom_text(data=spider.env.scores,
              aes(y=NMDS2*1.1, x=NMDS1*1.1, label=Label), color='blue')
g


## ----end
## ---- MDS macnally
## MUST READ IN THIS WAY..
macnally <- read.csv('../public/data/macnally_full.csv',strip.white=TRUE)
head(macnally)
macnally[1:5,1:5]


apply(macnally[,c(-1)],2,mean, na.rm=TRUE)
apply(macnally[,c(-1)],2,max)
apply(macnally[,c(-1)],2,sum)
apply(macnally[,c(-1)],2,var, na.rm=TRUE)


macnally.mds <- metaMDS(macnally[,-1], k=2,  plot=TRUE)
macnally.mds

library(vegan)
macnally.std <- wisconsin(macnally[,c(-1)]^0.25)
macnally.dist <- vegdist(macnally.std,"bray")

macnally.mds <- metaMDS(macnally.std, k=2, plot=TRUE)
macnally.mds <- metaMDS(macnally.dist, k=2, plot=TRUE)
macnally.mds <- metaMDS(macnally[,-1], k=2)


macnally.mds$stress

stressplot(macnally.mds)

macnally.mds.scores <- macnally.mds %>%
        fortify() %>%
    full_join(macnally %>% add_rownames(var='Label'))

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
g + ggforce::geom_mark_hull(data=macnally.mds.scores %>% filter(Score=='sites'),
                      aes(y=NMDS2, x=NMDS1, fill=HABITAT), expand=0) 
g + ggforce::geom_mark_hull(data=macnally.mds.scores %>% filter(Score=='sites'),
                      aes(y=NMDS2, x=NMDS1, fill=HABITAT), expand=0, concavity = 10) 


## Xmat <- model.matrix(~HABITAT, data=macnally)
## envfit <- envfit(macnally.mds, env=Xmat)
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
g


plot(envfit, col="gray")


bioenv(macnally.dist,
                      decostand(habitat,"standardize"))

adonis(macnally.dist ~ HABITAT, data=macnally)

simper(macnally.std, macnally$HABITAT)

macnally.disp <- betadisper(macnally.dist, macnally$HABITAT)
boxplot(macnally.disp)
plot(macnally.disp)
anova(macnally.disp)
permutest(macnally.disp, pairwise = TRUE)
TukeyHSD(macnally.disp)

macnally.disp <- betadisper(macnally.dist, macnally$HABITAT, type="median",bias.adjust = TRUE)
boxplot(macnally.disp)
plot(macnally.disp)
anova(macnally.disp)
permutest(macnally.disp, pairwise = TRUE)
TukeyHSD(macnally.disp)
## ----end


## ---- dune
dune <- read_csv('../public/data/dune.csv', trim_ws=TRUE)
dune <- dune %>% mutate(MANAGEMENT=factor(MANAGEMENT,  levels=c("NM","BF","HF","SF"))) %>%
  as.data.frame()
dune

#species means
apply(dune[,-1],2, mean, na.rm=TRUE)
#species maximums
apply(dune[,-1],2, max)
#species sums
apply(dune[,-1],2, sum, na.rm=TRUE)
#species variance
apply(dune[,-1],2, var, na.rm=TRUE)

dune.dist <- vegdist(wisconsin(dune[,-1]^0.25), "bray")

dune.mds = metaMDS(dune.dist, k=2)
## autoplot(dune.mds, geom="text")
plot(dune.mds, type="text", display="sites" )


dune.adonis<-adonis2(dune.dist~MANAGEMENT,  data=dune)
dune.adonis

management <-factor(dune$MANAGEMENT, levels=c("NM","BF","HF","SF"))
mm <- model.matrix(~management)
head(mm)
colnames(mm) <-gsub("management","",colnames(mm))
mm <- data.frame(mm)
dune.adonis<-adonis2(dune.dist~BF+HF+SF, data=mm,
                    perm=9999)
dune.adonis

dune.simper=simper(dune[,-1], dune[,1], permutations = 999)
summary(dune.simper)

#PERMDISP2 - multivariate homogeneity of group dispersions (variances)
dune.disp <- betadisper(dune.dist,  group=dune$MANAGEMENT)
boxplot(dune.disp)
plot(dune.disp)
anova(dune.disp)
permutest(dune.disp, pairwise = TRUE)
TukeyHSD(dune.disp)
## ----end


## ---- varveg
vareveg <- read.csv('../public/data/vareveg.csv')
head(vareveg)
vareenv <- read.csv('../public/data/vareenv.csv')
head(vareenv)


vareveg.dist <- vegdist(wisconsin(vareveg[,-1]^0.25),'bray')
#vareveg.dist <- vegdist(vareveg.std, "bray")
#environmental variables 
vareenv.std <- decostand(vareenv[,-1], "standardize")
vareenv.dist <- vegdist(vareenv.std, "euc")

#bioenv(vareveg.std, vareenv.std)

bioenv(vareveg.dist, vareenv.std)
#adonis(vareveg.std~Ca+Fe+Mn+Baresoil, data=vareenv)
adonis(vareveg.dist~Ca+Fe+Mn+Baresoil, data=vareenv.std)
adonis2(vareveg.dist~Ca+Fe+Mn+Baresoil, data=vareenv.std)

#meandist(vareveg.dist)
#vareveg.simper = simper(vareveg[,-1], vareenv.std$Ca, permutations=100)
#vareveg.simper
#summary(vareveg.simper)

## ----end

## ---- glmmTMB spider
dat.spider.1 <- spider.abund %>%
    as.data.frame() %>%
    mutate(Site = factor(1:n())) %>%
    pivot_longer(cols = -Site,
                 names_to = 'Species',
                 values_to = 'Abund')
library(glmmTMB)
spider.glmmTMB <- glmmTMB(Abund ~ 1 + rr(Species + 0|Site, d = 2),
                          family = nbinom2(),
                          dat = dat.spider.1
                          )
spider.loadings <- spider.glmmTMB$obj$env$report(spider.glmmTMB$fit$parfull)$fact_load[[1]] %>%
                                                                           as.data.frame() %>%
                                                                           mutate(Species = colnames(spider$abund))
fit <-
    ranef(spider.glmmTMB)[[1]]$Site %>%
    mutate(Site = 1:n())

ggplot(fit, aes(y = SpeciesAlopcune, x = SpeciesAlopacce)) +
    geom_text(aes(label = Site)) +
    geom_text(data = spider.loadings, aes(y = V2, x = V1, label = Species), color = 'blue')
## ----end

## ---- gllvm spiders
library(gllvm)
fitx <- gllvm(y = spider$abund, X=spider$x, family = "poisson", num.lv = 2)
fitx <- gllvm(y = spider$abund, X=spider$x, family = "negative.binomial", num.lv = 2)
fitx
par(mfrow = c(1,2))
plot(fitx, which = 1:2)
summary(fitx)
coefplot(fitx, mfrow = c(3,2), cex.ylab = 0.8)
crx <- getResidualCor(fitx)
corrplot(crx, diag = FALSE, type = "lower", method = "square", tl.srt = 25)

ordiplot(fitx, biplot = TRUE)
abline(h = 0, v = 0, lty=2)
## ----end

## ---- gllvm microbial

data(microbialdata)
X <- microbialdata$Xenv
y <- microbialdata$Y[, order(colMeans(microbialdata$Y > 0), 
                             decreasing = TRUE)[21:40]]
fit <- gllvm(y, X, formula = ~ pH + Phosp, family = poisson())
fit$logL
ordiplot(fit)
coefplot(fit)
Site<-data.frame(Site=X$Site)
Xsoils <- cbind(scale(X[, 1:3]),Site)
ftXph <- gllvm(y, Xsoils, formula = ~pH, family = "negative.binomial", 
               row.eff = ~(1|Site), num.lv = 2)
Xenv <- data.frame(X, Region = factor(X$Region),
                   Soiltype = factor(X$Soiltype))
ftXi <- gllvm(y, Xenv, formula = ~ SOM + pH + Phosp + Region, 
              family = "negative.binomial", row.eff = ~(1|Site), num.lv = 2,
              sd.errors = FALSE)
ftXi
ftXph

ph <- Xenv$pH
rbPal <- colorRampPalette(c('mediumspringgreen', 'blue'))
Colorsph <- rbPal(20)[as.numeric(cut(ph, breaks = 20))]
pchr = NULL
pchr[Xenv$Region == "Kil"] = 1
pchr[Xenv$Region == "NyA"] = 2
pchr[Xenv$Region == "Aus"] = 3
ordiplot(ftXi, main = "Ordination of sites",  
         symbols = TRUE, pch = pchr, s.colors = Colorsph)
legend("topleft", legend = c("Kil", "NyA", "Mayr"), pch = c(1, 2, 3), bty = "n")

ftNULL <- gllvm(y, X = data.frame(Site = X[,5]), 
              family = "negative.binomial", row.eff = ~(1|Site), num.lv = 2,
              sd.errors = FALSE)
1 - getResidualCov(ftXi)$trace/getResidualCov(ftNULL)$trace
## ----end
