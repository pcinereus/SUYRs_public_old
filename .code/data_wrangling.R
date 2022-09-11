## ----setup, include=FALSE, warnings=FALSE, message=FALSE----------------------
knitr::opts_chunk$set(echo = TRUE)


## ----libraries, results='markdown', eval=TRUE, message=FALSE, warning=FALSE----
library(tidyverse) #for data wrangling


## ----getData, results='markdown', eval=TRUE-----------------------------------
load(file='../data/manipulationDatasets.RData')


## ----getData1, results='markdown', eval=TRUE----------------------------------
dat.1 %>% head


## ----exploringdata, results='markdown', eval=TRUE-----------------------------
## replace this with code to explore imported data
## From now on I will not provide many code chunks.
## Instead you are incouraged to create your own chunks
## with discussed code
## In Rstudio, you can create a chunk with Cntr-Alt-I


## ----piping, results='markup'-------------------------------------------------
head(dat.1)
#OR
dat.1 %>% head()
#OR
dat.1 %>% head
#OR
dat.1 |> head()

