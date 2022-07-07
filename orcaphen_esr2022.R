######################################################################################
######## This R script contains code associated with the following reference: ########
## Ettinger, et al. 2022. Shifting phenology of an endangered apex predator mirrors ##
#######  changes in its favored prey. Endangered Species Reports. X(X):XX-XX. ########
########################## contact: ailene.ettinger@tnc.org ##########################
######################################################################################

#houskeeping
rm(list=ls()) 
options(stringsAsFactors = FALSE)

# Load libraries
library(dplyr)
library(mgcv)
library(scales)
library(RColorBrewer)
library(scales)
library(matrixStats)
library(plotfunctions)
library(igraph)
# 1. Read in the data
orcasum.days<-read.csv("data/skrwdays.csv")
# 2, Summarize whale days across the region
wdays<-as.data.frame(tapply(orcasum.days$AllSRpres,list(orcasum.days$year,orcasum.days$region),sum))
wdays.J<-as.data.frame(tapply(orcasum.days$Jpres,list(orcasum.days$year,orcasum.days$region),sum))
wdays.K<-as.data.frame(tapply(orcasum.days$Kpres,list(orcasum.days$year,orcasum.days$region),sum))
wdays.L<-as.data.frame(tapply(orcasum.days$Lpres,list(orcasum.days$year,orcasum.days$region),sum))

source("analyses/get_whaledays_lime.R")

# 3. Fit models to limekiln srkw data (this code takes a while)

source("analyses/fitsrkwmods_lime.R")

# 4. Fit models to Albion test fishery data

# 5. Fit regional models in jags