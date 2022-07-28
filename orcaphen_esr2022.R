######################################################################################
######## This R script contains code for lime kiln and albion test fishery datasets associated with the following reference: ########
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
library(brms)
library(rstan)

# 1. Read in the limekiln subset of data from Orca Master Dataset obtained from The Whale Museum https://whalemuseum.org/
limewdaysabs<-read.csv("data/limedat.csv")

# Albion Test Fishery data obtained and compiled from http://www.pac.dfo-mpo.gc.ca/fm-gp/species-especes/salmon-saumon/research-recherche/testfishery-pechedessai-eng.html

albiondat<-read.csv("data/albionddat.csv")

# CTC Escapement data (for comparison)

ctc<-read.csv("data/CTCEscapement.csv", header=TRUE)

# 3. Fit models to limekiln srkw data (this code takes a while)

source("analyses/fitsrkwmods_lime.R")

# 4. Fit models to Albion test fishery data

source("analyses/fitchinmods_albion.R")
