#########################################################################################
############### Orca phenology occupancy model from Ettinger, et al. 2022. ##############
# Shifting phenology of an endangered apex predator mirrors changes in its favored prey #
#######################  Endangered Species Reports. X(X):XX-XX. ########################
###################### contact: ailene.ettinger@tnc.org #################################
#########################################################################################
# based on the model presented in Strebel et al., 2014 ##################################
# Study of phenology by flexible estimation and modeling of seasonal detectability peaks# 
#########################################################################################

#housekeeping

rm(list=ls()) 
options(stringsAsFactorwants = FALSE)

# Load libraries
library(R2jags)
library(scales)

# Choose the data you want:
pod="J"#options= J,K,L
region="ps"#options=upper salish sea (uss) or puget sound (ps)

#Choose the credible intervals you want
lci<-0.125
uci<-0.875

prob<-0.5
# Read observation data from focal pod (created in orca_dataprep_occmodel.R)

  if(pod=="J"){dat<-read.csv("analyses/output/j_dat_assumeSRKW.csv",header=T)}
  if(pod=="K"){dat<-read.csv("analyses/output/k_dat_assumeSRKW.csv",header=T)}
  if(pod=="L"){dat<-read.csv("analyses/output/l_dat_assumeSRKW.csv",header=T)}

#choose region
dat<-dat[which(dat$region==region),]

#Add a column for "season" and restrict data to season that is appropriate to the region
#use may 1 for uss season, oct 1 for ps season as start dates
#use oct 31 for uss season, jan31 for ps season, as end dates

if(region == "ps"){
  season="1"#winter
  dat$season<-NA
  dat$season[dat$day>182]<-1#winter (July 1-Dec 31 = >182#should extend this to Jan 31
  }

if(region == "uss"){
  season="2"#summer
  dat$season<-NA
  dat$season[dat$day>91 & dat$day<304]<-2#summer (April 1-Oct 31)
}

dat<-dat[which(dat$season==season),]

dim(dat)


#-----------------------------------------------------------------
# Codes prepare data for jags and run the analysis
#-----------------------------------------------------------------

# Specify model in BUGS language
sink("analyses/splinesSiteOccS4psi.txt")
cat("
    model { 
    ### Define seasonal and annual patterns in occurrence probability
    for (m in 1:nyear) {  
    for (i in 1:n) {
    logit(psi[m,i]) <- lp[m,i]
    lp[m,i] <- mfe[m,i]+mre[m,i]
    mfe[m,i] <- a[m]*X[i,1]+b[m]*X[i,2]+c[m]*X[i,3]
    mre[m,i]<-sum(n.mre[m,i,1:nknots])
    for (k in 1:nknots) {
    n.mre[m,i,k]<-b.k[m,k]*Z[i,k]
    }
    }
    
    ### Random regression coefficients corresponding to the truncated polynomial functions
    for (k in 1:nknots) {
    b.k[m,k] ~ dnorm(0,taub)
    }
    
    ### Fixed regression coefficients corresponding to the 'plus' functions
    
    a[m] ~ dnorm(0,0.01)
    b[m] ~ dnorm(0,0.01)
    c[m] ~ dnorm(0,0.01)
    }
    
    ### precision for random regression coefficients corresponding to the truncated polynomial functions
    taub~dgamma(1.0E-6,1.0E-6)      
    
    # Specify priors for detection model
    for (i in 1:nsite){# 
    for (y in 1:nyear) {
    p[i,y] ~ dunif(0, 1)
    }
    }
    # Ecological submodel: Define state conditional on parameters
    
    for (y in 1:nyear) {  
    for (i in 1:n) {
    z[y,i] ~ dbern(psi[y,i])
    }    
    }
    
    # Observation model
    for (i in 1:nobs){
    muy[site[i],survey[i],year[i]] <- z[year[i],survey[i]]*p[site[i],year[i]]
    y[i] ~ dbin(muy[site[i],survey[i],year[i]], nrep[i])
    }
    
    }
    ",fill = TRUE)
    sink()

### The following procedure is based on the models presented in Crainiceanu et al. 2005 and in Gimenez et al. 2006 
# Degree of splines
degree <- 2

# covariate
covariate<-as.numeric(scale(range(dat$day)[1]:range(dat$day)[2]))

# covariate length
n <- length(covariate)

# location of knots
nk<-round((max(dat$day)-min(dat$day)+1)/4)
nknots<-ifelse(nk<35,nk,35)
knots<-quantile(unique(covariate),seq(0,1,length=(nknots+2))[-c(1,(nknots+2))])

#Note: the maximum number of knots is 35. thus, the annual model (for which nk=92 in many cases, but it is restricted to 35 by default) differs in flexibility than the seasonal model
#perhaps better to extract the seasonal peaks after fitting the whole year of data
# fixed effects matrix
 X<-NULL
for (l in 0:degree) {
  X<-cbind(X,covariate^l)  
}

# random coefficients matrix 
Z_K<-(abs(outer(covariate,knots,"-")))^3
OMEGA_all<-(abs(outer(knots,knots,"-")))^3
svd.OMEGA_all<-svd(OMEGA_all)
sqrt.OMEGA_all<-t(svd.OMEGA_all$v %*% (t(svd.OMEGA_all$u)*sqrt(svd.OMEGA_all$d)))
Z<-t(solve(sqrt.OMEGA_all,t(Z_K)))

# Input data
dat$site <- factor(dat$site)#
dat$site <- droplevels(dat$site)
fas<-sort(unique(dat$site))
dat$site <- as.integer(dat$site)

site <- dat$site
survey <- dat$day-min(dat$day)+1
nsurveys<-length(survey)
nobs <- length(unique(paste(dat$site,dat$day,dat$year)))
nrep <- dat$nrep
nsite <- length(unique(dat$site))
nyear <- length(unique(dat$year))
year <- as.numeric(factor(dat$year))
zst <- array(1, dim=c(nyear,n))  
y <- dat$ndet

# Simulation parameters

ni=20000; nc=2; nb=10000; nt=10

# List input data
jags.data <- list("site","survey","nobs","nrep","nsite","nyear","year","nknots","n","X","Z","nc", "nb", "ni", "nt","zst","y")

# Inits function
f.inits <- function(){list(a=rep(0,nyear), b=rep(0,nyear), c=rep(0,nyear), z=zst)}

# specify the parameters to be monitored
parameters <- c("a","b","c","lp","psi","taub","p")

### Run MCMC Analysis using jags

jags.out<-jags.parallel(jags.data,f.inits,parameters,"analyses/splinesSiteOccS4psi.txt",nc,ni,nb,nt)

#names(jags.out$BUGSoutput)
#diagnose the model
#plot(jags.out)
#traceplot(jags.out, dig=3)#


#Look at psi
out<-jags.out$BUGSoutput
jags.out$BUGSoutput$mean$psi#probability of presence (daily, across 40 years)

#Look at Rhat
range(out$summary[,8])
length(which(out$summary[,8]>1.1))/length(out$summary[,8])
rownames(out$summary[which(out$summary[,8]>1.1),])
  #look at problem psis & problem lps

dim(jags.out$BUGSoutput$mean$psi)
meanpsi<-rowMeans(jags.out$BUGSoutput$mean$psi)#meanannual mean prob occurrence
names(meanpsi)<-seq(min(dat$year),max(dat$year), by=1)

meanp<-jags.out$BUGSoutput$mean$p
rownames(meanp)<-fas
colnames(meanp)<-seq(min(dat$year),max(dat$year), by=1)

#-----------------------------------------------------------------
# Codes to summarize the output
#-----------------------------------------------------------------

### get estimated date of peak occurrence based on posterior distribution
# get date of peak probability of occurrence in each simulation

if(region == "uss"){color = "darkblue"
cols = c("darkblue","darkblue")}
if(region == "ps"){color = "goldenrod"
cols = c("goldenrod","goldenrod")}

prob<-c(0.2,0.3,0.4,0.5)
for(p in prob){
findmax.fn<-function(x) {
  mean(which(x==max(x)))
}
lpmax<-array(data=NA,dim=c(out$n.sims,nyear))
dimnames(lpmax)<-list(c(1:out$n.sims),c(sort(unique(dat$year))))
for (xj in sort(unique(as.numeric(factor(dat$year))))) { 
  lpmax[,xj]<-apply(out$sims.array[,,paste("lp[",xj[1],",",1:(max(dat$day)-min(dat$day)+1),"]",sep="")],MARGIN=c(if(out$n.chains>1) 1:2 else 1),findmax.fn)
}
lpmax<-lpmax+min(dat$day)-1
lpmax[lpmax==max(dat$day)]<-NA
lpmax[lpmax==min(dat$day)]<-NA
#would like to Extract and plot psi (probability of presence by day...)
dim(out$sims.list$psi)
doy<-seq(from=min(dat$day),to=max(dat$day),by=1)
}

# summarize estimates and look at change across the whole time series

ann.res<-array(NA, dim=c(max(dat$year)-min(dat$year)+1,3),dimnames=list(c(min(dat$year):max(dat$year)),c("mean","lci","uci")))
res<-apply(lpmax,c(2),mean,na.rm=T)
ann.res[names(res),"mean"]<-res
res<-apply(lpmax,c(2),quantile,probs=lci,na.rm=T)
ann.res[names(res),"lci"]<-res
res<-apply(lpmax,c(2),quantile,probs=uci,na.rm=T)
ann.res[names(res),"uci"]<-res


psi.ann<-array(NA, dim=c(max(dat$year)-min(dat$year)+1,3),dimnames=list(c(min(dat$year):max(dat$year)),c("mean","lci","uci")))
psi<-apply(out$sims.list$psi,c(2),mean,na.rm=T)
psi.ann[names(res),"mean"]<-psi
psi<-apply(out$sims.list$psi,c(2),quantile,probs=lci,na.rm=T)
psi.ann[names(res),"lci"]<-psi
psi<-apply(out$sims.list$psi,c(2),quantile,probs=uci,na.rm=T)
psi.ann[names(res),"uci"]<-psi


#plot mean psi across all years for the season
par(mfrow=c(1,2),mai=c(1,1,1,0.5))

#windows(height=6,width=8)
psi.med<-apply(out$sims.list$psi[,32:40,],c(3),median)
plot(doy,psi.med, type= "l", ylim=c(0,1.2), 
     ylab= "Probability of occurrence", xlab= "Day of Year", yaxt ="n",
     bty="l", lty=1,col=color, lwd=2)
axis(side = 2, at = seq(from = 0, to = 1.0, by = .2), labels = c("0","0.2","0.4","0.6","0.8", "1.0"))

psi.uci<-apply(out$sims.list$psi[,32:40,],c(3),quantile,probs=uci)
psi.lci<-apply(out$sims.list$psi[,32:40,],c(3),quantile,probs=lci)

polygon(c(rev(doy),doy),c(rev(psi.uci),psi.lci),col=alpha(color,0.05),lty=0)

#add point for mean peak day of year
pkocdoy<-mean(ann.res[,"mean"][32:40])
pkocdoy.lci<-quantile(ann.res[,"mean"][32:40],lci)
pkocdoy.uci<-quantile(ann.res[,"mean"][32:40],uci)

arrows(pkocdoy.lci,1.14,pkocdoy.uci,1.2,code = 0, col = cols[1], lwd = 3)
points(pkocdoy,1.14, pch = 21, bg=cols[1], cex = 2)

psi.med<-apply(out$sims.list$psi[,24:31,],c(3),median)

lines(doy,psi.med,col=cols[2], lwd=2, lty=2)
psi.uci<-apply(out$sims.list$psi[,24:31,],c(3),quantile,probs=uci)

psi.lci<-apply(out$sims.list$psi[,24:31,],c(3),quantile,probs=lci)

polygon(c(rev(doy),doy),c(rev(psi.uci),psi.lci),col=alpha(cols[2],0.05),lty=0)

#add point for mean peak day of year
pkocdoy<-mean(ann.res[,"mean"][24:31])
pkocdoy.lci<-quantile(ann.res[,"mean"][24:31],lci)
pkocdoy.uci<-quantile(ann.res[,"mean"][24:31],uci)

arrows(pkocdoy.lci,1.2,pkocdoy.uci,1.14,code = 0, col = cols[2], lty=2,lwd = 3)
points(pkocdoy,1.2, pch = 21, bg=cols[2], cex = 2)

if(season==1){mtext("C)",side = 3, line = 0, adj =0)}
if(season==2){  mtext("A)",side = 3, line = 0, adj =0)}

legend("topleft",legend=c(paste(unique(dat$year)[24],"-",unique(dat$year)[31],sep= ""),
                          paste(unique(dat$year)[32],"-",unique(dat$year)[40],sep= "")),lwd=c(2,2),lty=c(2,1),col=c(cols), bty="n")


# get estimate of trend in date of peak detectability over years
do.lm<-function(x) {
  lmres<-lm(x~as.numeric(names(x)))$coefficients
  return(lmres)
}
r<-matrix(NA,dim(lpmax)[1],2)
#just do lm from 1978-2017
for (o in 1:(dim(lpmax)[1])) {
  # if(!is.na(sum(lpmax[o,]))) {
  lm(lpmax[o,]~as.numeric(colnames(lpmax)))$coefficients->r[o,]
  #}    
}
slopevec<-as.vector(r[,2])
intercept<-mean(r[,1],na.rm=T)
slope<-mean(r[,2],na.rm=T)
intercept.lci<-quantile(r[,1],c(lci),na.rm=T)
intercept.uci<-quantile(r[,1],c(uci),na.rm=T)

slope.lci<-quantile(r[,2],c(lci),na.rm=T)
slope.uci<-quantile(r[,2],c(uci),na.rm=T)
