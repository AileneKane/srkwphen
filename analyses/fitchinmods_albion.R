
# Prep data for models 
albiondat$calDay<-as.integer(as.character(albiondat$calDay))
allyears<-unique(albiondat$year)
dat<-albiondat
season="allyear"#choices are "springsum" or "fall" or "allyear
head(dat)

dat<-dat[dat$year<2018,]
sumalb<-tapply(dat$cpue,list(dat$year),sum)
# 2. Compare albion data to CTC data
ctc<-ctc[ctc$Year<2018 & ctc$Yea>1979,]
#select only spring and summer
datsp<-dat[dat$calDay<213,]
cpuesptot<-aggregate(datsp$cpue,list(datsp$year),sum, na.rm= TRUE)
cpuetot<-aggregate(dat$cpue,list(dat$year),sum, na.rm= TRUE)

colnames(cpuetot)<-c("year","cpuetot")
colnames(cpuesptot)<-c("year","cpuesptot")

cpuetot<-cpuetot[cpuetot$year<2018,]
cpuesptot<-cpuesptot[cpuesptot$year<2018,]

ctc$tot<-as.numeric(ctc$SpringSummerTotalRun)
ctc$sp12<-ctc$SpringAge1.2ctc$total<-ctc$SpringAge1.2+ctc$SpringAge1.3+ctc$SummerAge0.3+ctc$SummerAge1.3+ctc$Harrison_Esc+ctc$LowerShuswap_Esc
pdf(file="analyses/orcaphen/figures/ctcalbion.pdf",height=15,width=8)
par(mfrow=c(4,2))
ctccols<-c(2,3,4,5)
for(i in ctccols){
  ctcnums<-ctc[ctc$Year>1990,]
  ctcnums<-ctcnums[,i]
  plot(ctcnums,cpuesptot$cpuesptot[as.numeric(cpuetot$year)>1990], pch=16, xlab = "Escapement estimates", ylab="Albion Test Fishery (CPUE, 1April-1Aug)",col="gray",bty="l", main = paste(colnames(ctc)[i]), 
       cex.lab=1.5, cex.axis=1.5,cex=2, cex.main=1.5)
mod<-lm(cpuesptot$cpuesptot[as.numeric(cpuesptot$year)>1990]~ctcnums)
if(summary(mod)$coef[2,4]<0.1){abline(mod, lwd=2)}
r<-cor.test(cpuesptot$cpuesptot[as.numeric(cpuesptot$year)>1990],ctcnums)
mtext(paste("cor = ",round(r$estimate,digits=2),", p = ",round (r$p.value, digits = 4)), side = 3,adj=0, line = -2, cex=1.1)
plot(ctcnums,cpuetot$cpuetot[as.numeric(cpuetot$year)>1990], pch=16, xlab = "Escapement estimates (#s)",  ylab="Albion Test Fishery (CPUE, 1April-20Oct)",col="gray",bty="l", main = paste(colnames(ctc)[i]), 
     cex.lab=1.5, cex.axis=1.5,cex=2, cex.main=1.5)
mod2<-lm(cpuetot$cpuetot[as.numeric(cpuetot$year)>1990]~ctcnums)
if(summary(mod)$coef[2,4]<0.1){abline(mod2, lwd=2)}
r<-cor.test(cpuetot$cpuetot[as.numeric(cpuetot$year)>1990],ctcnums)
mtext(paste("cor = ",round(r$estimate,digits=2),", p = ",round (r$p.value, digits = 4)), side = 3,adj=0, line = -2, cex=1.1)

}

dev.off()



#now fit splines suggested analysis
dat<-dat[dat$year>1993,]
dat$effort<-as.numeric(dat$effort)
dat$year2<-as.factor(dat$year)
dat$calDay<-as.numeric(dat$calDay)
dat$catch<-as.numeric(dat$catch)
dat$cpue1<-dat$cpue+.001
dat$logcpue<-log(dat$cpue1)

#Focal model is m2:
m2 <- brm(logcpue~ s(calDay) + (calDay|year2),
          data=dat, chains = 2,
          iter = 6000, warmup = 1000, thin = 10,
          control = list(adapt_delta = 0.99, max_treedepth=15))

plot(m2)

save(m2, file="analyses/output/albionchibrms.Rda")

#if already saved,
#load("analyses/output/albionchibrms.Rda")


albgam<-as.data.frame(cbind(dat$year,dat$calDay,fitted(m2),fitted(m2,probs=c(0.05,0.95)),fitted(m2,probs=c(0.25,0.75))))
colnames(albgam)[1:3]<-c("year","doy", "logcpue.est")
albgam$estcpue<-exp(as.numeric(albgam$logcpue.est))
colnames(albgam)[15]<-c("cpue.est")

#compare fitted to estimated cpue
sumest<-tapply(albgam$cpue.est,list(albgam$year),sum)

#write.csv(albgam,"analyses/output/albionbrmsests.csv", row.names = FALSE)


allyears<-unique(albgam$year)

seasons<-c("springsum","fall","allyear")
years<-c()
allseasons<-c()
firstobsdate<-c()
lastobsdate<-c()
midobsdate<-c()
peakobsdate<-c()
peakobsdate.sp<-c()
peakobsdate.fa<-c()
alltotal<-c()
alltotal.sp<-c()
alltotal.fa<-c()

for(y in allyears){
  datyr<-albgam[albgam$year==y,]
  if (dim(datyr)[1]<=1){
    first<-last<-mid<-peak<-NA
    total<-NA
  }
  if (dim(datyr)[1]>0){
    cpue<-as.numeric(datyr$cpue.est)
    cpuesp<-datyr$cpue.est[datyr$doy<213]#213= aug 1
    cpuefa<-datyr$cpue.est[datyr$doy>=213]
    #plot(datyr$doy,count, pch=21, bg="gray", main=paste(y))
    #if(y==min(allyears)){mtext(paste(sites[i], species[p]),side=3, line=3)}
    datdoy<-as.numeric(datyr$doy)
    datdoysp<-as.numeric(datyr$doy[datyr$doy<213])
    datdoyfa<-as.numeric(datyr$doy[datyr$doy>=213])
    first<-min(datdoy[which(cpue>0.005)])
    last<-max(datdoy[which(cpue>0.005)])
    total<-sum(cpue,na.rm=TRUE)
    totalsp<-sum(cpuesp,na.rm=TRUE)
    totalfa<-sum(cpuefa,na.rm=TRUE)
    
    mid<-datdoy[min(which(cumsum(cpue)>(total/2)))]#date at which half of fish have arrived
    peak<-min(datdoy[which(cpue==max(cpue, na.rm=TRUE))])#date of peak number of fish observed, if multiple dates with same number, choose first of these
    peaksp<-min(datdoysp[which(cpuesp==max(cpuesp, na.rm=TRUE))])#date of peak number of fish observed, if multiple dates with same number, choose first of these
    peakfa<-min(datdoyfa[which(cpuefa==max(cpuefa, na.rm=TRUE))])#date of peak number of fish observed, if multiple dates with same number, choose first of these
    #print(peak)
  }
  print(y);print(first);print(last);print(total); print(mid)
  years<-c(years,y)
  #allseasons<-c(allseasons,season[s])
  firstobsdate<-c(firstobsdate,first)
  lastobsdate<-c(lastobsdate,last)
  midobsdate<-c(midobsdate,mid)
  peakobsdate<-c(peakobsdate,peak)
  peakobsdate.fa<-c(peakobsdate.fa,peakfa)
  peakobsdate.sp<-c(peakobsdate.sp,peaksp)
  alltotal<-c(alltotal,total)
  alltotal.fa<-c(alltotal.fa,totalfa)
  alltotal.sp<-c(alltotal.sp,totalsp)
    }

#Save a file with these estimates in it
albchiphenest<-cbind("ck","albion",years,firstobsdate,lastobsdate,peakobsdate,peakobsdate.sp,peakobsdate.fa,midobsdate,alltotal,alltotal.sp,alltotal.fa)

colnames(albchiphenest)[1:3]<-c("sp","site","year")
#write.csv(albchiphenest,"analyses/output/albionchiphenestbrms.csv", row.names =FALSE)


#Now estimate trends in phenology
#restrict to time frame consistent with orcas
albchiphenest<-as.data.frame(albchiphenest)
albchiphenest<-albchiphenest[albchiphenest$year>1995,]
albchiphenest$year<-as.numeric(albchiphenest$year)
albchiphenest$firstobsdate<-as.numeric(albchiphenest$firstobsdate)
albchiphenest$peakobsdate<-as.numeric(albchiphenest$peakobsdate)
albchiphenest$lastobsdate
albchiphenest<-albchiphenest[albchiphenest$year<2018,]
firstmod<-lm(firstobsdate~year, data=albchiphenest)
firstcoefs<-coef(firstmod)
firstcoefs.50ci<-confint(firstmod,level = 0.50)
firstcoefs.75ci<-confint(firstmod,level = 0.75)
firstcoefs.95ci<-confint(firstmod,level = 0.95)

lastmod<-lm(lastobsdate~year, data=albchiphenest)
lastcoefs<-coef(lastmod)
lastcoefs.50ci<-confint(lastmod,level = 0.50)
lastcoefs.75ci<-confint(lastmod,level = 0.75)
lastcoefs.95ci<-confint(lastmod,level = 0.95)
peakmod<-lm(peakobsdate~year, data=albchiphenest)
peakcoefs<-coef(peakmod)
peakcoefs.50ci<-confint(peakmod,level = 0.50)
peakcoefs.75ci<-confint(peakmod,level = 0.75)
peakcoefs.95ci<-confint(peakmod,level = 0.95)
abundmod<-lm(alltotal~year, data=albchiphenest)
abundcoefs<-coef(abundmod)
abundcoefs.50ci<-confint(abundmod,level = 0.50)
abundcoefs.75ci<-confint(abundmod,level = 0.75)
abundcoefs.95ci<-confint(abundmod,level = 0.95)


allmodsums<-c(round(firstcoefs, digits=3),round(lastcoefs, digits=3),round(peakcoefs, digits=3))
allmodsums.50ci<-rbind(round(firstcoefs.50ci, digits=3),round(lastcoefs.50ci, digits=3),round(peakcoefs.50ci, digits=3))
allmodsums.75ci<-rbind(round(firstcoefs.75ci, digits=3),round(lastcoefs.75ci, digits=3),round(peakcoefs.75ci, digits=3))

allmodsums.95ci<-rbind(round(firstcoefs.95ci, digits=3),round(lastcoefs.95ci, digits=3),round(peakcoefs.95ci, digits=3))
phen<-c("first","first","last","last","peak","peak")
sums<-cbind("ck","albion",phen,allmodsums,allmodsums.50ci,allmodsums.75ci,allmodsums.95ci)
colnames(sums)<-c("sp","site","phen","est","ci25","ci75","ci12.5","ci87.5","ci2.5","ci97.5")

#write.csv(sums, "analyses/output/albionreturntrends_linmodyrsbrms.csv", row.names = TRUE)
