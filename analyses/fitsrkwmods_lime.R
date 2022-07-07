
#Fit multilevel bernouli gams with presence of SRKWs as the response, fit using brms package

 m2 <- brm(AllSRpres ~ s(day) + (day|year),
           data=limewdaysabs,
           family =bernoulli(), cores = 2,
           iter = 4000, warmup = 1000, thin = 10,
           control = list(adapt_delta = 0.99, max_treedepth=15))
 save(m2, file="analyses/output/sr.brms.Rda")
 #or, if already saved:
 #load("analyses/output/sr.brms.Rda")

j2 <- brm(Jpres ~ s(day) + (day|year),
          data=limewdaysabs,
          family =bernoulli(), cores = 2,
          iter = 4000, warmup = 1000, thin = 10,
control = list(adapt_delta = 0.99, max_treedepth=15))

save(j2, file="analyses/output/j.brms.Rda")
# or, if already saved:
# load("analyses/output/j.brms.Rda")
 
k2 <- brm(Kpres ~ s(day) + (day|year),
            data=limewdaysabs,
            family =bernoulli(), cores = 2,
            iter = 4000, warmup = 1000, thin = 10,
            control = list(adapt_delta = 0.99, max_treedepth=15))
save(k2, file="analyses/output/k.brms.Rda")
# or, if already saved:
# load("analyses/output/k.brms.Rda")

 l2 <- brm(Lpres ~ s(day) + (day|year),
            data=limewdaysabs,
            family =bernoulli(), chains = 2,
            iter = 4000, warmup = 1000, thin = 10,
            control = list(adapt_delta = 0.99, max_treedepth=15 ))
  save(l2, file="analyses/output/l.brms.Rda")
# or, if already saved:
# load("analyses/output/l.brms.Rda")

prob.occ.95<-cbind(limewdaysabs$year,limewdaysabs$day,fitted(m2),fitted(j2),fitted(k2),fitted(l2))

colnames(prob.occ.95)<-c("year","doy", paste("SRprob",colnames(fitted(m2)),sep="."),
                                    paste("Jprob",colnames(fitted(j2)),sep="."),
                                    paste("Kprob",colnames(fitted(k2)),sep="."),
                                    paste("Lprob",colnames(fitted(l2)),sep="."))
prob.occ.90<-cbind(limewdaysabs$year,limewdaysabs$day,fitted(m2,probs=c(0.05,0.95)),fitted(j2,probs=c(0.05,0.95)),fitted(k2,probs=c(0.05,0.95)),fitted(l2,probs=c(0.05,0.95)))

colnames(prob.occ.90)<-c("year","doy", paste("SRprob",colnames(fitted(m2,probs=c(0.05,0.95))),sep="."),
                         paste("Jprob",colnames(fitted(j2,probs=c(0.05,0.95))),sep="."),
                         paste("Kprob",colnames(fitted(k2,probs=c(0.05,0.95))),sep="."),
                         paste("Lprob",colnames(fitted(l2,probs=c(0.05,0.95))),sep="."))
prob.occ.50<-cbind(limewdaysabs$year,limewdaysabs$day,fitted(m2,probs=c(0.25,0.75)),fitted(j2,probs=c(0.25,0.75)),fitted(k2,probs=c(0.25,0.75)),fitted(l2,probs=c(0.25,0.75)))

colnames(prob.occ.50)<-c("year","doy", paste("SRprob",colnames(fitted(m2,probs=c(0.25,0.75))),sep="."),
                         paste("Jprob",colnames(fitted(j2,probs=c(0.25,0.75))),sep="."),
                         paste("Kprob",colnames(fitted(k2,probs=c(0.25,0.75))),sep="."),
                         paste("Lprob",colnames(fitted(l2,probs=c(0.25,0.75))),sep="."))

prob.occ.75<-cbind(limewdaysabs$year,limewdaysabs$day,fitted(m2,probs=c(0.125,0.875)),fitted(j2,probs=c(0.125,0.875)),fitted(k2,probs=c(0.125,0.875)),fitted(l2,probs=c(0.125,0.875)))

colnames(prob.occ.75)<-c("year","doy", paste("SRprob",colnames(fitted(m2,probs=c(0.125,0.875))),sep="."),
                         paste("Jprob",colnames(fitted(j2,probs=c(0.125,0.875))),sep="."),
                         paste("Kprob",colnames(fitted(k2,probs=c(0.125,0.875))),sep="."),
                         paste("Lprob",colnames(fitted(l2,probs=c(0.125,0.875))),sep="."))

#Save model results
write.csv(prob.occ.95,"analyses/output/lime_prob.occ.95.csv", row.names = FALSE)
write.csv(prob.occ.90,"analyses/output/lime_prob.occ.90.csv", row.names = FALSE)
write.csv(prob.occ.50,"analyses/output/lime_prob.occ.50.csv", row.names = FALSE)
write.csv(prob.occ.75,"analyses/output/lime_prob.occ.75.csv", row.names = FALSE)
