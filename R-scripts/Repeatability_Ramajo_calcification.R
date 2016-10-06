library(metafor)

meta.food.calcification.repeat<- read.csv(file.choose())
head(meta.food.calcification.repeat2)



#### Calculate effect sizes
meta.food.calcification.repeat2<-escalc(measure="ROM",m1i=Mean_Elevated, m2i=Mean_Ambient, sd1i=SD_Elevated,sd2i=SD_Ambient,n1i=N_Elevated,n2i=N_Ambient, data=meta.food.calcification.repeat,var.names=c("LRR","LRR_var"), digits=4)

summary(meta.food.calcification.repeat2)


#### BASIC INFO
#How many unqique observations? Number of papers?
length(unique(meta.food.calcification.repeat2$Observation_ID))
length(unique(meta.food.calcification.repeat2$Paper_no))

#####
####### PLOTTING THE DATA

######### FOREST PLOTS
# a few basic forest plots
forest(meta.food.calcification.repeat2$LRR,meta.food.calcification.repeat2$LRR_var,slab=meta.food.calcification.repeat2$Author,pch=19)

forest(meta.food.calcification.repeat2$LRR,meta.food.calcification.repeat2$LRR_var,slab=meta.food.calcification.repeat2$Paper_no,pch=19)

forest(meta.food.calcification.repeat2$LRR,meta.food.calcification.repeat2$LRR_var,slab=meta.food.calcification.repeat2$Food.supply,pch=19, order=order(meta.food.calcification.repeat2$Food.supply))


###### A more complicated, NICE forest plot
### decrease margins so the full space is used
par(mar=c(4,4,1,2))

### fit random-effects model (use slab argument to define study labels)
res <- rma.mv(yi=LRR,V=LRR_var, mods=~Food.supply,random= ~1|Paper_no,struct="CS",method="REML",digits=4,data=meta.food.calcification.repeat2, slab=paste(Author, Year, sep=", "))

### set up forest plot (with 2x2 table counts added; rows argument is used
### to specify exactly in which rows the outcomes will be plotted)
forest(res, xlim=c(-6, 3), at=log(c(.05, .25, 1, 4)), atransf=exp,
       ilab=cbind(meta.food.calcification.repeat2$Mean_Ambient, meta.food.calcification.repeat2$Mean_Elevated),
       ilab.xpos=c(-4,-3), cex=.70, ylim=c(-1,53),
       order=order(meta.food.calcification.repeat2$Food.supply), rows=c(2:16,22:36,42:47),
       xlab="LnRR", mlab="RE Model for All Studies", psize=1)

### set font expansion factor (as in forest() above) and use bold italic
### font and save original settings in object 'op'
op <- par(cex=.75, font=4)

### add text for the subgroups
text(-6, c(18,38,49), pos=4, c("High food",
                               "Low food",
                               "Med food"))

### switch to bold font
par(font=2)

### add column headings to the plot
text(c(-4,-3),52, c("Mean_Ambient", "Mean_Elevated"))
text(-6,                52, "Author(s) and Year",     pos=4)
text(3,                  52, "LnRR [95% CI]", pos=2)

### set par back to the original settings
par(op)


res <- rma.mv(yi=LRR,V=LRR_var, mods=~Food.supply,random= ~1|Paper_no,struct="CS",method="REML",digits=4,data=meta.food.calcification.repeat2, slab=paste(Author, Year, sep=", "))


### fit random-effects model in the three subgroups
res.m <-rma.mv(yi=LRR,V=LRR_var, subset=(Food.supply=="Med"),random= ~1|Paper_no,struct="CS",method="REML",digits=4,data=meta.food.calcification.repeat2)

res.l <-rma.mv(yi=LRR,V=LRR_var, subset=(Food.supply=="Low"),random= ~1|Paper_no,struct="CS",method="REML",digits=4,data=meta.food.calcification.repeat2)

res.h <-rma.mv(yi=LRR,V=LRR_var, subset=(Food.supply=="High"),random= ~1|Paper_no,struct="CS",method="REML",digits=4,data=meta.food.calcification.repeat2)

### add summary polygons for the three subgroups
addpoly(res.m, row=40.5, cex=.75, atransf=exp, mlab="RE Model for Subgroup")
addpoly(res.l, row= 20.5, cex=.75, atransf=exp, mlab="RE Model for Subgroup")
addpoly(res.h, row= 0.5, cex=.75, atransf=exp, mlab="RE Model for Subgroup")

addpoly(mod.med, row=40.5, cex=.75, atransf=exp, mlab="RE Model for Subgroup")
addpoly(mod.low, row= 20.5, cex=.75, atransf=exp, mlab="RE Model for Subgroup")
addpoly(mod.high, row= 0.5, cex=.75, atransf=exp, mlab="RE Model for Subgroup")
###This is the rma... can't have mods

#####

#### hierarchical model 
meta.food.calcification.repeat2<-escalc(measure="ROM",m1i=Mean_Elevated, m2i=Mean_Ambient, sd1i=SD_Elevated,sd2i=SD_Ambient,n1i=N_Elevated,n2i=N_Ambient, data=meta.food.calcification.repeat,var.names=c("LRR","LRR_var"), digits=4)

meta.food.calcification.repeat2$Ni <- unlist(lapply(split(meta.food.calcification.repeat2, meta.food.calcification.repeat2$study), function(x) rep(sum(x$N_Elevated) + x$N_Ambient[1], each=nrow(x))))
meta.food.calcification.repeat2$sdpi <- with(meta.food.calcification.repeat2, sqrt(((N_Elevated-1)*SD_Elevated^2 + (N_Ambient-1)*SD_Ambient^2)/(N_Elevated+N_Ambient-2)))
meta.food.calcification.repeat2$yi <- with(meta.food.calcification.repeat2, (Mean_Elevated-Mean_Ambient)/sdpi)
meta.food.calcification.repeat2$vi <- with(meta.food.calcification.repeat2, 1/N_Elevated + 1/N_Ambient + yi^2/(2*Ni))

calc.v <- function(x) {
  v <- matrix(1/x$N_Ambient[1] + outer(x$yi, x$yi, "*")/(2*x$Ni[1]), nrow=nrow(x), ncol=nrow(x))
  diag(v) <- x$vi
  v
}


V.calc.repeat <- bldiag(lapply(split(meta.food.calcification.repeat2, meta.food.calcification.repeat2$study), calc.v))
V.calc.repeat


###########Models

mod.overall.calc.repeat <- rma.mv(LRR, V.calc.repeat, mods = ~factor(Food.supply) - 1,random=list(~factor(DeltaCO2)|study, ~Food.supply|Paper_no), data=meta.food.calcification.repeat2)
profile(mod.overall.calc.repeat)
#over parameterized... 

# so simplify
mod.overall.calc.repeat.simple.3 <- rma.mv(LRR, V.calc.repeat, mods = ~factor(Food.supply)-1,random=~1|Paper_no, data=meta.food.calcification.repeat2)
profile(mod.overall.calc.repeat.simple.3)

logLik(mod.overall.calc.repeat)
logLik(mod.overall.calc.repeat.simple.3)
### same estimates though so it's okay. 

summary(mod.overall.calc.repeat)
anova(mod.overall.calc.repeat, L=c(1,-1,0)) #high vs. low
anova(mod.overall.calc.repeat, L=c(1,0,-1)) # high vs. med
anova(mod.overall.calc.repeat, L=c(0,1,-1)) # low vs. med


qqnorm(residuals(mod.overall.calc.repeat,type="pearson"),main="QQ plot: residuals")
qqline(residuals(mod.overall.calc.repeat,type="pearson"),col="red")
plot(fitted(mod.overall.calc.repeat) ~ residuals(mod.overall.calc.repeat))
# normal ish - enough

#Is UN the best?
mod.overall.calc.repeat.UN <- rma.mv(LRR, V.calc.repeat, mods = ~factor(Food.supply) - 1,random=list(~factor(DeltaCO2)|study, ~Food.supply|Paper_no), data=meta.food.calcification.repeat2, struct="UN")
summary(mod.overall.calc.repeat.UN)
logLik(mod.overall.calc.repeat)
logLik(mod.overall.calc.repeat.UN)
### slightly different logLik ... so different model 
anova(mod.overall.calc.repeat.UN, L=c(1,-1,0)) #high vs. low
anova(mod.overall.calc.repeat.UN, L=c(1,0,-1)) # high vs. med
anova(mod.overall.calc.repeat.UN, L=c(0,1,-1)) # low vs. med
# regardless, all are NS


#### Can do robust versions of these models
mod.overall.calc.repeat.UN.robust<-robust(mod.overall.calc.repeat.UN, cluster = meta.food.calcification.repeat2$study)
mod.overall.calc.repeat.UN.robust
anova(mod.overall.calc.repeat.UN.robust, L=c(1,-1, 0))


### robust.
mod.overall.calc.repeat.robust<-robust(mod.overall.calc.repeat.UN, cluster=meta.food.calcification.repeat2$study)
summary(mod.overall.calc.repeat.robust)
### unsure about these ... which is correct? 


### multiple models for each level
mod.low <- rma.mv(LRR, V.calc.repeat, mods = ~ 1,random = ~factor(DeltaCO2)|Paper_no, data=meta.food.calcification.repeat2, subset=(Food.supply=="Low"))
mod.high <- rma.mv(LRR, V.calc.repeat, mods = ~ 1,random = ~factor(DeltaCO2)|Paper_no, data=meta.food.calcification.repeat2, subset=(Food.supply=="High"))
mod.med <- rma.mv(LRR, V.calc.repeat, mods = ~1,random= ~factor(DeltaCO2)|Paper_no, data=meta.food.calcification.repeat2, subset=(Food.supply=="Med"))


###Adding in pH 
#overall by pH
mod.overall.calc.repeat.pH <- rma.mv(LRR, V.calc.repeat, mods = ~DeltaCO2:factor(Food.supply) - 1,random=list(~factor(DeltaCO2)|study, ~Food.supply|Paper_no), data=meta.food.calcification.repeat2)
summary(mod.overall.calc.repeat.pH)
profile(mod.overall.calc.repeat.pH)
# not very well estimated over parameterized
mod.overall.calc.repeat.simple.3.pH <- rma.mv(LRR, V.calc.repeat, mods = ~DeltaCO2:factor(Food.supply)-1,random=~1|Paper_no, data=meta.food.calcification.repeat2)
profile(mod.overall.calc.repeat.simple.3.pH)
logLik(mod.overall.calc.repeat.pH)
logLik(mod.overall.calc.repeat.simple.3.pH)
### same estimate
anova(mod.overall.calc.repeat.pH, L=c(1,-1,0)) #high vs. low


mod.overall.calc.repeat.pH.UN <- rma.mv(LRR, V.calc.repeat, mods = ~DeltaCO2:factor(Food.supply) - 1,random=list(~factor(DeltaCO2)|study, ~Food.supply|Paper_no), data=meta.food.calcification.repeat2, struct="UN")
summary(mod.overall.calc.repeat.pH.UN)
### not sure which is better
#same result though


## could try robust methods for by pH as well
mod.overall.calc.repeat.simple.3.pH.robust<-robust(mod.overall.calc.repeat.simple.3.pH, cluster = meta.food.calcification.repeat2$study)
print(mod.overall.calc.repeat.simple.3.pH.robust)
anova(mod.overall.calc.repeat.simple.3.pH.robust, L=c(1,-1, 0))






###### Checking biases
funnel(mod.overall.calc.repeat,ylim=c(1:6,by=2),yaxis="seinv",level=c(90, 95, 99),ylab="Precision (1/SE)" ,shade=c("white", "gray", "darkgray"), refline=0)
funnel(mod.1,ylim=c(1:6,by=2),yaxis="seinv",level=c(90, 95, 99),ylab="Precision (1/SE)" ,shade=c("white", "gray", "darkgray"), refline=0)

#Publication bias with more basic model (some can't handle ram.mv)
mod.basic.calc.repeat<-rma(LRR ~ factor(Food.supply),LRR_var,data=meta.food.calcification.repeat2)
funnel(mod.basic)
regtest(mod.basic.calc.repeat,model="rma",predictor="sei")
funnel(trimfill(mod.basic, side="right"))
##### The method can be used to estimate the number of studies missing from a meta-analysis due to the suppression of the most extreme results on one side of the funnel plot.
baujat(mod.basic.calc.repeat)
#### Baujat et al. (2002) proposed a diagnostic plot to detect sources of heterogeneity in meta-analytic data. 
#The plot shows the contribution of each study to the overall Q-test statistic for heterogeneity on the horizontal axis 
#versus the influence of each study (defined as the standardized squared difference between the overall estimate based 
#on a fixed-effects model with and without the ith study included in the model) on the vertical axis. An example of such a plot is shown below.

### calculate influence diagnostics
inf <- influence(mod.basic.calc)
### plot the influence diagnostics
plot(inf, layout=c(8,1))

########## Three studies have rather large residuals and may be considered outliers. 
#Only two of those studies actually have a strong influence on the results (as reflected, for example, 
#in their Cook's distances). Removal of these studies would reduce the amount of heterogeneity quite a bit 
#and increase the precision of the estimated average outcome (i.e., r-to-z transformed correlation) from the 
#random-effects model. However, instead of just removing those studies, one should examine them in detail to 
#determine what the reason may be for their unusual results. Outliers and influential cases can actually reveal 
#patterns that may lead to new insights about study characteristics that could be acting as potential moderators (Light & Pillemer, 1984).


###Sensitivity analysis ## analgous to drop1
# RULE:If rstandard >3 AND hatvalue >2 times average of hatvalues, 
# run analysis with those cases deleted to test for sensitivity.
rs.abund.wood <- rstandard(mod.overall.calc.repeat)
hat.abund.wood <- hatvalues(mod.overall.calc.repeat)/mean(hatvalues(mod.overall.calc.repeat))
plot(hat.abund.wood, rs.abund.wood$resid, ylim = c(-4.0,4))
text(hat.abund.wood, rs.abund.wood$resid, labels = meta.food.calcification.repeat2$Paper_no, cex= 1, pos = 2)
abline(h = -3)
abline(h = 3)
abline( v = 2)
#### should compare without Crook, Comeau, Hettinger outliers model fit







############# GGPLOT using summary data for each Food.supply level
library(ggplot2)

######## PLOTTING WITH A MODEL


summary(mod.1.calc.repeat)
head(mod.1.calc.repeat)

library(dplyr)
library(plyr)

cdata.LRR.calc.repeat<-ddply(meta.food.calcification.repeat2,~Food.supply,summarise,mean=mean(LRR))
head(cdata.LRR.calc.repeat)
cdata.LRR.calc.repeat_var<-ddply(meta.food.calcification.repeat2,~Food.supply,summarise,mean=mean(LRR_var))
head(cdata.LRR.calc.repeat_var)
cdata.LRR.calc.repeat$var<-cdata.LRR.calc.repeat_var$mean
head(cdata.LRR.calc.repeat)


plot.mod.overall.calc.repeat<-ggplot(cdata.LRR.calc.repeat,aes(x=Food.supply,y=mod.overall.calc.repeat$b,ymax=(mod.overall.calc.repeat$ci.ub),ymin=(mod.overall.calc.repeat$ci.lb),size=2))
plot.mod.overall.calc.repeat<-plot.mod.overall.calc.repeat+geom_pointrange(size=1)
plot.mod.overall.calc.repeat<-plot.mod.overall.calc.repeat #+coord_flip()
plot.mod.overall.calc.repeat<-plot.mod.overall.calc.repeat+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)+ylim(-.8, 0.8)
plot.mod.overall.calc.repeat<-plot.mod.overall.calc.repeat+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
plot.mod.overall.calc.repeat<-plot.mod.overall.calc.repeat+ xlab('Food.supply') +ylab ('Effect of CO2')
plot.mod.overall.calc.repeat<-plot.mod.overall.calc.repeat+theme(axis.text.x = element_text(size = 16, colour = 'black')) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.mod.overall.calc.repeat+theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))




plot.mod.overall.calc.repeat.robust<-ggplot(cdata.LRR.calc.repeat,aes(x=Food.supply,y=mod.overall.calc.repeat.robust$b,ymax=(mod.overall.calc.repeat.robust$ci.ub),ymin=(mod.overall.calc.repeat.robust$ci.lb),size=2))
plot.mod.overall.calc.repeat.robust<-plot.mod.overall.calc.repeat.robust+geom_pointrange(size=1)
plot.mod.overall.calc.repeat.robust<-plot.mod.overall.calc.repeat.robust+coord_flip()
plot.mod.overall.calc.repeat.robust<-plot.mod.overall.calc.repeat.robust+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.mod.overall.calc.repeat.robust<-plot.mod.overall.calc.repeat.robust+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
plot.mod.overall.calc.repeat.robust<-plot.mod.overall.calc.repeat.robust+ xlab('Food.supply') +ylab ('Calcification response to CO2')
plot.mod.overall.calc.repeat.robust<-plot.mod.overall.calc.repeat.robust+theme(axis.text.x = element_text(size = 16, colour = 'black')) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.mod.overall.calc.repeat.robust+theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))




plot.mod.1<-ggplot(cdata.LRR.calc.repeat,aes(x=Food.supply,y=mod.1$b,ymax=(mod.1$ci.ub),ymin=(mod.1$ci.lb),size=2))
plot.mod.1<-plot.mod.1+geom_pointrange(size=1)
plot.mod.1<-plot.mod.1+coord_flip()
plot.mod.1<-plot.mod.1+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.mod.1<-plot.mod.1+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
plot.mod.1<-plot.mod.1+ xlab('Food.supply') +ylab ('Effect of CO2')
plot.mod.1<-plot.mod.1+theme(axis.text.x = element_text(size = 16, colour = 'black')) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.mod.1+theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))


plot.mod.1.robust<-ggplot(cdata.LRR.calc.repeat,aes(x=Food.supply,y=mod.1.robust$b,ymax=(mod.1.robust$ci.ub),ymin=(mod.1.robust$ci.lb),size=2))
plot.mod.1.robust<-plot.mod.1.robust+geom_pointrange(size=1)
plot.mod.1.robust<-plot.mod.1.robust#+coord_flip()
plot.mod.1.robust<-plot.mod.1.robust+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)+ ylim(-0.5, 0.5)
plot.mod.1.robust<-plot.mod.1.robust+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
plot.mod.1.robust<-plot.mod.1.robust+ xlab('Food.supply') +ylab ('Calcification response to CO2')
plot.mod.1.robust<-plot.mod.1.robust+theme(axis.text.x = element_text(size = 16, colour = 'black')) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.mod.1.robust<-plot.mod.1.robust+theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))







########### PLOTTING WITHOUT A MODEL
###Can also plot this without a model ... just estimating from the LRR_variance
###Not sure which is better? has a point for each study


plot.not.mod<-ggplot(meta.food.calcification.repeat2,aes(x=Food.supply,y=LRR,ymax=(LRR + 1.96*sqrt(LRR_var)),ymin=(LRR - 1.96*sqrt(LRR_var)),size=2))
plot.not.mod<-plot.not.mod+geom_pointrange(size=1)
plot.not.mod<-plot.not.mod+coord_flip()
plot.not.mod<-plot.not.mod+geom_hline(aes(x=0, yintercept=0), lty=2,size=1)
plot.not.mod<-plot.not.mod+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
plot.not.mod<-plot.not.mod+ xlab('Food.supply') +ylab ('Effect of CO2')
plot.not.mod<-plot.not.mod+theme(axis.text.x = element_text(size = 16, colour = 'black')) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.not.mod+theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))


######## with bootstrapped CIS

meta.food.calcification.repeat.low<-subset(meta.food.calcification.repeat2, meta.food.calcification.repeat2$Food.supply=="Low")
meta.food.calcification.repeat.high<-subset(meta.food.calcification.repeat2, meta.food.calcification.repeat2$Food.supply=="High")
meta.food.calcification.repeat.med<-subset(meta.food.calcification.repeat2, meta.food.calcification.repeat2$Food.supply=="Med")

head(meta.food.calcification.repeat.med)

boot.func.2 <- function(dat, indices) {
  
  res <- try(rma(LRR, LRR_var, data=dat, subset=indices), silent=TRUE)
  
  if (is.element("try-error", class(res))) {
    NA
  } else {
    c(coef(res), vcov(res), res$tau2, res$se.tau2^2)
  }
  
}
res.boot.low.calcification.repeat.2 <- boot(meta.food.calcification.repeat.low, boot.func.2, R=10000)
res.boot.low.calcification.repeat.2.cis<-boot.ci(res.boot.low.calcification.repeat.2,index=1:2)
res.boot.high.calcification.repeat.2 <- boot(meta.food.calcification.repeat.high, boot.func.2, R=10000)
res.boot.high.calcification.repeat.2.cis<-boot.ci(res.boot.high.calcification.repeat.2,index=1:2)
res.boot.med.calcification.repeat.2 <- boot(meta.food.calcification.repeat.med, boot.func.2, R=10000)
res.boot.med.calcification.repeat.2.cis<-boot.ci(res.boot.med.calcification.repeat.2,index=1:2)



lower <- rbind((res.boot.high.calcification.repeat.2.cis[[8]][4]), (res.boot.low.calcification.repeat.2.cis[[8]][4]),  (res.boot.med.calcification.repeat.2.cis[[8]][4]))
upper <-rbind((res.boot.high.calcification.repeat.2.cis[[8]][5]), (res.boot.low.calcification.repeat.2.cis[[8]][5]), (res.boot.med.calcification.repeat.2.cis[[8]][5]))
food.supply <- factor(c('High','Low','Med'))
confidence_intervals <- as.data.frame(cbind(food.supply=levels(food.supply)[food.supply], lower, upper))

colnames(confidence_intervals)<- c("food.supply","lower","upper")

head(meta.food.calcification.repeat2)
cdata.LRR.calcification.repeat<-ddply(meta.food.calcification.repeat2,~Food.supply,summarise,mean=mean(LRR))
head(cdata.LRR.calcification.repeat)
cdata.LRR.calcification.repeat_var<-ddply(meta.food.calcification.repeat2,~Food.supply,summarise,mean=mean(LRR_var))
head(cdata.LRR.calcification.repeat_var)
cdata.LRR.calcification.repeat$var<-cdata.LRR.calcification.repeat_var$mean
head(cdata.LRR.calcification.repeat)

cdata.LRR.calcification.repeat<-cbind(cdata.LRR.calcification.repeat,lower, upper)
mean.scaled.calc<-rbind(res.boot.high.calcification.repeat.2.cis$t0, res.boot.low.calcification.repeat.2.cis$t0, res.boot.med.calcification.repeat.2.cis$t0)
cdata.LRR.calcification.repeat<-cbind(cdata.LRR.calcification.repeat,lower, upper, mean.scaled.calc)


##### With bootstrapped data

plot.not.mod.2.calcification.repeat<-ggplot(cdata.LRR.calcification.repeat, aes(x=Food.supply,y=mean,ymax=upper,ymin=lower))
plot.not.mod.2.calcification.repeat<-plot.not.mod.2.calcification.repeat+geom_pointrange(size=1)
plot.not.mod.2.calcification.repeat<-plot.not.mod.2.calcification.repeat#+coord_flip()
plot.not.mod.2.calcification.repeat<-plot.not.mod.2.calcification.repeat+geom_hline(aes(x=0, yintercept=0), lty=2,size=1) + ylim(-.5, .5)
plot.not.mod.2.calcification.repeat<-plot.not.mod.2.calcification.repeat+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
plot.not.mod.2.calcification.repeat<-plot.not.mod.2.calcification.repeat+ xlab('Food supply') +ylab ('LnRR Calcification')
plot.not.mod.2.calcification.repeat<-plot.not.mod.2.calcification.repeat+theme(axis.text.x = element_text(size = 16, colour = 'black')) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.not.mod.2.calcification.repeat<-plot.not.mod.2.calcification.repeat +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.not.mod.2.calcification.repeat

require(cowplot)
theme_set(theme_classic())


###with mean_scaled
plot.not.mod.2.calcification.repeat<-ggplot(cdata.LRR.calcification.repeat, aes(x=Food.supply,y=intrcpt,ymax=upper,ymin=lower))
plot.not.mod.2.calcification.repeat<-plot.not.mod.2.calcification.repeat+geom_pointrange(size=1)
plot.not.mod.2.calcification.repeat<-plot.not.mod.2.calcification.repeat#+coord_flip()
plot.not.mod.2.calcification.repeat<-plot.not.mod.2.calcification.repeat+geom_hline(aes(x=0, yintercept=0), lty=2,size=1) + ylim(-.5, .5)
plot.not.mod.2.calcification.repeat<-plot.not.mod.2.calcification.repeat+theme(panel.grid.minor=element_blank(), panel.grid.major=element_blank())+theme(axis.ticks = element_blank(), axis.text.x = element_blank())+ theme_bw()+theme(legend.position = "none")
plot.not.mod.2.calcification.repeat<-plot.not.mod.2.calcification.repeat+ xlab('Food supply') +ylab ('LnRR Calcification')
plot.not.mod.2.calcification.repeat<-plot.not.mod.2.calcification.repeat+theme(axis.text.x = element_text(size = 16, colour = 'black')) +theme(axis.text.y = element_text(size = 16, colour = 'black'))
plot.not.mod.2.calcification.repeat<-plot.not.mod.2.calcification.repeat +theme(axis.title.x = element_text(size = 20, colour = 'black'))+theme(axis.title.y = element_text(size = 20, colour = 'black'))
plot.not.mod.2.calcification.repeat





plot_grid(plot.not.mod.2.calcification.repeat,plot.not.mod.2.growth.repeat, ncol=2, align='h', labels=c('(a)', '(b)', label_size=12))




#### Plotting for publication

require(cowplot)
theme_set(theme_classic())



plot_grid(plot.not.mod.2.calcification.repeat,plot.not.mod.2.growth.repeat, ncol=2, align='h', labels=c('(a)', '(b)', label_size=12))



plot_grid(plot.mod.1.robust,plot.mod.overall.calc.repeat.growth.robust , ncol=2, align='h', labels=c('(a)', '(b)', label_size=12))


####################### PLots by pH, CO2, etc... 

######### By pH

mod.overall.calc.repeat.pH <- rma.mv(LRR, V.calc.repeat, mods = ~DeltaCO2*factor(Food.supply)-1,random=list(~factor(DeltaCO2)|study, ~factor(Food.supply)|Paper_no), data=meta.food.calcification.repeat2)
summary(mod.overall.calc.repeat.pH)

mod.3<-rma.mv(yi=LRR,V=V.calc.repeat, mods=~DeltapH,random= ~1|study,struct="CS",method="REML",digits=4,data=meta.food.calcification.repeat2)
mod.3

mod.overall.calc.repeat.pH <- rma.mv(LRR, V.calc.repeat, mods = ~DeltaCO2:factor(Food.supply),random=list(~factor(DeltaCO2)|study, ~factor(Food.supply)|Paper_no), data=meta.food.calcification.repeat2)
summary(mod.overall.calc.repeat.pH)
anova(mod.overall.calc.repeat.pH, L=c(1,-1,0))


mod.overall.calc.repeat.pH.robust<-robust(mod.overall.calc.repeat.pH, cluster=meta.food.calcification.repeat2$study)
summary(mod.overall.calc.repeat.pH.robust)
anova(mod.overall.calc.repeat.pH.robust, L=c(1,-1,0))


### calculate predicted relative risks for 0 to 60 degrees absolute latitude
preds <- predict(mod.3, newmods=c(0.2:.90), transf=exp)

### calculate point sizes by rescaling the standard errors
wi    <- 1/sqrt(meta.food.calcification.repeat2$LRR_var)
size  <- 0.5 + 3.0 * (wi - min(wi))/(max(wi) - min(wi))

### plot the relative risks against absolute latitude
plot(meta.food.calcification.repeat2$DeltapH, exp(meta.food.calcification.repeat2$LRR), pch=19, cex=size, xlab="Delta pH", ylab="LRR",las=1, bty="l", log="y",col=meta.food.calcification.repeat2$Food.supply)

### add predicted values (and corresponding CI bounds)
lines(c(0.2:.90), preds$pred)
lines(0.2:.90, preds$ci.lb, lty="dashed")
lines(0.2:.90, preds$ci.ub, lty="dashed")

### dotted line at RR=1 (no difference between groups)
abline(h=1, lty="dotted")



####### GGPLOT by CO2
head(meta.food.calcification.repeat2)

plot.LRR<- ggplot(meta.food.calcification.repeat2, aes(x=meta.food.calcification.repeat2$DeltapH, y=meta.food.calcification.repeat2$LRR, colour=Food.supply)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(fill="grey80", size=1.5, method='lm')

plot.LRR<- plot.LRR + theme_bw() + xlab(bquote('Delta pH')) + ylab(expression("LRR"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))

plot.LRR<- plot.LRR + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))
plot.LRR


plot.LRR.CO2<- ggplot(meta.food.calcification.repeat2, aes(x=meta.food.calcification.repeat2$DeltaCO2, y=meta.food.calcification.repeat2$LRR, colour=Food.supply)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(fill="grey80", size=1.5, method='lm')

plot.LRR.CO2<- plot.LRR.CO2 + theme_bw() + xlab(bquote('Delta CO2')) + ylab(expression("LRR"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))

plot.LRR.CO2<- plot.LRR.CO2 + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))
plot.LRR.CO2



plot.LRR.CO2.max<- ggplot(meta.food.calcification.repeat2, aes(x=meta.food.calcification.repeat2$CO2.treatment..elevated..uatm., y=meta.food.calcification.repeat2$LRR, colour=Food.supply)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(fill="grey80", size=1.5, method='lm')

plot.LRR.CO2.max<- plot.LRR.CO2.max + theme_bw() + xlab(bquote('Max CO2')) + ylab(expression("LRR"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))

plot.LRR.CO2.max<- plot.LRR.CO2.max + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))
plot.LRR.CO2.max


plot.LRR.CO2.min<- ggplot(meta.food.calcification.repeat2, aes(x=meta.food.calcification.repeat2$CO2.treatment..control..uatm., y=meta.food.calcification.repeat2$LRR, colour=Food.supply)) + geom_point(size=5) + guides(fill=FALSE) + scale_fill_manual(values=cbbPalette.all)+ geom_smooth(fill="grey80", size=1.5, method='lm')

plot.LRR.CO2.min<- plot.LRR.CO2.min + theme_bw() + xlab(bquote('Min CO2')) + ylab(expression("LRR"))  + theme(text = element_text(size=16), axis.text = element_text(size=16))+theme(axis.title.y = element_text(angle=90))

plot.LRR.CO2.min<- plot.LRR.CO2.min + theme(legend.text = element_text(colour="black", size = 16))+ theme(legend.title = element_text(colour="black", size=16))+theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank(), panel.grid.major.y = element_blank(), panel.grid.minor.y = element_blank(), axis.line=element_line(size=0.25), axis.ticks.length=unit(-0.25, "cm"), axis.ticks.margin=unit(0.5, "cm"))
plot.LRR.CO2.min



###  Visualizing taxon and Food supply 
ggplot(data = meta.food.calcification.repeat2, aes(y = LRR, x = Food.supply, col=Taxon.big)) + geom_point(size = 2) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = (LRR - 1.96*sqrt(LRR_var)), ymax = (LRR + 1.96*sqrt(LRR_var))), width = 0.1) + ggtitle("calcification responses") +
  ylab("log response ratio") + xlab("Food.supply")


##### Visiualizing taxon and Food.relative.score
ggplot(data = meta.food.calcification.repeat2, aes(y = LRR, x = Food.relative.score, col=Taxon.big)) + geom_point(size = 2) +
  geom_hline(yintercept = 0) +
  geom_errorbar(aes(ymin = (LRR - 1.96*sqrt(LRR_var)), ymax = (LRR + 1.96*sqrt(LRR_var))), width = 0.1) + ggtitle("calcification responses") +
  ylab("log response ratio") + xlab("Food.relative.score")







#### Bootstrapping method ######### Only for random effects model.... or independent models
meta.food.calcification.repeat.low<-subset(meta.food.calcification.repeat2, meta.food.calcification.repeat2$Food.supply=="Low")
meta.food.calcification.repeat.high<-subset(meta.food.calcification.repeat2, meta.food.calcification.repeat2$Food.supply=="High")

install.packages("boot")
library(boot)
boot.func <- function(data.boot) {
  
  res <- try(rma(yi, vi, data=data.boot), silent=TRUE)
  
  if (is.element("try-error", class(res))) {
    NA
  } else {
    c(coef(res), vcov(res), res$tau2, res$se.tau2^2)
  }
  
}

data.gen <- function(meta.food.calcification.repeat2, mle) {
  data.frame(yi=rnorm(nrow(meta.food.calcification.repeat2),
                      mle$mu, sqrt(mle$tau2 + meta.food.calcification.repeat2$vi)), 
             vi=meta.food.calcification.repeat2$vi)}

res.boot <- boot(meta.food.calcification.repeat2, boot.func, R=10000, sim="parametric", ran.gen=data.gen, mle=list(mu=coef(res), tau2=res$tau2))
res.boot

boot.ci(res.boot, type=c("norm", "basic", "stud", "perc"), index=1:2)

#low parametric
data.gen.low <- function(meta.food.calcification.repeat.low, mle) {
  data.frame(yi=rnorm(nrow(meta.food.calcification.repeat.low),
                      mle$mu, sqrt(mle$tau2 + meta.food.calcification.repeat.low$vi)), 
             vi=meta.food.calcification.repeat.low$vi)}

res.boot.low <- boot(meta.food.calcification.repeat.low, boot.func, R=10000, sim="parametric", ran.gen=data.gen.low, mle=list(mu=coef(res), tau2=res$tau2))
res.boot.low
boot.ci(res.boot.low,index=1:2)


###OR non parametric for low
boot.func.2 <- function(dat, indices) {
  
  res <- try(rma(LRR, LRR_var, data=dat, subset=indices), silent=TRUE)
  
  if (is.element("try-error", class(res))) {
    NA
  } else {
    c(coef(res), vcov(res), res$tau2, res$se.tau2^2)
  }
  
}
res.boot.low.2 <- boot(meta.food.calcification.repeat.low, boot.func.2, R=10000)
boot.ci(res.boot.low.2,index=1:2)

###be careful bc vi was already in the datafram as something else
#### Can probably estimate bootstrapped CI2 from full model - just needs to be in the function ... when would we do that?? 




#high
data.gen.high <- function(meta.food.calcification.repeat.high, mle) {
  data.frame(yi=rnorm(nrow(meta.food.calcification.repeat.high),
                      mle$mu, sqrt(mle$tau2 + meta.food.calcification.repeat.high$vi)), 
             vi=meta.food.calcification.repeat.high$vi)}

res.boot.high <- boot(meta.food.calcification.repeat.high, boot.func, R=10000, sim="parametric", ran.gen=data.gen.high, mle=list(mu=coef(res), tau2=res$tau2))
boot.ci(res.boot.high,index=1:2)

######### Not necessarily recommended... 

# how do I do without a model?? is it just random effects? 
random.mod <- rma(yi, vi,data=meta.food.calcification.repeat.low)
confint(random.mod)

head(meta.food.calcification.repeat.low)

########### Bootstrapping WITHOUT an rma object:






##############################################
#Jarrett's page ... a bit confusing

#Bootstrap 95% Confidence Intervals

#Set number of bootstraps

B=999

#Create dataframe to store predicted values

boot.mat=matrix(NA,nrow=nrow(meta.food.calcification.repeat2),ncol=B)

#Get bootstrapped values and store in a matrix (include progress bar)

pb=txtProgressBar(min=0,max=B,style=3)

for(j in 1:B) \{
  
  df2=df[sample(nrow(df),replace=T),]
  
  df2=df2[order(df2[,"richness"]),]
  
  df2=groupedData(value~richness|Study.id.marine,data=df2)
  
  bestmod2=try(update(bestmod,data=df2),TRUE)
  
  ifelse(isTRUE(class(bestmod2)=="try-error"),next,bestmod2)
  
  if(bestmod.name=="Linear") \{
    
    boot.mat[,j]=fixef(bestmod2)[1]+fixef(bestmod2)\hich\af4\dbch\af31505\loch\f4 [2]*pred$richness 
    
    \} else if(bestmod.name=="Logarithmic") \{
      
      boot.mat[,j]=fixef(bestmod2)[1]+fixef(bestmod2)[2]*log(pred$richness)
      
      \} else if(bestmod.name=="Power") \{
        
        boot.mat[,j]=fixef(bestmod2)[1]*pred$richness^fixef(bestmod2)[2]
        
        \} \hich\af4\dbch\af31505\loch\f4 else if(bestmod.name=="Exponential") \{
          
          boot.mat[,j]=exp(fixef(bestmod2)[1]+fixef(bestmod2)[2]*pred$richness)
          
          \} else \{
            
            boot.mat[,j]=(max(df$value)*pred$richness)/(fixef(bestmod2)[1]+pred$richness)
            
            \} 
  
  setTxtProgressBar(pb,j) \}

close(pb)

#Extract 95% CIs on bootstrapped values

CIboot.mat=adply(boot.mat,1,function(x) quantile(x,prob=c(0.025,0.5,0.975),na.rm=T) )

colnames(CIboot.mat)=c("iter","LL","M","UL")

CIboot.mat$richness=pred[,1]

#Prep dataframe for plotting

boo\hich\af4\dbch\af31505\loch\f4 t.mat2=melt(boot.mat)

colnames(boot.mat2)=c("index","B","value")


