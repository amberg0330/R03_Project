setwd("~/Desktop/John Hughes/Report")
hpv=read.table("hpv.adj.txt",header=T,sep="\t") #cases w/ round!!

###############################statewise IR#################################
setwd("~/Desktop/John Hughes/data")
state=read.csv("statecruderate.csv",header=T)
library(ggplot2)
qplot(year,rate,data=state,color=site,group=site,ylab="Crude IR per 100,000 people",ylim=c(0,17))+
  geom_line()+geom_point()+
  facet_wrap(~gender)+
  labs(title="Texas Crude IR of HPV related cancer 2006-2012")

state.aj=read.csv("stateaj.csv",header=T)
qplot(year,rate,data=state.aj,color=site,group=site,ylab="Crude IR per 100,000 people",ylim=c(0,17))+
  geom_line()+geom_point()+
  facet_wrap(~gender)+
  labs(title="Texas age-adjusted IR of HPV related cancer 2006-2012")


###############################exploratory analysis#####################################
setwd("~/Desktop/John Hughes/data")
dat=read.table("dat.txt",sep="\t",header=T)
#delete singularity
#-17:per.blkAA
#-14:AGE20.29
#delete correlated covariates
x=dat[,c(11:34)]
cor(dat[,c(30,31)]) #-30: mm.hh.inc 
cor(dat[,c(30,31,32,33)]) #-33:perc.under18inpov
x.tol=x[,-c(4,7)] #-22:per.allinpov???
X=scale(x,center=T,scale=F)
dat.tol=cbind(dat[,c(1:10)],x);dim(dat.tol) #original data without highly correlated covariates
dat.cen=cbind(dat[,c(1:10)],X);dim(dat.cen) #centered data

dat.tol=na.omit(dat.tol);dim(dat.tol) #only 213 left 
dat.cen=na.omit(dat.cen);dim(dat.cen)
################################Linear model###########################
##check the correlation of covariates
View(cor(x.tol))

#correlation?
poo = as.vector(cor(x.tol))
poo[poo == 1] = 0
range(poo)
plot(poo, pch = 20)

library(plyr)
colwise(class)(dat.tol)
################calculate type II error###############
#Anova from car package--use type II sum of squares 
#drop1 
library(car)

##################
###linear model###
##################

fit.lm=lm(CIR.tol~per.m+AGE0.14+AGE15.19+AGE.30+per.wht+per.AIAN+per.asian+
            per.oth+per.hl+Chlam.IR+per.smk.tol+perc.bechelor+perc.high+HI_18+HI18_64+
            percap.inc+perc.allinpov+perc.under18inpov+perc.urb,data=dat.tol)
summary(step(fit.lm))

fit.lm1=lm(CIR.tol~per.m+AGE0.14+AGE15.19+per.wht+per.AIAN+
            per.oth+per.hl+per.smk.tol+perc.bechelor+HI_18+
            perc.allinpov+perc.urb,data=dat.tol)
summary(fit.lm1)

#by significance 
fit.lm2=lm(CIR.tol~per.m+AGE0.14+AGE15.19+per.wht+per.AIAN+
            per.oth+per.hl+per.smk.tol+perc.bechelor+HI_18+
            perc.allinpov+perc.urb,data=dat.tol)
summary(fit.lm2)

#using Type II to select model
#for linear model-F test
Anova(fit.lm1,type='II')

#backwards
drop1(fit.lm,test='F')
new.lm=update(fit.lm1,.~.-perc.high)
drop1(new.lm,test='F')

#forwards
basic.lm=lm(CIR.tol~1,data=dat.tol)
add1(basic.lm,test='F',
     scope= ~per.m+AGE0.14+AGE15.19+AGE.30+per.wht+per.AIAN+per.asian+
       per.oth+per.hl+Chlam.IR+per.smk.tol+perc.bechelor+perc.high+HI_18+HI18_64+
       percap.inc+perc.allinpov+perc.under18inpov+perc.urb)
new.lm2=update(basic.lm,.~.+AGE.30)
add1(new.lm2,test='F',
     scope=~per.m+AGE0.14+AGE15.19+AGE.30+per.wht+per.AIAN+per.asian+
       per.oth+per.hl+Chlam.IR+per.smk.tol+perc.bechelor+perc.high+HI_18+HI18_64+
       percap.inc+perc.allinpov+perc.under18inpov+perc.urb)
new.lm2=update(new.lm2,.~.+per.smk.tol)
add1(new.lm2,test='F',
     scope=~per.m+AGE0.14+AGE15.19+AGE.30+per.wht+per.AIAN+per.asian+
       per.oth+per.hl+Chlam.IR+per.smk.tol+perc.bechelor+perc.high+HI_18+HI18_64+
       percap.inc+perc.allinpov+perc.under18inpov+perc.urb)

#stepwise selection
library(MASS)
best.lm1=stepAIC(fit.lm)
best.lm2=stepAIC(basic.lm,test='F',
                 scope=~per.m+AGE0.14+AGE15.19+AGE.30+per.wht+per.AIAN+per.asian+
                   per.oth+per.hl+Chlam.IR+per.smk.tol+perc.bechelor+perc.high+HI_18+HI18_64+
                   percap.inc+perc.allinpov+perc.under18inpov+perc.urb)


##################
##Poisson model###
##################
fit.pos<-glm(w.can~per.m+AGE0.14+AGE15.19+AGE.30+per.wht+per.AIAN+per.asian+
              per.oth+per.hl+Chlam.IR+per.smk.tol+perc.bechelor+perc.high+HI_18+HI18_64+
              percap.inc+perc.allinpov+perc.under18inpov+perc.urb-1+offset(log(pop.tol)),family=poisson(link="log"),data=dat.tol)
summary(step(fit.pos))

fit.pos1=glm(w.can~offset(log(pop.tol))+per.m+AGE0.14+AGE15.19+AGE.30+
               per.hl+Chlam.IR+perc.bechelor+perc.high+
               percap.inc-1,family=poisson(link="log"),data=dat.tol)
summary(fit.pos1)

#by significance
fit.pos2=glm(w.can~per.m+AGE0.14+AGE15.19+AGE.30+
               per.hl+Chlam.IR+perc.bechelor+perc.high+
               percap.inc-1+offset(log(pop.tol)),family=poisson(link="log"),data=dat.tol)
summary(fit.pos2)

Anova(fit.pos1,type='II')

######################
##quasipoisson model##
######################
fit.qpos=glm(w.can~per.m+AGE0.14+AGE.30+per.wht+
               per.oth+per.hl+Chlam.IR+perc.bechelor+perc.high+HI18_64+
               percap.inc+perc.allinpov+perc.under18inpov+perc.urb-1+offset(log(pop.tol)),family=quasipoisson(link="log"),data=dat.tol)
summary(fit.qpos)
#summary(step(fit.qpos)) #AIC is not defined in this model-have to use hand 

#use negative binomial model for overdispersion
#####################
##Negative Binomial##
#####################
library(MASS)
library(lmtest)
fit.nb=glm.nb(w.cancer~offset(log(pop.tol))+per.m+AGE0.14+AGE.30+per.wht+
                per.oth+per.hl+perc.bechelor+perc.high+HI18_64+
                percap.inc+perc.allinpov+perc.urb,data=dat.tol)
summary(fit.nb)

######################################quadratic forms######################################

#################################
##Poisson model-Quadtratic form##
#################################
##original data##
##step()##
fit.pos.qf=glm(w.cancer~offset(log(pop.tol))+per.m+AGE0.14+AGE15.19+AGE.30+per.wht+per.AIAN+per.asian+
              per.oth+per.hl+Chlam.IR+per.smk.tol+perc.bechelor+perc.high+HI_18+HI18_64+
              percap.inc+perc.allinpov+perc.urb+
              I(per.m^2)+I(AGE0.14^2)+I(AGE15.19^2)+I(AGE.30^2)+I(per.wht^2)+I(per.AIAN^2)+I(per.asian^2)+
              I(per.oth^2)+I(per.hl^2)+I(Chlam.IR^2)+I(per.smk.tol^2)+I(perc.bechelor^2)+I(perc.high^2)+I(HI_18^2)+I(HI18_64^2)+
              I(percap.inc^2)+I(perc.allinpov^2)+I(perc.urb^2),family=poisson(link="log"),data=dat.cen)
summary(step(fit.pos.qf))

fit.pos.qf1=glm(w.cancer~offset(log(pop.tol))+per.m+AGE0.14+AGE.30+
                 per.oth+per.smk.tol+perc.bechelor+perc.high+HI18_64+
                 percap.inc+perc.allinpov+perc.urb+
                 I(perc.urb^2),family=poisson(link="log"),data=dat.cen)
summary(fit.pos.qf1)

#original data
fit.pos.qf2=glm(w.cancer~offset(log(pop.tol))+per.m+AGE0.14+AGE.30+
                 per.hl+per.smk.tol+perc.bechelor+
                 percap.inc+perc.allinpov+perc.urb+
                 I(per.hl^2)+I(per.smk.tol^2)+
                 I(perc.urb^2),family=poisson(link="log"),data=dat.tol)
summary(fit.pos.qf2)

######################################
##Quasipoisson model-Quadtratic form##
######################################
##centered data
fit.qpos.qf1=glm(w.cancer~offset(log(pop.tol))+per.m+AGE0.14+AGE.30+per.wht+
                 per.oth+per.hl+per.smk.tol+perc.bechelor+perc.high+HI_18+HI18_64+
                 percap.inc+perc.allinpov+perc.urb+
                 I(AGE.30^2)+
                 I(per.hl^2)+I(per.smk.tol^2)+I(HI_18^2)+
                 I(perc.urb^2),family=quasipoisson(link="log"),data=dat.cen)
summary(fit.qpos.qf1)

##original data
fit.qpos.qf2=glm(w.cancer~offset(log(pop.tol))+per.m+AGE0.14+AGE.30+
                   per.hl+per.smk.tol+perc.bechelor+perc.high+HI18_64+
                   percap.inc+perc.allinpov+perc.urb+
                   I(AGE.30^2)+
                   I(per.hl^2)+I(per.smk.tol^2)+I(perc.bechelor^2)+
                   I(perc.urb^2),family=quasipoisson(link="log"),data=dat.tol)
summary(fit.qpos.qf2)




