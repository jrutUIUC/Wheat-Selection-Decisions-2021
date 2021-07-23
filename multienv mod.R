setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2021")
load('TrainingSetJul13.2021.RData')

#Model for yield,test weight,height
mod<- asreml(fixed=predicted.value~1, random=~loc+year+site+study+germplasmName+germplasmName:site, 
             weights=wt, data=Yield, na.action = na.method(y='omit', x='omit'), family=asr_gaussian(dispersion = 1), workspace=64e6)
blups<- predict(mod, classify='germplasmName', ignore=c('(Intercept)'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(mod)$varcomp['germplasmName','component']
rel<- 1-(pev/Vg) 
mean(rel) #0.6133151

#Model for other traits
mod<- asreml(fixed=predicted.value~1, random=~year+study+germplasmName+germplasmName:year, 
             weights=wt, data=Sev, na.action = na.method(y='omit', x='omit'), family=asr_gaussian(dispersion = 1), workspace=64e6)

blups<- predict(mod, classify='germplasmName', ignore=c('(Intercept)'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(mod)$varcomp['germplasmName','component']
rel<- 1-(pev/Vg) 
mean(rel) #0.4229114
