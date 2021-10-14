setwd("/Users/jrut/Documents/GitHub/Wheat-Selection-Decisions-2021")
library(ASExtras4)
library(reshape)
data2020to21<- All[which(as.numeric(as.character(All$year))>19),]
#exclude<- c(intersect(which("Pr[0-9]_Car_19"==as.character(data2019to21$study)), which("Aliased"==as.character(data2019to21$status))), intersect(which("Pr[0-9]_Stj_19"==as.character(data2019to21$study)), which("Aliased"==as.character(data2019to21$status))))
#data2019to21<- data2019to21[-exclude,]

#make converge
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >2, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

#yield only model
yld2020to21<- data2020to21[which(data2020to21$trait=='Yield'),]
famod<- asreml(fixed=predicted.value~1, random=~study+fa(site, 2):germplasmName, 
               family=asr_gaussian(dispersion = 1), na.action = na.method(y='omit', x='omit'),
               workspace=64e6, data=yld2020to21, weights=wt)
famod<- mkConv(famod)
smry<- fa.asreml(famod)
smry
cormat<- cov2cor(smry$gammas$`fa(site, 2):germplasmName`$Gmat)
p<- predict(famod, classify= 'site:germplasmName')$pvals

#summarize yield results
yldSmry<- cast(p, germplasmName~site, value='predicted.value')
yldSouth<- rowMeans(yldSmry[,c("Car_20", "Neo_20","Neo_21","Stj_20", "Stj_21")])
yldNorth<- rowMeans(yldSmry[,c("Urb_20", "Urb_21")])
yldSmry<- cbind(yldSmry, yldSouth, yldNorth)
locwts<- c(0.16, 0.16, 0.16, 0.16, 0.16, 0.1, 0.1)
yldIndex<- as.matrix(yldSmry[,c("Car_20", "Neo_20","Neo_21","Stj_20", "Stj_21","Urb_20", "Urb_21")]) %*%  matrix(locwts)
yldSmry<- cbind(yldSmry, yldIndex)

## test weight only model
#yield only model
tw2020to21<- data2020to21[which(data2020to21$trait=='TestWeight'),]
mod<- asreml(fixed=predicted.value~1, random=~study+site+ germplasmName + site:germplasmName, 
               family=asr_gaussian(dispersion = 1), na.action = na.method(y='omit', x='omit'),
               workspace=64e6, data=tw2020to21, weights=wt)
p<- predict(mod, classify= 'germplasmName')$pvals
twtSmry<- p
twtSmry<- twtSmry[,c(1:2)]
colnames(twtSmry)[2]<- 'TestWeight'

## days to heading only model
dth2020to21<- data2020to21[which(data2020to21$trait=='Heading'),]
mod<- asreml(fixed=predicted.value~1, random=~study+site+ germplasmName + site:germplasmName, 
             family=asr_gaussian(dispersion = 1), na.action = na.method(y='omit', x='omit'),
             workspace=64e6, data=dth2020to21, weights=wt)
p<- predict(mod, classify= 'germplasmName')$pvals
dthSmry<- p
dthSmry<- dthSmry[,c(1:2)]
colnames(dthSmry)[2]<- 'Heading'

## height only model
ht2020to21<- data2020to21[which(data2020to21$trait=='Height'),]
mod<- asreml(fixed=predicted.value~1, random=~study+site+ germplasmName + site:germplasmName, 
             family=asr_gaussian(dispersion = 1), na.action = na.method(y='omit', x='omit'),
             workspace=64e6, data=ht2020to21, weights=wt)
p<- predict(mod, classify= 'germplasmName')$pvals
htSmry<- p
htSmry<- htSmry[,c(1:2)]
colnames(htSmry)[2]<- 'Height'


## scab model ## wont fit in two stage

#single stage scab model
library(plyr)
data<- read.csv('All scb data 2021.csv', as.is=TRUE)
dataB<- read.csv('All trials 2015-2020.csv', as.is=TRUE)
data<-rbind.fill(data, dataB)
data<- data[which(as.numeric(as.character(data$studyYear))>2019),]
data<- droplevels.data.frame(data[c(grep('Pr', data$studyName), grep('Adv', data$studyName)),])

#Fit multivariate model for blups
data$studyName<- as.factor(data$studyName)
data$blockNumber<- as.factor(data$blockNumber)
data$germplasmName<- as.factor(data$germplasmName)
data$rowNumber<- as.factor(data$rowNumber)
data$colNumber<- as.factor(data$colNumber)
amod1<- asreml(fixed=cbind(FHB.DON.content...ppm.CO_321.0001154, FHB.grain.incidence.....CO_321.0001155, FHB.incidence.....CO_321.0001149,
                           FHB.severity.....CO_321.0001440, Heading.time...Julian.date..JD..CO_321.0001233)~1+trait,
               random= ~at(trait):studyName+at(trait):studyName:blockNumber+us(trait):germplasmName,
               residual = ~id(units):us(trait),
               data=data,  na.action = na.method(y='include', x='include'))
p<- predict(amod1, classify= 'trait:germplasmName')$pvals
scbSmry<- cast(p, germplasmName~trait, value='predicted.value')[,c(1:3)]
colnames(scbSmry)<- c('germplasmName', 'DON', 'FDK')

#merge all data
all<- merge(yldSmry, twtSmry, by='germplasmName', all=TRUE)
all<- merge(all, dthSmry, all=TRUE)
all<- merge(all, htSmry, all=TRUE)
all<- merge(all, scbSmry)

#calculte net merit
#Net merit function
#starting price of wheat and soybean, five year average based on macrotrends.net
wheat_price0<- mean(c(5.4621, 4.9414, 4.9757, 4.4014, 4.3945))
soybean_price<- mean(c(9.3785, 8.9298, 9.3456, 9.7820, 9.8753))
#wheat price fcn
wheatPrice<- function(fdk, don, twt, wheat_price0){
  if(don==0){
    donDiscount<- 0
  }else{
    donDiscount<- sqrt(don)*-0.2
  }
  if(fdk==0){
    fdkDiscount<- 0
  }else{
    fdkDiscount<- sqrt(fdk)*-0.04
  }
  twtDiscount<- c(58-twt)*-.2
  twtDiscount[which(twtDiscount>0)]<- 0
  wheat_price<- wheat_price0+donDiscount+fdkDiscount+twtDiscount
  return(wheat_price)
}
#net merit function
netMerit<- function(headings, yields, dons, fdks, twt, wheat_price0, soybean_price){
  wheat_price1<- wheatPrice(fdks, dons, twt, wheat_price0)
  soy_yld_gain<- 0.5* (138.4524-headings)
  soy_profit_gain<- soy_yld_gain*soybean_price
  wheat_profit<- yields*wheat_price1
  total_profit<- wheat_profit + soy_profit_gain
  return(total_profit)
}

nets<- c()
for(i in 1:nrow(all)){
  net<- netMerit(headings=all[i,'Heading'], yields=all[i,'yldIndex'], dons=all[i,'DON'],
                 fdks=all[i,'FDK'], twt=all[i,'TestWeight'],wheat_price0, soybean_price)
  nets<- append(nets,  net)
}
all<- cbind(all, nets)

#add the testing cohort
trials2021<- All[which(as.numeric(as.character(All$year))==21),]
trials2020<- All[which(as.numeric(as.character(All$year))==20),]
Pr2021<- all$germplasmName %in% trials2021[grep('Pr', trials2021$study),'germplasmName']
Adv2021<- all$germplasmName %in% trials2021[grep('Adv', trials2021$study),'germplasmName']
Adv2020<- all$germplasmName %in% trials2020[grep('Adv', trials2020$study),'germplasmName']
Cohort<- rep(NA, nrow(all))
Cohort[which(Adv2020=='TRUE' & Adv2021=='TRUE')]<- 'Advanced 2'
Cohort[which(Adv2020=='FALSE' & Adv2021=='TRUE')]<- 'Advanced 1'
Cohort[which(Pr2021=='TRUE')]<- 'Preliminary'
all<- cbind(all, Cohort)
for(i in 2:17){
  all[,i]<- round(all[,i], 2)
}
Preliminary<- all[which(all$Cohort=='Preliminary'),]
Adv1<- all[which(all$Cohort=='Advanced 1'),]
Adv2<- all[which(all$Cohort=='Advanced 2'),]


excellist<- list(combined=all, Pr=Preliminary, Adv1=Adv1, Adv2=Adv2)
WriteXLS::WriteXLS(excellist, ExcelFileName=paste('SelectionFile2021', ".xls", sep=""), 
                   SheetNames=c('combined', 'Pr', 'Adv1', 'Adv2'))

#cor(all[,-1], use='pairwise.complete.obs')[,'nets']
