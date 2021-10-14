setwd("/Users/jrut/Documents/GitHub/Wheat-Selection-Decisions-2021")
library(asreml)
library(ASExtras4)
library(reshape)
load('TrainingSetJul24.2021.RData')
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

#########Get across env mean for traits
#days to heading only model
dth2020to21<- data2020to21[which(data2020to21$trait=='Heading'),]
mod<- asreml(fixed=predicted.value~1+ germplasmName, random=~study+site + site:germplasmName, 
             family=asr_gaussian(dispersion = 1), na.action = na.method(y='omit', x='omit'),
             workspace=64e6, data=dth2020to21, weights=wt)
p<- predict(mod, classify= 'germplasmName')$pvals
dthSmry1<- p
dthSmry<- dthSmry1[,c(1:2)]
colnames(dthSmry)[2]<- 'Heading'

#test weight only model
tw2020to21<- data2020to21[which(data2020to21$trait=='TestWeight'),]
mod<- asreml(fixed=predicted.value~1+ germplasmName, random=~study+site + site:germplasmName, 
             family=asr_gaussian(dispersion = 1), na.action = na.method(y='omit', x='omit'),
             workspace=64e6, data=tw2020to21, weights=wt)
p<- predict(mod, classify= 'germplasmName')$pvals
twSmry1<- p
twSmry<- twSmry1[,c(1:2)]
colnames(twSmry)[2]<- 'TestWeight'

## height only model
ht2020to21<- data2020to21[which(data2020to21$trait=='Height'),]
mod<- asreml(fixed=predicted.value~1+ germplasmName, random=~study+site + site:germplasmName, 
             family=asr_gaussian(dispersion = 1), na.action = na.method(y='omit', x='omit'),
             workspace=64e6, data=ht2020to21, weights=wt)
p<- predict(mod, classify= 'germplasmName')$pvals
htSmry1<- p
htSmry<- htSmry1[,c(1:2)]
colnames(htSmry)[2]<- 'Height'


#yield across environments model
yld2020to21<- data2020to21[which(data2020to21$trait=='Yield'),]
mod<- asreml(fixed=predicted.value~1+germplasmName, random=~study+site + site:germplasmName, 
               family=asr_gaussian(dispersion = 1), na.action = na.method(y='omit', x='omit'),
               workspace=64e6, data=yld2020to21, weights=wt)
mod<- mkConv(mod)
p<- predict(mod, classify= 'germplasmName')$pvals
yldSmry1<- p
yldSmry<- yldSmry1[,c(1:2)]
colnames(yldSmry)[2]<- 'Yield'

#merge traits
yldSmry1<- data.frame(trait='Yield', yldSmry1)
dthSmry1<- data.frame(trait='Heading', dthSmry1)
htSmry1<- data.frame(trait='Height', htSmry1)
twSmry1<- data.frame(trait='TestWeight', twSmry1)
allTrt<- rbind(yldSmry1, dthSmry1, htSmry1, twSmry1)
yldSmry1$predicted.value<- yldSmry1$predicted.value- mean(yldSmry1$predicted.value)
dthSmry1$predicted.value<- dthSmry1$predicted.value- mean(dthSmry1$predicted.value)
twSmry1$predicted.value<- twSmry1$predicted.value- mean(twSmry1$predicted.value)
allTrt$wt<- 1/(allTrt$std.error^2)
allTrt<- allTrt[-which(allTrt$germplasmName== 'Pio 25R47'),]

#relationship matrix
load("/Users/jrut/Documents/Wheat/2022/Crossing/K_19_20_21.RData")

name2vec<-c()
for(i in 1:nrow(allTrt)){
  a<- match(as.character(allTrt[i,'germplasmName']), row.names(K_19_20_21))
  b<- match(as.character(allTrt[i,'germplasmName']), gsub("IL", "", row.names(K_19_20_21)))
  row<- na.omit(c(a, b))
  if(length(row)>0){
    name2<- row.names(K_19_20_21)[row] 
  }else{
    name2=NA
  }
  name2vec<- append(name2vec, name2)
}
allTrt<- data.frame(allTrt, name2vec, stringsAsFactors = FALSE)
allTrt[which(allTrt$germplasmName=='Kaskaskia'),'name2vec']<- 'KASKASKIA'
allTrt[which(allTrt$germplasmName=='Pio25R74'),'name2vec']<- 'PIO-25R74'
allTrt[which(allTrt$germplasmName=='Pio 25R74'),'name2vec']<- 'PIO-25R74'
allTrt<- droplevels.data.frame(allTrt[-which(is.na(allTrt$name2vec)),])
key<- unique(allTrt[,c('germplasmName', 'name2vec')])
genoSub<- geno_df_19_20_21[which(row.names(geno_df_19_20_21) %in% allTrt$name2vec),]
library(rrBLUP)
Gmat<- A.mat(genoSub)
diag(Gmat)<- diag(Gmat) + 0.00002
allTrt$name2vec<- as.factor(allTrt$name2vec)
cor(cast(allTrt, germplasmName~trait, value='predicted.value')[,-1])

##Use the previous stage data
name2vec<- key[match(data2020to21$germplasmName, key$germplasmName),'name2vec']
data2020to21<- data.frame(data2020to21, name2= name2vec)
data2020to21<- data2020to21[-which(is.na(data2020to21$name2)),]
data2020to21<- data2020to21[-which(data2020to21$loc=='Scb'),]
data2020to21<- droplevels.data.frame(data2020to21)
data2020to21$trait<- as.factor(data2020to21$trait)
data2020to21<- data2020to21[order(data2020to21$trait),]
data2020to21<- data2020to21[-which(is.na(data2020to21$wt)),]

mod<- asreml(fixed=predicted.value~trait, random= ~study:trait+site:trait+ us(trait):vm(name2, Gmat), 
             family=asr_gaussian(dispersion = 1), residual= ~id(units):id(trait), 
             na.action = na.method(y='include', x='include'),
             workspace=64e6, asmv=trait, data=data2020to21, weights=wt)
mod<- update(mod)
pAg<- predict(mod, classify= 'name2:trait', ignore=c('site'), pworkspace= 64e6)$pvals
save(pAg, file='AgronomicTraitGBLUPS.RData')

## scab model ## wont fit in two stage

#single stage scab model
library(plyr)
data<- read.csv('All scb data 2021.csv', as.is=TRUE)
dataB<- read.csv('All trials 2015-2020.csv', as.is=TRUE)
data<-rbind.fill(data, dataB)
data<- data[which(as.numeric(as.character(data$studyYear))>2019),]
data<- droplevels.data.frame(data[c(grep('Pr', data$studyName), grep('Adv', data$studyName)),])

##match up data with the marker data
name2<- key[match(data$germplasmName, key$germplasmName),'name2vec']
data<- data.frame(data, name2)
data<- data[-which(is.na(data$name2)),]

#Fit multivariate model for blups
data$studyName<- as.factor(data$studyName)
data$blockNumber<- as.factor(data$blockNumber)
data$germplasmName<- as.factor(data$germplasmName)
data$rowNumber<- as.factor(data$rowNumber)
data$colNumber<- as.factor(data$colNumber)
amod1<- asreml(fixed=cbind(FHB.DON.content...ppm.CO_321.0001154, FHB.grain.incidence.....CO_321.0001155)~trait,
               random= ~trait:studyName+trait:studyName:blockNumber+us(trait):vm(name2, Gmat),
               residual = ~id(units):us(trait),workspace=64e6, 
               data=data,  na.action = na.method(y='include', x='include'))
p<- predict(amod1, classify= 'trait:name2',pworkspace=64e6)$pvals
save(p, file='ScabTraitGBLUPS.RData')

#### Calculate net merit of all combinations in the crossing block
setwd("/Users/jrut/Documents/Wheat/2022/crossing")
prnts<- read.csv('2021_ParentsCross_list.xlsx - Sheet1.csv')
genoPrnts<- key[match(prnts$Name., key$germplasmName),'name2vec']
df<- data.frame(Numb= c(1:44), name2=  genoPrnts)

#make data tables
ag<- cast(pAg, name2~trait, value='predicted.value')
scb<- cast(p, name2~trait, value='predicted.value')
df2<- merge(df, ag, by='name2', all.x=TRUE, all.y = FALSE)
df3<- merge(df2, scb, by='name2', all.x=TRUE, all.y = FALSE)
colnames(df3)[c(7,8)]<- c('DON', 'FDK')
df4<- na.omit(df3)
#get net merit of all combinations
cn<- t(combn(df4$Numb, 2)) 
#cnr<- cn[,c(2,1)]
#cns<- rbind(cn, cnr)

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
    donDiscount<- don*-0.2
  }
  if(fdk==0){
    fdkDiscount<- 0
  }else{
    fdkDiscount<- fdk*(-0.04/5)
  }
  twtDiscount<- c(58-twt)*-.2
  twtDiscount[which(twtDiscount>0)]<- 0
  wheat_price<- wheat_price0+donDiscount+fdkDiscount+twtDiscount
  return(wheat_price)
}

#net merit function
netMerit<- function(headings, yields, dons, fdks, twt, wheat_price0, soybean_price){
  wheat_price1<- wheatPrice(fdks, dons, twt, wheat_price0)
  soy_yld_gain<- 0.5* (137.0577-headings)
  soy_profit_gain<- soy_yld_gain*soybean_price
  wheat_profit<- yields*wheat_price1
  total_profit<- wheat_profit + soy_profit_gain
  return(total_profit)
}

for(i in 1:nrow(cn)){
  numbs<- cn[i, ]
  pair<- df3[match(numbs, df3$Numb),]
  F1trt<- pair[,c('Heading', 'Height', 'TestWeight', 'Yield', 'DON', 'FDK')]
  F10<- colMeans(F1trt)
  F1<- data.frame(F10)
  merit<- round(netMerit(F1['Heading',], F1['Yield',], F1['DON',], F1['FDK',], F1['TestWeight',], wheat_price0, soybean_price),2)
  crossA= paste(pair$name2, collapse= ' x ')
  crossB= paste(pair$name2[2:1], collapse= ' x ')
  crossA1= paste(pair$Numb, collapse= ' x ')
  crossB1= paste(pair$Numb[2:1], collapse= ' x ')
  row<- c(crossA, crossB, crossA1, crossB1, F10, merit)
  names(row)[1:4]<- c('cross', 'reciprocal', 'cross_numb', 'reciprocal_numb')
  names(row)[11]<- c('merit')
  if(i==1){
    rows<- row
  }else{
    rows<- rbind(rows, row)
  }
}
crossRank<- as.data.frame(rows)
row.names(crossRank)<- as.character(c(1:nrow(crossRank)))
ranking<- rank(-as.numeric(as.character(crossRank$merit)))
crossRank<- data.frame(crossRank, ranking)
crossRank2<- crossRank
crossRank2$cross<- crossRank$reciprocal
crossRank2$reciprocal<- crossRank$cross
crossRank2$cross_numb<- crossRank$reciprocal_numb
crossRank2$reciprocal_numb<- crossRank$cross_numb
crossRankFinal<- rbind(crossRank, crossRank2)
crossRankFinal$merit<- as.numeric(as.character(crossRankFinal$merit))
relativeProfit<- round(crossRankFinal$merit- mean(crossRankFinal$merit), 2)
crossRankFinal$relativeProfit<- relativeProfit

for(i in 5:10){
  crossRankFinal[,i]<- round(as.numeric(as.character(crossRankFinal[,i])),2)
}

p1<- matrix(unlist(strsplit(as.character(crossRankFinal$cross_numb), split= " x ")), nrow=2)[1,]
p2<- matrix(unlist(strsplit(as.character(crossRankFinal$cross_numb), split= " x ")), nrow=2)[2,]
crossRankFinal$p1<- p1
crossRankFinal$p2<- p2
crossRankFinal<- crossRankFinal[,c(14:15, 1:13)]

write.csv(crossRankFinal, file='crossRankings2021_22.csv')
crossRank21_22<- crossRankFinal
save(crossRank21_22, file='crossRank21_22.RData')
