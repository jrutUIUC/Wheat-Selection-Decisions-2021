setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2021")
a<- read.csv('Adv_Neo_21blues.csv', row.names=1)
b<- read.csv('Adv_Scb_21blues.csv', row.names=1)
c<- read.csv('Adv_Stj_21blues.csv', row.names=1)
d<- read.csv('Adv_Urb_21blues.csv', row.names=1)
a2<- read.csv('Pr_Neo_21blues.csv', row.names=1)
b2<- read.csv('Pr_Scb_21blues.csv', row.names=1)
c2<- read.csv('Pr_Stj_21blues.csv', row.names=1)
d2<- read.csv('Pr_Urb_21blues.csv', row.names=1)
e<- read.csv('predicted values table.csv', row.names=1)
e<- e[,c('studyName', 'name', 'trait','predicted.value','std.error','status')]
colnames(e)<- colnames(d)
all<- rbind(a, b, c, d, a2, b2, c2, d2,  e)

#make converge
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >2, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

#subset the studies
std<- unique(all$study)
std2<- c(as.character(std[grep('Pr', std)]), as.character(std[grep('Adv', std)]))
std3<- std2[-grep('Sbmv', std2)]
all2<- all[which(all$study %in% std3),]

#subset the traits
traitsel<- c('Grain.yield...bu.ac','Heading.time...Julian.date..JD..CO_321.0001233',
             'Grain.test.weight...lbs.bu', 'Plant.height.inches',
             'FHB.incidence.....CO_321.0001149', 'FHB.severity.....CO_321.0001440',
             'FHB.grain.incidence.....', 'FHB.incidence.....','FHB.severity.....',
             'FHB.DON.content...ppm.', 'Heading.time...Julian.date..JD..')
traitsel2<- c('Yield','Heading',
             'TestWeight', 'Height',
             'Inc', 'Sev',
             'Fdk', 'Inc','Sev',
             'DON', 'Heading')
all3<- all2[which(all2$trait %in% traitsel),]

#rename the traits
all3$trait<- as.character(all3$trait)
for(i in 1:length(traitsel)){
  all3[which(all3$trait==traitsel[i]),'trait']<- traitsel2[i] 
}

#add weight to the file
var<- all3$std.error^2
wt<- 1/var
all3<- data.frame(all3, wt)

#add more factors
loc<- matrix(unlist(strsplit(as.character(all3$study), split="_")), nrow=3)[2,]
loc[which(loc=='StJ')]<- 'Stj'
year<- matrix(unlist(strsplit(as.character(all3$study), split="_")), nrow=3)[3,]
site<- paste(loc, year, sep="_")
all3<- data.frame(site, loc, year, all3)

#individual trait data only
All<- all3
Yield<- all3[which(all3$trait=='Yield'),]
Heading<- all3[which(all3$trait=='Heading'),]
TestWeight<- all3[which(all3$trait=='TestWeight'),]
Height<- all3[which(all3$trait=='Height'),]
Sev<- all3[which(all3$trait=='Sev'),]
Inc<- all3[which(all3$trait=='Inc'),]
Fdk<- all3[which(all3$trait=='Fdk'),]
DON<- all3[which(all3$trait=='DON'),]
save(All, Yield, Heading, TestWeight, Height, Sev, Inc, Fdk, DON, file='TrainingSetJul13.2021.RData')

#fitting a model
mod<- asreml(fixed=predicted.value~1, random=~loc+year+site+study+germplasmName+germplasmName:site, 
             weights=wt, data=Yield, na.action = na.method(y='omit', x='omit'), family=asr_gaussian(dispersion = 1), workspace=64e6)
blups<- predict(mod, classify='germplasmName', ignore=c('(Intercept)'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(mod)$varcomp['germplasmName','component']
rel<- 1-(pev/Vg) 
mean(rel) #0.6133151

#fitting a model SCB
mod<- asreml(fixed=predicted.value~1, random=~site+study+germplasmName+germplasmName:site, 
             weights=wt, data=Sev, na.action = na.method(y='omit', x='omit'), family=asr_gaussian(dispersion = 1), workspace=64e6)
mod<- mkConv(mod)
blups<- predict(mod, classify='germplasmName', ignore=c('(Intercept)'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(mod)$varcomp['germplasmName','component']
rel<- 1-(pev/Vg) 
mean(rel) #0.4229114

#fitting a model M

mod<- asreml(fixed=predicted.value~1+trait, random=~at(trait):loc+at(trait):year+at(trait):site+at(trait):study+at(trait):germplasmName+at(trait):germplasmName:site, 
             weights=wt, asmv='trait', data=yldtw, 
             residual = ~id(units):us(trait), na.action = na.method(y='include', x='include'), family=asr_gaussian(dispersion = 1), workspace=64e6)
mod<- mkConv(mod)
blups<- predict(mod, classify='germplasmName', ignore=c('(Intercept)'))$pvals
pev<- blups[,'std.error']^2
Vg<- summary(mod)$varcomp['germplasmName','component']
rel<- 1-(pev/Vg) 
mean(rel) #0.6133151




usite<- unique(all3$site)
for(i in 1:length(usite)){
  sub<- all3[which(all3$site == usite[i]),]
  
  #identify checks
  wide<- cast(sub, germplasmName~study)
  wide0<- wide[,-1]
  wide0<-apply(wide0, 2, function(x){
    if(0 %in% x)
    x[which(x==0)]<- NA
    return(x)
  })
  wide<- data.frame(germplasmName=wide[,1], wide0)
  checks<- as.character(na.omit(wide)$germplasmName)
  sub$type<- 'test'
  sub[sub$germplasmName %in% checks,'type']<- 'check'
  sub$type<- as.factor(sub$type)
  
  sub2<- sub[which(sub$trait=='Yield'),]
  mod<- asreml(fixed=predicted.value~1+at(type, 'check'):study, random=~germplasmName, 
               weights=wt, data=sub2, family=asr_gaussian(dispersion = 1), workspace=64e6)
  coefficients(mod)$fixed
  wald(mod)
  
}

