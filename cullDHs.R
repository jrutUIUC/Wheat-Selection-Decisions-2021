setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2021")
dhdata<- read.csv('DHphenotype_download.csv')

#combine with data from flatfile
setwd("~/Documents/Wheat/2021/Data Upload Files/check urbana dh data")
dhdata1<- read.csv('dhfield_urb_21-upload.csv', row.names=1)
colnames(dhdata1)[1]<- 'observationUnitName'
dhdata<- merge(dhdata, dhdata1, by='observationUnitName')

#add the row column information
setwd("~/Documents/Wheat/2021/HarvestMaps")
mirus<- read.csv('AugDHfield2021_mirusfile.csv', row.names=1)[,c(1,2,4)]
colnames(mirus)[3]<-'observationUnitName'
dhdata<- merge(dhdata, mirus, by='observationUnitName', all.x = TRUE, all.y=FALSE)
gids<- levels(dhdata$germplasmName)

#add the neoga data
setwd("~/Documents/Wheat/2021/Data Upload Files/check neoga yield data")
neoga<- read.csv('2021-06-22T113715phenotype_download.csv')
neoga2<-read.csv('Neoga2021_mirisfile_curated_bothdays.csv', row.names=1)
colnames(neoga2)[5]<- 'observationUnitName'
colnames(neoga2)[4]<- 'row'
neoga<- merge(neoga, neoga2, by='observationUnitName')
neogaDH<- neoga[grep('DH', neoga$observationUnitName),]
library(breedbase)
neogaDH$Test.Weight<- convert_lbsbu_gL(neogaDH$Test.Weight)
neogaDH$kgpha<- convert_buac_kgHa(neogaDH$buperacre, "wheat")

#combine data
commoncols<- intersect(colnames(dhdata), colnames(neogaDH))
neogaDH<- neogaDH[,commoncols]
dhdata<- dhdata[,commoncols]
dhdata<- droplevels.data.frame(rbind(dhdata, neogaDH))


#blocking factor
dhdata$blockNumber<- paste(dhdata$blockNumber, dhdata$studyName, sep="_")
dhdata$row<- paste(dhdata$row, dhdata$studyName, sep="_")
dhdata$range<- paste(dhdata$range, dhdata$studyName, sep="_")

#change to factors
dhdata$blockNumber<- as.factor(dhdata$blockNumber)
dhdata$studyName<- as.factor(dhdata$studyName)
dhdata$germplasmName<- as.factor(dhdata$germplasmName)
dhdata$entryType<- as.factor(dhdata$entryType)
dhdata$row<- as.factor(dhdata$row)
dhdata$range<- as.factor(dhdata$range)
dhdata$locationName<- as.factor(dhdata$locationName)

#make sure checks are recorded as checks
dhdata[which(dhdata$germplasmName=='Kaskaskia'),'entryType']<- 'check'
dhdata[which(dhdata$germplasmName=='07-4415'),'entryType']<- 'check'

#fit models, neoga only
library(asreml)
dhdataNeo<- dhdata[which(dhdata$locationName=='Neoga, IL'),]
dhdataNeo0<- dhdataNeo
traits<- c("Heading.time...Julian.date..JD..CO_321.0001233", "Plant.height...cm.CO_321.0001301","kgpha","Test.Weight")
for(i in 1:length(traits)){
  dhdataNeo<- dhdataNeo0
  colnames(dhdataNeo)[match(traits[i], colnames(dhdataNeo))]<- 'Y'
  asreml.options(extra=100, maxit=500)
  mod<- asreml(fixed= Y~1+at(entryType, 'check'):germplasmName,
               random= ~at(entryType, 'check'):blockNumber+at(entryType, 'test'):germplasmName+at(entryType, 'check'):row, data=dhdataNeo)
  mod<- update(mod)
  blups<- na.omit(predict(mod, classify='at(entryType, test):germplasmName', ignore=c('(Intercept)'))$pvals)
  pev<- blups[,'std.error']^2
  Vg<- summary(mod)$varcomp['at(entryType, test):germplasmName','component']
  rel<- 1-(pev/Vg)
  blups<- data.frame(blups, rel, trait=traits[i])
  if(i==1){
    blupsAll<- blups
  }else{
    blupsAll<- rbind(blupsAll, blups)
  }
  mod2<- asreml(fixed= Y~1+at(entryType, 'check'):germplasmName+at(entryType, 'test'):germplasmName,
               random= ~at(entryType, 'check'):blockNumber+at(entryType, 'check'):row,data=dhdataNeo)
  mod2<- update(mod2)
  blues<- na.omit(predict(mod2, classify='at(entryType, test):germplasmName')$pvals)
  blues<- data.frame(blues, trait=traits[i])
  if(i==1){
    bluesAll<- blues
  }else{
    bluesAll<- rbind(bluesAll, blues)
  }
}
blupsAll<- blupsAll[-c(grep('07-4415', blupsAll$germplasmName), grep('Kaskaskia', blupsAll$germplasmName)),]
unique(blupsAll[,c('trait', 'rel')])
bluesNeoga<- data.frame(loc='Neoga', bluesAll)

#fit models, urbana only
dhdataUrb<- droplevels.data.frame(dhdata[which(dhdata$locationName=='Urbana, IL'),])
dhdataUrb0<- dhdataUrb
library(asreml)
traits<- c("Heading.time...Julian.date..JD..CO_321.0001233", "Plant.height...cm.CO_321.0001301","kgpha","Test.Weight")
for(i in 1:length(traits)){
  dhdataUrb<- dhdataUrb0
  colnames(dhdataUrb)[match(traits[i], colnames(dhdataUrb))]<- 'Y'
  asreml.options(extra=100, maxit=500)
  mod<- asreml(fixed= Y~1+at(entryType, 'check'):germplasmName, 
             random= ~at(entryType, 'check'):blockNumber+
              at(entryType, 'test'):germplasmName+at(entryType, 'check'):row, data=dhdataUrb)
  mod<- update(mod)
  blups<- na.omit(predict(mod, classify='at(entryType, test):germplasmName',ignore=c('(Intercept)'))$pvals)
  pev<- blups[,'std.error']^2
  Vg<- summary(mod)$varcomp['at(entryType, test):germplasmName','component']
  rel<- 1-(pev/Vg)
  blups<- data.frame(blups, rel, trait=traits[i])
  if(i==1){
    blupsAll<- blups
  }else{
    blupsAll<- rbind(blupsAll, blups)
  }
  mod2<- asreml(fixed= Y~1+at(entryType, 'check'):germplasmName+at(entryType, 'test'):germplasmName,
                random= ~at(entryType, 'check'):blockNumber+at(entryType, 'check'):row, data=dhdataUrb)
  mod2<- update(mod2)
  blues<- na.omit(predict(mod2, classify='at(entryType, test):germplasmName')$pvals)
  blues<- data.frame(blues, trait=traits[i])
  if(i==1){
    bluesAll<- blues
  }else{
    bluesAll<- rbind(bluesAll, blues)
  }
}
blupsAll<- blupsAll[-c(grep('07-4415', blupsAll$germplasmName), grep('Kaskaskia', blupsAll$germplasmName)),]
unique(blupsAll[,c('trait', 'rel')])
bluesUrbana<- data.frame(loc='Urbana', bluesAll)

#all blues
dhmeans<- rbind(bluesNeoga, bluesUrbana)
dhmeans<- data.frame(dhmeans, wt= 1/(dhmeans$std.error^2))

for(i in 1:length(traits)){
  sub<- dhmeans[which(dhmeans$trait==traits[i]),]
  modME<- asreml(fixed=predicted.value~1+loc, random= ~germplasmName, weights=wt,family = asr_gaussian(dispersion = 1),
                 data=sub)
  blups<- predict(modME, classify='germplasmName',ignore=c('(Intercept)', 'loc'))$pvals
  pev<- blups[,'std.error']^2
  Vg<- summary(modME)$varcomp['germplasmName','component']
  rel<- 1-(pev/Vg)
  blups<- data.frame(blups, trait=traits[i], rel)
  if(i==1){
    blupsAll<- blups
  }else{
    blupsAll<- rbind(blupsAll, blups)
  }
}
wide<- cast(blupsAll, germplasmName~trait, value='predicted.value')

##Select based on net merit
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
  soy_yld_gain<- 0.5* (max(headings)-headings)
  soy_profit_gain<- soy_yld_gain*soybean_price
  wheat_profit<- yields*wheat_price1
  total_profit<- wheat_profit + soy_profit_gain
  return(total_profit)
}

wide$Test.Weight_imperial<- wide$Test.Weight/convert_lbsbu_gL(1)
wide$Yield_imperial<- wide$kgpha/convert_buac_kgHa(1, "wheat")
convert_buac_kgHa(wide$Yield_imperial, "wheat")[1:10]
wide$Test.Weight_imperial<- wide$Test.Weight_imperial+58
wide$Yield_imperial<- wide$Yield_imperial+80
wide$Heading.time...Julian.date..JD..CO_321.0001233<- wide$Heading.time...Julian.date..JD..CO_321.0001233+136

nets<- c()
for(i in 1:nrow(wide)){
  net<- netMerit(wide[i,'Heading.time...Julian.date..JD..CO_321.0001233'], wide[i,'Yield_imperial'], 0,
                 0, wide[i,'Test.Weight_imperial'],wheat_price0, soybean_price)
  nets<- append(nets,  net)
}
wide<- data.frame(wide, nets)
wide<- wide[-which(wide$germplasmName %in% c('07-4415', 'Kaskaskia')),]
notes<- dhdata[,c('notes','germplasmName','locationName')]
notes<- notes[-which(notes$germplasmName %in% c('07-4415', 'Kaskaskia')),]
notes<- cast(notes, germplasmName~locationName, value='notes')
notes<- notes[,c(1,3)]
colnames(notes)<- c("germplasmName", "notes")

wide<- merge(wide, notes, by='germplasmName')
wide$decision<- 'select'
wide[grep('exclude_plot:true', wide$notes),'decision']<- 'discard'
wide<- wide[order(-wide$nets),]
wide[122:541, 'decision']<- 'discard'
wide[which(wide$Heading.time...Julian.date..JD..CO_321.0001233>139),'decision']<- 'discard' 
wide[which(wide$Plant.height...cm.CO_321.0001301>10),'decision']<- 'discard' 
table(wide$decision)

#make file for inventory
selectedDH<- as.character(wide[which(wide$decision =='select'),'germplasmName'])
dhdata$decision<- 'discard'
dhdata[match(selectedDH, dhdata$germplasmName),'decision']<- 'select'
inventoryfile<- dhdata[,c('observationUnitDbId', 'studyDbId','plotNumber', 'germplasmName','observationUnitName', 'decision')]

#combine with aug_urb file for inventory
setwd("~/Documents/Wheat/2021/HarvestMaps")
meta<- read.csv('AugUrbphenotypedownload.csv')
aug<- read.csv('AugUrb2021_mirusfile.csv', row.names=1)
#colnames(aug)[4]<- 'observationUnitName'
aug<- merge(meta, aug, by='observationUnitName')
aug<- aug[,c('observationUnitDbId','studyDbId','plotNumber.x', 'germplasmName.x', 'observationUnitName', 'bag')]
augInventory<- aug[which(aug$bag=='yes'),]
colnames(augInventory)<- colnames(inventoryfile)
augInventory$decision='select'

setwd("~/Documents/Wheat/2021/Seed inventory and cleaning")
inventoryfile<- rbind(inventoryfile,augInventory)
inventoryfile<- inventoryfile[-grep('Neo', inventoryfile$observationUnitName),]
inventoryfile<- inventoryfile[order(as.numeric(as.character(inventoryfile$plotNumber))),]
inventoryfile<- inventoryfile[order(as.numeric(as.character(inventoryfile$studyDbId))),]

write.csv(inventoryfile, file='Stg2inventory_fb.csv', row.names=FALSE)
setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2021")