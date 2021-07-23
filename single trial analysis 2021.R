setwd("~/Documents/GitHub/Wheat-Selection-Decisions-2021")
library(asreml)
library(asremlPlus)
library(reshape)
data<- read.csv('TrialsForAnalysis_July12.2021.csv', as.is=TRUE)

#suppress bad yield trial plots on the Maxwell field
pltex<- read.csv('maxwell_excludedplots.csv')
data[match(pltex[,1], data$observationUnitName),'Grain.yield...kg.ha.CO_321.0001218']<- NA
data[match(pltex[,1], data$observationUnitName),'Grain.test.weight...g.l.CO_321.0001210']<- NA
data[match(pltex[,1], data$observationUnitName),'Grain.moisture.content.....CO_321.0001198']<- NA
row.names(data)<- data$observationUnitName

#suppress suspicious values
data[match('Pr_Urb_21-plot902', data$observationUnitName),'Grain.test.weight...g.l.CO_321.0001210']<- NA
data[match('Adv_Urb_21-plot454', data$observationUnitName),'Grain.test.weight...g.l.CO_321.0001210']<- NA

#exclude the moisture data from analysis
data<- data[,-c(33)]

#make vector of ids that were submitted to cooperatives
coopgids<-c('07-19334','16-23941','16-23972','16-36752',
'16-36206','16-8605','16-22039','13-1960','15-2639',
'16-8048','17-17739','07-4415', '02-18228',
'07-19334','17-23904','17-25205','17-8626','15-4957','17-23874','17-29544')



############################
##Functions to be used later
############################

#mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

#function to check model convergence and update until converged (tolerate a 1.5% change in components)
mkConv<- function(mod){
  pctchg<- summary(mod)$varcomp[,c('%ch')]
  while(any(pctchg >2, na.rm=TRUE)){
    mod<-suppressWarnings(update(mod))
    pctchg<- summary(mod)$varcomp[,c('%ch')]
  }
  return(mod)
}

#create vector of traits recorded in the trial
selectTraitcols<- function(trl, trtnms, thresh=0.2){
  ppres<- c()
  for(j in 1:length(trtnms)){
    vec<-trl[,trtnms[j]]
    ppres<- append(ppres, length(na.omit(vec))/length(vec))
  }
  return(trtnms[which(ppres>thresh)])
}

#function to add means, MSE, LSD, and CV
addRows<- function(df, varNm, label, vec){
  lenvec<- length(vec)
  df[nrow(df)+1,]<-df[nrow(df),]
  df[nrow(df), c(1:varNm)]<- rep("", varNm)
  df[nrow(df), varNm]<- label
  df[nrow(df), -c(1:varNm)]<- round(vec, 5)
  return(df)
}

#convert yld to bu/acre
convYld<- function(y){
  x<- y/c(60 * 0.453592 * 2.47105)
  return(x)
}

#convert tw to lbs/bu
convTwt<- function(y){
  x<- y/1000 *2.2046 *35.2391
  return(x)
}  

############################
## Data curation
############################
#convert yield and test weight to common units
data[,'Grain.yield...kg.ha.CO_321.0001218']<- convYld(data[,'Grain.yield...kg.ha.CO_321.0001218'])
data[,'Grain.test.weight...g.l.CO_321.0001210']<- convTwt(data[,'Grain.test.weight...g.l.CO_321.0001210'])
data[,"Plant.height...cm.CO_321.0001301"]<- data[,"Plant.height...cm.CO_321.0001301"] *0.393701
colnames(data)<- gsub('kg.ha.CO_321.0001218', "bu.ac", colnames(data))
colnames(data)<- gsub('g.l.CO_321.0001210', "lbs.bu", colnames(data))
colnames(data)<- gsub('..cm.CO_321.0001301', "inches", colnames(data))

#shorten the trait names to avoid errors in model fitting
ixCOs<- grep('CO_', colnames(data))
colnames(data)[ixCOs]<- matrix(unlist(strsplit(colnames(data)[ixCOs], split="CO_")), nrow=2)[1,]

############################
##Subset each trial
############################

#get unique study names
stdnms<- unique(data$studyName)

#indicate which require a single trial analysis summary file
stdnmsCoop<- stdnms[grep(c("NU|PNU|A5S|P5S|A6S|P6S|SU|UE|VT|CC|NN"), stdnms)]

#empty vector trial names with no data will be added
nodata<-c()

#empty vector of trait and trial combinations with only one replication
trtnonrep<-c()

#loop begins
for(i in 1:length(stdnms)){
  cat(i, '\n')
  #subset single trial
  trl<- data[data$studyName==stdnms[i],]
  trl<- trl[order(trl$plotNumber), ]
  ############################
  ##Extract design information
  ############################
  
  #get vector of all possible traits
  trtnms<- colnames(data)[31:36]
  #vector of traits used in the selected trial
  ttrt<- selectTraitcols(trl, trtnms)
  
  if(length(ttrt)==0){
    nodata<- append(nodata, stdnms[i])
  }else{
    
    #blocking for the traits measured
    blkfac<- 'blockNumber'
    minBlknos<- c()
    for(a in 1:length(ttrt)){
      minBlkno<- length(unique(na.omit(trl[,c(blkfac, ttrt[a])])[,blkfac]))
      minBlknos<- append(minBlknos, minBlkno)
    }
    minBlkno<- max(minBlknos)
    
    #remove traits with no replication
    repTF<- minBlknos==1
    if(any(repTF)){
      trtnonrep<- append(trtnonrep, paste(stdnms[i], ttrt[which(repTF)], sep="-"))
      ttrt<- ttrt[-which(repTF)]
    }
    
    #single-trait or multitrait model?
    if(length(ttrt)==1){
      uvvmv<- "UV"
      clasfy<- 'germplasmDbId'
    }else{
      uvvmv<- "MV"
      clasfy<- 'germplasmDbId:trait'
    }

    #############################
    ## Create the fixed formula
    #############################
    if(uvvmv== "UV"){
      fxform<- paste(ttrt, "~1+germplasmDbId", sep="")
    }
    if(uvvmv== "MV"){
      fxform<- paste('cbind(', paste(ttrt, collapse=", "), ")~1+trait+us(trait):germplasmDbId", sep="")
    
    #add the blocking factor if any
    if(uvvmv== "MV"){
      if(minBlkno>1 & minBlkno<4){
        fxform<- paste(fxform, "at(trait):blockNumber", sep="+")
      }
    }else{
      fxform<- paste(fxform, "blockNumber", sep="+")
    }
    
    #convert to formula
    fxform<- as.formula(fxform)
    
    #############################
    ## Create the random formula
    #############################
    
    if(uvvmv== "MV"){
      rform<-"~at(trait):rowNumber+at(trait):colNumber"
    }else{
      #add the blocking factor if any
      if(minBlkno>=4){
        rform<- paste(rform, 'at(trait):blockNumber', sep="+")
      }
      rform<- as.formula(rform)
    }
    
    if(uvvmv== "UV"){
      rform<-"~rowNumber+colNumber"
    }else{
      #add the blocking factor if any
      if(minBlkno>=4){
        rform<- paste(rform, 'blockNumber', sep="+")
      }
      rform<- as.formula(rform)
    }

    
    #############################
    ## Convert variables to factors
    #############################
    trl$germplasmDbId<- as.factor(as.character(trl$germplasmDbId))
    trl$blockNumber<- as.factor(as.character(trl$blockNumber))
    trl$entryType<- as.factor(as.character(trl$entryType))
    trl$replicate<- as.factor(as.character(trl$replicate))
    trl$rowNumber<- as.factor(as.character(trl$rowNumber))
    trl$colNumber<- as.factor(as.character(trl$colNumber))
    
    #################################
    ## Fit model and extract results
    #################################
    
    if(uvvmv=='MV'){
      if(class(rform)=='logical'){
        mod<- suppressWarnings(asreml(fixed=fxform, residual=~id(units):us(trait), data=trl, trace=FALSE, aom=T, workspace=64e6))
      }else{
        mod<- suppressWarnings(asreml(fixed=fxform, random=rform, residual= ~id(units):us(trait), data=trl, trace=FALSE, aom=T, workspace=64e6))
      }
      mod<- mkConv(mod)
      p<- suppressWarnings(predictPlus(mod, classify = clasfy, meanLSD.type='factor.combination', LSDby = 'trait', pworkspace=64e7))
    }
    if(uvvmv=='UV'){
      if(class(rform)=='logical'){
        mod<- suppressWarnings(asreml(fixed=fxform, data=trl, trace=FALSE, aom=T,workspace=64e6))
      }else{
        mod<- suppressWarnings(asreml(fixed=fxform, random=rform, data=trl, trace=FALSE, aom=T,workspace=64e6))
      }
      mod<- mkConv(mod)
      p<- suppressWarnings(predictPlus(mod, classify = clasfy, pworkspace=64e6))
      
    }
    blues<- p$predictions
    
    #compute the LSD
    LSDs<- p$LSD[,'meanLSD']
    names(LSDs)<- row.names(p$LSD)
    
    #add study name to the blues table
    if(uvvmv=='UV'){
      df<- data.frame(studyName=stdnms[i], trait=ttrt, blues)
      df<- df[,c(1,3,2, 4:ncol(df))]
    }
    if(uvvmv=='MV'){
      df<- data.frame(studyName=stdnms[i], blues)
    }
    
    #get residuals
    jpeg(file=paste(stdnms[i], "-residuals.jpeg", sep=""))
    resids<- resid(mod,type="stdCond")
    plot(resid(mod, type="stdCond"), main=stdnms[i])
    dev.off()
    
    #make potential outlier table
    mltTrl<- melt(trl, id.vars=c('observationUnitName','germplasmDbId'), measure.vars=ttrt)
    mltTrl<-mltTrl[order(mltTrl$variable),]
    mltTrl<-mltTrl[order(mltTrl$observationUnitName),]
    residsTab<- cbind(mltTrl, resids)
    outTab<- residsTab[which(sqrt(resids^2)>3),]
    
    if(stdnms[i] %in% stdnmsCoop){
      #################################
      ## create basic single trial analysis summary table 
      ## Not very useful for selection, but some people like to see it
      #################################
      #means
      smryA<- cast(df, studyName+germplasmDbId~trait, value='predicted.value')
      Means<- colMeans(smryA[,-c(1:2)], na.rm=TRUE)
      #standard errors
      smryB<- cast(df, studyName+germplasmDbId~trait, value='standard.error')
      SE<- colMeans(smryB[,-c(1:2)], na.rm=T)
      #trial metadata
      metaCols<- c("studyName", "studyDescription", "studyYear", "studyDesign", "locationName", "germplasmName","germplasmDbId")
      meta<- trl[match(smryA$germplasmDbId, trl$germplasmDbId), metaCols]
      #merge results with meta
      df2<- merge(meta, smryA[,-1], by='germplasmDbId')
      #convert factors to characters
      df2[,c(1:7)] <- lapply(df2[,c(1:7)], as.character)
      #add the IL to Illinois names
      df2$germplasmName[which(df2$germplasmName %in% coopgids)]<- paste("IL", df2$germplasmName[which(df2$germplasmName %in% coopgids)], sep="")
      trl$germplasmName[which(trl$germplasmName %in% coopgids)]<- paste("IL", trl$germplasmName[which(trl$germplasmName %in% coopgids)], sep="")
      #re-order the rows according to the raw data
      trl$replicate<- as.numeric(as.character(trl$replicate))
      trl1<- trl[which(trl$replicate==min(trl$replicate, na.rm=T)),]
      df2<- df2[match(trl1$germplasmName,df2$germplasmName),]
      
      #count number of rows
      nrowdf2<- nrow(df2)
      
      #add means, SE, LSD, and CV
      df2<- addRows(df2, 7, "MEAN", Means)
      df2<- addRows(df2, 7, "SE", SE)
      df2<- addRows(df2, 7, "LSD", LSDs)
      df2[,-c(1:7)]<- round(df2[,-c(1:7)],2)
      
      #add the number of replicates
      nrep<-c()
      for(k in 1:length(ttrt)){
        tb<- table(na.omit(trl[,c('germplasmDbId', ttrt[k])])[,1])
        tb[which(tb==0)]<- NA
        nrep<- append(nrep, mean(tb, na.rm=T))
      }
      df2<- addRows(df2, 7, "No. of Reps", nrep)
      df2[,-c(1:7)]<- round(df2[,-c(1:7)],2)
      
      #add CV
      df2<- addRows(df2, 7, "CV", round(sqrt(nrep)*SE / Means *100,1))
      
      #add rankings
      df2rnk<- df2[, ttrt]
      for(k in 1:ncol(df2rnk)){
        const<- 1
        if(ttrt[k] %in% c("Grain.test.weight...lbs.bu", "Grain.yield...bu.ac") ){
          const<- -1
        }
        rnk<- rank(df2rnk[1:nrowdf2,k]*const, ties.method='min')
        rnk[which(is.na(df2rnk[1:nrowdf2,k]))]<- NA
        df2rnk[,k]<- c(rnk, rep(NA, nrow(df2)-nrowdf2)) 
      }
      colnames(df2rnk)<- paste(colnames(df2rnk), "Rank", sep="-")
      df3<- cbind(df2, df2rnk)
      df3<- df3[,c(colnames(df3)[1:7], sort(colnames(df3[,-c(1:7)])))]
      
      #rename rank cols
      colnames(df3)[grep("Rank", colnames(df3))]<- 'rank'
      
      #add entno to UE
      if(length(grep('UE', df3$studyName[1]))==1){
        ueent<- read.csv('UEentno.csv')
        df3<- merge(ueent, df3, by='germplasmName')
        df3<- df3[order(df3$ent),]
      }
      
      #write cc csv file
      if(length(grep('CC', df3$studyName[1]))==1){
        write.csv(df3, file='CC_Urb_21.csv')
      }
      
      #write to an excel sheet
      excellist<- list(results=df3, rawdata=trl)
      WriteXLS::WriteXLS(excellist, ExcelFileName=paste(stdnms[i], ".xls", sep=""), 
                         SheetNames=c('results', 'rawdata'))
    }  
    
    
    #################################
    ## combine all analysis results into one table
    #################################    
    if(exists('dfall')){
      dfall<- rbind(dfall, df)
      outTabs<- rbind(outTabs, outTab)
    }else{
      dfall<- df
      outTabs<- outTab
    }
    
  }
  }
}

name<- data[match(dfall$germplasmDbId, data$germplasmDbId),'germplasmName']
dfall<- data.frame(name, dfall)
write.csv(dfall, file='predicted.values.table_MVsta.csv')
write.csv(outTabs, file='possible.outliers.csv')

##subset the CC data
cc<- read.csv('CC_Urb_21.csv')
ccraw<- read.csv('CC_Urb_21raw.csv')
A5S<- as.character(unique(dfall[grep('A5S', dfall$studyName),'name']))
A5S[na.omit(match(coopgids, A5S))]<- paste("IL", "", A5S[na.omit(match(coopgids, A5S))], sep="")
A5S<- sort(A5S)

P5S<- as.character(unique(dfall[grep('P5S', dfall$studyName),'name']))
P5S[na.omit(match(coopgids, P5S))]<- paste("IL", "", P5S[na.omit(match(coopgids, P5S))], sep="")
P5S<- sort(P5S)

UE<- as.character(unique(dfall[grep('UE', dfall$studyName),'name']))
UE[na.omit(match(coopgids, UE))]<- paste("IL", "", UE[na.omit(match(coopgids, UE))], sep="")
UE<- sort(UE)


A5S_Urb_21<- cc[c(na.omit(match(A5S, cc$germplasmName)),151:155), ]
A5S_Urb_21$rank[1:(nrow(A5S_Urb_21)-5)]<- rank(A5S_Urb_21[1:(nrow(A5S_Urb_21)-5), 'Heading.time...Julian.date..JD..'])
A5S_Urb_21$rank.1[1:(nrow(A5S_Urb_21)-5)]<- rank(A5S_Urb_21[1:(nrow(A5S_Urb_21)-5), 'Plant.height.inches'])

P5S_Urb_21<- cc[c(na.omit(match(P5S, cc$germplasmName)),151:155), ]
P5S_Urb_21$rank[1:(nrow(P5S_Urb_21)-5)]<- rank(P5S_Urb_21[1:(nrow(P5S_Urb_21)-5), 'Heading.time...Julian.date..JD..'])
P5S_Urb_21$rank.1[1:(nrow(P5S_Urb_21)-5)]<- rank(P5S_Urb_21[1:(nrow(P5S_Urb_21)-5), 'Plant.height.inches'])

UE_Urb_21<- cc[c(na.omit(match(UE, cc$germplasmName)),151:155), ]
UE_Urb_21$rank[1:(nrow(UE_Urb_21)-5)]<- rank(UE_Urb_21[1:(nrow(UE_Urb_21)-5), 'Heading.time...Julian.date..JD..'])
UE_Urb_21$rank.1[1:(nrow(UE_Urb_21)-5)]<- rank(UE_Urb_21[1:(nrow(UE_Urb_21)-5), 'Plant.height.inches'])
ueent<- read.csv('UEentno.csv')
UE_Urb_21<- merge(ueent, UE_Urb_21, by='germplasmName')
UE_Urb_21<- UE_Urb_21[order(UE_Urb_21$ent),]


excellist<- list(results=A5S_Urb_21, rawdata=ccraw)
WriteXLS::WriteXLS(excellist, ExcelFileName=paste('A5S_Urb_21', ".xls", sep=""), 
                   SheetNames=c('results', 'rawdata'))

excellist<- list(results=P5S_Urb_21, rawdata=ccraw)
WriteXLS::WriteXLS(excellist, ExcelFileName=paste('P5S_Urb_21', ".xls", sep=""), 
                   SheetNames=c('results', 'rawdata'))

ccraw2<- ccraw[which(ccraw$germplasmName %in% UE_Urb_21$germplasmName),]
excellist<- list(results=UE_Urb_21, rawdata=ccraw2)
WriteXLS::WriteXLS(excellist, ExcelFileName=paste('UE_Urb_21', ".xls", sep=""), 
                   SheetNames=c('results', 'rawdata'))

