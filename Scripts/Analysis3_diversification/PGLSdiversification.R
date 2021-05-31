############
## DIVERSIFICATION ANALYSIS
############

# This script aims to analyse whether diversification rate correlates with brain size measures. Diversification was assessed using the CLADS2 algorithm developed within Helene Morlon's lab by . Maliet.
# It simply corresponds to correlation analyses (while accounting for phylogeny). The same procedure than in previous scripts is done thus, it is decomposed into three steps:
# 1) Cumpute the correlation for each brain area
# 2) Compute the plots
# 3) Carry on the sensitivity/stability etc... analyses

##--------------

# This script allows to perform PGLS analysis to assess the correlation of brain size with spatial co-occurrence data, specifically in areas for which competitive models were shown to fit the best.

# THIS IS STILL IN PROGRESS; FURTHER CLEANING/COMPLETION NEEDED

###Set working directory
setwd("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates")

#Import environment
rm(list=ls())
load("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/REnvironments/Data_spatial_primate.RData")
load("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/REnvironments/geography_traits.RData")

##--------
#Home made functions
#To source all phylogenetics functions (biogeobears + models of evolution)
source("~/PhD/Meta_analysis/Meta_analysis_cognition_primates/Scripts/Functions.R")

#My toolkit
source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/diagnostic_fcns.R")

##--------

##--------
#Libraries

library(ape)
library(phylolm)
library(phytools)
library(caper)
library(nlme)
library(geiger)

##--------

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

################
# Loading the data


################
# Fitting linear models

repetitionBootstrap=200

#####---------------
#Fixing data
#####---------------

summaryData$Family<- NA
summaryData$DietaryGuild <- NA
summaryData$FrugivoryPercent <- NA
summaryData$FolivoryPercent <- NA
summaryData$Brain <- NA
summaryData$Neocortex <- NA 
summaryData$Cerebellum <- NA
summaryData$Hippocampus <- NA
summaryData$Striatum <- NA
summaryData$MOB <- NA 
summaryData$Optic.Tract <- NA 
summaryData$Bodymass <- NA

for(i in 1:nrow(summaryData)){
  #Frugivory 
  value <- max(c(summaryData$Diet_frug_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                 summaryData$Diet_frug_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]]), na.rm=TRUE)
  #NOTE: #Here, I add the max (compared to other script) so as to keep the species as "frugivorous" if threshold falls in between values
  
  value <- value[!is.na(value)]
  if(!is.finite(value)){
    value <- NULL
  }
  #Because there are issues if of length 1
  if(length(value)==1){
    value <- c(value, value)
  }
  
  if(length(value)>0){
    summaryData$FrugivoryPercent[i] <- sample(value, 1)
  }  else{
    summaryData$FrugivoryPercent[i] <- NA
  }
  
  #Folivory 
  value <- max(c(summaryData$Diet_leaves_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                 summaryData$Diet_leaves_willems[summaryData$Species_abbrv==summaryData$Species_abbrv[i]]), na.rm=TRUE)
  #NOTE: #Here, I add the max (compared to other script) so as to keep the species as "folivorous" if threshold falls in between values
  
  value <- value[!is.na(value)]
  #Because there are issues if of length 1
  if(!is.finite(value)){
    value <- NULL
  }
  if(length(value)==1){
    value <- c(value, value)
  }
  
  if(length(value)>0){
    summaryData$FolivoryPercent[i] <- sample(value, 1)
  }  else{
    summaryData$FolivoryPercent[i] <- NA
  }
  
  ##--------------
  #Associate the dietary guild
  if(!is.na(as.numeric(as.character(summaryData$FrugivoryPercent[i])))&as.numeric(as.character(summaryData$FrugivoryPercent[i]))>=frugivoryThreshold){
    summaryData$DietaryGuild[i] <- "Fruit"
  } else if(!is.na(as.numeric(as.character(summaryData$FolivoryPercent[i])))&as.numeric(as.character(summaryData$FolivoryPercent[i]))>=folivoryThreshold){
    summaryData$DietaryGuild[i] <- "Leaf"
  }  else if(!is.na(summaryData$Guild_init_decasien[i])){
    summaryData$DietaryGuild[i] <- summaryData$Guild_init_decasien[i]
    if(summaryData$DietaryGuild[i] =="Om"|summaryData$DietaryGuild[i]=="Frug/Fol"){
      summaryData$DietaryGuild[i] <- "Fruit"
    }
  }  else {
    summaryData$DietaryGuild[i] <- "Other"
  }
  ##--------------
  
  #Brain volume
  value <- mean(c(
    summaryData$Brain_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
    summaryData$Brain_volume_mm3_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
    summaryData$Brain_volume_mm3_todorov[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
    summaryData$Brain_volume_mm3_grueter[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
    summaryData$Brain_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]]), na.rm=TRUE)
  #NOTE: I use here the mean and will take into account variability later. This can be biased by redundancy between datasets but I had no better solution in mind.
  
  value <- value[!is.na(value)]
  #Because there are issues if of length 1
  if(length(value)==1){
    value <- c(value, value)
  }
  
  if(length(value)>0){
    summaryData$Brain[i] <- sample(value, 1)
  }  else{
    summaryData$Brain[i] <- NA
  }
  
  #Striatum
  value <- mean(c(
    summaryData$Striatum_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
    summaryData$Striatum_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]]), na.rm=TRUE)
  #NOTE: I use here the mean and will take into account variability later. This can be biased by redundancy between datasets but I had no better solution in mind.
  
  value <- value[!is.na(value)]
  #Because there are issues if of length 1
  if(length(value)==1){
    value <- c(value, value)
  }
  
  if(length(value)>0){
    summaryData$Striatum[i] <- sample(value, 1)
  }  else{
    summaryData$Striatum[i] <- NA
  }
  
  #Neocortex
  value <- mean(c(summaryData$Neocortex_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                  summaryData$Neocortex_volume_mm3_powell_mosaic[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                  summaryData$Neocortex_volume_mm3_todorov[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                  summaryData$Neocortex_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]]), na.rm=TRUE)
  #NOTE: I use here the mean and will take into account variability later. This can be biased by redundancy between datasets but I had no better solution in mind.
  
  value <- value[!is.na(value)]
  #Because there are issues if of length 1
  if(length(value)==1){
    value <- c(value, value)
  }
  
  if(length(value)>0){
    summaryData$Neocortex[i] <- sample(value, 1)
  }  else{
    summaryData$Neocortex[i] <- NA
  }
  
  #Cerebellum
  value <- mean(c(summaryData$Cerebellum_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                  summaryData$Cerebellum_volume_mm3_powell_mosaic[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                  summaryData$Cerebellum_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]]), na.rm=TRUE)
  #NOTE: I use here the mean and will take into account variability later. This can be biased by redundancy between datasets but I had no better solution in mind.
  
  value <- value[!is.na(value)]
  #Because there are issues if of length 1
  if(length(value)==1){
    value <- c(value, value)
  }
  
  if(length(value)>0){
    summaryData$Cerebellum[i] <- sample(value, 1)
  }  else{
    summaryData$Cerebellum[i] <- NA
  }
  
  #Hippocampus
  value <- mean(c(summaryData$Hippocampus_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                  summaryData$Hippocampus_volume_mm3_todorov[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                  summaryData$Hippocampus_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]]), na.rm=TRUE)
  #NOTE: I use here the mean and will take into account variability later. This can be biased by redundancy between datasets but I had no better solution in mind.
  
  value <- value[!is.na(value)]
  #Because there are issues if of length 1
  if(length(value)==1){
    value <- c(value, value)
  }
  
  if(length(value)>0){
    summaryData$Hippocampus[i] <- sample(value, 1)
  }  else{
    summaryData$Hippocampus[i] <- NA
  }
  
  #MOB
  value <- mean(c(
    summaryData$MOB_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]]), na.rm=TRUE)
  #NOTE: I use here the mean and will take into account variability later. This can be biased by redundancy between datasets but I had no better solution in mind.
  
  value <- value[!is.na(value)]
  #Because there are issues if of length 1
  if(length(value)==1){
    value <- c(value, value)
  }
  
  if(length(value)>0){
    summaryData$MOB[i] <- sample(value, 1)
  }  else{
    summaryData$MOB[i] <- NA
  }
  
  # #Optic.Tract #Too low sample
  # value <- c(
  #   summaryData$Optic.Tract_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
  # 
  # value <- value[!is.na(value)]
  # #Because there are issues if of length 1
  # if(length(value)==1){
  #   value <- c(value, value)
  # }
  # 
  # if(length(value)>0){
  #   summaryData$Optic.Tract[i] <- sample(value, 1)
  # }
  # else{
  #   summaryData$Optic.Tract[i] <- NA
  # }
  
  #Bodymass
  value <- mean(c(summaryData$Body_mass_g_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                  summaryData$Body_mass_g_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                  summaryData$Body_mass_g_pearce[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                  summaryData$Body_mass_g_grueter[summaryData$Species_abbrv==summaryData$Species_abbrv[i]]), na.rm=TRUE)
  #NOTE: I use here the mean and will take into account variability later. This can be biased by redundancy between datasets but I had no better solution in mind.
  
  value <- value[!is.na(value)]
  #Because there are issues if of length 1
  if(length(value)==1){
    value <- c(value, value)
  }
  
  if(length(value)>0){
    summaryData$Bodymass[i] <- sample(value, 1)
  }  else{
    summaryData$Bodymass[i] <- NA
  }
}

summaryData$Family<- summaryData$Family[match(summaryData$Species_abbrv,summaryData$Species_abbrv)]

#All the potential output predictors:

summaryData$ratioBrain <- summaryData$Brain*1.036*(10**-3)/summaryData$Bodymass #Following decasien for multiplication by 1.036
summaryData$EQ <- summaryData$Brain*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775) #Following decasien, according to #Jerison, H. J. Evolution of the Brain and Intelligence (Academic, 1973).

summaryData$Brain.log <- log(summaryData$Brain)
summaryData$EQ.log <- log(summaryData$EQ)

summaryData$ratioNeocortex <- summaryData$Neocortex/summaryData$Brain
summaryData$ratioNeocortex.log  <- log(summaryData$ratioNeocortex)

summaryData$ratioHippocampus <- summaryData$Hippocampus/summaryData$Brain
summaryData$ratioHippocampus.log <- log(summaryData$ratioHippocampus)

summaryData$ratioCerebellum <- summaryData$Cerebellum/summaryData$Brain

summaryData$ratioStriatum <- summaryData$Striatum/summaryData$Brain

summaryData$ratioMOB <- summaryData$MOB/summaryData$Brain
summaryData$ratioMOB.log <- log(summaryData$ratioMOB)

##---
#Load consensus tree

phylo_all <-read.nexus("Raw_data/Tree/consensusTree_10kTrees_Primates_Version3.nex")
phylo_init <- phylo_all

phylo <- force.ultrametric(tree=phylo_init, method="extend")#method="nnls")
is.ultrametric(phylo)

#check if this differ (was done outside the loop and was ok)
plot(phylo_init$edge.length, phylo$edge.length, xlab="Original tree", ylab="Post-chronos tree", pch=20, bty="n")
abline(a=0, b=1, col="gray", lwd=0.5)

phylo <- multi2di(phylo)
is.binary(phylo)

# Diversification/speciation rates
f=67 # Assumed availability based on current knowledge

table_MAPS_rates <- read.table(paste0("Scripts/Analysis3_diversification/diversification/MAPS_speciation_rates_tips_ClaDS2_tree_primate_complete_f",f,".csv"), sep=";",header=T)
head(table_MAPS_rates)
str(table_MAPS_rates)
#Check, should be perfectly linear
plot(table_MAPS_rates$Diversification_rate, table_MAPS_rates$Diversification_rate)

##Give correct names for phylogeny: Match with substr

toSeparate <- as.data.frame(table_MAPS_rates$Species)
colnames(toSeparate) <- "Species"
library(tidyr)
toSeparate <- separate(toSeparate, col="Species", into=c("Name1", "Name2"), sep="_")
toSeparate$Name2 <- substr(toSeparate$Name2, 1, 4)
toSeparate$Species_abbrv <- paste(toSeparate$Name1, toSeparate$Name2, sep="_")

table_MAPS_rates$Species_abbrv <- toSeparate$Species_abbrv
table_MAPS_rates$SpeciesForPhylogeny <- summaryData$SpeciesForPhylogeny[match(table_MAPS_rates$Species_abbrv, 
                                                                              summaryData$Species_abbrv)]

check <- as.data.frame(table(table_MAPS_rates$Species_abbrv))
nrow(check[check$Freq>1,])
check[check$Freq>1,]

table_MAPS_rates[table_MAPS_rates$Species_abbrv %in% check$Var1[check$Freq>1],]

#Correct: those with same name because classification too fine
table_MAPS_rates$SpeciesForPhylogeny[table_MAPS_rates$Species_abbrv %in% check$Var1[check$Freq>1] & table_MAPS_rates$SpeciesForPhylogeny != table_MAPS_rates$Species] <- NA
                   
table_MAPS_rates[table_MAPS_rates$Species_abbrv %in% check$Var1[check$Freq>1],]
#Ok

#Check
summaryData$SpeciesForPhylogeny
table_MAPS_rates$Species[is.na(table_MAPS_rates$SpeciesForPhylogeny)]
#ok

summaryData$Diversification_rate <- table_MAPS_rates$Diversification_rate[match(summaryData$SpeciesForPhylogeny,table_MAPS_rates$SpeciesForPhylogeny)]
hist(summaryData$Diversification_rate)#ok

##----
## Run analysis for each trait

traitToStudy=c("EQ.log", "ratioBrain", "ratioHippocampus.log", "ratioNeocortex.log", "ratioCerebellum", "ratioStriatum", "ratioMOB.log")
traitName=c("EQ (log)", "Brain", "Hippocampus (log)", "Neocortex (log)", "Cerebellum", "Striatum", "MOB (log)")
vifVector <- rep(NA, times=length(traitName))  

results.df <- as.data.frame(matrix(NA, nrow=4*length(traitToStudy), ncol=7))
colnames(results.df) <- c("", "Est.", "CI2.5%", "CI97.5%", "Sd", "t", "p-value")

for(a in 1:length(traitToStudy)){
  #Matching the brain trait to the dataset with predictors
  dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[a])]
  dataRangePrimate$Diversification_rate <-  summaryData$Diversification_rate[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny)]
  dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,4])&!is.na(dataRangePrimate[,5]),]
  
  #Readjust phylo tree
  phyloConsensus <- drop.tip(phylo,
                             phylo$tip.label[
                               which(phylo$tip.label
                                     %nin%dataRangePrimate_rdc$Species)])
  
  #Studied trait + sample size
  firstRowWrite <- which(is.na(results.df[,1]))[1]
  results.df[firstRowWrite, 1] <- paste("Diversification ", traitName[a], " (N=", nrow(dataRangePrimate_rdc), ")", sep="") 
  results.df[(firstRowWrite+1):(firstRowWrite+3), 1] <- c("Intercept", traitName[a], "Lambda") 
  
  ##--------
  ## Main model: based on consensus tree
  ##--------
  
  rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species
  modelBrain <- phylolm(formula = Diversification_rate ~ Trait, 
                        data=dataRangePrimate_rdc ,phy=phyloConsensus,model="lambda",measurement_error=FALSE,boot=repetitionBootstrap)
  
  # modelBrain <- gls(Trait ~ Overlap_average.sqrt  correlation = corPagel(1, phyloConsensus, form = ~Species),
  #                   data = dataRangePrimate_rdc, method = "ML")
  # CI <- intervals(modelBrain)
  # results <- anova(modelBrain)
  
  #Complete df results
  
  #commented is when I used pgls + phylolm
  results.df[(firstRowWrite+1):(firstRowWrite+2), 2] <- roundIntelligent(summary(modelBrain)$coefficients[,1], digit=2) #roundIntelligent(summary(modelBrain)$tTable[,1], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+2), 3] <- roundIntelligent(summary(modelBrain)$coefficients[,4], digit=2) #roundIntelligent(CI$coef[,1], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+2), 4] <- roundIntelligent(summary(modelBrain)$coefficients[,5], digit=2) #roundIntelligent(CI$coef[,3], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+2), 5] <- roundIntelligent(summary(modelBrain)$coefficients[,3], digit=2) #roundIntelligent(summary(modelBrain)$tTable[,2], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+2), 6] <- c("-", roundIntelligent(summary(modelBrain)$coefficients[2,3], digit=2))
  results.df[(firstRowWrite+1):(firstRowWrite+2), 7] <- c("-", roundIntelligent(summary(modelBrain)$coefficients[2,6], digit=2))
  
  # val1 <- ifelse(is.na(summary(modelBrainPgls)$param.CI$lambda$ci.val[1]), 0, round(summary(modelBrainPgls)$param.CI$lambda$ci.val[1], digit=2))
  # val2 <- ifelse(is.na(summary(modelBrainPgls)$param.CI$lambda$ci.val[2]), 1, round(summary(modelBrainPgls)$param.CI$lambda$ci.val[2], digit=2))
  # 
  # lambda <- paste(round(summary(modelBrainPgls)$param.CI$lambda$opt, digit=2), " [",val1,",",val2,"]", sep="")
  #         
  results.df[firstRowWrite+3,2] <- roundIntelligent(summary(modelBrain)$bootmean[4])
  results.df[firstRowWrite+3,3] <- roundIntelligent(summary(modelBrain)$bootconfint95[1,4])#CI$corStruct[1,1])
  results.df[firstRowWrite+3,4] <- roundIntelligent(summary(modelBrain)$bootconfint95[2,4])#CI$corStruct[1,3])
  
  #Check assumptions
  assign(paste("modelBrainDiversification", traitName[a], sep="_"), modelBrain)
  #assign(paste(modelBrainPgls, traitName[a], sep="_"), modelBrainPgls)
  
  #vifVectorDiversification[a] <- max(vif(modelBrain))

}

results.df[2:ncol(results.df)] <- replaceNAtable(results.df[2:ncol(results.df)], "")
results.df_diversification <- results.df


repetitionTrees=50
repetitionModels=30

sensitivityEstimate <- as.data.frame(matrix(NA, ncol=2, nrow=repetitionTrees*repetitionModels*length(traitToStudy)))
sensitivityPvalue <- as.data.frame(matrix(NA, ncol=2, nrow=repetitionTrees*repetitionModels*length(traitToStudy)))
sensitivityLambda <- rep(NA, times=repetitionTrees*repetitionModels*length(traitToStudy))
traitVectorSensitivity <- rep(NA, times=repetitionTrees*repetitionModels*length(traitToStudy))

for(d in 1:repetitionModels){#data
    
    summaryData$Family<- NA
    summaryData$DietaryGuild <- NA
    summaryData$FrugivoryPercent <- NA
    summaryData$FolivoryPercent <- NA
    summaryData$Brain <- NA
    summaryData$Neocortex <- NA 
    summaryData$Cerebellum <- NA
    summaryData$Hippocampus <- NA
    summaryData$Striatum <- NA
    summaryData$MOB <- NA 
    summaryData$Optic.Tract <- NA 
    summaryData$Bodymass <- NA
    
    for(i in 1:nrow(summaryData)){
      #Frugivory 
      value <- max(c(summaryData$Diet_frug_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                     summaryData$Diet_frug_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]]), na.rm=TRUE)
      #NOTE: #Here, I add the max (compared to other script) so as to keep the species as "frugivorous" if threshold falls in between values
      
      value <- value[!is.na(value)]
      if(!is.finite(value)){
        value <- NULL
      }
      #Because there are issues if of length 1
      if(length(value)==1){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryData$FrugivoryPercent[i] <- sample(value, 1)
      }  else{
        summaryData$FrugivoryPercent[i] <- NA
      }
      
      #Folivory 
      value <- max(c(summaryData$Diet_leaves_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                     summaryData$Diet_leaves_willems[summaryData$Species_abbrv==summaryData$Species_abbrv[i]]), na.rm=TRUE)
      #NOTE: #Here, I add the max (compared to other script) so as to keep the species as "folivorous" if threshold falls in between values
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(!is.finite(value)){
        value <- NULL
      }
      if(length(value)==1){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryData$FolivoryPercent[i] <- sample(value, 1)
      }  else{
        summaryData$FolivoryPercent[i] <- NA
      }
      
      ##--------------
      #Associate the dietary guild
      if(!is.na(as.numeric(as.character(summaryData$FrugivoryPercent[i])))&as.numeric(as.character(summaryData$FrugivoryPercent[i]))>=frugivoryThreshold){
        summaryData$DietaryGuild[i] <- "Fruit"
      } else if(!is.na(as.numeric(as.character(summaryData$FolivoryPercent[i])))&as.numeric(as.character(summaryData$FolivoryPercent[i]))>=folivoryThreshold){
        summaryData$DietaryGuild[i] <- "Leaf"
      }  else if(!is.na(summaryData$Guild_init_decasien[i])){
        summaryData$DietaryGuild[i] <- summaryData$Guild_init_decasien[i]
        if(summaryData$DietaryGuild[i] =="Om"|summaryData$DietaryGuild[i]=="Frug/Fol"){
          summaryData$DietaryGuild[i] <- "Fruit"
        }
      }  else{
        summaryData$DietaryGuild[i] <- "Other"
      }
      ##--------------
      
      #Brain volume
      value <- c(
        summaryData$Brain_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
        summaryData$Brain_volume_mm3_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
        summaryData$Brain_volume_mm3_todorov[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
        summaryData$Brain_volume_mm3_grueter[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
        summaryData$Brain_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value)==1){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryData$Brain[i] <- sample(value, 1)
      }  else{
        summaryData$Brain[i] <- NA
      }
      
      #Striatum
      value <- c(
        summaryData$Striatum_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
        summaryData$Striatum_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value)==1){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryData$Striatum[i] <- sample(value, 1)
      }  else{
        summaryData$Striatum[i] <- NA
      }
      
      #Neocortex
      value <- c(summaryData$Neocortex_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                 summaryData$Neocortex_volume_mm3_powell_mosaic[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                 summaryData$Neocortex_volume_mm3_todorov[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                 summaryData$Neocortex_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value)==1){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryData$Neocortex[i] <- sample(value, 1)
      }  else{
        summaryData$Neocortex[i] <- NA
      }
      
      #Cerebellum
      value <- c(summaryData$Cerebellum_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                 summaryData$Cerebellum_volume_mm3_powell_mosaic[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                 summaryData$Cerebellum_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value)==1){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryData$Cerebellum[i] <- sample(value, 1)
      }  else{
        summaryData$Cerebellum[i] <- NA
      }
      
      #Hippocampus
      value <- c(summaryData$Hippocampus_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                 summaryData$Hippocampus_volume_mm3_todorov[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                 summaryData$Hippocampus_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value)==1){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryData$Hippocampus[i] <- sample(value, 1)
      }  else{
        summaryData$Hippocampus[i] <- NA
      }
      
      #MOB
      value <- c(
        summaryData$MOB_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value)==1){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryData$MOB[i] <- sample(value, 1)
      }  else{
        summaryData$MOB[i] <- NA
      }
      
      # #Optic.Tract #Too low sample
      # value <- c(
      #   summaryData$Optic.Tract_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
      # 
      # value <- value[!is.na(value)]
      # #Because there are issues if of length 1
      # if(length(value)==1){
      #   value <- c(value, value)
      # }
      # 
      # if(length(value)>0){
      #   summaryData$Optic.Tract[i] <- sample(value, 1)
      # }
      # else{
      #   summaryData$Optic.Tract[i] <- NA
      # }
      
      #Bodymass
      value <- c(summaryData$Body_mass_g_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                 summaryData$Body_mass_g_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                 summaryData$Body_mass_g_pearce[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                 summaryData$Body_mass_g_grueter[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value)==1){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryData$Bodymass[i] <- sample(value, 1)
      } else{
        summaryData$Bodymass[i] <- NA
      }
    }
  
    summaryData$Family<- summaryData$Family[match(summaryData$Species_abbrv,summaryData$Species_abbrv)]
    
    #All the potential output predictors:
    
    summaryData$ratioBrain <- summaryData$Brain*1.036*(10**-3)/summaryData$Bodymass #Following decasien for multiplication by 1.036
    summaryData$EQ <- summaryData$Brain*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775) #Following decasien, according to #Jerison, H. J. Evolution of the Brain and Intelligence (Academic, 1973).
    
    summaryData$Brain.log <- log(summaryData$Brain)
    summaryData$EQ.log <- log(summaryData$EQ)
    
    summaryData$ratioNeocortex <- summaryData$Neocortex/summaryData$Brain
    summaryData$ratioNeocortex.log  <- log(summaryData$ratioNeocortex)
    
    summaryData$ratioHippocampus <- summaryData$Hippocampus/summaryData$Brain
    summaryData$ratioHippocampus.log <- log(summaryData$ratioHippocampus)
    
    summaryData$ratioCerebellum <- summaryData$Cerebellum/summaryData$Brain
    
    summaryData$ratioStriatum <- summaryData$Striatum/summaryData$Brain
    
    summaryData$ratioMOB <- summaryData$MOB/summaryData$Brain
    summaryData$ratioMOB.log <- log(summaryData$ratioMOB)
    
    for(t in 1:repetitionTrees){
      
      #Dowload tree
      phylo_all <-read.nexus("Raw_data/Tree/TreeBlock_10kTrees_Primates_Version3.nex")
      phylo_init <- phylo_all[[t]]
      
      phylo <- force.ultrametric(tree=phylo_init, method="extend")#method="nnls")
      is.ultrametric(phylo)
      
      #check if this differ (was done outside the loop and was ok)
      plot(phylo_init$edge.length, phylo$edge.length, xlab="Original tree", ylab="Post-chronos tree", pch=20, bty="n")
      abline(a=0, b=1, col="gray", lwd=0.5)
      
      phylo <- multi2di(phylo)
      is.binary(phylo)
    
    ##----
    ## RUN MODEL FOR SENSITIVITY
    
    for(a in 1:length(traitToStudy)){
      
      rowToAdd=which(is.na(sensitivityLambda))[1]
      traitVectorSensitivity[rowToAdd] <- traitToStudy[a]
      
      #Matching the brain trait to the dataset with predictors
      dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[a])]
      dataRangePrimate$Diversification_rate <-  summaryData$Diversification_rate[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny)]
      dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,4])&!is.na(dataRangePrimate[,5]),]
      
      #Readjust phylo tree
      phyloConsensus <- drop.tip(phylo,
                                 phylo$tip.label[
                                   which(phylo$tip.label
                                         %nin%dataRangePrimate_rdc$Species)])
      
      #Studied trait + sample size
      firstRowWrite <- which(is.na(results.df[,1]))[1]
      results.df[firstRowWrite, 1] <- paste("Diversification ", traitName[a], " (N=", nrow(dataRangePrimate_rdc), ")", sep="") 
      results.df[(firstRowWrite+1):(firstRowWrite+2), 1] <- c("Intercept",traitName[a], "Lambda") 
      
      ##--------
      ## Main model: based on consensus tree
      ##--------

      rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species
      modelBrain <- phylolm(formula = Diversification_rate ~ Trait, 
                            data=dataRangePrimate_rdc ,phy=phyloConsensus,model="lambda",measurement_error=FALSE,boot=repetitionBootstrap)
      results <- anova(modelBrain)

      #Complete df results
      sensitivityEstimate[rowToAdd,] <- summary(modelBrain)$coefficients[,1]#round(summary(modelBrain)$coefficients[,1], digit=4)
      sensitivityPvalue[rowToAdd,] <- summary(modelBrain)$coefficients[,4]#c("-", round(results[1:2,5], digit=3))
      sensitivityLambda[rowToAdd] <- summary(modelBrain)$optpar#CI$corStruct[1,2])#round(summary(modelBrain)$param.CI$lambda$opt, digit=2)
      traitVectorSensitivity[rowToAdd] <- traitToStudy[a]
    }
  }
}

#Calculate Dfbetas

dfBetasEstimate <- as.data.frame(matrix(NA, ncol=2, nrow=nrow(dataRangePrimate)*length(traitToStudy)))
dfBetasPvalue <- as.data.frame(matrix(NA, ncol=2, nrow=nrow(dataRangePrimate)*length(traitToStudy)))
dfBetasLambda <- rep(NA, times=nrow(dataRangePrimate)*length(traitToStudy))
traitVectordfBetas <- rep(NA, times=nrow(dataRangePrimate)*length(traitToStudy))

for(a in 1:length(traitToStudy)){
  
  #Matching the brain trait to the dataset with predictors
  dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[a])]
  dataRangePrimate$Diversification_rate <-  summaryData$Diversification_rate[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny)]
  dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,4])&!is.na(dataRangePrimate[,5]),]
  
  for(r in 1:nrow(dataRangePrimate_rdc)){
      
    dataRangePrimate_rdc2 <- dataRangePrimate_rdc [-r,] 
    #Readjust phylo tree
    phyloConsensus <- drop.tip(phylo,
                               phylo$tip.label[
                                 which(phylo$tip.label
                                       %nin%dataRangePrimate_rdc2$Species)])
    
    ##--------
    ## Main model: based on consensus tree
    ##--------
    
    rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species
    modelBrain <- phylolm(formula = Diversification_rate ~ Trait, 
                          data=dataRangePrimate_rdc ,phy=phyloConsensus,model="lambda",measurement_error=FALSE,boot=repetitionBootstrap)
    
    rowToAdd=which(is.na(dfBetasLambda))[1]
    traitVectordfBetas[rowToAdd] <- traitToStudy[a]
    
    dfBetasEstimate[rowToAdd,] <- summary(modelBrain)$coefficients[,1]
    dfBetasPvalue[rowToAdd,] <- summary(modelBrain)$coefficients[,4]
    dfBetasLambda[rowToAdd] <- summary(modelBrain)$optpar
    traitVectordfBetas[rowToAdd] <- traitToStudy[a]
  }
}



fraction.v <- c(60, 70, 80, 90, 95)

sensitivityFractionTaxonEstimate <- as.data.frame(matrix(NA, ncol=2, nrow=length(fraction.v)*length(traitToStudy)))
sensitivityFractionTaxonPvalue <- as.data.frame(matrix(NA, ncol=2, nrow=length(fraction.v)*length(traitToStudy)))
sensitivityFractionTaxonLambda <- rep(NA, times=length(fraction.v)*length(traitToStudy))
traitVectorSensitivityFractionTaxon <- rep(NA, times=length(fraction.v)*length(traitToStudy))

for (f in fraction.v){
  
  table_MAPS_rates <- read.table(paste0("Scripts/Analysis3_diversification/diversification/MAPS_speciation_rates_tips_ClaDS2_tree_primate_complete_f",f,".csv"), sep=";",header=T)
  head(table_MAPS_rates)
  str(table_MAPS_rates)
  #Check, should be perfectly linear
  plot(table_MAPS_rates$Diversification_rate, table_MAPS_rates$Diversification_rate)
  
  ##Give correct names for phylogeny: Match with substr
  
  toSeparate <- as.data.frame(table_MAPS_rates$Species)
  colnames(toSeparate) <- "Species"
  library(tidyr)
  toSeparate <- separate(toSeparate, col="Species", into=c("Name1", "Name2"), sep="_")
  toSeparate$Name2 <- substr(toSeparate$Name2, 1, 4)
  toSeparate$Species_abbrv <- paste(toSeparate$Name1, toSeparate$Name2, sep="_")
  
  table_MAPS_rates$Species_abbrv <- toSeparate$Species_abbrv
  table_MAPS_rates$SpeciesForPhylogeny <- summaryData$SpeciesForPhylogeny[match(table_MAPS_rates$Species_abbrv, 
                                                                                summaryData$Species_abbrv)]
  
  check <- as.data.frame(table(table_MAPS_rates$Species_abbrv))
  nrow(check[check$Freq>1,])
  check[check$Freq>1,]
  
  table_MAPS_rates[table_MAPS_rates$Species_abbrv %in% check$Var1[check$Freq>1],]
  
  #Correct: those with same name because classification too fine
  table_MAPS_rates$SpeciesForPhylogeny[table_MAPS_rates$Species_abbrv %in% check$Var1[check$Freq>1] & table_MAPS_rates$SpeciesForPhylogeny != table_MAPS_rates$Species] <- NA
  
  table_MAPS_rates[table_MAPS_rates$Species_abbrv %in% check$Var1[check$Freq>1],]
  #Ok
  
  #Check
  summaryData$SpeciesForPhylogeny
  table_MAPS_rates$Species[is.na(table_MAPS_rates$SpeciesForPhylogeny)]
  #ok
  
  summaryData$Diversification_rate <- table_MAPS_rates$Diversification_rate[match(summaryData$SpeciesForPhylogeny,table_MAPS_rates$SpeciesForPhylogeny)]
  hist( summaryData$Diversification_rate)#ok
  
  ##----
  ## Run analysis for each trait
  
  traitToStudy=c("EQ.log", "ratioBrain", "ratioHippocampus.log", "ratioNeocortex.log", "ratioCerebellum", "ratioStriatum", "ratioMOB.log")
  traitName=c("EQ (log)", "Brain", "Hippocampus (log)", "Neocortex (log)", "Cerebellum", "Striatum", "MOB (log)")
  vifVector <- rep(NA, times=length(traitName))  
  
  results.df <- as.data.frame(matrix(NA, nrow=4*length(traitToStudy), ncol=8))
  colnames(results.df) <- c("", "Est.", "CI2.5%", "CI97.5%", "Sd", "Df", "F", "p-value")
  
  for(a in 1:length(traitToStudy)){
    
    rowToAdd=which(is.na(sensitivityFractionTaxonLambda))[1]
    traitVectorSensitivityFractionTaxon[rowToAdd] <- traitToStudy[a]
    
    #Matching the brain trait to the dataset with predictors
    dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[a])]
    dataRangePrimate$Diversification_rate <-  summaryData$Diversification_rate[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny)]
    dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,4])&!is.na(dataRangePrimate[,5]),]
    
    #Readjust phylo tree
    phyloConsensus <- drop.tip(phylo,
                               phylo$tip.label[
                                 which(phylo$tip.label
                                       %nin%dataRangePrimate_rdc$Species)])
    
    ##--------
    ## Main model: based on consensus tree
    ##--------
    
    rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species
    modelBrain <- phylolm(formula = Diversification_rate ~ Trait, 
                          data=dataRangePrimate_rdc ,phy=phyloConsensus,model="lambda",measurement_error=FALSE,boot=repetitionBootstrap)
    
    #Complete df results
    sensitivityFractionTaxonEstimate[rowToAdd,] <- summary(modelBrain)$coefficients[,1]
    sensitivityFractionTaxonPvalue[rowToAdd,] <- summary(modelBrain)$coefficients[,4]
    sensitivityFractionTaxonLambda[rowToAdd] <- summary(modelBrain)$optpar
    traitVectorSensitivityFractionTaxon[rowToAdd] <- traitToStudy[a]
    
  }
}

###----
# Create sensitivity tables

#Summary sensitivity

#Dfbetas
library(tidyr)
minEst <- tibble(cbind(aggregate(dfBetasEstimate, by=list(traitVectordfBetas), min), aggregate(dfBetasLambda, by=list(traitVectordfBetas), min)[,2]))
colnames(minEst) <- c("Trait", "Intercept", "Overlap", "N co-occurrence", "Lambda")
minEst <- pivot_longer(minEst, col=2:5, names_to="Variable", values_to="Estimate")

maxEst <- tibble(cbind(aggregate(dfBetasEstimate, by=list(traitVectordfBetas), max), aggregate(dfBetasLambda, by=list(traitVectordfBetas), max)[,2]))
colnames(maxEst) <- c("Trait", "Intercept", "Overlap", "N co-occurrence", "Lambda")
maxEst <- pivot_longer(maxEst, col=2:5, names_to="Variable", values_to="Estimate")

estimate <- results.df[results.df[,2]!="",2]
whatIs <- rep(traitToStudy, each=4)
numberForFinishingMatch <- rep(1:4, times=length(traitToStudy))

#medianEst <- aggregate(dfBetasEstimate, by=list(traitVectordfBetas), median)
minEst$EstimateTrue <- estimate[match(paste(minEst$Trait,numberForFinishingMatch),paste(whatIs,numberForFinishingMatch))]

dfBetasDiversification <- cbind(minEst,
                         maxEst[,3])
colnames(dfBetasDiversification) <- c("Trait", "Variable", "Est. min.", "Est.", "Est. max.")


dfBetasDiversification[,3] <- sapply(dfBetasDiversification[,3], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)

dfBetasDiversification[,5] <- sapply(dfBetasDiversification[,5], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)

# pid <- aggregate(dfBetasPvalue, by=list(traitVectordfBetas), FUN= function(x)
#   length(x>0.05)
# )

#Sensitivity to phylo/data
minEst <- tibble(cbind(aggregate(sensitivityEstimate, by=list(traitVectorSensitivity), min), aggregate(sensitivityLambda, by=list(traitVectorSensitivity), min)[,2]))
colnames(minEst) <- c("Trait", "Intercept", "Overlap", "N co-occurrence", "Lambda")
minEst <- pivot_longer(minEst, col=2:5, names_to="Variable", values_to="Estimate")

maxEst <- tibble(cbind(aggregate(sensitivityEstimate, by=list(traitVectorSensitivity), max), aggregate(sensitivityLambda, by=list(traitVectorSensitivity), max)[,2]))
colnames(maxEst) <- c("Trait", "Intercept", "Overlap", "N co-occurrence", "Lambda")
maxEst <- pivot_longer(maxEst, col=2:5, names_to="Variable", values_to="Estimate")

estimate <- results.df[results.df[,2]!="",2]
whatIs <- rep(traitToStudy, each=4)
numberForFinishingMatch <- rep(1:4, times=length(traitToStudy))

#medianEst <- aggregate(sensitivityEstimate, by=list(traitVectorSensitivity), median)
minEst$EstimateTrue <- estimate[match(paste(minEst$Trait,numberForFinishingMatch),paste(whatIs,numberForFinishingMatch))]

sensitivityDiversification <- cbind(minEst,
                             maxEst[,3])
colnames(sensitivityDiversification) <- c("Trait", "Variable", "Est. min.", "Est.", "Est. max.")

sensitivityDiversification[,3] <- sapply(sensitivityDiversification[,3], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)

sensitivityDiversification[,5] <- sapply(sensitivityDiversification[,5], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)

#Sensitivity to sampling frac
minEst <- tibble(cbind(aggregate(sensitivityFractionTaxonEstimate, by=list(traitVectorsensitivityFractionTaxon), min), aggregate(sensitivityFractionTaxonLambda, by=list(traitVectorsensitivityFractionTaxon), min)[,2]))
colnames(minEst) <- c("Trait", "Intercept", "Overlap", "N co-occurrence", "Lambda")
minEst <- pivot_longer(minEst, col=2:5, names_to="Variable", values_to="Estimate")

maxEst <- tibble(cbind(aggregate(sensitivityFractionTaxonEstimate, by=list(traitVectorsensitivityFractionTaxon), max), aggregate(sensitivityFractionTaxonLambda, by=list(traitVectorsensitivityFractionTaxon), max)[,2]))
colnames(maxEst) <- c("Trait", "Intercept", "Overlap", "N co-occurrence", "Lambda")
maxEst <- pivot_longer(maxEst, col=2:5, names_to="Variable", values_to="Estimate")

estimate <- results.df[results.df[,2]!="",2]
whatIs <- rep(traitToStudy, each=4)
numberForFinishingMatch <- rep(1:4, times=length(traitToStudy))

#medianEst <- aggregate(sensitivityFractionTaxonEstimate, by=list(traitVectorsensitivityFractionTaxon), median)
minEst$EstimateTrue <- estimate[match(paste(minEst$Trait,numberForFinishingMatch),paste(whatIs,numberForFinishingMatch))]

sensitivityFractionTaxonDiversification <- cbind(minEst,
                                    maxEst[,3])
colnames(sensitivityFractionTaxonDiversification) <- c("Trait", "Variable", "Est. min.", "Est.", "Est. max.")

sensitivityFractionTaxonDiversification[,3] <- sapply(sensitivityFractionTaxonDiversification[,3], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)

sensitivityFractionTaxonDiversification[,5] <- sapply(sensitivityFractionTaxonDiversification[,5], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)

summarySensitivityDiversification <- cbind(dfBetasDiversification, sensitivityDiversification[,3:5], sensitivityFractionTaxonDiversification[,3:5])

knitr::kable(summarySensitivityGradient, escape=TRUE, booktabs = TRUE,
             caption = "Sensitivity analysis of phylogenetic regressions to detect the assess the diversification pattern | Depicted is the minimum and maximum of estimates when one observation was removed at a time (DfBetas), when varying the used phylogenetic tree and the data sampling (Phylogeny/Data), or when the sampling fraction varied (Sampling fraction)") %>% 
  kableExtra::kable_styling(latex_options = "striped") %>%
  kableExtra::kable_styling(latex_options="scale_down") %>%
  kableExtra::add_header_above(c("Model:" = 2, "DfBetas" = 3, "Phylogeny/Data" = 3, "Sampling fraction" = 3))

###---
## Create main table

whichToBold <- which(as.numeric(results.df_diversification$"p-value")<=0.05)
toPlotBold <- rep(FALSE, times=nrow(results.df_diversification))
toPlotBold[whichToBold] <- TRUE

replace <- sapply(as.numeric(results.df_diversification$"p-value"), function(x) if(!is.na(x)){pvalueRound(x, text=FALSE)}else{x})
results.df_diversification$"p-value"[!is.na(replace)] <- replace

knitr::kable(results.df_diversification, escape=TRUE, booktabs = TRUE,
             caption = "Model estimates and significance of phylogenetic regressions to assess the diversification pattern | Est.=Estimate, CI2.5%=Lower border of the CI95%, CI97.5%=Upper border of the CI95%, Sd= Standard deviation, t= Statitsitics t-vale. The brain area (as well as the associated sample size)
             are indicate prior each list of estimates. the transformation (logarithm or square-root) if indicated in parenthese by the abbreviation (log or sqrt).") %>%
  kableExtra::column_spec(2:ncol(results.df_diversification), bold = toPlotBold) %>%
  kableExtra::kable_styling(latex_options = "striped") %>%
  kableExtra::kable_styling(latex_options="scale_down")
