####---------------------------------
#### PGLSdirectionSelection - BAYESIAN UNCOMPLETE JUST CHECK FIRST PART (MAIN MODELS)
####---------------------------------

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
# Fitting linear models

repetitionBootstrap=1000

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

##--
# Ratio body mass
summaryData$ratioBrain <- summaryData$Brain/summaryData$Bodymass
summaryData$ratioBrain.log <- log(summaryData$ratioBrain)

summaryData$EQ <- summaryData$Brain*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775) #Following decasien, according to #Jerison, H. J. Evolution of the Brain and Intelligence (Academic, 1973).

summaryData$Brain.log <- log(summaryData$Brain)
summaryData$EQ.log <- log(summaryData$EQ)

summaryData$ratioNeocortex <- summaryData$Neocortex/summaryData$Bodymass
summaryData$ratioNeocortex.log  <- log(summaryData$ratioNeocortex)

summaryData$ratioHippocampus <- summaryData$Hippocampus/summaryData$Bodymass
summaryData$ratioHippocampus.log <- log(summaryData$ratioHippocampus)

summaryData$ratioCerebellum <- summaryData$Cerebellum/summaryData$Bodymass
summaryData$ratioCerebellum.log <- log(summaryData$ratioCerebellum )

summaryData$ratioStriatum <- summaryData$Striatum/summaryData$Bodymass
summaryData$ratioStriatum.log <- log(summaryData$Striatum/summaryData$Bodymass)

summaryData$ratioMOB <- summaryData$MOB/summaryData$Bodymass
summaryData$ratioMOB.log <- log(summaryData$ratioMOB)

##--
# Ratio EQ

# summaryData$ratioBrain <- summaryData$Brain/summaryData$Bodymass
# summaryData$ratioBrain.log <- log(summaryData$ratioBrain)
# summaryData$EQ <- summaryData$Brain*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775) #Following decasien, according to #Jerison, H. J. Evolution of the Brain and Intelligence (Academic, 1973).
# 
# summaryData$Brain.log <- log(summaryData$Brain)
# summaryData$EQ.log <- log(summaryData$EQ)
# 
# summaryData$ratioNeocortex <- summaryData$Neocortex*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775)#summaryData$Neocortex/summaryData$Brain
# summaryData$ratioNeocortex.log  <- log(summaryData$ratioNeocortex)
# 
# summaryData$ratioHippocampus <- summaryData$Hippocampus*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775)#summaryData$Hippocampus/summaryData$Brain
# summaryData$ratioHippocampus.log <- log(summaryData$ratioHippocampus)
# 
# summaryData$ratioCerebellum <- summaryData$Cerebellum*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775)#summaryData$Cerebellum/summaryData$Brain
# 
# summaryData$ratioStriatum <- summaryData$Striatum*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775)#summaryData$Striatum/summaryData$Brain
# 
# summaryData$ratioMOB <- summaryData$MOB*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775)#summaryData$MOB/summaryData$Brain
# summaryData$ratioMOB.log <- log(summaryData$ratioMOB)


##--
# Ratio brain

# summaryData$ratioBrain <- summaryData$Brain/summaryData$Bodymass
# summaryData$ratioBrain.log <- log(summaryData$ratioBrain)
# summaryData$EQ <- summaryData$Brain*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775) #Following decasien, according to #Jerison, H. J. Evolution of the Brain and Intelligence (Academic, 1973).
# 
# summaryData$Brain.log <- log(summaryData$Brain)
# summaryData$EQ.log <- log(summaryData$EQ)
# 
# summaryData$ratioNeocortex <- summaryData$Neocortex/summaryData$Brain
# summaryData$ratioNeocortex.log  <- log(summaryData$ratioNeocortex)
# 
# summaryData$ratioHippocampus <- summaryData$Hippocampus/summaryData$Brain
# summaryData$ratioHippocampus.log <- log(summaryData$ratioHippocampus)
# 
# summaryData$ratioCerebellum <- summaryData$Cerebellum/summaryData$Brain
# 
# summaryData$ratioStriatum <- summaryData$Striatum/summaryData$Brain
# 
# summaryData$ratioMOB <- summaryData$MOB/summaryData$Brain
# summaryData$ratioMOB.log <- log(summaryData$ratioMOB)

##--
# RAW

# summaryData$ratioBrain <- summaryData$Brain*1.036*(10**-3)/summaryData$Bodymass #Following decasien for multiplication by 1.036
# summaryData$EQ <- summaryData$Brain*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775) #Following decasien, according to #Jerison, H. J. Evolution of the Brain and Intelligence (Academic, 1973).
# 
# summaryData$Brain.log <- log(summaryData$Brain)
# summaryData$EQ.log <- log(summaryData$EQ)
# 
# summaryData$ratioNeocortex <- summaryData$Neocortex
# summaryData$ratioNeocortex.log  <- log(summaryData$ratioNeocortex)
# 
# summaryData$ratioHippocampus <- summaryData$Hippocampus
# summaryData$ratioHippocampus.log <- log(summaryData$ratioHippocampus)
# 
# summaryData$ratioCerebellum <- summaryData$Cerebellum
# 
# summaryData$ratioStriatum <- summaryData$Striatum

# summaryData$ratioMOB <- summaryData$MOB
# summaryData$ratioMOB.log <- log(summaryData$ratioMOB)

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

traitToStudy=c("EQ.log", "ratioBrain.log", "ratioHippocampus.log", "ratioNeocortex.log", "ratioCerebellum.log", "ratioStriatum.log", "ratioMOB.log")
traitName=c("EQ (log)", "Brain (/bodymass, log)", "Hippocampus (/bodymass, log)", "Neocortex (/bodymass, log)", "Cerebellum (/bodymass, log)", "Striatum (/bodymass, log)", "MOB (/bodymass, log)")
vifVector <- rep(NA, times=length(traitName))  

results.df <- as.data.frame(matrix(NA, nrow=5*length(traitToStudy), ncol=7))
colnames(results.df) <- c("", "Est.", "CI2.5%", "CI97.5%", "Sd", "t", "p-value")


##----


for(a in 1:length(traitToStudy)){
  #Matching the brain trait to the dataset with predictors
  dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[a])]
  dataRangePrimate$Bodymass <- summaryData$Bodymass[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny)]
  dataRangePrimate$Bodymass.log <- log(dataRangePrimate$Bodymass)
  
  dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,2])&!is.na(dataRangePrimate[,3])&!is.na(dataRangePrimate[,4]),]
  
  #Readjust phylo tree
  phyloConsensus <- drop.tip(phylo,
                             phylo$tip.label[
                               which(phylo$tip.label
                                     %nin%dataRangePrimate_rdc$Species)])
  
  #Studied trait + sample size
  firstRowWrite <- which(is.na(results.df[,1]))[1]
  results.df[firstRowWrite, 1] <- paste(traitName[a], " (N=", nrow(dataRangePrimate_rdc), ")", sep="") 
  results.df[(firstRowWrite+1):(firstRowWrite+4), 1] <- c("Intercept", "% of overlapped home range", "Number of sympatric frugivorous (sqrt)", "Lambda") 
  
  
  ##--------
  ## Main model: based on consensus tree
  ##--------
  dataRangePrimate_rdc$Number_species_cooccurrence.sqrt <- sqrt(dataRangePrimate_rdc$Number_species_cooccurrence)
  dataRangePrimate_rdc$Overlap_average.sqrt <- sqrt(dataRangePrimate_rdc$Overlap_average)
  rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species
  
  #Readjust phylo tree
  phyloConsensus <- drop.tip(phylo,
                             phylo$tip.label[
                               which(phylo$tip.label
                                     %nin%dataRangePrimate_rdc$Species)])
  
  #Clean tree to avoid problem, save and reload
  source("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/Scripts/Functions.R")
  library(BioGeoBEARS)
  phyloConsensus <- cleanTree(tr=phyloConsensus, trfn="C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/Raw_data/Tree/phyloConsensus2.nex")
  phyloConsensus <- force.ultrametric(tree=phyloConsensus, method="extend")#method="nnls")
  
  ##MCMC glmm
  set.seed(a)
  
  library(MCMCglmm)
  
  #Create the necessary var/cov for mcmcglmm to be fitted
  inv.phylo <-inverseA(phyloConsensus,nodes="TIPS",scale=TRUE)
  
  #Compute model with 3 chains
  
  # qqplot(x, y,
  #        xlab=deparse(substitute(x)), ylab=deparse(substitute(y)),
  #        common.scale=TRUE, preserveLabels=FALSE, lwd=1,
  #        lcol="red", ...)
  # 
  
  prior<-list(G=list(G1=list(V=V,nu=nu)),R=list(V=V,nu=nu))#Neutral prior
  
  model1<-MCMCglmm(Trait ~ Overlap_average.sqrt + Number_species_cooccurrence.sqrt, , random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
                   data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
  # model2<-MCMCglmm(Trait ~ Overlap_average.sqrt + Number_species_cooccurrence.sqrt, , random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
  #                  data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
  # model3<-MCMCglmm(Trait ~ Overlap_average.sqrt + Number_species_cooccurrence.sqrt, , random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
  #                  data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
  # 
  #Check convergence with Gelman-Rubin
  #gelmanRubinValues[a] <- gelman.diag(list(model1$Sol, model2$Sol, model3$Sol))$mpsrf[1]
  #heidel.diag(model1)
  
  #HPDinterval(model1$Sol)#ok the summary corresponds to the HDP
  
  #Calculate Lambda, from: http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm
  
  lambda <- model1$VCV[,'Species']/
    (model1$VCV[,'Species']+model1$VCV[,'units'])
  meanLambda <- mean(lambda)
  modeLambda <- posterior.mode(lambda)
  CIlambda <- HPDinterval(lambda)
  
  #Complete df results
  
  #commented is when I used pgls
  results.df[(firstRowWrite+1):(firstRowWrite+3), 2] <- roundIntelligent(summary(model1)$solutions[,1], digit=2) #roundIntelligent(summary(modelBrain)$tTable[,1], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+3), 3] <- roundIntelligent(summary(model1)$solutions[,2], digit=2) #roundIntelligent(CI$coef[,1], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+3), 4] <- roundIntelligent(summary(model1)$solutions[,3], digit=2) #roundIntelligent(CI$coef[,3], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+3), 5] <- roundIntelligent(summary(model1)$solutions[,4], digit=2)#sd: useless in Bayesian// put the effective sample size instead #roundIntelligent(summary(modelBrain)$coefficients[,2], digit=2) #roundIntelligent(summary(modelBrain)$tTable[,2], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+3), 6] <- NA#t: useless in bayesian c("-", roundIntelligent(summary(modelBrain)$coefficients[2,3], digit=2))
  results.df[(firstRowWrite+1):(firstRowWrite+3), 7] <- c("-", roundIntelligent(summary(model1)$solutions[2,5], digit=2),  roundIntelligent(summary(model1)$solutions[3,5], digit=2))
  
  # val1 <- ifelse(is.na(summary(modelBrainPgls)$param.CI$lambda$ci.val[1]), 0, round(summary(modelBrainPgls)$param.CI$lambda$ci.val[1], digit=2))
  # val2 <- ifelse(is.na(summary(modelBrainPgls)$param.CI$lambda$ci.val[2]), 1, round(summary(modelBrainPgls)$param.CI$lambda$ci.val[2], digit=2))
  # 
  # lambda <- paste(round(summary(modelBrainPgls)$param.CI$lambda$opt, digit=2), " [",val1,",",val2,"]", sep="")
  #         
  results.df[firstRowWrite+4,2] <- meanLambda 
  results.df[firstRowWrite+4,3] <- CIlambda[1,1]#CI$corStruct[1,1])
  results.df[firstRowWrite+4,4] <- CIlambda[1,2]#CI$corStruct[1,3])
  
  #Check assumptions
  #assign(paste("modelBrainDiversification", traitName[a], sep="_"), model1)
  #assign(paste(modelBrainPgls, traitName[a], sep="_"), modelBrainPgls)
  
  #vifVectorDiversification[a] <- max(vif(modelBrain))
  
}

save.image("REnvironments/PGLSdirectionSelection_bayesian.RData")
