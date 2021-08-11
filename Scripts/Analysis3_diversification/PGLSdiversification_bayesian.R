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

# The frequentist approach suffered from heteroskedasticity. Thus, we favoured Bayesian computation for this analysis. #We thank V. lauret and M. Queroue for discussion on these models.

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

summaryData$ratioBrain <- summaryData$Brain/summaryData$Bodymass #Following decasien for multiplication by 1.036
summaryData$EQ <- summaryData$Brain*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775) #Following decasien, according to #Jerison, H. J. Evolution of the Brain and Intelligence (Academic, 1973).

summaryData$Brain.log <- log(summaryData$Brain)
summaryData$EQ.log <- log(summaryData$EQ)

summaryData$ratioNeocortex <- summaryData$Neocortex/summaryData$Bodymass
summaryData$ratioNeocortex.log  <- log(summaryData$ratioNeocortex)

summaryData$ratioHippocampus <- summaryData$Hippocampus/summaryData$Bodymass
summaryData$ratioHippocampus.log <- log(summaryData$ratioHippocampus)

summaryData$ratioCerebellum <- summaryData$Cerebellum/summaryData$Bodymass
summaryData$ratioCerebellum.log <- log(summaryData$ratioCerebellum)

summaryData$ratioStriatum <- summaryData$Striatum/summaryData$Bodymass
summaryData$ratioStriatum.log <- log(summaryData$ratioStriatum)

summaryData$ratioMOB <- summaryData$MOB/summaryData$Bodymass
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
f=73 # Assumed availability based on current knowledge

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
summary(summaryData$Diversification_rate)
hist(summaryData$Diversification_rate)#ok

##
# Remove Homo sapiens

dataRangePrimate <- dataRangePrimate[-grep("Homo", dataRangePrimate$Species), ] 

##----
## Run analysis for each trait
set.seed(26)
traitToStudy=c("EQ.log", "ratioBrain", "ratioHippocampus.log", "ratioNeocortex.log", "ratioCerebellum.log", "ratioStriatum.log", "ratioMOB.log")
traitName=c("EQ (log)", "Brain (/bodymass, log)", "Hippocampus (/bodymass, log)", "Neocortex (/bodymass, log)", "Cerebellum (/bodymass, log)", "Striatum (/bodymass, log)", "MOB (/bodymass, log)")
vifVector <- rep(NA, times=length(traitName))  

results.df <- as.data.frame(matrix(NA, nrow=4*length(traitToStudy), ncol=7))
colnames(results.df) <- c("", "Est.", "HDP2.5%", "HDP97.5%", "Eff. samp", "NA", "pMCMC")

gelmanRubinValues <- rep(NA, times=length(traitToStudy))

thin=50
burnin=5000
nitt=50000
V=1
nu=0.02

for(a in 1:length(traitToStudy)){
  #Matching the brain trait to the dataset with predictors
  dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[a])]
  dataRangePrimate$Diversification_rate <-  summaryData$Diversification_rate[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny)]
  dataRangePrimate$Bodymass <-  summaryData$Bodymass[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny)]
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
  #Readjust phylo tree
  phyloConsensus <- drop.tip(phylo,
                             phylo$tip.label[
                               which(phylo$tip.label
                                     %nin%dataRangePrimate_rdc$Species)])
  
  #Clean tree to avoid problem, save and reload
  source("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/Scripts/Functions.R")
  library(BioGeoBEARS)
  phyloConsensus <- cleanTree(tr=phyloConsensus, trfn="C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/Raw_data/Tree/phyloConsensus.nex")
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
  
  model1<-MCMCglmm(Diversification_rate~Trait, random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
                   data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
  model2<-MCMCglmm(Diversification_rate~Trait, random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
                   data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
  model3<-MCMCglmm(Diversification_rate~Trait, random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
                   data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
  
  #Check convergence with Gelman-Rubin
  gelmanRubinValues[a] <- gelman.diag(list(model1$Sol, model2$Sol, model3$Sol))$mpsrf[1]
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
  results.df[(firstRowWrite+1):(firstRowWrite+2), 2] <- roundIntelligent(summary(model1)$solutions[,1], digit=2) #roundIntelligent(summary(modelBrain)$tTable[,1], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+2), 3] <- roundIntelligent(summary(model1)$solutions[,2], digit=2) #roundIntelligent(CI$coef[,1], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+2), 4] <- roundIntelligent(summary(model1)$solutions[,3], digit=2) #roundIntelligent(CI$coef[,3], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+2), 5] <- roundIntelligent(summary(model1)$solutions[,4], digit=2)#sd: useless in Bayesian// put the effective sample size instead #roundIntelligent(summary(modelBrain)$coefficients[,2], digit=2) #roundIntelligent(summary(modelBrain)$tTable[,2], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+2), 6] <- NA#t: useless in bayesian c("-", roundIntelligent(summary(modelBrain)$coefficients[2,3], digit=2))
  results.df[(firstRowWrite+1):(firstRowWrite+2), 7] <- c("-", roundIntelligent(summary(model1)$solutions[2,5], digit=2))
  
  # val1 <- ifelse(is.na(summary(modelBrainPgls)$param.CI$lambda$ci.val[1]), 0, round(summary(modelBrainPgls)$param.CI$lambda$ci.val[1], digit=2))
  # val2 <- ifelse(is.na(summary(modelBrainPgls)$param.CI$lambda$ci.val[2]), 1, round(summary(modelBrainPgls)$param.CI$lambda$ci.val[2], digit=2))
  # 
  # lambda <- paste(round(summary(modelBrainPgls)$param.CI$lambda$opt, digit=2), " [",val1,",",val2,"]", sep="")
  #         
  results.df[firstRowWrite+3,2] <- meanLambda 
  results.df[firstRowWrite+3,3] <- CIlambda[1,1]#CI$corStruct[1,1])
  results.df[firstRowWrite+3,4] <- CIlambda[1,2]#CI$corStruct[1,3])
  
  #Check assumptions
  assign(paste("modelBrainDiversification", traitName[a], sep="_"), model1)
  #assign(paste(modelBrainPgls, traitName[a], sep="_"), modelBrainPgls)
  
  #vifVectorDiversification[a] <- max(vif(modelBrain))
  
}

results.df_diversification_init <- results.df
results.df[2:ncol(results.df)] <- replaceNAtable(results.df[2:ncol(results.df)], "")
results.df <- results.df[, -6] #delete column that was created for frequentist and is useless
results.df_diversification <- results.df

save.image("REnvironments/PGLSdiversification.RData")

##############
## PLOTTING RESULTS
##############

#Remaining to do in plotting:

#1) create layout + matrix for all plots

#finish the plot part:write for overlap + adjust axis and color + add regression
pdf("Plots/diversificationPGLS.pdf", height=25, width=18)

numberPlots = ifelse((length(traitToStudy)-1)/2==floor((length(traitToStudy)-1)/2), (length(traitToStudy)-1), (length(traitToStudy)-1) + 1)#if odd number, readjust

layout(mat=matrix(1:numberPlots, ncol=2),
       widths=c(5,5), heights=rep(5, times=numberPlots/2))
par(mar=c(3.5, 3.5, 1, 1), mgp=c(2.5, 0.5, 0), xpd=TRUE, cex=1.2)

countLabel=0
for(i in 1:length(traitToStudy)){ 
  if(i!=2){
    countLabel=countLabel+1
    #Matching the brain trait to the dataset with predictors
    dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[i])]
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
    #Readjust phylo tree
    phyloConsensus <- drop.tip(phylo,
                               phylo$tip.label[
                                 which(phylo$tip.label
                                       %nin%dataRangePrimate_rdc$Species)])
    
    #Clean tree to avoid problem, save and reload
    source("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/Scripts/Functions.R")
    library(BioGeoBEARS)
    phyloConsensus <- cleanTree(tr=phyloConsensus, trfn="C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/Raw_data/Tree/phyloConsensus.nex")
    phyloConsensus <- force.ultrametric(tree=phyloConsensus, method="extend")#method="nnls")
    
    ##MCMC glmm
    set.seed(a)
    
    library(MCMCglmm)
    
    #Create the necessary var/cov for mcmcglmm to be fitted
    inv.phylo <-inverseA(phyloConsensus,nodes="TIPS",scale=TRUE)
    
    #Compute model with 3 chains
    
    prior<-list(G=list(G1=list(V=V,nu=nu)),R=list(V=V,nu=nu))#Neutral prior
    
    model1<-MCMCglmm(Diversification_rate~Trait, random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
                     data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
    
    
    CI <- as.data.frame(summary(model1)$solutions[,1:3])
    
    ##----
    #Plot against N co-occ
    
    xmin=min(round(dataRangePrimate_rdc[,c(4)], digit=2))
    xmax=max(round(dataRangePrimate_rdc[,c(4)], digit=2))
    ymin=0#min(round(dataRangePrimate_rdc[,c(5)], digit=2))
    ymax=max(round(dataRangePrimate_rdc[,c(5)], digit=2))
    
    ##With number of co-occurring species
    plot(dataRangePrimate_rdc[,c(4)], dataRangePrimate_rdc[,c(5)], xlab=traitName[i], ylab="Diversification rate",
         font.lab=2, cex.lab=1.25,
         xlim=c(xmin, xmax), ylim=c(ymin, ymax),
         las=1, type="n", tcl=-0.25, bty="n",
         xaxt="n",xaxs="i",yaxs="i", yaxt="n", xpd=TRUE)
    
    #Add grid
    addGrid(
      cexAxisX=1.15, cexAxisY=1.15,
      xmin=xmin, xmax=xmax, xintsmall=(xmax-xmin)/20, xintbig=(xmax-xmin)/5,
      ymin=ymin, ymax=ymax, yintsmall=(ymax-ymin)/20, yintbig=(ymax-ymin)/5,
      axisPlot=TRUE, round=TRUE, digit=c(2,2))
    axis(side=1, at=round(seq(from=xmin, to=xmax, by=(xmax-xmin)/5), digit=1), labels=round(seq(from=xmin, to=xmax, by=(xmax-xmin)/5), digit=1), las=1, tcl=-0.25)
    addLabel(xfrac=0.05, yfrac=0.05, label=countLabel, circle=TRUE, radiuscircle=(xmax-xmin)/35, circle.bg="black", font.col="white")
    
    #Add result model
    ymean <- c(CI[1,1], CI[1,1] + CI[2,1]*max(dataRangePrimate_rdc[,c(4)]))
    ylower <- c(CI[1,2], CI[1,2] + (CI[2,2])*max(dataRangePrimate_rdc[,c(4)]))
    yupper <- c(CI[1,3], CI[1,3] + (CI[2,3])*max(dataRangePrimate_rdc[,c(4)]))
    lines(
      x=c(min(dataRangePrimate_rdc[,c(4)]),max(dataRangePrimate_rdc[,c(4)])),
      y=ymean
    )
    polygon(
      x=c(min(dataRangePrimate_rdc[,c(4)]), max(dataRangePrimate_rdc[,c(4)]), max(dataRangePrimate_rdc[,c(4)]), min(dataRangePrimate_rdc[,c(4)]), min(dataRangePrimate_rdc[,c(4)])),
      y=c(yupper, rev(ylower), ylower[1]),
      col=adjustcolor("black", alpha.f=0.1),
      border=NA#,
      #lty=2
    )
    
    #Add background tree
    col=list(col.edge=setNames(rep("darkgrey",nrow(phyloConsensus$edge)),as.character(phyloConsensus$edge[,2])),
             col.node=setNames(rep("black",max(phyloConsensus$edge)),as.character(1:max(phyloConsensus$edge))))
    
    phylomorphospace(phyloConsensus,dataRangePrimate_rdc[,c(4,5)], add=TRUE, label="true", lty=3,
                     control=col, xpd=TRUE)
  }
}

dev.off()

save.image("REnvironments/PGLSdiversification.RData")

#Calculate Dfbetas

dfBetasEstimate <- as.data.frame(matrix(NA, ncol=2, nrow=nrow(dataRangePrimate)*length(traitToStudy)))
dfBetasPvalue <- as.data.frame(matrix(NA, ncol=2, nrow=nrow(dataRangePrimate)*length(traitToStudy)))
dfBetasLambda <- rep(NA, times=nrow(dataRangePrimate)*length(traitToStudy))
traitVectordfBetas <- rep(NA, times=nrow(dataRangePrimate)*length(traitToStudy))
gelmanRubinValuesdfBetas <- rep(NA, times=nrow(dataRangePrimate)*length(traitToStudy))

for(a in 1:length(traitToStudy)){
  
  #Matching the brain trait to the dataset with predictors
  dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[a])]
  dataRangePrimate$Diversification_rate <-  summaryData$Diversification_rate[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny)]
  dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,4])&!is.na(dataRangePrimate[,5]),]
  
  for(r in 1:nrow(dataRangePrimate_rdc)){
    rowToAdd=which(is.na(dfBetasLambda))[1]
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
    #Readjust phylo tree
    phyloConsensus <- drop.tip(phylo,
                               phylo$tip.label[
                                 which(phylo$tip.label
                                       %nin%dataRangePrimate_rdc$Species)])
    
    #Clean tree to avoid problem, save and reload
    source("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/Scripts/Functions.R")
    library(BioGeoBEARS)
    phyloConsensus <- cleanTree(tr=phyloConsensus, trfn="C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/Raw_data/Tree/phyloConsensus.nex")
    phyloConsensus <- force.ultrametric(tree=phyloConsensus, method="extend")#method="nnls")
    
    ##MCMC glmm
    set.seed(a)
    
    library(MCMCglmm)
    
    #Create the necessary var/cov for mcmcglmm to be fitted
    inv.phylo <-inverseA(phyloConsensus,nodes="TIPS",scale=TRUE)
    
    #Compute model with 3 chains
    
    prior<-list(G=list(G1=list(V=V,nu=nu)),R=list(V=V,nu=nu))#Neutral prior
    
    model1<-MCMCglmm(Diversification_rate~Trait, random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
                     data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
    # model2<-MCMCglmm(Diversification_rate~Trait, random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
    #                  data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
    # model3<-MCMCglmm(Diversification_rate~Trait, random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
    #                  data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
    
    #Check convergence with Gelman-Rubin
    #gelmanRubinValuesdfBetas[rowToAdd] <- gelman.diag(list(model1$Sol, model2$Sol, model3$Sol))$mpsrf[1]
    #heidel.diag(model1)
    
    #HPDinterval(model1$Sol)#ok the summary corresponds to the HDP
    
    #Calculate Lambda, from: http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm
    
    lambda <- model1$VCV[,'Species']/
      (model1$VCV[,'Species']+model1$VCV[,'units'])
    meanLambda <- mean(lambda)
    modeLambda <- posterior.mode(lambda)
    CIlambda <- HPDinterval(lambda)
    
    
    #Complete df results
    dfBetasEstimate[rowToAdd,] <- summary(model1)$solutions[,1]
    dfBetasPvalue[rowToAdd,] <- summary(model1)$solutions[,5]
    dfBetasLambda[rowToAdd] <- meanLambda
    traitVectordfBetas[rowToAdd] <- traitName[a]
  }
}

save.image("REnvironments/PGLSdiversification.RData")


repetitionTrees=50
repetitionModels=50

sensitivityEstimate <- as.data.frame(matrix(NA, ncol=2, nrow=repetitionTrees*repetitionModels*length(traitToStudy)))
sensitivityPvalue <- as.data.frame(matrix(NA, ncol=2, nrow=repetitionTrees*repetitionModels*length(traitToStudy)))
sensitivityLambda <- rep(NA, times=repetitionTrees*repetitionModels*length(traitToStudy))
traitVectorSensitivity <- rep(NA, times=repetitionTrees*repetitionModels*length(traitToStudy))
gelmanRubinValuesSensitivity <- rep(NA, times=repetitionTrees*repetitionModels*length(traitToStudy))

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
  
  summaryData$ratioBrain <- summaryData$Brain/summaryData$Bodymass #Following decasien for multiplication by 1.036
  summaryData$EQ <- summaryData$Brain*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775) #Following decasien, according to #Jerison, H. J. Evolution of the Brain and Intelligence (Academic, 1973).
  
  summaryData$Brain.log <- log(summaryData$Brain)
  summaryData$EQ.log <- log(summaryData$EQ)
  
  summaryData$ratioNeocortex <- summaryData$Neocortex/summaryData$Bodymass
  summaryData$ratioNeocortex.log  <- log(summaryData$ratioNeocortex)
  
  summaryData$ratioHippocampus <- summaryData$Hippocampus/summaryData$Bodymass
  summaryData$ratioHippocampus.log <- log(summaryData$ratioHippocampus)
  
  summaryData$ratioCerebellum <- summaryData$Cerebellum/summaryData$Bodymass
  
  summaryData$ratioStriatum <- summaryData$Striatum/summaryData$Bodymass
  
  summaryData$ratioMOB <- summaryData$MOB/summaryData$Bodymass
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
      #Readjust phylo tree
      phyloConsensus <- drop.tip(phylo,
                                 phylo$tip.label[
                                   which(phylo$tip.label
                                         %nin%dataRangePrimate_rdc$Species)])
      
      #Clean tree to avoid problem, save and reload
      source("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/Scripts/Functions.R")
      library(BioGeoBEARS)
      phyloConsensus <- cleanTree(tr=phyloConsensus, trfn="C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/Raw_data/Tree/phyloConsensus.nex")
      phyloConsensus <- force.ultrametric(tree=phyloConsensus, method="extend")#method="nnls")
      
      ##MCMC glmm
      set.seed(a)
      
      library(MCMCglmm)
      
      #Create the necessary var/cov for mcmcglmm to be fitted
      inv.phylo <-inverseA(phyloConsensus,nodes="TIPS",scale=TRUE)
      
      #Compute model with 3 chains
      
      prior<-list(G=list(G1=list(V=V,nu=nu)),R=list(V=V,nu=nu))#Neutral prior
      
      model1<-MCMCglmm(Diversification_rate~Trait, random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
                       data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
      # model2<-MCMCglmm(Diversification_rate~Trait, random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
      #                  data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
      # model3<-MCMCglmm(Diversification_rate~Trait, random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
      #                  data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
      
      #Check convergence with Gelman-Rubin
      #gelmanRubinValuesSensitivity[rowToAdd] <- gelman.diag(list(model1$Sol, model2$Sol, model3$Sol))$mpsrf[1]
      #heidel.diag(model1)
      
      #HPDinterval(model1$Sol)#ok the summary corresponds to the HDP
      
      #Calculate Lambda, from: http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm
      
      lambda <- model1$VCV[,'Species']/
        (model1$VCV[,'Species']+model1$VCV[,'units'])
      meanLambda <- mean(lambda)
      modeLambda <- posterior.mode(lambda)
      CIlambda <- HPDinterval(lambda)
      
      
      #Complete df results
      sensitivityEstimate[rowToAdd,] <- summary(model1)$solutions[,1]
      sensitivityPvalue[rowToAdd,] <- summary(model1)$solutions[,5]
      sensitivityLambda[rowToAdd] <- meanLambda
      traitVectorSensitivity[rowToAdd] <- traitName[a]
    }
  }
  
  #progress(d/repetitionModels*100)
}

save.image("REnvironments/PGLSdiversification.RData")


fraction.v <- c(60, 70, 80, 90, 95)

sensitivityFractionTaxonEstimate <- as.data.frame(matrix(NA, ncol=2, nrow=length(fraction.v)*length(traitToStudy)))
sensitivityFractionTaxonPvalue <- as.data.frame(matrix(NA, ncol=2, nrow=length(fraction.v)*length(traitToStudy)))
sensitivityFractionTaxonLambda <- rep(NA, times=length(fraction.v)*length(traitToStudy))
traitVectorSensitivityFractionTaxon <- rep(NA, times=length(fraction.v)*length(traitToStudy))
gelmanRubinValuesFractionTaxon <- rep(NA, times=length(fraction.v)*length(traitToStudy))

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
    #Readjust phylo tree
    phyloConsensus <- drop.tip(phylo,
                               phylo$tip.label[
                                 which(phylo$tip.label
                                       %nin%dataRangePrimate_rdc$Species)])
    
    #Clean tree to avoid problem, save and reload
    source("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/Scripts/Functions.R")
    library(BioGeoBEARS)
    phyloConsensus <- cleanTree(tr=phyloConsensus, trfn="C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/Raw_data/Tree/phyloConsensus.nex")
    phyloConsensus <- force.ultrametric(tree=phyloConsensus, method="extend")#method="nnls")
    
    ##MCMC glmm
    set.seed(a)
    
    library(MCMCglmm)
    
    #Create the necessary var/cov for mcmcglmm to be fitted
    inv.phylo <-inverseA(phyloConsensus,nodes="TIPS",scale=TRUE)
    
    #Compute model with 3 chains
    
    prior<-list(G=list(G1=list(V=V,nu=nu)),R=list(V=V,nu=nu))#Neutral prior
    
    model1<-MCMCglmm(Diversification_rate~Trait, random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
                     data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
    # model2<-MCMCglmm(Diversification_rate~Trait, random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
    #                  data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
    # model3<-MCMCglmm(Diversification_rate~Trait, random=~Species, ginverse=list(Species=inv.phylo$Ainv), family="gaussian",
    #                  data=dataRangePrimate_rdc, nitt=nitt,burnin=burnin,thin=thin, prior=prior)
    
    #Check convergence with Gelman-Rubin
    #gelmanRubinValuesFractionTaxon[rowToAdd] <- gelman.diag(list(model1$Sol, model2$Sol, model3$Sol))$mpsrf[1]
    #heidel.diag(model1)
    
    #HPDinterval(model1$Sol)#ok the summary corresponds to the HDP
    
    #Calculate Lambda, from: http://www.mpcm-evolution.com/practice/online-practical-material-chapter-11/chapter-11-1-simple-model-mcmcglmm
    
    lambda <- model1$VCV[,'Species']/
      (model1$VCV[,'Species']+model1$VCV[,'units'])
    meanLambda <- mean(lambda)
    modeLambda <- posterior.mode(lambda)
    CIlambda <- HPDinterval(lambda)
    
    
    #Complete df results
    sensitivityFractionTaxonEstimate[rowToAdd,] <- summary(model1)$solutions[,1]
    sensitivityFractionTaxonPvalue[rowToAdd,] <- summary(model1)$solutions[,5]
    sensitivityFractionTaxonLambda[rowToAdd] <- meanLambda
    traitVectorSensitivityFractionTaxon[rowToAdd] <- traitName[a]
    
  }
}

save.image("REnvironments/PGLSdiversification.RData")

###----
# Create sensitivity tables

#Summary sensitivity

#Dfbetas
library(tidyr)
minEst <- tibble(cbind(aggregate(dfBetasEstimate, by=list(traitVectordfBetas), min), aggregate(dfBetasLambda, by=list(traitVectordfBetas), min)[,2]))
colnames(minEst) <- c("Model", "Intercept", "Trait", "Lambda")
minEst <- pivot_longer(minEst, col=2:4, names_to="Variable", values_to="Estimate")

maxEst <- tibble(cbind(aggregate(dfBetasEstimate, by=list(traitVectordfBetas), max), aggregate(dfBetasLambda, by=list(traitVectordfBetas), max)[,2]))
colnames(maxEst) <- c("Model", "Intercept", "Trait", "Lambda")
maxEst <- pivot_longer(maxEst, col=2:4, names_to="Variable", values_to="Estimate")

estimate <- results.df[results.df[,2]!="",2]
whatIs <- rep(traitName, each=3)
numberForFinishingMatch <- rep(1:3, times=length(traitName))

#medianEst <- aggregate(dfBetasEstimate, by=list(traitVectordfBetas), median)
minEst$EstimateTrue <- estimate[match(paste(minEst$Model,numberForFinishingMatch),paste(whatIs,numberForFinishingMatch))]

dfBetasDiversification <- cbind(minEst,
                                maxEst[,3])
colnames(dfBetasDiversification) <- c("Model", "Variable", "Est. min.", "Est.", "Est. max.")


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
colnames(minEst) <- c("Model", "Intercept", "Trait", "Lambda")
minEst <- pivot_longer(minEst, col=2:4, names_to="Variable", values_to="Estimate")

maxEst <- tibble(cbind(aggregate(sensitivityEstimate, by=list(traitVectorSensitivity), max), aggregate(sensitivityLambda, by=list(traitVectorSensitivity), max)[,2]))
colnames(maxEst) <- c("Model", "Intercept", "Trait", "Lambda")
maxEst <- pivot_longer(maxEst, col=2:4, names_to="Variable", values_to="Estimate")

estimate <- results.df[results.df[,2]!="",2]
whatIs <- rep(traitName, each=3)
numberForFinishingMatch <- rep(1:3, times=length(traitName))

#medianEst <- aggregate(sensitivityEstimate, by=list(traitVectorSensitivity), median)
minEst$EstimateTrue <- estimate[match(paste(minEst$Model,numberForFinishingMatch),paste(whatIs,numberForFinishingMatch))]

sensitivityDiversification <- cbind(minEst,
                                    maxEst[,3])
colnames(sensitivityDiversification) <- c("Model", "Variable", "Est. min.", "Est.", "Est. max.")

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
minEst <- tibble(cbind(aggregate(sensitivityFractionTaxonEstimate, by=list(traitVectorSensitivityFractionTaxon), min), aggregate(sensitivityFractionTaxonLambda, by=list(traitVectorSensitivityFractionTaxon), min)[,2]))
colnames(minEst) <- c("Model", "Intercept", "Trait", "Lambda")
minEst <- pivot_longer(minEst, col=2:4, names_to="Variable", values_to="Estimate")

maxEst <- tibble(cbind(aggregate(sensitivityFractionTaxonEstimate, by=list(traitVectorSensitivityFractionTaxon), max), aggregate(sensitivityFractionTaxonLambda, by=list(traitVectorSensitivityFractionTaxon), max)[,2]))
colnames(maxEst) <- c("Model", "Intercept", "Trait", "Lambda")
maxEst <- pivot_longer(maxEst, col=2:4, names_to="Variable", values_to="Estimate")

estimate <- results.df[results.df[,2]!="",2]
whatIs <- rep(traitName, each=3)
numberForFinishingMatch <- rep(1:3, times=length(traitName))

#medianEst <- aggregate(sensitivityFractionTaxonEstimate, by=list(traitVectorsensitivityFractionTaxon), median)
minEst$EstimateTrue <- estimate[match(paste(minEst$Model,numberForFinishingMatch),paste(whatIs,numberForFinishingMatch))]

sensitivityFractionTaxonDiversification <- cbind(minEst,
                                                 maxEst[,3])
colnames(sensitivityFractionTaxonDiversification) <- c("Model", "Variable", "Est. min.", "Est.", "Est. max.")

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
summarySensitivityDiversification$Model[-seq(from=1, to=nrow(summarySensitivityDiversification), by=3)] <- ""

knitr::kable(summarySensitivityDiversification, escape=TRUE, booktabs = TRUE,
             caption = "Sensitivity analysis of phylogenetic regressions to detect the assess the diversification pattern | Depicted is the minimum and maximum of estimates when one observation was removed at a time (DfBetas), when varying the used phylogenetic tree and the data sampling (Phylogeny/Data), or when the sampling fraction varied (Sampling fraction)") %>% 
  kableExtra::kable_styling(latex_options = "striped") %>%
  kableExtra::kable_styling(latex_options="scale_down") %>%
  kableExtra::add_header_above(c("Model:" = 2, "DfBetas" = 3, "Phylogeny/Data" = 3, "Sampling fraction" = 3))

###---
## Create main table

whichToBold <- which(as.numeric(results.df_diversification$"pMCMC")<=0.05)
toPlotBoldDiversification <- rep(FALSE, times=nrow(results.df_diversification))
toPlotBoldDiversification[whichToBold] <- TRUE

replace <- sapply(as.numeric(results.df_diversification$"pMCMC"), function(x) if(!is.na(x)&x!=""&x!="-"){pvalueRound(x, text=FALSE)}else{x})
results.df_diversification$"pMCMC"[!is.na(replace)] <- replace[!is.na(replace)]

knitr::kable(results.df_diversification, escape=TRUE, booktabs = TRUE,
             caption = "Model estimates and significance of Bayesian phylogenetic regressions to assess the diversification pattern | Est.=Estimate, HDP2.5%=Lower border of the 95% Highest Posterior Density, HDP97.5%=Upper border of the 95% Highest Posterior Density, Eff. samp.=Effective sample (adjusted for autocorrelation). The brain area (as well as the associated sample size)
             are indicate prior each list of estimates. the transformation (logarithm or square-root) if indicated in parenthese by the abbreviation (log or sqrt).") %>%
  kableExtra::column_spec(2:ncol(results.df_diversification), bold = toPlotBoldDiversification) %>%
  kableExtra::kable_styling(latex_options = "striped") %>%
  kableExtra::kable_styling(latex_options="scale_down")

save.image("REnvironments/PGLSdiversification.RData")

load("REnvironments/PGLSdiversification.RData")
library(coda)
valueAutoCorr <- rep(NA, times=length(traitToStudy))
for(a in 1:length(traitToStudy)){
modeltest <- get(paste("modelBrainDiversification", traitName[a], sep="_"))
corr <- coda::autocorr.diag(modeltest$Sol)
corr <- abs(corr)
valueAutoCorr[a] <- max(apply(corr[-1,], 2, max))
}

save.image("REnvironments/PGLSdiversification_withautocorr.RData")