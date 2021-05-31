####---------------------------------
#### PGLSdirectionSelection
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

summaryData$ratioNeocortex <- summaryData$Neocortex*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775)#summaryData$Neocortex/summaryData$Brain
summaryData$ratioNeocortex.log  <- log(summaryData$ratioNeocortex)

summaryData$ratioHippocampus <- summaryData$Hippocampus*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775)#summaryData$Hippocampus/summaryData$Brain
summaryData$ratioHippocampus.log <- log(summaryData$ratioHippocampus)

summaryData$ratioCerebellum <- summaryData$Cerebellum*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775)#summaryData$Cerebellum/summaryData$Brain

summaryData$ratioStriatum <- summaryData$Striatum*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775)#summaryData$Striatum/summaryData$Brain

summaryData$ratioMOB <- summaryData$MOB*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775)#summaryData$MOB/summaryData$Brain
summaryData$ratioMOB.log <- log(summaryData$ratioMOB)


# summaryData$ratioBrain <- summaryData$Brain*1.036*(10**-3)/summaryData$Bodymass #Following decasien for multiplication by 1.036
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

summaryData$ratioMOB <- summaryData$MOB
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

traitToStudy=c("EQ.log", "ratioBrain", "ratioHippocampus.log", "ratioNeocortex.log", "ratioCerebellum", "ratioStriatum", "ratioMOB.log")
traitName=c("EQ (log)", "Brain", "Hippocampus (log)", "Neocortex (log)", "Cerebellum", "Striatum", "MOB (log)")
vifVector <- rep(NA, times=length(traitName))  
  
results.df <- as.data.frame(matrix(NA, nrow=5*length(traitToStudy), ncol=7))
colnames(results.df) <- c("", "Est.", "CI2.5%", "CI97.5%", "Sd", "t", "p-value")

##----
## Run analysis for each trait

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


#Create data for phylogenetic regression
# comp_data <- comparative.data(phy = phyloConsensus, data= dataRangePrimate_rdc,
#                               names.col = Species, vcv = TRUE)

# #Output distribution
# hist(comp_data$data$Trait)
# 
# #Predictor distribution
# hist(sqrt(comp_data$data$Number_species_cooccurrence))
# hist(sqrt(comp_data$data$Overlap_average))

# #Make predictor more symmetrical
# comp_data$data$Number_species_cooccurrence.sqrt <- sqrt(comp_data$data$Number_species_cooccurrence)
# comp_data$data$Overlap_average.sqrt <- sqrt(comp_data$data$Overlap_average)
# 
# modelBrainPgls <- pgls(formula = Trait ~ Overlap_average.sqrt + Number_species_cooccurrence.sqrt, 
#                        data=comp_data, lambda="ML",param.CI=0.05)

dataRangePrimate_rdc$Number_species_cooccurrence.sqrt <- sqrt(dataRangePrimate_rdc$Number_species_cooccurrence)
dataRangePrimate_rdc$Overlap_average.sqrt <- sqrt(dataRangePrimate_rdc$Overlap_average)
rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species

# plot(dataRangePrimate_rdc$Number_species_cooccurrence.sqrt, dataRangePrimate_rdc$Overlap_average.sqrt)
     
modelBrain <- phylolm(formula = Trait ~ Overlap_average.sqrt + Number_species_cooccurrence.sqrt, 
                       data=dataRangePrimate_rdc ,phy=phyloConsensus,model="lambda",measurement_error=FALSE,boot=repetitionBootstrap)

# modelBrain <- gls(Trait ~ Overlap_average.sqrt  correlation = corPagel(1, phyloConsensus, form = ~Species),
#                   data = dataRangePrimate_rdc, method = "ML")
# CI <- intervals(modelBrain)
# results <- anova(modelBrain)

#Complete df results

#commented is when I used pgls + phylolm
results.df[(firstRowWrite+1):(firstRowWrite+3), 2] <- roundIntelligent(summary(modelBrain)$coefficients[,1], digit=2) #roundIntelligent(summary(modelBrain)$tTable[,1], digit=2)
results.df[(firstRowWrite+1):(firstRowWrite+3), 3] <- roundIntelligent(summary(modelBrain)$coefficients[,4], digit=2) #roundIntelligent(CI$coef[,1], digit=2)
results.df[(firstRowWrite+1):(firstRowWrite+3), 4] <- roundIntelligent(summary(modelBrain)$coefficients[,5], digit=2) #roundIntelligent(CI$coef[,3], digit=2)
results.df[(firstRowWrite+1):(firstRowWrite+3), 5] <- roundIntelligent(summary(modelBrain)$coefficients[,3], digit=2) #roundIntelligent(summary(modelBrain)$tTable[,2], digit=2)
results.df[(firstRowWrite+1):(firstRowWrite+3), 6] <- c("-", roundIntelligent(summary(modelBrain)$coefficients[2:3,3], digit=2))
results.df[(firstRowWrite+1):(firstRowWrite+3), 7] <- c("-", roundIntelligent(summary(modelBrain)$coefficients[2:3,6], digit=2))

# val1 <- ifelse(is.na(summary(modelBrainPgls)$param.CI$lambda$ci.val[1]), 0, round(summary(modelBrainPgls)$param.CI$lambda$ci.val[1], digit=2))
# val2 <- ifelse(is.na(summary(modelBrainPgls)$param.CI$lambda$ci.val[2]), 1, round(summary(modelBrainPgls)$param.CI$lambda$ci.val[2], digit=2))
# 
# lambda <- paste(round(summary(modelBrainPgls)$param.CI$lambda$opt, digit=2), " [",val1,",",val2,"]", sep="")
#         
results.df[firstRowWrite+4,2] <- roundIntelligent(summary(modelBrain)$bootmean[5])
results.df[firstRowWrite+4,3] <- roundIntelligent(summary(modelBrain)$bootconfint95[1,5])#CI$corStruct[1,1])
results.df[firstRowWrite+4,4] <- roundIntelligent(summary(modelBrain)$bootconfint95[2,5])#CI$corStruct[1,3])

#Check assumptions
# pdf(file=paste("Figure/Diagnostics", traitName[a], ".pdf", sep=""))
assign(paste("modelBrain", traitName[a], sep="_"), modelBrain)
#assign(paste(modelBrainPgls, traitName[a], sep="_"), modelBrainPgls)

#since only two variables
cortest <- cor.test(dataRangePrimate_rdc$Number_species_cooccurrence.sqrt, dataRangePrimate_rdc$Overlap_average.sqrt)

vifVector[a] <- 1/(1-cortest$estimate**2)#max(vif(modelBrain))
#diagnostics.plot()
# dev.off()
}

results.df[2:ncol(results.df)] <- replaceNAtable(results.df[2:ncol(results.df)], "")
results.df_gradient <- results.df

######
## DFbetas

#Calculate Dfbetas

dfBetasEstimate <- as.data.frame(matrix(NA, ncol=3, nrow=nrow(dataRangePrimate)*length(traitToStudy)))
dfBetasPvalue <- as.data.frame(matrix(NA, ncol=3, nrow=nrow(dataRangePrimate)*length(traitToStudy)))
dfBetasLambda <- rep(NA, times=nrow(dataRangePrimate)*length(traitToStudy))
traitVectordfBetas <- rep(NA, times=nrow(dataRangePrimate)*length(traitToStudy))

library(svMisc)
for(a in 1:length(traitToStudy)){
  #Matching the brain trait to the dataset with predictors
  dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[a])]
  dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,2])&!is.na(dataRangePrimate[,3])&!is.na(dataRangePrimate[,4]),]
  
  #Readjust phylo tree
  phyloConsensus <- drop.tip(phylo,
                             phylo$tip.label[
                               which(phylo$tip.label
                                     %nin%dataRangePrimate_rdc$Species)])
  
  for(r in 1:nrow(dataRangePrimate_rdc)){
    dataRangePrimate_rdc$Number_species_cooccurrence.sqrt <- sqrt(dataRangePrimate_rdc$Number_species_cooccurrence)
    dataRangePrimate_rdc$Overlap_average.sqrt <- sqrt(dataRangePrimate_rdc$Overlap_average)
    rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species

    dataRangePrimate_rdc2 <- dataRangePrimate_rdc [-r,] 
    #Readjust phylo tree
    phyloConsensus <- drop.tip(phylo,
                               phylo$tip.label[
                                 which(phylo$tip.label
                                       %nin%dataRangePrimate_rdc2$Species)])

    ##--------
    ## Main model: based on consensus tree
    ##--------
    
    rowToAdd=which(is.na(dfBetasLambda))[1]

    rownames(dataRangePrimate_rdc2) <- dataRangePrimate_rdc2$Species

    modelBrain <- phylolm(Trait ~ Overlap_average.sqrt + Number_species_cooccurrence.sqrt, phy=phyloConsensus, data = dataRangePrimate_rdc2, model = "lambda", measurement_error=FALSE)
   
   dfBetasEstimate[rowToAdd,] <- summary(modelBrain)$coefficients[,1]
   dfBetasPvalue[rowToAdd,] <- summary(modelBrain)$coefficients[,4]
   dfBetasLambda[rowToAdd] <- summary(modelBrain)$optpar
   traitVectordfBetas[rowToAdd] <- traitToStudy[a]

  progress(r/nrow(dataRangePrimate_rdc)*100)
  }
}

##############
## PLOTTING RESULTS
##############

#Remaining to do in plotting:

#1) create layout + matrix for all plots

#finish the plot part:write for overlap + adjust axis and color + add regression

layout(mat=rbind(1:length(traitToStudy), 
                 (length(traitToStudy)+1):(length(traitToStudy)+length(traitToStudy))),
       widths=c(5,5), heights=rep(5, times=length(traitToStudy)))
par(mar=c(3, 3, 1, 1), mgp=c(2, 0.5, 0), xpd=TRUE)

for(i in 1:length(traitToStudy)){ 
  
  #Matching the brain trait to the dataset with predictors
  dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[i])]
  dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,2])&!is.na(dataRangePrimate[,3])&!is.na(dataRangePrimate[,4]),]
  rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species
  
  #Readjust phylo tree
  phyloConsensus <- drop.tip(phylo,
                             phylo$tip.label[
                               which(phylo$tip.label
                                     %nin%dataRangePrimate_rdc$Species)])
  
  # #Create data for phylogenetic regression
  # comp_data <- comparative.data(phy = phyloConsensus, data= dataRangePrimate_rdc,
  #                               names.col = Species, vcv = TRUE)
  # 
  # #Make predictor more symmetrical
  # comp_data$data$Number_species_cooccurrence.sqrt <- sqrt(comp_data$data$Number_species_cooccurrence)
  # comp_data$data$Overlap_average.sqrt <- sqrt(comp_data$data$Overlap_average)
  # 
  # modelBrain <- pgls(formula = Trait ~ Overlap_average.sqrt + Number_species_cooccurrence.sqrt , data = comp_data, lambda = "ML")
  # 
  dataRangePrimate_rdc$Number_species_cooccurrence.sqrt <- sqrt(dataRangePrimate_rdc$Number_species_cooccurrence)
  dataRangePrimate_rdc$Overlap_average.sqrt <- sqrt(dataRangePrimate_rdc$Overlap_average)
  
  dataRangePrimate_rdc$Number_species_cooccurrence.sqrt.z <- scale(dataRangePrimate_rdc$Number_species_cooccurrence.sqrt)
  dataRangePrimate_rdc$Overlap_average.sqrt.z <- scale(dataRangePrimate_rdc$Overlap_average.sqrt)
  
  modelBrain <- phylolm(formula = Trait ~ Overlap_average.sqrt + Number_species_cooccurrence.sqrt, 
                        data=dataRangePrimate_rdc, phy=phyloConsensus, model="lambda", measurement_error=FALSE, boot=repetitionBootstrap)
  
  CI <- cbind(modelBrain$bootmean, t(summary(modelBrain)$bootconfint95))
    
  ##----
  #Plot against N co-occ
  
  xmin=min(round(dataRangePrimate_rdc[,c(2)]), digit=2)
  xmax=max(round(dataRangePrimate_rdc[,c(2)]), digit=2)
  ymin=min(round(dataRangePrimate_rdc[,c(4)]), digit=2)
  ymax=max(round(dataRangePrimate_rdc[,c(4)]), digit=2)
  
  ##With number of co-occurring species
  plot(dataRangePrimate_rdc[,c(2)], dataRangePrimate_rdc[,c(4)], xlab="Number of sympatric\nfrugivorous species", ylab=traitName[i],
       las=1, type="n", tcl=-0.25, bty="n",
       xaxt="n",xaxs="i",yaxs="i", yaxt="n", xpd=TRUE)
  
  #Add grid
  addGrid(
    xmin=xmin, xmax=xmax, xintsmall=(xmax-xmin)/20, xintbig=(xmax-xmin)/5,
    ymin=ymin, ymax=ymax, yintsmall=(ymax-ymin)/20, yintbig=(ymax-ymin)/5,
    axisPlot=TRUE, round=TRUE, digit=c(2,2))
  addLabel(xfrac=0.05, yfrac=0.05, label=paste(i, "a", sep=""), circle=TRUE, radiuscircle=(xmax-xmin)/25, circle.bg="black", font.col="white")
  
  #Add result model
  ymean <- c(CI[1,1], CI[1,1] + CI[3,1]*max(dataRangePrimate_rdc[,c(4)]))
  ylower <- c(CI[1,2], CI[1,2] + (CI[3,2])*max(dataRangePrimate_rdc[,c(4)]))
  yupper <- c(CI[1,3], CI[1,3] + (CI[3,3])*max(dataRangePrimate_rdc[,c(4)]))
  lines(
    x=c(0,max(dataRangePrimate_rdc[,c(2)])),
    y=ymean
  )
  polygon(
    x=c(0, max(dataRangePrimate_rdc[,c(2)]), max(dataRangePrimate_rdc[,c(2)]), 0, 0),
    y=c(yupper, rev(ylower), ylower[1]),
    col=adjustcolor("black", alpha.f=0.1),
    border=NA#,
    #lty=2
  )
  
  #Add background tree
  col=list(col.edge=setNames(rep("darkgrey",nrow(phyloConsensus$edge)),as.character(phyloConsensus$edge[,2])),
           col.node=setNames(rep("black",max(phyloConsensus$edge)),as.character(1:max(phyloConsensus$edge))))
  
  phylomorphospace(phyloConsensus,dataRangePrimate_rdc[,c(2,4)], add=TRUE, label="true", lty=3,
                   control=col, xpd=TRUE)

  ##----
  #Plot against overlap
  
  modelBrain <- phylolm(formula = Trait ~ Overlap_average.sqrt + Number_species_cooccurrence.sqrt, 
                        data=dataRangePrimate_rdc ,phy=phyloConsensus,model="lambda",measurement_error=FALSE,boot=repetitionBootstrap)
  
  CI <- cbind(modelBrain$bootmean, t(summary(modelBrain)$bootconfint95))
  
  xmin=min(round(dataRangePrimate_rdc[,c(3)]), digit=2)
  xmax=max(round(dataRangePrimate_rdc[,c(3)]), digit=2)
  ymin=min(round(dataRangePrimate_rdc[,c(4)]), digit=2)
  ymax=max(round(dataRangePrimate_rdc[,c(4)]), digit=2)
  
  ##With overlap
  plot(dataRangePrimate_rdc[,c(3)], dataRangePrimate_rdc[,c(4)], xlab="Average surfacic overlap\nwith sympatric frugivorous species", ylab=traitName[i],
       las=1, type="n", tcl=-0.25, bty="n",
       xaxt="n",xaxs="i",yaxs="i", yaxt="n", xpd=TRUE)
  
  #Add grid
  addGrid(
    xmin=xmin, xmax=xmax, xintsmall=(xmax-xmin)/20, xintbig=(xmax-xmin)/5,
    ymin=ymin, ymax=ymax, yintsmall=(ymax-ymin)/20, yintbig=(ymax-ymin)/5,
    axisPlot=TRUE, round=TRUE, digit=c(2,2))
  addLabel(xfrac=0.05, yfrac=0.05, label=paste(i, "b", sep=""), circle=TRUE, radiuscircle=(xmax-xmin)/25, circle.bg="black", font.col="white")
  
  #Add background tree
  col=list(col.edge=setNames(rep("darkgrey",nrow(phyloConsensus$edge)),as.character(phyloConsensus$edge[,2])),
           col.node=setNames(rep("black",max(phyloConsensus$edge)),as.character(1:max(phyloConsensus$edge))))
  
  phylomorphospace(phyloConsensus,dataRangePrimate_rdc[,c(2,4)], add=TRUE, label="true", lty=3,
                   control=col)
  
  #Add result model
  ymean <- c(CI[1,1], CI[1,1] + CI[3,1]*max(dataRangePrimate_rdc[,c(4)]))
  ylower <- c(CI[1,2], CI[1,2] + (CI[3,2])*max(dataRangePrimate_rdc[,c(4)]))
  yupper <- c(CI[1,3], CI[1,3] + (CI[3,3])*max(dataRangePrimate_rdc[,c(4)]))
  lines(
    x=c(0,max(dataRangePrimate_rdc[,c(3)])),
    y=ymean
  )
  polygon(
    x=c(0, max(dataRangePrimate_rdc[,c(3)]), max(dataRangePrimate_rdc[,c(3)]), 0, 0),
    y=c(yupper, rev(ylower), ylower[1]),
    col=adjustcolor("black", alpha.f=0.1),
    border=NA#,
    #lty=2
  )
  
  #Overlay points
  points(dataRangePrimate_rdc[,c(3)], dataRangePrimate_rdc[,c(4)], pch=19, col="red",xpd=TRUE)
  
}

dev.off()

#Apparently you can access CI following the method in https://www.pnas.org/content/pnas/110/22/9001.full.pdf?with-ds=yes but to read

##------
## Sensitivity to phylogeny: using a block of trees
##------

##Remains to do: -> summarizing the output + check if it works

repetitionTrees=50
repetitionModels=30
  
sensitivityEstimate <- as.data.frame(matrix(NA, ncol=3, nrow=repetitionTrees*repetitionModels*length(traitToStudy)))
sensitivityPvalue <- as.data.frame(matrix(NA, ncol=3, nrow=repetitionTrees*repetitionModels*length(traitToStudy)))
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
  
  
  for(t in 1:repetitionTrees){
    
    #Dowload tree
    phylo_all <-read.nexus("Raw_data/Tree/TreeBlock_10kTrees_Primates_Version3.nex")
    phylo_init <- phylo_all[[t]]
    
    phylo <- force.ultrametric(tree=phylo_init, method="extend")#method="nnls")
    is.ultrametric(phylo)
    
    #check if this differ (was done outside the loop and was ok)
    # plot(phylo_init$edge.length, phylo$edge.length, xlab="Original tree", ylab="Post-chronos tree", pch=20, bty="n")
    # abline(a=0, b=1, col="gray", lwd=0.5)
    
    phylo <- multi2di(phylo)
    is.binary(phylo)
  ##----
  ## RUN MODEL FOR SENSITIVITY
  
  for(a in 1:length(traitToStudy)){
    
    rowToAdd=which(is.na(sensitivityLambda))[1]
    traitVectorSensitivity[rowToAdd] <- traitToStudy[a]
    #Matching the brain trait to the dataset with predictors
    
    dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[a])]
    dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,2])&!is.na(dataRangePrimate[,3])&!is.na(dataRangePrimate[,4]),]
    rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species
    
    #Readjust phylo tree
    phyloTree <- drop.tip(phylo,
                               phylo$tip.label[
                                 which(phylo$tip.label
                                       %nin%dataRangePrimate_rdc$Species)])
    
    ##--------
    ## Main model: based on consensus tree
    ##--------
    
    
    #Create data for phylogenetic regression
    # comp_data <- comparative.data(phy = phyloTree, data= dataRangePrimate_rdc,
    #                               names.col = Species, vcv = TRUE)
    # 
    # #Output distribution
    # hist(comp_data$data$Trait)
    # 
    # #Predictor distribution
    # hist(sqrt(comp_data$data$Number_species_cooccurrence))
    # hist(sqrt(comp_data$data$Overlap_average))
    # 
    # #Make predictor more symmetrical
    # comp_data$data$Number_species_cooccurrence.sqrt <- sqrt(comp_data$data$Number_species_cooccurrence)
    # # comp_data$data$Overlap_average.sqrt <- sqrt(comp_data$data$Overlap_average)
    # 
    # comp_data$data$Number_species_cooccurrence.sqrt.z <- scale(comp_data$data$Number_species_cooccurrence.sqrt)
    # comp_data$data$Overlap_average.sqrt.z <- scale(comp_data$data$Overlap_average.sqrt)
    
    dataRangePrimate_rdc$Overlap_average.sqrt <- sqrt(dataRangePrimate_rdc$Overlap_average)
    dataRangePrimate_rdc$Number_species_cooccurrence.sqrt <- sqrt(dataRangePrimate_rdc$Number_species_cooccurrence)
      
    modelBrain <- phylolm(formula = Trait ~ Overlap_average.sqrt + Number_species_cooccurrence.sqrt, 
                          data=dataRangePrimate_rdc ,phy=phyloConsensus,model="lambda",measurement_error=FALSE)
    
    # modelBrain <- gls(Trait ~ Overlap_average.sqrt  correlation = corPagel(1, phyloConsensus, form = ~Species),
    #                   data = dataRangePrimate_rdc, method = "ML")
    # CI <- intervals(modelBrain)
    # results <- anova(modelBrain)
    
    #Complete df results
    sensitivityEstimate[rowToAdd,] <- summary(modelBrain)$coefficients[,1]#round(summary(modelBrain)$coefficients[,1], digit=4)
    sensitivityPvalue[rowToAdd,] <- summary(modelBrain)$coefficients[,4]#c("-", round(results[1:2,5], digit=3))
    sensitivityLambda[rowToAdd] <- summary(modelBrain)$optpar#CI$corStruct[1,2])#round(summary(modelBrain)$param.CI$lambda$opt, digit=2)
    traitVectorSensitivity[rowToAdd] <- traitToStudy[a]
    
    }
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

dfBetasGradient <- cbind(minEst,
                         maxEst[,3])
colnames(dfBetasGradient) <- c("Trait", "Variable", "Est. min.", "Est.", "Est. max.")


dfBetasGradient[,3] <- sapply(dfBetasGradient[,3], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)

dfBetasGradient[,5] <- sapply(dfBetasGradient[,5], function(x) {
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

sensitivityGradient <- cbind(minEst,
                             maxEst[,3])
colnames(sensitivityGradient) <- c("Trait", "Variable", "Est. min.", "Est.", "Est. max.")

sensitivityGradient[,3] <- sapply(sensitivityGradient[,3], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)
  
sensitivityGradient[,5] <- sapply(sensitivityGradient[,5], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)

summarySensitivityGradient <- cbind(dfBetasGradient, sensitivityGradient[,3:5])
summarySensitivityGradient$Trait[-seq(from=1, to=nrow(summarySensitivityGradient), by=4)] <- ""

knitr::kable(summarySensitivityGradient, escape=TRUE, booktabs = TRUE,
             caption = "Sensitivity analysis of phylogenetic regressions to assess the selection gradient direction | Depicted is the minimum and maximum of estimates when one observation was removed at a time (DfBetas) or when varying the used phylogenetic tree and the data sampling (Phylogeny/Data)") %>% 
  kableExtra::kable_styling(latex_options = "striped") %>%
  kableExtra::kable_styling(latex_options="scale_down") %>%
  kableExtra::add_header_above(c("Model:" = 2, "DfBetas" = 3, "Phylogeny/Data" = 3))

###---
## Create main table

whichToBold <- which(as.numeric(results.df_gradient$"p-value")<=0.05)
toPlotBold <- rep(FALSE, times=nrow(results.df_gradient))
toPlotBold[whichToBold] <- TRUE

replace <- sapply(as.numeric(results.df_gradient$"p-value"), function(x) if(!is.na(x)){pvalueRound(x, text=FALSE)}else{x})
results.df_gradient$"p-value"[!is.na(replace)] <- replace
  
knitr::kable(results.df_gradient, escape=TRUE, booktabs = TRUE,
             caption = "Model estimates and significance of phylogenetic regressions to assess the selection gradient direction | Est.=Estimate, CI2.5%=Lower border of the CI95%, CI97.5%=Upper border of the CI95%, Sd= Standard deviation, t= Statitsitics t-vale. The brain area (as well as the associated sample size)
             are indicate prior each list of estimates. the transformation (logarithm or square-root) if indicated in parenthese by the abbreviation (log or sqrt).") %>%
  kableExtra::column_spec(2:ncol(results.df_gradient), bold = toPlotBold) %>%
  kableExtra::kable_styling(latex_options = "striped") %>%
  kableExtra::kable_styling(latex_options="scale_down")


save.image("REnvironments/PGLSdirectionSelection.RData")


# #Sensitivity to phylo/data
# minEst <- aggregate(sensitivityEstimate, by=list(traitVectorSensitivity), min)
# maxEst <- aggregate(sensitivityEstimate, by=list(traitVectorSensitivity), max)
# medianEst <- aggregate(sensitivityEstimate, by=list(traitVectorSensitivity), median)
# pid <- aggregate(sensitivityPvalue, by=list(traitVectorSensitivity), FUN= function(x)
#   length(x>0.05)
#   )
# lambdamin <- aggregate(sensitivityLambda, by=list(traitVectorSensitivity), min)
# lambdamax <- aggregate(sensitivityLambda, by=list(traitVectorSensitivity), max)
# lambdamedian <- aggregate(sensitivityLambda, by=list(traitVectorSensitivity), median)
# 

#Plot results of sensitivity
# 
# ##With number of co-occurring species
# xmin=-3
# xmax=3
# ymin=0
# ymax=length(unique(traitVectorSensitivity))+1
# plot(0, 0, xlab="Estimate (scaled variables)", ylab="",
#      las=1, type="n", tcl=-0.25, bty="n",
#      xaxt="n",xaxs="i",yaxs="i", yaxt="n", xpd=TRUE,
#      xlim=c(xmin, xmax),
#      ylim=c(ymin, ymax))
# 
# #Add grid
# addGrid(
#   xmin=xmin, xmax=xmax, xintsmall=(xmax-xmin)/20, xintbig=(xmax-xmin)/5,
#   ymin=ymin, ymax=ymax, yintsmall=(ymax-ymin)/20, yintbig=(ymax-ymin)/5,
#   axisPlot=FALSE)
# axis(side=1, at=seq(from=ymin, to=ymax, by=(ymax-ymin)/5), 
#      labels=round(seq(from=ymin, to=ymax, by=(ymax-ymin)/5), digit=2), tcl=-0.25)
# 
# #Add model names
# text(x=-1, y=1:(ymax-1), labels=unique(traitVectorSensitivity))
# 
# #Add results     
# 
# colour.v <- brewer.pal(n = 4, name = "Pastel1")
# for(i in 1:length(unique(traitVectorSensitivity))){
#   #Estimate model
#   errorBars(location=c(i-0.3, i-0.1, i+0.1), 
#             meanPt=medianEst[%NAME VARIABLE%==unique(traitVectorSensitivity)[i]],
#             barValue=rep(100, times=4),
#             upperBarValue=maxEst[%NAME VARIABLE%==unique(traitVectorSensitivity)[i]],
#             lowerBarValue=minEst[%NAME VARIABLE%==unique(traitVectorSensitivity)[i]],
#             col=colour.v[1:3],
#             horiz=TRUE,
#             symmetrical=FALSE)
#   
#   #Lambda model
#   errorBars(location=c(i+0.2), 
#             meanPt=lambdamedian[%NAME VARIABLE%==unique(traitVectorSensitivity)[i]],
#             barValue=rep(100, times=4),
#             upperBarValue=lambdamax[%NAME VARIABLE%==unique(traitVectorSensitivity)[i]],
#             lowerBarValue=lambdamin[%NAME VARIABLE%==unique(traitVectorSensitivity)[i]],
#             col=colour.v[4],
#             horiz=TRUE,
#             symmetrical=FALSE)
# }
# 
# legend("bottomright", legend=c("Intercept", "Overlap", "Nspecies", "Lambda"), col=colour.v[1:4], lty=c(1,1,1,1), bty="n")
