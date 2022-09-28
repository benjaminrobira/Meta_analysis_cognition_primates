####---------------------------------
#### PGLSdirectionSelection
####---------------------------------

# This script allows to perform PGLS analysis to assess the correlation of bodymass with spatial co-occurrence data, specifically in areas for which competitive models were shown to fit the best.

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

summaryData$Bodymass.log <- log(summaryData$Bodymass)
traitToStudy=c("Bodymass.log")
traitName=c("Bodymass (log)")
vifVector <- rep(NA, times=length(traitName))  

results.df <- as.data.frame(matrix(NA, nrow=5*length(traitToStudy), ncol=7))
colnames(results.df) <- c("", "Est.", "CI2.5%", "CI97.5%", "Sd", "t", "p-value")

##----
## Run analysis for each trait

for(a in 1:length(traitToStudy)){
  #Matching the brain trait to the dataset with predictors
  dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[a])]
  dataRangePrimate$Bodymass <- summaryData$Bodymass[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny)]
  # dataRangePrimate$Bodymass.log <- log(dataRangePrimate$Bodymass)
  
  dataRangePrimate_rdc <- dataRangePrimate[which(!is.na(dataRangePrimate[,2])&!is.na(dataRangePrimate[,3])&!is.na(dataRangePrimate[,4])),]
  #keep only frugivorous
  dataRangePrimate_rdc$Diet <- summaryData$DietaryGuild[match(dataRangePrimate_rdc$Species,summaryData$SpeciesForPhylogeny)]
  dataRangePrimate_rdc <- dataRangePrimate_rdc[dataRangePrimate_rdc$Diet=="Fruit",]
  
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
  # modelBodymassPgls <- pgls(formula = Trait ~ Overlap_average.sqrt + Number_species_cooccurrence.sqrt, 
  #                        data=comp_data, lambda="ML",param.CI=0.05)
  
  dataRangePrimate_rdc$Number_species_cooccurrence.sqrt <- sqrt(dataRangePrimate_rdc$Number_species_cooccurrence)
  dataRangePrimate_rdc$Overlap_average.sqrt <- sqrt(dataRangePrimate_rdc$Overlap_average)
  rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species
  
  # plot(dataRangePrimate_rdc$Number_species_cooccurrence.sqrt, dataRangePrimate_rdc$Overlap_average.sqrt)
  
  modelBodymass <- phylolm(formula = Trait ~ Overlap_average + Number_species_cooccurrence.sqrt, 
                        data=dataRangePrimate_rdc ,phy=phyloConsensus,model="lambda", measurement_error=FALSE,boot=repetitionBootstrap)
  
  # modelBodymass <- gls(Trait ~ Overlap_average.sqrt  correlation = corPagel(1, phyloConsensus, form = ~Species),
  #                   data = dataRangePrimate_rdc, method = "ML")
  # CI <- intervals(modelBodymass)
  # results <- anova(modelBodymass)
  
  #Complete df results
  
  #commented is when I used pgls + phylolm
  results.df[(firstRowWrite+1):(firstRowWrite+3), 2] <- roundIntelligent(summary(modelBodymass)$coefficients[,1], digit=2) #roundIntelligent(summary(modelBodymass)$tTable[,1], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+3), 3] <- roundIntelligent(summary(modelBodymass)$coefficients[,4], digit=2) #roundIntelligent(CI$coef[,1], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+3), 4] <- roundIntelligent(summary(modelBodymass)$coefficients[,5], digit=2) #roundIntelligent(CI$coef[,3], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+3), 5] <- roundIntelligent(summary(modelBodymass)$coefficients[,2], digit=2) #roundIntelligent(summary(modelBodymass)$tTable[,2], digit=2)
  results.df[(firstRowWrite+1):(firstRowWrite+3), 6] <- c("-", roundIntelligent(summary(modelBodymass)$coefficients[2:3,3], digit=2))
  results.df[(firstRowWrite+1):(firstRowWrite+3), 7] <- c("-", roundIntelligent(summary(modelBodymass)$coefficients[2:3,6], digit=2))
  
  # val1 <- ifelse(is.na(summary(modelBodymassPgls)$param.CI$lambda$ci.val[1]), 0, round(summary(modelBodymassPgls)$param.CI$lambda$ci.val[1], digit=2))
  # val2 <- ifelse(is.na(summary(modelBodymassPgls)$param.CI$lambda$ci.val[2]), 1, round(summary(modelBodymassPgls)$param.CI$lambda$ci.val[2], digit=2))
  # 
  # lambda <- paste(round(summary(modelBodymassPgls)$param.CI$lambda$opt, digit=2), " [",val1,",",val2,"]", sep="")
  #         
  results.df[firstRowWrite+4,2] <- roundIntelligent(summary(modelBodymass)$optpar)#bootmean[5])
  results.df[firstRowWrite+4,3] <- roundIntelligent(summary(modelBodymass)$bootconfint95[1,5])#CI$corStruct[1,1])
  results.df[firstRowWrite+4,4] <- roundIntelligent(summary(modelBodymass)$bootconfint95[2,5])#CI$corStruct[1,3])
  
  #Check assumptions
  # pdf(file=paste("Figure/Diagnostics", traitName[a], ".pdf", sep=""))
  assign(paste("modelBodymass", traitName[a], sep="_"), modelBodymass)
  #assign(paste(modelBodymassPgls, traitName[a], sep="_"), modelBodymassPgls)
  
  #since only two variables
  cortest <- cor.test(dataRangePrimate_rdc$Number_species_cooccurrence.sqrt, dataRangePrimate_rdc$Overlap_average)
  
  vifVector[a] <- 1/(1-cortest$estimate**2)#max(vif(modelBodymass))
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
  
  #keep only frugivorous
  dataRangePrimate_rdc$Diet <- summaryData$DietaryGuild[match(dataRangePrimate_rdc$Species,summaryData$SpeciesForPhylogeny)]
  dataRangePrimate_rdc <- dataRangePrimate_rdc[dataRangePrimate_rdc$Diet=="Fruit",]
  
  #Readjust phylo tree
  phyloConsensus <- drop.tip(phylo,
                             phylo$tip.label[
                               which(phylo$tip.label
                                     %nin%dataRangePrimate_rdc$Species)])
  
  for(r in 1:nrow(dataRangePrimate_rdc)){
    dataRangePrimate_rdc$Number_species_cooccurrence.sqrt <- sqrt(dataRangePrimate_rdc$Number_species_cooccurrence)
    dataRangePrimate_rdc$Overlap_average.sqrt <- sqrt(dataRangePrimate_rdc$Overlap_average)
    rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species
    
    dataRangePrimate_rdc2 <- dataRangePrimate_rdc[-r,] 
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
    
    modelBodymass <- phylolm(Trait ~ Overlap_average + Number_species_cooccurrence.sqrt, phy=phyloConsensus, data = dataRangePrimate_rdc2, model = "lambda", measurement_error=FALSE)
    
    dfBetasEstimate[rowToAdd,] <- summary(modelBodymass)$coefficients[,1]
    dfBetasPvalue[rowToAdd,] <- summary(modelBodymass)$coefficients[,4]
    dfBetasLambda[rowToAdd] <- summary(modelBodymass)$optpar
    traitVectordfBetas[rowToAdd] <- traitName[a]
    
    progress(r/nrow(dataRangePrimate_rdc)*100)
  }
}

##------
## Sensitivity to phylogeny: using a block of trees
##------

##Remains to do: -> summarizing the output + check if it works

repetitionTrees=50
repetitionModels=50

sensitivityEstimate <- as.data.frame(matrix(NA, ncol=3, nrow=repetitionTrees*repetitionModels*length(traitToStudy)))
sensitivityPvalue <- as.data.frame(matrix(NA, ncol=3, nrow=repetitionTrees*repetitionModels*length(traitToStudy)))
sensitivityLambda <- rep(NA, times=repetitionTrees*repetitionModels*length(traitToStudy))
traitVectorSensitivity <- rep(NA, times=repetitionTrees*repetitionModels*length(traitToStudy))

for(d in 1:repetitionModels){#data
  
  summaryData$Family<- NA
  summaryData$DietaryGuild <- NA
  summaryData$FrugivoryPercent <- NA
  summaryData$FolivoryPercent <- NA
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
      traitVectorSensitivity[rowToAdd] <- traitName[a]
      #Matching the brain trait to the dataset with predictors
      
      dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[a])]
      dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,2])&!is.na(dataRangePrimate[,3])&!is.na(dataRangePrimate[,4]),]
      rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species
      
      #keep only frugivorous
      dataRangePrimate_rdc$Diet <- summaryData$DietaryGuild[match(dataRangePrimate_rdc$Species,summaryData$SpeciesForPhylogeny)]
      dataRangePrimate_rdc <- dataRangePrimate_rdc[dataRangePrimate_rdc$Diet=="Fruit",]
      
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
      
      modelBodymass <- phylolm(formula = Trait ~ Overlap_average + Number_species_cooccurrence.sqrt, 
                            data=dataRangePrimate_rdc, phy=phyloConsensus,model="lambda",measurement_error=FALSE)
      
      # modelBodymass <- gls(Trait ~ Overlap_average.sqrt  correlation = corPagel(1, phyloConsensus, form = ~Species),
      #                   data = dataRangePrimate_rdc, method = "ML")
      # CI <- intervals(modelBodymass)
      # results <- anova(modelBodymass)
      
      #Complete df results
      sensitivityEstimate[rowToAdd,] <- summary(modelBodymass)$coefficients[,1]#round(summary(modelBodymass)$coefficients[,1], digit=4)
      sensitivityPvalue[rowToAdd,] <- summary(modelBodymass)$coefficients[,4]#c("-", round(results[1:2,5], digit=3))
      sensitivityLambda[rowToAdd] <- summary(modelBodymass)$optpar#CI$corStruct[1,2])#round(summary(modelBodymass)$param.CI$lambda$opt, digit=2)
      traitVectorSensitivity[rowToAdd] <- traitName[a]
      
    }
  }
  
  progress(d/repetitionModels*100)
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
whatIs <- rep(traitName, each=4)
numberForFinishingMatch <- rep(1:4, times=length(traitName))

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
whatIs <- rep(traitName, each=4)
numberForFinishingMatch <- rep(1:4, times=length(traitName))

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

results.df_gradient_init <- results.df_gradient
results.df_gradient <- results.df_gradient_init

whichToBold <- which(as.numeric(results.df_gradient$"p-value")<=0.05)
toPlotBold <- rep(FALSE, times=nrow(results.df_gradient))
toPlotBold[whichToBold] <- TRUE

replace <- sapply(as.numeric(results.df_gradient$"p-value"), function(x) if(!is.na(x)&x!=""&x!="-"){pvalueRound(x, text=FALSE)}else{x})
results.df_gradient$"p-value"[!is.na(replace)] <- replace[!is.na(replace)]

knitr::kable(results.df_gradient, escape=TRUE, booktabs = TRUE,
             caption = "Model estimates and significance of phylogenetic regressions to assess the selection gradient direction | Est.=Estimate, CI2.5%=Lower border of the CI95%, CI97.5%=Upper border of the CI95%, Sd= Standard deviation, t= Statistics t-value. The brain areas (as well as the associated sample size)
             are indicated prior to each list of estimates. the transformation (logarithm or square-root) if indicated in parentheses by the abbreviation (log or sqrt).") %>%
  kableExtra::column_spec(2:ncol(results.df_gradient), bold = toPlotBold) %>%
  kableExtra::kable_styling(latex_options = "striped") %>%
  kableExtra::kable_styling(latex_options="scale_down")

save.image("REnvironments/PGLSbodymass.RData")

results.df_gradient_Bodymass <- results.df_gradient
summarySensitivityGradient_Bodymass <- summarySensitivityGradient

rm(list=setdiff(ls(), c("modelBodymass", "results.df_gradient_Bodymass", "summarySensitivityGradient_Bodymass")))
save.image("REnvironments/PGLSbodymass_rdc.RData")

