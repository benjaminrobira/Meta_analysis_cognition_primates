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

# traitToStudy=c("EQ.log", "ratioBrain.log", "ratioHippocampus.log", "ratioNeocortex.log", "ratioCerebellum.log", "ratioStriatum.log", "ratioMOB.log")
# traitName=c("EQ (log)", "Brain (/bodymass, log)", "Hippocampus (/bodymass, log)", "Neocortex (/bodymass, log)", "Cerebellum (/bodymass, log)", "Striatum (/bodymass, log)", "MOB (/bodymass, log)")
vifVector <- NA#rep(NA, times=length(traitName))  

results.df <- as.data.frame(matrix(NA, nrow=5, ncol=7))
colnames(results.df) <- c("", "Est.", "CI2.5%", "CI97.5%", "Sd", "t", "p-value")

##----
## Run analysis for each trait

#Matching the brain trait to the dataset with predictors
dataRangePrimate$Diversification_rate <-  summaryData$Diversification_rate[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny)]
dataRangePrimate$Bodymass <- summaryData$Bodymass[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny)]
dataRangePrimate$Bodymass.log <- log(dataRangePrimate$Bodymass)

dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,2])&!is.na(dataRangePrimate[,3])&!is.na(dataRangePrimate[,4]),]
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
results.df[firstRowWrite, 1] <- paste("Diversification", " (N=", nrow(dataRangePrimate_rdc), ")", sep="") 
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
# modelBrainPgls <- pgls(formula = Diversification_rate ~Overlap_average.sqrt + Number_species_cooccurrence.sqrt, 
#                        data=comp_data, lambda="ML",param.CI=0.05)

dataRangePrimate_rdc$Number_species_cooccurrence.sqrt <- sqrt(dataRangePrimate_rdc$Number_species_cooccurrence)
dataRangePrimate_rdc$Overlap_average.sqrt <- sqrt(dataRangePrimate_rdc$Overlap_average)
rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species

# plot(dataRangePrimate_rdc$Number_species_cooccurrence.sqrt, dataRangePrimate_rdc$Overlap_average.sqrt)

modelBrain <- phylolm(formula = Diversification_rate ~Overlap_average + Number_species_cooccurrence.sqrt, 
                      data=dataRangePrimate_rdc ,phy=phyloConsensus,model="lambda", measurement_error=FALSE,boot=repetitionBootstrap)

# modelBrain <- gls(Diversification_rate ~Overlap_average.sqrt  correlation = corPagel(1, phyloConsensus, form = ~Species),
#                   data = dataRangePrimate_rdc, method = "ML")
# CI <- intervals(modelBrain)
# results <- anova(modelBrain)

#Complete df results

#commented is when I used pgls + phylolm
results.df[(firstRowWrite+1):(firstRowWrite+3), 2] <- roundIntelligent(summary(modelBrain)$coefficients[,1], digit=2) #roundIntelligent(summary(modelBrain)$tTable[,1], digit=2)
results.df[(firstRowWrite+1):(firstRowWrite+3), 3] <- roundIntelligent(summary(modelBrain)$coefficients[,4], digit=2) #roundIntelligent(CI$coef[,1], digit=2)
results.df[(firstRowWrite+1):(firstRowWrite+3), 4] <- roundIntelligent(summary(modelBrain)$coefficients[,5], digit=2) #roundIntelligent(CI$coef[,3], digit=2)
results.df[(firstRowWrite+1):(firstRowWrite+3), 5] <- roundIntelligent(summary(modelBrain)$coefficients[,2], digit=2) #roundIntelligent(summary(modelBrain)$tTable[,2], digit=2)
results.df[(firstRowWrite+1):(firstRowWrite+3), 6] <- c("-", roundIntelligent(summary(modelBrain)$coefficients[2:3,3], digit=2))
results.df[(firstRowWrite+1):(firstRowWrite+3), 7] <- c("-", roundIntelligent(summary(modelBrain)$coefficients[2:3,6], digit=2))

# val1 <- ifelse(is.na(summary(modelBrainPgls)$param.CI$lambda$ci.val[1]), 0, round(summary(modelBrainPgls)$param.CI$lambda$ci.val[1], digit=2))
# val2 <- ifelse(is.na(summary(modelBrainPgls)$param.CI$lambda$ci.val[2]), 1, round(summary(modelBrainPgls)$param.CI$lambda$ci.val[2], digit=2))
# 
# lambda <- paste(round(summary(modelBrainPgls)$param.CI$lambda$opt, digit=2), " [",val1,",",val2,"]", sep="")
#         
results.df[firstRowWrite+4,2] <- roundIntelligent(summary(modelBrain)$optpar)#bootmean[5])
results.df[firstRowWrite+4,3] <- roundIntelligent(summary(modelBrain)$bootconfint95[1,5])#CI$corStruct[1,1])
results.df[firstRowWrite+4,4] <- roundIntelligent(summary(modelBrain)$bootconfint95[2,5])#CI$corStruct[1,3])

#Check assumptions
# pdf(file=paste("Figure/Diagnostics", traitName[a], ".pdf", sep=""))
assign(paste("modelBrainDiversificationAndSympatry", sep="_"), modelBrain)
#assign(paste(modelBrainPgls, traitName[a], sep="_"), modelBrainPgls)

#since only two variables
cortest <- cor.test(dataRangePrimate_rdc$Number_species_cooccurrence.sqrt, dataRangePrimate_rdc$Overlap_average)

vifVector <- 1/(1-cortest$estimate**2)#max(vif(modelBrain))
#diagnostics.plot()
# dev.off()

results.df[2:ncol(results.df)] <- replaceNAtable(results.df[2:ncol(results.df)], "")
results.df_diversificationAndSympatry <- results.df


##############
## PLOTTING RESULTS
##############

#Remaining to do in plotting:

#1) create layout + matrix for all plots

#finish the plot part:write for overlap + adjust axis and color + add regression
pdf("Plots/diversificationAndSympatryPGLS.pdf", height=10, width=20)
source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
layout(mat=t(c(1,2)),
       widths=c(5,5), heights=rep(5))
par(mar=c(3.5, 3.5, 1, 1), mgp=c(2.5, 0.5, 0), xpd=TRUE, cex=1.2)

    countLabel=0
    countLabel=countLabel+1
    #Matching the brain trait to the dataset with predictors
    dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,2])&!is.na(dataRangePrimate[,3])&!is.na(dataRangePrimate[,4]),]
    rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species
    
    #keep only frugivorous
    dataRangePrimate_rdc$Diet <- summaryData$DietaryGuild[match(dataRangePrimate_rdc$Species,summaryData$SpeciesForPhylogeny)]
    dataRangePrimate_rdc <- dataRangePrimate_rdc[dataRangePrimate_rdc$Diet=="Fruit",]
    
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
    # modelBrain <- pgls(formula = Diversification_rate ~Overlap_average.sqrt + Number_species_cooccurrence.sqrt , data = comp_data, lambda = "ML")
    # 
    dataRangePrimate_rdc$Number_species_cooccurrence.sqrt <- sqrt(dataRangePrimate_rdc$Number_species_cooccurrence)
    dataRangePrimate_rdc$Overlap_average.sqrt <- sqrt(dataRangePrimate_rdc$Overlap_average)
    
    dataRangePrimate_rdc$Number_species_cooccurrence.sqrt.z <- scale(dataRangePrimate_rdc$Number_species_cooccurrence.sqrt)
    dataRangePrimate_rdc$Overlap_average.z <- scale(dataRangePrimate_rdc$Overlap_average)
    
    modelBrain <- phylolm(formula = Diversification_rate ~ Overlap_average.z + Number_species_cooccurrence.sqrt, 
                          data=dataRangePrimate_rdc, phy=phyloConsensus, model="lambda", measurement_error=FALSE, boot=repetitionBootstrap)
    
    CI <- cbind(modelBrain$bootmean, t(summary(modelBrain)$bootconfint95))
    
    ##----
    #Plot against N co-occ
    
    
    xmin=0#min(round(dataRangePrimate_rdc[,c(2)], digit=2))
    xmax=max(round(dataRangePrimate_rdc$Number_species_cooccurrence.sqrt, digit=2))
    ymin=min(round(dataRangePrimate_rdc[,c(4)], digit=2))
    ymax=max(round(dataRangePrimate_rdc[,c(4)], digit=2))
    
    par(mar=c(3.5, 3.5, 1, 1), mgp=c(2.5, 0.5, 0), xpd=TRUE, cex=1.2)
    ##With number of co-occurring species
    plot(dataRangePrimate_rdc$Number_species_cooccurrence.sqrt, dataRangePrimate_rdc[,c(4)], xlab="Number of sympatric frugivorous species (sqrt)", ylab="Diversification",
         font.lab=2, cex.lab=1.25,xlim=c(xmin, xmax), 
         las=1, type="n", tcl=-0.25, bty="n",
         xaxt="n",xaxs="i",yaxs="i", yaxt="n", xpd=TRUE)
    
    #Add grid
    addGrid(
      cexAxisX=1.15, cexAxisY=1.15,
      xmin=xmin, xmax=xmax, xintsmall=(xmax-xmin)/20, xintbig=(xmax-xmin)/5,
      ymin=ymin, ymax=ymax, yintsmall=(ymax-ymin)/20, yintbig=(ymax-ymin)/5,
      axisPlot=TRUE, round=TRUE, digit=c(2,2))
    axis(side=1, at=round(seq(from=xmin, to=xmax, by=(xmax-xmin)/5), digit=1), labels=round(seq(from=xmin, to=xmax, by=(xmax-xmin)/5), digit=1), las=1, tcl=-0.25)
    addLabel(xfrac=0.05, yfrac=0.05, label="A", circle=TRUE, radiuscircle=(xmax-xmin)/35, circle.bg="black", font.col="white")
    
    #Add result model
    ymean <- c(CI[1,1], CI[1,1] + CI[3,1]*max(dataRangePrimate_rdc[,c(4)]))
    ylower <- c(CI[1,2], CI[1,2] + (CI[3,2])*max(dataRangePrimate_rdc[,c(4)]))
    yupper <- c(CI[1,3], CI[1,3] + (CI[3,3])*max(dataRangePrimate_rdc[,c(4)]))
    lines(
      x=c(0,max(dataRangePrimate_rdc$Number_species_cooccurrence.sqrt)),
      y=ymean
    )
    polygon(
      x=c(0, max(dataRangePrimate_rdc$Number_species_cooccurrence.sqrt), max(dataRangePrimate_rdc$Number_species_cooccurrence.sqrt), 0, 0),
      y=c(yupper, rev(ylower), ylower[1]),
      col=adjustcolor("black", alpha.f=0.1),
      border=NA#,
      #lty=2
    )
    
    #Add background tree
    col=list(col.edge=setNames(rep("darkgrey",nrow(phyloConsensus$edge)),as.character(phyloConsensus$edge[,2])),
             col.node=setNames(rep("black",max(phyloConsensus$edge)),as.character(1:max(phyloConsensus$edge))))
    
    toPlot <- cbind(dataRangePrimate_rdc$Number_species_cooccurrence.sqrt, dataRangePrimate_rdc[,c(4)])
    rownames(toPlot) <- rownames(dataRangePrimate_rdc)
                              
    phylomorphospace(phyloConsensus,toPlot, add=TRUE, label="true", lty=3,
                     control=col, xpd=TRUE)
    
    ##----
    #Plot against overlap
    
    modelBrain <- phylolm(formula = Diversification_rate ~ Overlap_average + Number_species_cooccurrence.sqrt.z, 
                          data=dataRangePrimate_rdc ,phy=phyloConsensus,model="lambda",measurement_error=FALSE,boot=repetitionBootstrap)
    
    CI <- cbind(modelBrain$bootmean, t(summary(modelBrain)$bootconfint95))
    
    xmin=0#min(round(dataRangePrimate_rdc[,c(3)], digit=2))
    xmax=1#max(round(dataRangePrimate_rdc[,c(3)], digit=2))
    ymin=min(round(dataRangePrimate_rdc[,c(4)], digit=2))
    ymax=max(round(dataRangePrimate_rdc[,c(4)], digit=2))
    
    par(mar=c(3.5, 1, 1, 3), mgp=c(2.5, 0.5, 0), xpd=TRUE)
    
    ##With overlap
    plot(dataRangePrimate_rdc[,c(3)], dataRangePrimate_rdc[,c(4)], xlab="Average surfacic overlap with sympatric frugivorous species", ylab="",
         font.lab=2, cex.lab=1.25,
         las=1, type="n", tcl=-0.25, bty="n",
         xaxt="n",xaxs="i",yaxs="i", yaxt="n", xpd=TRUE)
    
    #Add grid
    addGrid(
      cexAxisX=1.15, cexAxisY=1.15,
      xmin=xmin, xmax=xmax, xintsmall=(xmax-xmin)/20, xintbig=(xmax-xmin)/5,
      ymin=ymin, ymax=ymax, yintsmall=(ymax-ymin)/20, yintbig=(ymax-ymin)/5,
      axisPlot=FALSE, round=TRUE, digit=c(2,2))
    axis(side=1, at=round(seq(from=xmin, to=xmax, by=(xmax-xmin)/5), digit=1), labels=round(seq(from=xmin, to=xmax, by=(xmax-xmin)/5), digit=1), las=1, tcl=-0.25)
    
    addLabel(xfrac=0.05, yfrac=0.05, label="B", circle=TRUE, radiuscircle=(xmax-xmin)/35, circle.bg="black", font.col="white")
    
    #Add background tree
    col=list(col.edge=setNames(rep("darkgrey",nrow(phyloConsensus$edge)),as.character(phyloConsensus$edge[,2])),
             col.node=setNames(rep("black",max(phyloConsensus$edge)),as.character(1:max(phyloConsensus$edge))))
    
    phylomorphospace(phyloConsensus,dataRangePrimate_rdc[,c(3,4)], add=TRUE, label="true", lty=3,
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

dev.off()

######
## DFbetas

#Calculate Dfbetas

dfBetasEstimate <- as.data.frame(matrix(NA, ncol=3, nrow=nrow(dataRangePrimate)))
dfBetasPvalue <- as.data.frame(matrix(NA, ncol=3, nrow=nrow(dataRangePrimate)))
dfBetasLambda <- rep(NA, times=nrow(dataRangePrimate))

library(svMisc)
#Matching the brain trait to the dataset with predictors
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
  
  modelBrain <- phylolm(Diversification_rate ~ Overlap_average + Number_species_cooccurrence.sqrt, phy=phyloConsensus, data = dataRangePrimate_rdc2, model = "lambda", measurement_error=FALSE)
  
  dfBetasEstimate[rowToAdd,] <- summary(modelBrain)$coefficients[,1]
  dfBetasPvalue[rowToAdd,] <- summary(modelBrain)$coefficients[,4]
  dfBetasLambda[rowToAdd] <- summary(modelBrain)$optpar
  
  progress(r/nrow(dataRangePrimate_rdc)*100)
}

#Apparently you can access CI following the method in https://www.pnas.org/content/pnas/110/22/9001.full.pdf?with-ds=yes but to read

##------
## Sensitivity to phylogeny: using a block of trees
##------

##Remains to do: -> summarizing the output + check if it works

repetitionTrees=50
repetitionModels=50

sensitivityEstimate <- as.data.frame(matrix(NA, ncol=3, nrow=repetitionTrees*repetitionModels))
sensitivityPvalue <- as.data.frame(matrix(NA, ncol=3, nrow=repetitionTrees*repetitionModels))
sensitivityLambda <- rep(NA, times=repetitionTrees*repetitionModels)

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
      
      rowToAdd=which(is.na(sensitivityLambda))[1]
      #Matching the brain trait to the dataset with predictors
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
      
      modelBrain <- phylolm(formula = Diversification_rate ~Overlap_average + Number_species_cooccurrence.sqrt, 
                            data=dataRangePrimate_rdc, phy=phyloConsensus,model="lambda",measurement_error=FALSE)
      
      # modelBrain <- gls(Diversification_rate ~Overlap_average.sqrt  correlation = corPagel(1, phyloConsensus, form = ~Species),
      #                   data = dataRangePrimate_rdc, method = "ML")
      # CI <- intervals(modelBrain)
      # results <- anova(modelBrain)
      
      #Complete df results
      sensitivityEstimate[rowToAdd,] <- summary(modelBrain)$coefficients[,1]#round(summary(modelBrain)$coefficients[,1], digit=4)
      sensitivityPvalue[rowToAdd,] <- summary(modelBrain)$coefficients[,4]#c("-", round(results[1:2,5], digit=3))
      sensitivityLambda[rowToAdd] <- summary(modelBrain)$optpar#CI$corStruct[1,2])#round(summary(modelBrain)$param.CI$lambda$opt, digit=2)
      
    }
  
  progress(d/repetitionModels*100)
}

###---------
### Sensitivity to fraction of taxon

fraction.v <- c(60, 70, 80, 90, 95)

sensitivityFractionTaxonEstimate <- as.data.frame(matrix(NA, ncol=3, nrow=length(fraction.v)))
sensitivityFractionTaxonPvalue <- as.data.frame(matrix(NA, ncol=3, nrow=length(fraction.v)))
sensitivityFractionTaxonLambda <- rep(NA, times=length(fraction.v))

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
  ## Run model
  dataRangePrimate$Diversification_rate <-  summaryData$Diversification_rate[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny)]
  dataRangePrimate$Bodymass <- summaryData$Bodymass[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny)]
  dataRangePrimate$Bodymass.log <- log(dataRangePrimate$Bodymass)
  
  dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,2])&!is.na(dataRangePrimate[,3])&!is.na(dataRangePrimate[,4]),]
  #keep only frugivorous
  dataRangePrimate_rdc$Diet <- summaryData$DietaryGuild[match(dataRangePrimate_rdc$Species,summaryData$SpeciesForPhylogeny)]
  dataRangePrimate_rdc <- dataRangePrimate_rdc[dataRangePrimate_rdc$Diet=="Fruit",]
  
  
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
  
  rowToAdd=which(is.na(sensitivityFractionTaxonLambda))[1]
  
  rownames(dataRangePrimate_rdc2) <- dataRangePrimate_rdc2$Species
  
  modelBrain <- phylolm(Diversification_rate ~ Overlap_average + Number_species_cooccurrence.sqrt, phy=phyloConsensus, data = dataRangePrimate_rdc2, model = "lambda", measurement_error=FALSE)
  
  sensitivityFractionTaxonEstimate[rowToAdd,] <- summary(modelBrain)$coefficients[,1]
  sensitivityFractionTaxonPvalue[rowToAdd,] <- summary(modelBrain)$coefficients[,4]
  sensitivityFractionTaxonLambda[rowToAdd] <- summary(modelBrain)$optpar
  
}

###----
# Create sensitivity tables

dfBetasEstimate <- dfBetasEstimate[!is.na(dfBetasEstimate[,1]),]
dfBetasPvalue <- dfBetasPvalue[!is.na(dfBetasPvalue)]
dfBetasLambda <- dfBetasLambda[!is.na(dfBetasLambda)]

sensitivityEstimate <- sensitivityEstimate[!is.na(sensitivityEstimate[,1]),]
sensitivityPvalue <- sensitivityPvalue[!is.na(sensitivityPvalue)]
sensitivityLambda <- sensitivityLambda[!is.na(sensitivityLambda)]

sensitivityFractionTaxonEstimate <- sensitivityFractionTaxonEstimate[!is.na(sensitivityFractionTaxonEstimate[,1]),]
sensitivityFractionTaxonPvalue <- sensitivityFractionTaxonPvalue[!is.na(sensitivityFractionTaxonPvalue)]
sensitivityFractionTaxonLambda <- sensitivityFractionTaxonLambda[!is.na(sensitivityFractionTaxonLambda)]

#Summary sensitivity

#Dfbetas
library(tidyr)
minEst <- tibble(cbind(aggregate(dfBetasEstimate, by=list(rep(1, times=nrow(dfBetasEstimate))), FUN=min), aggregate(dfBetasLambda, by=list(rep(1, times=nrow(dfBetasEstimate))), min)[,2]))
minEst <- minEst[,-1]
colnames(minEst) <- c("Intercept", "Overlap", "N co-occurrence", "Lambda")
minEst <- pivot_longer(minEst, col=1:4, names_to="Variable", values_to="Estimate")

maxEst <- tibble(cbind(aggregate(dfBetasEstimate, by=list(rep(1, times=nrow(dfBetasEstimate))), max), aggregate(dfBetasLambda, by=list(rep(1, times=nrow(dfBetasEstimate))), max)[,2]))
maxEst <- maxEst[,-1]
colnames(maxEst) <- c("Intercept", "Overlap", "N co-occurrence", "Lambda")
maxEst <- pivot_longer(maxEst, col=1:4, names_to="Variable", values_to="Estimate")

estimate <- results.df[results.df[,2]!="",2]
#medianEst <- aggregate(dfBetasEstimate, by=list(traitVectordfBetas), median)
minEst$EstimateTrue <- estimate

dfBetasdiversificationAndSympatry <- cbind(minEst,
                         maxEst[,2])
colnames(dfBetasdiversificationAndSympatry) <- c("Variable", "Est. min.", "Est.", "Est. max.")


dfBetasdiversificationAndSympatry[,2] <- sapply(dfBetasdiversificationAndSympatry[,2], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)

dfBetasdiversificationAndSympatry[,4] <- sapply(dfBetasdiversificationAndSympatry[,4], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)


#sensitivityFractionTaxon
library(tidyr)
minEst <- tibble(cbind(aggregate(sensitivityFractionTaxonEstimate, by=list(rep(1, times=nrow(sensitivityFractionTaxonEstimate))), FUN=min), aggregate(sensitivityFractionTaxonLambda, by=list(rep(1, times=nrow(sensitivityFractionTaxonEstimate))), min)[,2]))
minEst <- minEst[,-1]
colnames(minEst) <- c("Intercept", "Overlap", "N co-occurrence", "Lambda")
minEst <- pivot_longer(minEst, col=1:4, names_to="Variable", values_to="Estimate")

maxEst <- tibble(cbind(aggregate(sensitivityFractionTaxonEstimate, by=list(rep(1, times=nrow(sensitivityFractionTaxonEstimate))), max), aggregate(sensitivityFractionTaxonLambda, by=list(rep(1, times=nrow(sensitivityFractionTaxonEstimate))), max)[,2]))
maxEst <- maxEst[,-1]
colnames(maxEst) <- c("Intercept", "Overlap", "N co-occurrence", "Lambda")
maxEst <- pivot_longer(maxEst, col=1:4, names_to="Variable", values_to="Estimate")

estimate <- results.df[results.df[,2]!="",2]
#medianEst <- aggregate(sensitivityFractionTaxonEstimate, by=list(traitVectorsensitivityFractionTaxon), median)
minEst$EstimateTrue <- estimate

sensitivityFractionTaxondiversificationAndSympatry <- cbind(minEst,
                         maxEst[,2])
colnames(sensitivityFractionTaxondiversificationAndSympatry) <- c("Variable", "Est. min.", "Est.", "Est. max.")


sensitivityFractionTaxondiversificationAndSympatry[,2] <- sapply(sensitivityFractionTaxondiversificationAndSympatry[,2], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)

sensitivityFractionTaxondiversificationAndSympatry[,4] <- sapply(sensitivityFractionTaxondiversificationAndSympatry[,4], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)


#sensitivity
library(tidyr)
minEst <- tibble(cbind(aggregate(sensitivityEstimate, by=list(rep(1, times=nrow(sensitivityEstimate))), FUN=min), aggregate(sensitivityLambda, by=list(rep(1, times=nrow(sensitivityEstimate))), min)[,2]))
minEst <- minEst[,-1]
colnames(minEst) <- c("Intercept", "Overlap", "N co-occurrence", "Lambda")
minEst <- pivot_longer(minEst, col=1:4, names_to="Variable", values_to="Estimate")

maxEst <- tibble(cbind(aggregate(sensitivityEstimate, by=list(rep(1, times=nrow(sensitivityEstimate))), max), aggregate(sensitivityLambda, by=list(rep(1, times=nrow(sensitivityEstimate))), max)[,2]))
maxEst <- maxEst[,-1]
colnames(maxEst) <- c("Intercept", "Overlap", "N co-occurrence", "Lambda")
maxEst <- pivot_longer(maxEst, col=1:4, names_to="Variable", values_to="Estimate")

estimate <- results.df[results.df[,2]!="",2]
#medianEst <- aggregate(sensitivityEstimate, by=list(traitVectorsensitivity), median)
minEst$EstimateTrue <- estimate

sensitivitydiversificationAndSympatry <- cbind(minEst,
                         maxEst[,2])
colnames(sensitivitydiversificationAndSympatry) <- c("Variable", "Est. min.", "Est.", "Est. max.")


sensitivitydiversificationAndSympatry[,2] <- sapply(sensitivitydiversificationAndSympatry[,2], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)

sensitivitydiversificationAndSympatry[,4] <- sapply(sensitivitydiversificationAndSympatry[,4], function(x) {
  x <- roundIntelligent(as.numcharac(x))
  return(x)
}
)

summarySensitivitydiversificationAndSympatry <- cbind(dfBetasdiversificationAndSympatry, sensitivitydiversificationAndSympatry[,2:4], sensitivityFractionTaxondiversificationAndSympatry[,2:4])

knitr::kable(summarySensitivitydiversificationAndSympatry, escape=TRUE, booktabs = TRUE,
             caption = "Sensitivity analysis of phylogenetic regressions to assess the selection diversificationAndSympatry direction | Depicted is the minimum and maximum of estimates when one observation was removed at a time (DfBetas) or when varying the used phylogenetic tree and the data sampling (Phylogeny/Data)") %>% 
  kableExtra::kable_styling(latex_options = "striped") %>%
  kableExtra::kable_styling(latex_options="scale_down") %>%
  kableExtra::add_header_above(c("Model:" = 2, "DfBetas" = 3, "Phylogeny/Data" = 3, "Sampling fraction"))

###---
## Create main table

results.df_diversificationAndSympatry_init <- results.df_diversificationAndSympatry
results.df_diversificationAndSympatry <- results.df_diversificationAndSympatry_init

whichToBold <- which(as.numeric(results.df_diversificationAndSympatry$"p-value")<=0.05)
toPlotBoldSympatry <- rep(FALSE, times=nrow(results.df_diversificationAndSympatry))
toPlotBoldSympatry[whichToBold] <- TRUE

replace <- sapply(as.numeric(results.df_diversificationAndSympatry$"p-value"), function(x) if(!is.na(x)&x!=""&x!="-"){pvalueRound(x, text=FALSE)}else{x})
results.df_diversificationAndSympatry$"p-value"[!is.na(replace)] <- replace[!is.na(replace)]

knitr::kable(results.df_diversificationAndSympatry, escape=TRUE, booktabs = TRUE,
             caption = "Model estimates and significance of phylogenetic regressions to assess the correlation between diversification and sympatry | Est.=Estimate, CI2.5%=Lower border of the CI95%, CI97.5%=Upper border of the CI95%, Sd= Standard deviation, t= Statitsitics t-vale. The brain area (as well as the associated sample size)
             are indicate prior each list of estimates. the transformation (logarithm or square-root) if indicated in parenthese by the abbreviation (log or sqrt).") %>%
  kableExtra::column_spec(2:ncol(results.df_diversificationAndSympatry), bold = toPlotBoldSympatry) %>%
  kableExtra::kable_styling(latex_options = "striped") %>%
  kableExtra::kable_styling(latex_options="scale_down")


save.image("REnvironments/PGLSdiversificationAndSympatry.RData")