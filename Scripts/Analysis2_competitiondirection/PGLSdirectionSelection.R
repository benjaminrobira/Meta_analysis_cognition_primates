####---------------------------------
#### PGLSdirectionSelection
####---------------------------------

# This script allows to perform PGLS analysis to assess the correlation of brain size with spatial co-occurrence data, specifically in areas for which competitive models were shown to fit the best.

# THIS IS STILL IN PROGRESS; FURTHER CLEANING/COMPLETION NEEDED

###Set working directory
setwd("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates")

#Import environment
rm(list=ls())
load("Data_spatial_primate.RData")
load("geography_traits.RData")

##--------
#Home made functions
#To source all phylogenetics functions (biogeobears + models of evolution)
source("~/PhD/Meta_analysis/Cognition_metaanalysis/Functions.R")

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

##--------

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


################
# Fitting linear models

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
summaryData$ratioNeocortex <- summaryData$Neocortex/summaryData$Brain
summaryData$ratioHippocampus <- summaryData$Hippocampus/summaryData$Brain
summaryData$ratioCerebellum <- summaryData$Cerebellum/summaryData$Brain
summaryData$ratioStriatum <- summaryData$Striatum/summaryData$Brain
summaryData$ratioMOB <- summaryData$MOB/summaryData$Brain
summaryData$ratioMOB.log <- log(summaryData$ratioMOB)
  
##---
#Load consensus tree

phylo_all <-read.nexus("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3.nex")
phylo_init <- phylo_all

phylo <- force.ultrametric(tree=phylo_init, method="extend")#method="nnls")
is.ultrametric(phylo)

#check if this differ (was done outside the loop and was ok)
plot(phylo_init$edge.length, phylo$edge.length, xlab="Original tree", ylab="Post-chronos tree", pch=20, bty="n")
abline(a=0, b=1, col="gray", lwd=0.5)

phylo <- multi2di(phylo)
is.binary(phylo)


traitToStudy=c("EQ", "ratioBrain", "ratioHippocampus", "ratioCerebellum", "ratioStriatum", "ratioMOB.log")
traitName=c("EQ", "Brain", "Hippocampus", "Cerebellum", "Striatum", "MOB")

results.df <- as.data.frame(matrix(NA, nrow=5*length(traitToStudy), ncol=8))
colnames(results.df) <- c("", "Est.", "CI2.5%", "CI97.5%", "Sd", "Df", "F", "p-value")

##----
## Run analysis for each trait

for(a in 1:length(traitToStudy)){
#Matching the brain trait to the dataset with predictors
dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[a])]
dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,2])&!is.na(dataRangePrimate[,3])&!is.na(dataRangePrimate[,4]),]

#Readjust phylo tree
phyloConsensus <- drop.tip(phylo,
                  phylo$tip.label[
                    which(phylo$tip.label
                          %nin%dataRangePrimate_rdc$Species)])

#Studied trait + sample size
firstRowWrite <- which(is.na(results.df[,1]))[1]
results.df[firstRowWrite, 1] <- paste(traitName[a], " (N=", nrow(dataRangePrimate_rdc), ")", sep="") 
results.df[(firstRowWrite+1):(firstRowWrite+4), 1] <- c("Intercept", "% of overlapped home range", "Number of sympatric frugivorous (sqrt)", "λ") 

##--------
## Main model: based on consensus tree
##--------


# #Create data for phylogenetic regression
# comp_data <- comparative.data(phy = phyloConsensus, data= dataRangePrimate_rdc,
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
# comp_data$data$Overlap_average.sqrt <- sqrt(comp_data$data$Overlap_average)

dataRangePrimate_rdc$Number_species_cooccurrence.sqrt <- sqrt(dataRangePrimate_rdc$Number_species_cooccurrence)
dataRangePrimate_rdc$Overlap_average.sqrt <- sqrt(dataRangePrimate_rdc$Overlap_average)
rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species

modelBrain <- phylolm(formula = Trait ~ Overlap_average.sqrt + Number_species_cooccurrence.sqrt, 
                      data=dataRangePrimate_rdc ,phy=phyloConsensus,model="lambda",measurement_error=FALSE,boot=repetitionBootstrap)

summary(modelBrain)$coefficients
results <- anova.pgls(modelBrain)

#Complete df results
results.df[(firstRowWrite+1):(firstRowWrite+3), 1] <- round(summary(modelBrain)$coefficients[,1], digit=4)
results.df[(firstRowWrite+1):(firstRowWrite+3), 2] <- round(summary(modelBrain)$coefficients[,2], digit=4)
results.df[(firstRowWrite+1):(firstRowWrite+3), 3] <- round(summary(modelBrain)$coefficients[,5], digit=4)
results.df[(firstRowWrite+1):(firstRowWrite+3), 4] <- round(summary(modelBrain)$coefficients[,6], digit=4)
results.df[(firstRowWrite+1):(firstRowWrite+3), 5] <- c("ns", results[1:2,1])
results.df[(firstRowWrite+1):(firstRowWrite+3), 6] <- c("ns", round(results[1:2,4], digit=4))
results.df[(firstRowWrite+1):(firstRowWrite+3), 7] <- c("ns", round(results[1:2,5], digit=3))

val1 <- ifelse(is.na(summary(modelBrain)$param.CI$lambda$ci.val[1]), 0, round(summary(modelBrain)$param.CI$lambda$ci.val[1], digit=2))
val2 <- ifelse(is.na(summary(modelBrain)$param.CI$lambda$ci.val[2]), 1, round(summary(modelBrain)$param.CI$lambda$ci.val[2], digit=2))

lambda <- paste(round(summary(modelBrain)$param.CI$lambda$opt, digit=2), " [",val1,",",val2,"]", sep="")
        
results.df[firstRowWrite+4,1] <- lambda

#Check assumptions
# pdf(file=paste("Figure/Diagnostics", traitName[a], ".pdf", sep=""))
diagnostics.plot(modelBrain)
# dev.off()
}

######
## REMAINS TO DO: VIF, DFBETAS ETC...


##############
## PLOTTING RESULTS
##############

#Remaining to do in plotting:

#1) create layout + matrix for all plots

#finish the plot part:write for overlap + adjust axis and color + add regression


layout(mat=rbind(1:length(traitToStudy), 
                 (length(traitToStudy)+1):(length(traitToStudy)+1+length(traitToStudy))),
       widths=c(5,5), heights=rep(5, times=length(traitToStudy)))
par(mar=c(3, 3, 1, 1), mgp=c(2, 0.5, 0), xpd=TRUE)
  
repetitionBootstrap=2000

for(i in 1:length(traitToStudy)){ 
  
  #Matching the brain trait to the dataset with predictors
  dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[i])]
  dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,2])&!is.na(dataRangePrimate[,3])&!is.na(dataRangePrimate[,4]),]
  
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
  rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species
  modelBrain <- phylolm(formula = Trait ~ Overlap_average.sqrt + Number_species_cooccurrence.sqrt, 
          data=dataRangePrimate_rdc ,phy=phyloConsensus,model="lambda",measurement_error=FALSE,boot=repetitionBootstrap)
  
  ##----
  #Plot against N co-occ

  rownames(dataRangePrimate_rdc) <- dataRangePrimate_rdc$Species
  
  xmin=min(round(dataRangePrimate_rdc[,c(2)]), digit=2)
  xmax=max(round(dataRangePrimate_rdc[,c(2)]), digit=2)
  ymin=min(round(dataRangePrimate_rdc[,c(4)]), digit=2)
  ymax=max(round(dataRangePrimate_rdc[,c(4)]), digit=2)
  
  ##With number of co-occurring species
  plot(dataRangePrimate_rdc[,c(2)], dataRangePrimate_rdc[,c(4)], xlab="Trait 1", ylab="Trait 2",
       las=1, type="n", tcl=-0.25, bty="n",
       xaxt="n",xaxs="i",yaxs="i", yaxt="n", xpd=TRUE)
  
  #Add grid
  addGrid(
    xmin=xmin, xmax=xmax, xintsmall=(xmax-xmin)/20, xintbig=(xmax-xmin)/5,
    ymin=ymin, ymax=ymax, yintsmall=(ymax-ymin)/20, yintbig=(ymax-ymin)/5,
    axisPlot=TRUE, round=TRUE, digit=c(2,2))
  addLabel(xfrac=0.05, yfrac=0.05, label=paste(i, "a", sep=""), circle=TRUE, radiuscircle=0.5, circle.bg="black", font.col="white")
  
  #Add background tree
  col=list(col.edge=setNames(rep("darkgrey",nrow(phyloConsensus$edge)),as.character(phyloConsensus$edge[,2])),
           col.node=setNames(rep("black",max(phyloConsensus$edge)),as.character(1:max(phyloConsensus$edge))))
  
  phylomorphospace(phyloConsensus,dataRangePrimate_rdc[,c(2,4)], add=TRUE, label="true", lty=3,
                   control=col)
  
  #Add result model
  ymean <- c(summary(modelBrain)$coefficients[1,1], summary(modelBrain)$coefficients[1,1] + summary(modelBrain)$coefficients[3,1]*max(dataRangePrimate_rdc[,c(4)]))
  ylower <- c(summary(modelBrain)$coefficients[1,4], summary(modelBrain)$coefficients[1,4] + (summary(modelBrain)$coefficients[3,4])*max(dataRangePrimate_rdc[,c(4)]))
  yupper <- c(summary(modelBrain)$coefficients[1,5], summary(modelBrain)$coefficients[1,5] + (summary(modelBrain)$coefficients[3,5])*max(dataRangePrimate_rdc[,c(4)]))
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
    
  #Overlay points
  points(dataRangePrimate_rdc[,c(2)], dataRangePrimate_rdc[,c(4)], pch=19, col="red",xpd=TRUE)
  
  ##----
  #Plot against overlap
  
  xmin=min(round(dataRangePrimate_rdc[,c(3)]), digit=2)
  xmax=max(round(dataRangePrimate_rdc[,c(3)]), digit=2)
  ymin=min(round(dataRangePrimate_rdc[,c(4)]), digit=2)
  ymax=max(round(dataRangePrimate_rdc[,c(4)]), digit=2)
  
  ##With number of co-occurring species
  plot(dataRangePrimate_rdc[,c(3)], dataRangePrimate_rdc[,c(4)], xlab="Trait 1", ylab="Trait 2",
       las=1, type="n", tcl=-0.25, bty="n",
       xaxt="n",xaxs="i",yaxs="i", yaxt="n", xpd=TRUE)
  
  #Add grid
  addGrid(
    xmin=xmin, xmax=xmax, xintsmall=(xmax-xmin)/20, xintbig=(xmax-xmin)/5,
    ymin=ymin, ymax=ymax, yintsmall=(ymax-ymin)/20, yintbig=(ymax-ymin)/5,
    axisPlot=TRUE, round=TRUE, digit=c(2,2))
  addLabel(xfrac=0.05, yfrac=0.05, label=paste(i, "b", sep=""), circle=TRUE, radiuscircle=0.5, circle.bg="black", font.col="white")
  
  #Add background tree
  col=list(col.edge=setNames(rep("darkgrey",nrow(phyloConsensus$edge)),as.character(phyloConsensus$edge[,2])),
           col.node=setNames(rep("black",max(phyloConsensus$edge)),as.character(1:max(phyloConsensus$edge))))
  
  phylomorphospace(phyloConsensus,dataRangePrimate_rdc[,c(2,4)], add=TRUE, label="true", lty=3,
                   control=col)
  
  #Add result model
  ymean <- c(summary(modelBrain)$coefficients[1,1], summary(modelBrain)$coefficients[1,1] + summary(modelBrain)$coefficients[3,1]*max(dataRangePrimate_rdc[,c(4)]))
  ylower <- c(summary(modelBrain)$coefficients[1,4], summary(modelBrain)$coefficients[1,4] + (summary(modelBrain)$coefficients[3,4])*max(dataRangePrimate_rdc[,c(4)]))
  yupper <- c(summary(modelBrain)$coefficients[1,5], summary(modelBrain)$coefficients[1,5] + (summary(modelBrain)$coefficients[3,5])*max(dataRangePrimate_rdc[,c(4)]))
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

#Apparently you can access CI following the method in https://www.pnas.org/content/pnas/110/22/9001.full.pdf?with-ds=yes but to read

##------
## Sensitivity to phylogeny: using a block of trees
##------

##Remains to do: -> summarizing the output + check if it works
  
sensitivityEstimate <- as.data.frame(matrix(NA, ncol=3, nrow=50*30*length(traitToStudy)))
sensitivityPvalue <- as.data.frame(matrix(NA, ncol=3, nrow=50*30*length(traitToStudy)))
sensitivityLambda <- rep(NA, times=50*30*length(traitToStudy))
traitVectorSensitivity <- rep(NA, times=50*30*length(traitToStudy))

repetitionTrees=50
repetitionModels=30

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
      }
      summaryData$Family<- summaryData$Family[match(summaryData$Species_abbrv,summaryData$Species_abbrv)]
      
      
      ##----
      ## RUN MODEL FOR SENSITIVITY
      
      for(a in 1:length(traitToStudy)){
        
        rowToAdd=which(is.na(sensitivityLambda))[1]
        traitVectorSensitivity[rowToAdd] <- traitToStudy[a]
        #Matching the brain trait to the dataset with predictors
        dataRangePrimate$Trait <- summaryData[match(dataRangePrimate$Species,summaryData$SpeciesForPhylogeny), which(colnames(summaryData)==traitToStudy[a])]
        dataRangePrimate_rdc <- dataRangePrimate[!is.na(dataRangePrimate[,2])&!is.na(dataRangePrimate[,3])&!is.na(dataRangePrimate[,4]),]
        
        #Readjust phylo tree
        phyloTree <- drop.tip(phylo,
                                   phylo$tip.label[
                                     which(phylo$tip.label
                                           %nin%dataRangePrimate_rdc$Species)])
        
        ##--------
        ## Main model: based on consensus tree
        ##--------
        
        
        #Create data for phylogenetic regression
        comp_data <- comparative.data(phy = phyloTree, data= dataRangePrimate_rdc,
                                      names.col = Species, vcv = TRUE)
        
        #Output distribution
        hist(comp_data$data$Trait)
        
        #Predictor distribution
        hist(sqrt(comp_data$data$Number_species_cooccurrence))
        hist(sqrt(comp_data$data$Overlap_average))
        
        #Make predictor more symmetrical
        comp_data$data$Number_species_cooccurrence.sqrt <- sqrt(comp_data$data$Number_species_cooccurrence)
        comp_data$data$Overlap_average.sqrt <- sqrt(comp_data$data$Overlap_average)
        
        comp_data$data$Number_species_cooccurrence.sqrt.z <- scale(comp_data$data$Number_species_cooccurrence.sqrt)
        comp_data$data$Overlap_average.sqrt.z <- scale(comp_data$data$Overlap_average.sqrt)
        
        modelBrain <- pgls(formula = Trait ~ Overlap_average.sqrt.z + Number_species_cooccurrence.sqrt.z , data = comp_data, lambda = "ML")
        summary(modelBrain)$coefficients
        results <- anova.pgls(modelBrain)
        
        #Complete df results
        sensitivityEstimate[rowToAdd,] <- round(summary(modelBrain)$coefficients[,1], digit=4)
        sensitivityPvalue[rowToAdd,] <- c("ns", round(results[1:2,5], digit=3))
        sensitivityLambda[rowToAdd] <- round(summary(modelBrain)$param.CI$lambda$opt, digit=2)
        traitVectorSensitivity[rowToAdd] <- traitToStudy[a]
        
    }
  }
  
#Summary sensitivity

minEst <- aggregate(sensitivityEstimate, by=list(traitVectorSensitivity), min)
maxEst <- aggregate(sensitivityEstimate, by=list(traitVectorSensitivity), max)
medianEst <- aggregate(sensitivityEstimate, by=list(traitVectorSensitivity), median)

pid <- aggregate(sensitivityPvalue, by=list(traitVectorSensitivity), FUN= function(x)
  length(x>0.05)
  )
  
lambdamin <- aggregate(sensitivityLambda, by=list(traitVectorSensitivity), min)
lambdamax <- aggregate(sensitivityLambda, by=list(traitVectorSensitivity), max)
lambdamedian <- aggregate(sensitivityLambda, by=list(traitVectorSensitivity), median)

#Plot results of sensitivity


##With number of co-occurring species
xmin=-3
xmax=3
ymin=0
ymax=length(unique(traitVectorSensitivity))+1
plot(0, 0, xlab="Estimate (scaled variables)", ylab="",
     las=1, type="n", tcl=-0.25, bty="n",
     xaxt="n",xaxs="i",yaxs="i", yaxt="n", xpd=TRUE,
     xlim=c(xmin, xmax),
     ylim=c(ymin, ymax))

#Add grid
addGrid(
  xmin=xmin, xmax=xmax, xintsmall=(xmax-xmin)/20, xintbig=(xmax-xmin)/5,
  ymin=ymin, ymax=ymax, yintsmall=(ymax-ymin)/20, yintbig=(ymax-ymin)/5,
  axisPlot=FALSE)
axis(side=1, at=seq(from=ymin, to=ymax, by=(ymax-ymin)/5), 
     labels=round(seq(from=ymin, to=ymax, by=(ymax-ymin)/5), digit=2), tcl=-0.25)

#Add model names
text(x=-1, y=1:(ymax-1), labels=unique(traitVectorSensitivity))

#Add results     

colour.v <- brewer.pal(n = 4, name = "Pastel1")
for(i in 1:length(unique(traitVectorSensitivity))){
  #Estimate model
  errorBars(location=c(i-0.3, i-0.1, i+0.1), 
            meanPt=medianEst[%NAME VARIABLE%==unique(traitVectorSensitivity)[i]],
            barValue=rep(100, times=4),
            upperBarValue=maxEst[%NAME VARIABLE%==unique(traitVectorSensitivity)[i]],
            lowerBarValue=minEst[%NAME VARIABLE%==unique(traitVectorSensitivity)[i]],
            col=colour.v[1:3],
            horiz=TRUE,
            symmetrical=FALSE)
  
  #Lambda model
  errorBars(location=c(i+0.2), 
            meanPt=lambdamedian[%NAME VARIABLE%==unique(traitVectorSensitivity)[i]],
            barValue=rep(100, times=4),
            upperBarValue=lambdamax[%NAME VARIABLE%==unique(traitVectorSensitivity)[i]],
            lowerBarValue=lambdamin[%NAME VARIABLE%==unique(traitVectorSensitivity)[i]],
            col=colour.v[4],
            horiz=TRUE,
            symmetrical=FALSE)
}

legend("bottomright", legend=c("Intercept", "Overlap", "Nspecies", "Lambda"), col=colour.v[1:4], lty=c(1,1,1,1), bty="n")
