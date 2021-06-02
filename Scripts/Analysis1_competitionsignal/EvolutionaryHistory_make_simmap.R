####---------------------------------
#### EvolutionaryHistoryClusterRun
####---------------------------------

# This script allows running one the cluster the different scenario (non-competitive vs. competitive) of brain size evolution in frugivorous primates. 
# This is the example script for low frugivory and low folivory threshold i.e.:

# a=1
# frugivoryThreshold=frugivoryThresholdVector[a]

# b=1
# folivoryThreshold=folivoryThresholdVector[b] 


##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###Set working directory
setwd("/Users/bperez/Documents/GitHub/Meta_analysis_cognition_primates/Processed_data/")

dir.create(file.path("Sample_size"), showWarnings = FALSE)
dir.create(file.path("extdata"), showWarnings = FALSE)
dir.create(file.path("OutputEvolModel"), showWarnings = FALSE)
dir.create(file.path("Dataplot"), showWarnings = FALSE)

#Import environment
rm(list=ls())

load("../BioGeoBEARS/geography_traits_biogeobears.RData")


#Libraries

# require(devtools)
# install_version("phytools", version = "0.6.99", repos = "http://cran.us.r-project.org", lib="/users/biodiv/bperez/packages/")

# Phylogenetics
# library(caper)
# library(ape)
library(phytools)
# library(geiger)
# library(MCMCglmm, lib.loc = "/users/biodiv/bperez/packages/")
# library(ellipsis, lib.loc = "/users/biodiv/bperez/packages/")
library(RPANDA)
# library(BioGeoBEARS)
# library(optimx)
# library(svMisc, lib.loc = "/users/biodiv/bperez/packages/")

#Parallelizing
# library(parallel)

##--------
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

frugivoryThresholdVector <- seq(from=20, to=40, by=20)
folivoryThresholdVector <- seq(from=40, to=60, by=20)
geographicThresholdVector <- c(10,30)/100
randomSampling=10
numberSimulations=10
numberTrees=1

totModels=randomSampling*numberSimulations*numberTrees*length(frugivoryThresholdVector)*length(folivoryThresholdVector)*length(geographicThresholdVector)


progression=0

##--------

# Set the frugivory and folivory threshold to consider:
a=1
for (a in c(1,2)){
  print(c("a: ", a))
frugivoryThreshold=frugivoryThresholdVector[a]

b=1
for (b in c(1,2)){
  print(c("b: ", b))
folivoryThreshold=folivoryThresholdVector[b] 

##------------
## Initialise output files:
##------------

#Data files
summaryBrainFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))
summaryEQFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))
summaryNeocortexFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))
summaryHippocampusFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))
summaryCerebellumFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))
summaryStriatumFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))
summaryMOBFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))
summaryOpticTractFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))

#Sample size files:
repetition=length(frugivoryThresholdVector)*length(folivoryThresholdVector)*length(geographicThresholdVector)*randomSampling
checkSampleFruit <- rep(NA, times=repetition)
checkSampleLeaf <- rep(NA, times=repetition)
checkSampleRange <- rep(NA, times=repetition)
checkSampleBrain <-  rep(NA, times=repetition)
checkSampleEQ <-  rep(NA, times=repetition)
checkSampleNeocortex <-  rep(NA, times=repetition)
checkSampleHippocampus <- rep(NA, times=repetition)
checkSampleCerebellum <- rep(NA, times=repetition)
checkSampleStriatum <- rep(NA, times=repetition)
checkSampleMOB <- rep(NA, times=repetition)
checkSampleRange <- rep(NA, times=repetition)

### Run the different evolutionary models  ###

for(c in 1:length(geographicThresholdVector)){
  
  print(c("c: ", c))

  parrallel_run <- function(d){

    
    ## Adding the co-occurence
    
    summaryData$geographicCode <- matrixRangingSensitivity[match(summaryData$SpeciesForPhylogeny,matrixRangingSensitivity$SpeciesForPhylogeny),which(thresholdPresenceRange==geographicThresholdVector[c])]   
    
    load(paste("../BioGeoBEARS/BSM_output_file", c, ".Rdata", sep=""))
    
    summaryData_init <- summaryData
    
    base::set.seed(d)
    base::set.seed(d)
    rm(.Random.seed)
    base::set.seed(d)
    base::set.seed(d)
    
    
    #####  To run with different d values 
    
    summaryData <- summaryData_init 
    
    #Source functions & packages
    
    #Capitalize first letter
    firstup <- function(x) {
      substr(x, 1, 1) <- toupper(substr(x, 1, 1))
      x
    }
    
    #####   Random sampling for covariate in case of multiple sources  ####
    
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
    
    for (i in 1:nrow(summaryData)){
      #Frugivory 
      value <- c(summaryData$Diet_frug_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                 summaryData$Diet_frug_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value)==1){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryData$FrugivoryPercent[i] <- sample(value, 1)
      }
      else{
        summaryData$FrugivoryPercent[i] <- NA
      }
      
      #Folivory 
      value <- c(summaryData$Diet_leaves_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                 summaryData$Diet_leaves_willems[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value)==1){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryData$FolivoryPercent[i] <- sample(value, 1)
      }
      else{
        summaryData$FolivoryPercent[i] <- NA
      }
      
      ##### Associate the dietary guild
      if(!is.na(as.numeric(as.character(summaryData$FrugivoryPercent[i])))&as.numeric(as.character(summaryData$FrugivoryPercent[i]))>=frugivoryThreshold){
        summaryData$DietaryGuild[i] <- "Fruit"
      } else if(!is.na(as.numeric(as.character(summaryData$FolivoryPercent[i])))&as.numeric(as.character(summaryData$FolivoryPercent[i]))>=folivoryThreshold){
        summaryData$DietaryGuild[i] <- "Leaf"
      }
      else if(!is.na(summaryData$Guild_init_decasien[i])){
        summaryData$DietaryGuild[i] <- summaryData$Guild_init_decasien[i]
        if(summaryData$DietaryGuild[i] =="Om"|summaryData$DietaryGuild[i]=="Frug/Fol"){
          summaryData$DietaryGuild[i] <- "Fruit"
        }
      }
      else {
        summaryData$DietaryGuild[i] <- "Other"
      }
      
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
      }
      else{
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
      }
      else{
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
      }
      else{
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
      }
      else{
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
      }
      else{
        summaryData$Hippocampus[i] <- NA
      }
      
      #MOB
      
      value <- c(summaryData$MOB_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value)==1){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryData$MOB[i] <- sample(value, 1)
      }
      else{
        summaryData$MOB[i] <- NA
      }
      
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
      }
      else{
        summaryData$Bodymass[i] <- NA
      }
      
    }
    
    summaryData$Family<- summaryData$Family[match(summaryData$Species_abbrv,summaryData$Species_abbrv)]
    
    
    #save data for plotting
    if(a==1&b==1&c==1&d==1){
      summaryDataForPlot <- summaryData
      write.table(summaryData, "Dataplot.txt", row.names=FALSE, col.names=TRUE, sep="\t")
    }
    
    #SaveDataGeneral
    write.table(summaryData, paste("Dataplot/Dataplot", a, "_", b, "_", c, "_", d,".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    

    ##--------
    # Evolutionary history of diet
    ##--------
    
    
    vectorDiet <-  summaryData$DietaryGuild
    names(vectorDiet) <-  summaryData$SpeciesForPhylogeny
    vectorDiet <- vectorDiet[vectorDiet!="Other"&!is.na(summaryData$geographicCode)]
    
  
    #Load and save tree corresponding to species with diet
    phylo <- read.tree("../Raw_data/Tree/Tree_biogeobears.nex")
    
    phylo <- drop.tip(phylo,
                      phylo$tip.label[
                        which(phylo$tip.label
                              %nin%names(vectorDiet))]) 
    simmapdiet1 <- make.simmap(tree=phylo, vectorDiet, model="ARD", pi="estimated", nsim=numberSimulations) #inequal and not symmetrical rate of transition from folivory to frugivory etc...
    
    save(list = "simmapdiet1", file=paste("Simmap/Output_simmap_transition", a, "_", b, "_", c, "_", d, ".Rdata", sep=""), version=2)
    
  }
  
  # Run parallel
  #parallel::mclapply(1:randomSampling, parrallel_run, mc.cores=10, mc.preschedule = T)

  for (d in 1:10){
    print(c("d: ", d))
    parrallel_run(d)
  }
}

}
}
