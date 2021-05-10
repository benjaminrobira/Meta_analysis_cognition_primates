##############################
## BRAIN EVOLUTIONARY HISTORY
##############################

#---------------------------------
##This part only runs the (non-)competitive models of brain size evolution in frugivorous primates.
#---------------------------------

setwd("//10.8.4.232/shared-2/robira/Metaanalysescognition/Evolutionary_history")

###----
##Load data

load("geography_traits_biogeobears.RData")

##-----
## Run the different evolutionary models

for(a in 1:length(frugivoryThresholdVector)){
  
  frugivoryThreshold=frugivoryThresholdVector[a]
  
  for(b in 1:length(folivoryThresholdVector)){
    
    folivoryThreshold=folivoryThresholdVector[b] 
    
    for(c in 1:length(geographicThresholdVector)){
      
      ###------------------
      ##Adding the co-occurence
      ###------------------
      
      summaryData$geographicCode <- matrixRangingSensitivity[match(summaryData$SpeciesForPhylogeny,matrixRangingSensitivity$SpeciesForPhylogeny),which(thresholdPresenceRange==geographicThresholdVector[c])]   
      
      load(paste("BioGeoBEARS/BSM_output_file", c, ".Rdata", sep=""))
      
      ##--------
      #for(e in 1:numberTrees){
      
      #Parallelizing loop
      
      #for(d in 1:randomSampling){
      cores=detectCores()
      cl <- makeCluster(cores[1]-2) #not to overload your computer
      registerDoParallel(cl)
      foreach(d=1:randomSampling#, .packages=c('caper',
              #'RPANDA',
              #'BioGeoBEARS',
              #'phytools',
              #'ape',
              #'phytools')
              ) %dopar% {
                
                setwd("//10.8.4.232/shared-2/robira/Metaanalysescognition/Evolutionary_history")
                number= (a-1)*randomSampling + (b-1)*randomSampling + (c-1)*randomSampling + d-1
                
                #Source functions & packages
                
                ###
                #Functions
                #To source all phylogenetics functions (biogeobears + models of evolution)
                source("//10.8.4.232/shared-2/robira/Metaanalysescognition/Evolutionary_history/Functions.R")
                
                #My toolkit
                source("//10.8.4.232/shared-2/robira/Metaanalysescognition/Evolutionary_history/toolbox.R")
                
                #Capitalize first letter
                firstup <- function(x) {
                  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
                  x
                }
                ###
                
                #Phylogenetics
                library(caper)
                library(MCMCglmm)
                library(RPANDA)
                library(BioGeoBEARS)
                library(phytools)
                library(ape)
                library(phytools)
                library(geiger)
                library(optimx)
                
                #####---------------
                #Random sampling for covariate in case of multiple sources
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
                  
                  ##--------------
                  #Associate the dietary guild
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
                  value <- c(
                    summaryData$MOB_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
                  
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
                
                #write.table(summaryData, paste("OutputEvolModel/Data", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
                #save data for plotting
                
                if(a==1&b==1&c==1&d==1){
                  summaryDataForPlot <- summaryData
                  write.table(summaryData, "OutputEvolModel/Dataplot.txt", row.names=FALSE, col.names=TRUE, sep="\t")
                }
                #Sample size
                counter=counter+1
                checkSampleFruit[counter] <- length(summaryData$DietaryGuild[summaryData$DietaryGuild=="Fruit"])
                checkSampleLeaf[counter]  <- length(summaryData$DietaryGuild[summaryData$DietaryGuild=="Leaf"])
                checkSampleRange[counter]  <- length(summaryData$geographicCode[!is.na(summaryData$geographicCode)])
                checkSampleBrain[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&summaryData$DietaryGuild=="Fruit",])
                checkSampleEQ[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Bodymass)&summaryData$DietaryGuild=="Fruit",])
                checkSampleNeocortex[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Neocortex)&summaryData$DietaryGuild=="Fruit",]) 
                checkSampleHippocampus[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Hippocampus)&summaryData$DietaryGuild=="Fruit",])
                checkSampleCerebellum[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Cerebellum)&summaryData$DietaryGuild=="Fruit",])
                checkSampleStriatum[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Striatum)&summaryData$DietaryGuild=="Fruit",])
                checkSampleMOB[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$MOB)&summaryData$DietaryGuild=="Fruit",])
                #checkSampleOptic_tract[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Optic.Tract)&summaryData$DietaryGuild=="Fruit",])#too low sample
                
                ##--------
                # Evolutionary history of diet
                ##--------
                
                #library(RPANDA)
                #write.table(summaryData, paste("OutputEvolModel/Data", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
                
                vectorDiet <-  summaryData$DietaryGuild
                names(vectorDiet) <-  summaryData$SpeciesForPhylogeny
                vectorDiet <- vectorDiet[vectorDiet!="Other"&!is.na(summaryData$geographicCode)]
                #Load and save tree corresponding to species with diet
                options(warn=1)
                
                phylo <- read.tree("//10.8.4.232/shared-2/robira/Metaanalysescognition/Evolutionary_history/Tree/Tree_biogeobears.nex")
                
                phylo <- drop.tip(phylo,
                                  phylo$tip.label[
                                    which(phylo$tip.label
                                          %nin%names(vectorDiet))]) 
                simmapdiet1 <- make.simmap(tree=phylo, vectorDiet, model="ARD", pi="estimated", nsim=numberSimulations)#inequal and not symmetrical rate of transition from folivory to frugivory etc...
                #simmapdiet2 <- make.simmap(tree=phylo, vectorDiet, model="ARD", pi="equal", nsim=numberSimulations)#inequal and not symmetrical rate of transition from folivory to frugivory etc...
                #simmapdiet3 <- make.simmap(tree=phylo, vectorDiet, model="ARD", pi=c(0,1), nsim=numberSimulations)#inequal and not symmetrical rate of transition from folivory to frugivory etc...
                
                #QtransitionRateEst[which(!is.na(QtransitionRateEst[,1]))[1],] <- as.vector(simmapdiet1[[1]]$Q[,1])
                write.table(as.vector(simmapdiet1[[1]]$Q[,1]), paste("OutputEvolModel/Output_simmap_transition", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
                #QtransitionRateEqual[which(!is.na(QtransitionRateEqual[,1]))[1],] <- as.vector(simmapdiet2[[1]]$Q[,1])
                #QtransitionRatePrior[which(!is.na(QtransitionRatePrior[,1]))[1],] <- as.vector(simmapdiet3[[1]]$Q[,1])
                
                # Evolutionary history of traits (~brain size) with and without competition
                #following Drury et al.'s approach to see whether one trait (here one brain area size) phylogenetical history is better described if considering competition
                #https://academic.oup.com/sysbio/article/65/4/700/1753588
                #https://journals.plos.org/plosbiology/article?rev=1&id=10.1371/journal.pbio.2003563#pbio.2003563.ref028
                ##--------
                
                #Create variable of rinterest
                summaryData$ratioBrain <- summaryData$Brain*1.036*(10**-3)/summaryData$Bodymass #Following decasien for multiplication by 1.036
                summaryData$EQ <- summaryData$Brain*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775) #Following decasien, according to #Jerison, H. J. Evolution of the Brain and Intelligence (Academic, 1973).
                
                #Make it having symmetrical (and if possible gaussian distribution, since it seems to be a prerequisite of the analysis)
                hist(summaryData$ratioBrain)
                hist(log(summaryData$ratioBrain))#good enough normal
                #summaryData$ratioBrain.log <- log(summaryData$ratioBrain)
                
                hist(summaryData$EQ)
                #summaryData$EQ.log <- log(summaryData$EQ)
                #hist(summaryData$EQ.log)
                
                hist(summaryData$Brain)
                summaryData$Brain.log <- log(summaryData$Brain)
                hist(summaryData$Brain.log)
                
                # cor.test(summaryData$EQ, summaryData$Brain)
                # cor.test(summaryData$EQ, summaryData$ratioBrain)
                # cor.test(summaryData$Brain, summaryData$ratioBrain)
                
                #Reload tree to have same than used for biogeobears
                phylo <- read.tree("//10.8.4.232/shared-2/robira/Metaanalysescognition/Evolutionary_history/Tree/Tree_biogeobears.nex")
                
                colnames(summaryData)[colnames(summaryData)=="DietaryGuild"] <- "Guild"
                # tableDrury <- cbind(summaryData$SpeciesForPhylogeny,summaryData$EQ,summaryData$Guild,summaryData$geographicCode) 
                # colnames(tableDrury) <- c("SpeciesForPhylogeny", "Trait", "Guild", "geographicCode")
                # 
                # write.table(tableDrury , "table_data_rdc_jonathan_drury.txt", sep="\t")
                # 
                resultBrainFrugivory <- runComparisonModelsCompetition(
                  simmap=simmapdiet1,
                  data=summaryData[!is.na(summaryData$Brain.log)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
                  subgroup="Fruit",
                  numberMaps=numberSimulations,
                  trait="Brain.log",
                  tree=phylo,
                  ana_events_tables=BSM_output$RES_ana_events_tables,
                  clado_events_tables=BSM_output$RES_clado_events_tables
                )
                
                resultEQFrugivory <- runComparisonModelsCompetition(
                  simmap=simmapdiet1,
                  data=summaryData[!is.na(summaryData$EQ)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
                  subgroup="Fruit",
                  numberMaps=numberSimulations,
                  trait="EQ",
                  tree=phylo,
                  ana_events_tables=BSM_output$RES_ana_events_tables,
                  clado_events_tables=BSM_output$RES_clado_events_tables
                )
                
                #Neocortex
                summaryData$ratioNeocortex <- summaryData$Neocortex/ summaryData$Brain
                hist(summaryData$ratioNeocortex )
                #summaryData$ratioNeocortex.log <- log(summaryData$ratioNeocortex)
                #hist(summaryData$ratioNeocortex.log)
                
                resultNeocortexFrugivory <- runComparisonModelsCompetition(
                  simmap=simmapdiet1,
                  data=summaryData[!is.na(summaryData$ratioNeocortex)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
                  subgroup="Fruit",
                  numberMaps=numberSimulations,
                  trait="ratioNeocortex",
                  tree=phylo,
                  ana_events_tables=BSM_output$RES_ana_events_tables,
                  clado_events_tables=BSM_output$RES_clado_events_tables
                )
                
                #Hippocampus
                summaryData$ratioHippocampus <- summaryData$Hippocampus/ summaryData$Brain
                hist(summaryData$ratioHippocampus )
                
                resultHippocampusFrugivory <- runComparisonModelsCompetition(
                  simmap=simmapdiet1,
                  data=summaryData[!is.na(summaryData$ratioHippocampus)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
                  subgroup="Fruit",
                  numberMaps=numberSimulations,
                  trait="ratioHippocampus",
                  tree=tree,
                  ana_events_tables=BSM_output$RES_ana_events_tables,
                  clado_events_tables=BSM_output$RES_clado_events_tables
                )
                
                #Cerebellum
                summaryData$ratioCerebellum <- summaryData$Cerebellum/ summaryData$Brain
                hist(summaryData$ratioCerebellum )
                #summaryData$ratioCerebellum.log <- log(summaryData$ratioCerebellum)
                #hist(summaryData$ratioCerebellum.log)
                
                resultCerebellumFrugivory <- runComparisonModelsCompetition(
                  simmap=simmapdiet1,
                  data=summaryData[!is.na(summaryData$ratioCerebellum)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
                  subgroup="Fruit",
                  numberMaps=numberSimulations,
                  trait="ratioCerebellum",
                  tree=phylo,
                  ana_events_tables=BSM_output$RES_ana_events_tables,
                  clado_events_tables=BSM_output$RES_clado_events_tables
                )
                
                #Striatum
                summaryData$ratioStriatum <- summaryData$Striatum/ summaryData$Brain
                hist(summaryData$ratioStriatum)
                
                resultStriatumFrugivory <- runComparisonModelsCompetition(
                  simmap=simmapdiet1,
                  data=summaryData[!is.na(summaryData$ratioStriatum)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
                  subgroup="Fruit",
                  numberMaps=numberSimulations,
                  trait="ratioStriatum",
                  tree=phylo,
                  ana_events_tables=BSM_output$RES_ana_events_tables,
                  clado_events_tables=BSM_output$RES_clado_events_tables
                )
                
                #MOB
                summaryData$ratioMOB <- summaryData$MOB/ summaryData$Brain
                hist(summaryData$ratioMOB)
                summaryData$ratioMOB.log <- log(summaryData$ratioMOB)
                hist(summaryData$ratioMOB.log)
                
                resultMOBFrugivory <- runComparisonModelsCompetition(
                  simmap=simmapdiet1,
                  data=summaryData[!is.na(summaryData$ratioMOB)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
                  subgroup="Fruit",
                  numberMaps=numberSimulations,
                  trait="ratioMOB.log",
                  tree=phylo,
                  ana_events_tables=BSM_output$RES_ana_events_tables,
                  clado_events_tables=BSM_output$RES_clado_events_tables
                )
                
                start=which(is.na(summaryBrainFrugivory[,1]))[1]
                end=which(is.na(summaryBrainFrugivory[,1]))[1] + numberSimulations - 1
                
                summaryBrainFrugivory[start:end,1:5] <-
                  resultBrainFrugivory[,(ncol(resultBrainFrugivory)-4):ncol(resultBrainFrugivory)]
                summaryEQFrugivory[start:end,1:5] <-
                  resultEQFrugivory[,(ncol(resultEQFrugivory)-4):ncol(resultEQFrugivory)]
                summaryNeocortexFrugivory[start:end,1:5] <-
                  resultNeocortexFrugivory[, (ncol(resultNeocortexFrugivory)-4):ncol(resultNeocortexFrugivory)]
                summaryHippocampusFrugivory[start:end,1:5] <-
                  resultHippocampusFrugivory[, (ncol(resultHippocampusFrugivory)-4):ncol(resultHippocampusFrugivory)]
                summaryCerebellumFrugivory[start:end,1:5] <-
                  resultCerebellumFrugivory[, (ncol(resultCerebellumFrugivory)-4):ncol(resultCerebellumFrugivory)]
                summaryStriatumFrugivory[start:end,1:5] <-
                  resultStriatumFrugivory[, (ncol(resultStriatumFrugivory)-4):ncol(resultStriatumFrugivory)]
                summaryMOBFrugivory[start:end,1:5] <-
                  resultMOBFrugivory[, (ncol(resultMOBFrugivory)-4):ncol(resultMOBFrugivory)]
                summaryOpticTractFrugivory[start:end,1:5] <-
                  resultOptic.TractFrugivory[, (ncol(resultOptic.TractFrugivory)-4):ncol(resultOptic.TractFrugivory)]
                
                
                #Part added because now parallelizing code // to avoid to many changes
                summaryBrainFrugivory <- summaryBrainFrugivory[start:end,]
                summaryEQFrugivory <- summaryEQFrugivory[start:end,]
                summaryNeocortexFrugivory <- summaryNeocortexFrugivory[start:end,]
                summaryHippocampusFrugivory <- summaryHippocampusFrugivory[start:end,]
                summaryCerebellumFrugivory <- summaryCerebellumFrugivory[start:end,]
                summaryStriatumFrugivory <- summaryStriatumFrugivory[start:end,]
                summaryMOBFrugivory <- summaryMOBFrugivory[start:end,]
                
                write.table(summaryBrainFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Brain", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
                write.table(summaryEQFrugivory, paste("OutputEvolModel/Output_evolutionary_history_EQ", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
                write.table(summaryNeocortexFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Neocortex", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
                write.table(summaryHippocampusFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Hippocampus", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
                write.table(summaryCerebellumFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Cerebellum", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
                write.table(summaryStriatumFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Striatum", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
                write.table(summaryMOBFrugivory, paste("OutputEvolModel/Output_evolutionary_history_MOB", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")

              }
      #stop cluster
      stopCluster(cl)
      #}
    }
  }
}

save.image("geography_trait_models.RData")
