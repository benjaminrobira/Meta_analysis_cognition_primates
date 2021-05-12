####---------------------------------
#### EvolutionaryHistoryClusterRun
####---------------------------------

# This script allows running one the ENS cluster the different scenario (non-competitive vs. competitive) of brain size evolution in frugivorous primates. 

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###Set working directory
# scp -r /Users/bperez/ownCloud/Recherche/These/ENS/Autres/Benjamin/Evolutionary_history/* bperez@jord.biologie.ens.fr:/users/biodiv/bperez/data/others/Benji/Evolutionary_history/
# setwd("/Users/bperez/ownCloud/Recherche/These/ENS/Autres/Benjamin/Evolutionary_history/")
setwd("/users/biodiv/bperez/data/others/Benji/Evolutionary_history_2/")

dir.create(file.path("Sample_size"), showWarnings = FALSE)
dir.create(file.path("extdata"), showWarnings = FALSE)
dir.create(file.path("OutputEvolModel"), showWarnings = FALSE)
dir.create(file.path("Dataplot"), showWarnings = FALSE)

#Import environment
rm(list=ls())
load("geography_traits_biogeobears.RData")

##--------
#Home made functions
#To source all phylogenetics functions (biogeobears + models of evolution)
source("Functions.R")

##--------

##--------
#Libraries

#Phylogenetics
library(caper)
library(ape)
library(phytools)
library(geiger)
library(MCMCglmm, lib.loc = "/users/biodiv/bperez/packages/")
library(ellipsis, lib.loc = "/users/biodiv/bperez/packages/")
library(RPANDA, lib.loc = "/users/biodiv/bperez/packages/")
library(BioGeoBEARS)
library(optimx)
library(svMisc, lib.loc = "/users/biodiv/bperez/packages/")

#??
library(ellipsis)

#Parallelizing
# library(snow)
# library(foreach)
# library(doParallel)
library(parallel)

##--------
##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

frugivoryThresholdVector <- seq(from=20, to=40, by=20)
folivoryThresholdVector <- seq(from=40, to=60, by=20)
geographicThresholdVector <- c(10,30)/100
randomSampling=10
numberSimulations=10
numberTrees=1

totModels=randomSampling*numberSimulations*numberTrees*length(frugivoryThresholdVector)*length(folivoryThresholdVector)*length(geographicThresholdVector)

# QtransitionRateEst <- as.data.frame(matrix(NA, ncol=2, nrow=totModels))
#QtransitionRateEqual <- as.data.frame(matrix(NA, ncol=2, nrow=totModels))
#QtransitionRatePrior <- as.data.frame(matrix(NA, ncol=2, nrow=totModels))

progression=0

##--------

# Set the frugivory and folivory threshold to consider:
a=1
frugivoryThreshold=frugivoryThresholdVector[a]

b=2
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
#checkSampleOptic_tract <- rep(NA, times=)
checkSampleRange <- rep(NA, times=repetition)

### Run the different evolutionary models  ###

for(c in 1:length(geographicThresholdVector)){
  
  ## Adding the co-occurence
  
  summaryData$geographicCode <- matrixRangingSensitivity[match(summaryData$SpeciesForPhylogeny,matrixRangingSensitivity$SpeciesForPhylogeny),which(thresholdPresenceRange==geographicThresholdVector[c])]   
  
  load(paste("BioGeoBEARS/BSM_output_file", c, ".Rdata", sep=""))
  
  summaryData_init <- summaryData
  source("Functions.R")
  source("toolbox.R")
  
  parrallel_run <- function(d){
    
    
    runComparisonModelsCompetition <- function(
      
      #######################
      ## Run the phylogenetic comparison between random evolution (BM, OU) or competitive scenario with MC= exclusion from competitive
      # Taxa = mutual exclusion, or DDlin and DDexp, which are positive or negative density dependance to taxa of same "group" (here feeding group)
      # of the evolutionary rate, either linearily (lin), or exponentially (exp)
      ######################
      
      simmap, #The simmap object for grouping
      numberMaps, #Number of simmap and BSM 
      tree, #The phylo object of the tree
      data, #Data including: species column name as inTree. This column should be name "SpeciesForPhylogeny". a "Guild" column has to specify the grouping variable.
      #Other variables can be associated, namely the trait and subgroup variable
      subgroup, #Character string indicating the group of te "Guild" column that if of interest for the analysis 
      trait, #Column name of variable of interest. String
      ana_events_tables,
      clado_events_tables
    ){ 
      
      #isolate data from subgroup of interest 
      group.map<-simmap
      data.grouped<-subset(data, Guild==subgroup)
      mass<-data.grouped[,which(colnames(data.grouped)=="SpeciesForPhylogeny")]
      names(mass)<-data.grouped[,which(colnames(data.grouped)=="SpeciesForPhylogeny")]
      nc<-name.check(tree,mass)
      data.grouped.tree<-drop.tip(tree,nc$tree_not_data)
      subdata<-data.grouped[which(data.grouped$SpeciesForPhylogeny%in%data.grouped.tree$tip.label),]
      
      #set up analyses
      
      N=numberMaps #number of stochastic maps to analyze
      
      T=trait #which trait to analyze
      res.mat<-matrix(nrow=N,ncol=53)
      colnames(res.mat)<-c("subgroup","trait","N","BM.lnL","BM.sig2","BM.z0","BM.AICc","BM.conv",
                           "OU.lnL","OU.sig2","OU.alpha","OU.z0","OU.AICc","OU.conv",
                           "EB.lnL","EB.sig2","EB.alpha","EB.z0","EB.AICc","EB.conv",
                           "MCgeo.dietmap","MCgeo.lnL","MCgeo.sig2","MCgeo.S","MCgeo.z0","MCgeo.AICc","MCgeo.conv",
                           "DDlingeo.dietmap","DDlingeo.lnL","DDlingeo.sig2","DDlingeo.b","DDlingeo.z0","DDlingeo.AICc","DDlingeo.conv",
                           "DDexpgeo.dietmap","DDexpgeo.lnL","DDexpgeo.sig2","DDexpgeo.r","DDexpgeo.z0","DDexpgeo.AICc","DDexpgeo.conv",
                           "BM.delaic","OU.delaic","EB.delaic","MCgeo.delaic","DDlingeo.delaic","DDexpgeo.delaic",
                           "BM.wi","OU.wi","EB.wi","MCgeo.wi","DDlingeo.wi","DDexpgeo.wi")
      
      #fit BM, OU, EB, MC, DDexp, and DDlin models to subgroup trait data
      #write.table("worked", "/OutputEvolModel/workedinitFunction.txt", row.names=FALSE, col.names=TRUE, sep="\t")
      
      for(i in 1:N){
        
        group.map2 <-drop.tip.simmap(group.map[[i]],
                                     group.map[[i]]$tip.label[which(!group.map[[i]]$tip.label%in%tree$tip.label)])
        j=which(colnames(subdata)==T)
        M<-subdata[,j]
        subtree<-data.grouped.tree
        names(M)<-subdata$SpeciesForPhylogeny
        M<-subset(M,M!='NA')
        
        nc<-name.check(subtree,M)
        if(is.list(nc)){
          subtree<-drop.tip(subtree,nc$tree_not_data)
        }
        
        tryCatch(
          {o2<-fitContinuous(subtree,M,model="BM")#, ncores=1)
          BM.log_lik<-o2$opt$lnL
          BM.sig2<-o2$opt$sigsq
          BM.z0<-o2$opt$z0
          BM.aicc<-o2$opt$aicc
          BM.conv<-as.numeric(tail(o2$res[,length(o2$res[1,])],n=1))
          }, error=function(e){
            cat("ERROR :",conditionMessage(e), "STEP : ", i, "MODEL : BM", "\n")
            BM.log_lik<-NA
            BM.sig2<-NA
            BM.z0<-NA
            BM.aicc<-NA
            BM.conv<-NA
          }
        )
        
        tryCatch(
          {
            o3<-fitContinuous(subtree,M,model="OU")#, ncores=1)
            OU.log_lik<-o3$opt$lnL
            OU.sig2<-o3$opt$sigsq
            OU.alpha<-o3$opt$alpha
            OU.z0<-o3$opt$z0
            OU.aicc<-o3$opt$aicc
            OU.conv<-as.numeric(tail(o3$res[,length(o3$res[1,])],n=1))
          }, error=function(e){
            cat("ERROR :",conditionMessage(e), "STEP : ", i,"MODEL : OU", "\n")
            OU.log_lik<-NA
            OU.sig2<-NA
            OU.alpha<-NA
            OU.z0<-NA
            OU.aicc<-NA
            OU.conv<-NA
          }
        )
        
        tryCatch(
          {
            o32<-fitContinuous(subtree,M,model="EB")#, ncores=1)
            EB.log_lik<-o32$opt$lnL
            EB.sig2<-o32$opt$sigsq
            EB.alpha<-o32$opt$a
            EB.z0<-o32$opt$z0
            EB.aicc<-o32$opt$aicc
            EB.conv<-as.numeric(tail(o32$res[,length(o32$res[1,])],n=1))
          }, error=function(e){
            cat("ERROR :",conditionMessage(e), "STEP : ", i,"MODEL : EB", "\n")
            EB.log_lik<-NA
            EB.sig2<-NA
            EB.alpha<-NA
            EB.z0<-NA
            EB.aicc<-NA
            EB.conv<-NA
          }
        )
        
        tryCatch(
          {
            o4<-fit_t_comp_subgroup(full.phylo=tree,
                                    ana.events=ana_events_tables[[i]],
                                    clado.events=clado_events_tables[[i]],
                                    stratified=FALSE,map=group.map2,data=M,
                                    trim.class=subgroup,model="MC",par=NULL,method="Nelder-Mead",bounds=NULL)
            MCgeo.lnL<-o4$LH
            MCgeo.sig2<-o4$sig2
            MCgeo.S<-o4$S
            MCgeo.z0<-o4$z0
            MCgeo.aicc<-o4$aicc
            MCgeo.conv<-o4$convergence
          }, error=function(e){
            cat("ERROR :",conditionMessage(e), "STEP : ", i, "MODEL : MC", "\n")
            MCgeo.lnL<-NA
            MCgeo.sig2<-NA
            MCgeo.S<-NA
            MCgeo.z0<-NA 
            MCgeo.aicc<-NA
            MCgeo.conv<-NA
          }
        )
        
        tryCatch(
          {
            o5<-fit_t_comp_subgroup(full.phylo=tree,
                                    ana.events=ana_events_tables[[i]],
                                    clado.events=clado_events_tables[[i]],
                                    stratified=FALSE,map=group.map2,data=M,trim.class=subgroup,
                                    model="DDexp",par=NULL,method="Nelder-Mead",bounds=NULL)
            DDexpgeo.lnL<-o5$LH
            DDexpgeo.sig2<-o5$sig2
            DDexpgeo.r<-o5$r
            DDexpgeo.z0<-o5$z0
            DDexpgeo.aicc<-o5$aicc
            DDexpgeo.conv<-o5$convergence
          }, error=function(e){
            cat("ERROR :",conditionMessage(e), "STEP : ", i, "MODEL : DDexp", "\n")
            DDexpgeo.lnL<-NA
            DDexpgeo.sig2<-NA
            DDexpgeo.r<-NA
            DDexpgeo.z0<-NA
            DDexpgeo.aicc<-NA
            DDexpgeo.conv<-NA
          }
        )
        
        tryCatch(
          {
            o6<-fit_t_comp_subgroup(full.phylo=tree,
                                    ana.events=ana_events_tables[[i]],
                                    clado.events=clado_events_tables[[i]],
                                    stratified=FALSE,map=group.map2,data=M,trim.class=subgroup,
                                    model="DDlin",par=NULL,method="Nelder-Mead",bounds=NULL)
            DDlingeo.lnL<-o6$LH
            DDlingeo.sig2<-o6$sig2
            DDlingeo.b<-o6$b
            DDlingeo.z0<-o6$z0
            DDlingeo.aicc<-o6$aicc
            DDlingeo.conv<-o6$convergence
          }, error=function(e){
            cat("ERROR :",conditionMessage(e), "STEP : ", i, "MODEL : DDlin", "\n")
            DDlingeo.lnL<-NA
            DDlingeo.sig2<-NA
            DDlingeo.b<-NA
            DDlingeo.z0<-NA
            DDlingeo.aicc<-NA
            DDlingeo.conv<-NA
          }
        )
        
        if(length(which(is.na(c(BM.aicc,OU.aicc,EB.aicc,MCgeo.aicc,DDlingeo.aicc,DDexpgeo.aicc))))==0){
          BM.delaic<-BM.aicc-min(BM.aicc,OU.aicc,EB.aicc,MCgeo.aicc,DDlingeo.aicc,DDexpgeo.aicc)
          OU.delaic<-OU.aicc-min(BM.aicc,OU.aicc,EB.aicc,MCgeo.aicc,DDlingeo.aicc,DDexpgeo.aicc)
          EB.delaic<-EB.aicc-min(BM.aicc,OU.aicc,EB.aicc,MCgeo.aicc,DDlingeo.aicc,DDexpgeo.aicc)
          MCgeo.delaic<-MCgeo.aicc-min(BM.aicc,OU.aicc,EB.aicc,MCgeo.aicc,DDlingeo.aicc,DDexpgeo.aicc)
          DDlingeo.delaic<-DDlingeo.aicc-min(BM.aicc,OU.aicc,EB.aicc,MCgeo.aicc,DDlingeo.aicc,DDexpgeo.aicc)
          DDexpgeo.delaic<-DDexpgeo.aicc-min(BM.aicc,OU.aicc,EB.aicc,MCgeo.aicc,DDlingeo.aicc,DDexpgeo.aicc)
          all=sum(exp(-0.5*BM.delaic),exp(-0.5*OU.delaic),exp(-0.5*EB.delaic),exp(-0.5*MCgeo.delaic),exp(-0.5*DDlingeo.delaic),exp(-0.5*DDexpgeo.delaic))
          BM.wi<-exp(-0.5*BM.delaic)/all
          OU.wi<-exp(-0.5*OU.delaic)/all
          EB.wi<-exp(-0.5*EB.delaic)/all
          MCgeo.wi<-exp(-0.5*MCgeo.delaic)/all
          DDlingeo.wi<-exp(-0.5*DDlingeo.delaic)/all
          DDexpgeo.wi<-exp(-0.5*DDexpgeo.delaic)/all
          
          MCgeo.iter <- NA
          DDlingeo.iter <- NA
          DDexpgeo.iter <- NA
          int<-c(subgroup,names(subdata)[j],length(subtree$tip.label),BM.log_lik,BM.sig2,BM.z0,BM.aicc,BM.conv,
                 OU.log_lik,OU.sig2,OU.alpha,OU.z0,OU.aicc,OU.conv,
                 EB.log_lik,EB.sig2,EB.alpha,EB.z0,EB.aicc,EB.conv,
                 MCgeo.iter,MCgeo.lnL,MCgeo.sig2,MCgeo.S,MCgeo.z0,MCgeo.aicc,MCgeo.conv,
                 DDlingeo.iter,DDlingeo.lnL,DDlingeo.sig2,DDlingeo.b,DDlingeo.z0,DDlingeo.aicc,DDlingeo.conv,
                 DDexpgeo.iter,DDexpgeo.lnL,DDexpgeo.sig2,DDexpgeo.r,DDexpgeo.z0,DDexpgeo.aicc,DDexpgeo.conv,
                 BM.delaic,OU.delaic,EB.delaic,MCgeo.delaic,DDlingeo.delaic,DDexpgeo.delaic,
                 BM.wi,OU.wi,EB.wi,MCgeo.wi,DDlingeo.wi,DDexpgeo.wi)
          
          res.mat[i,]<-int
        }
        else{
          res.mat[i,] <- rep(NA, times=ncol(res.mat))
        }
      }
      
      res.mat <- as.data.frame(res.mat)
      resm <- data.frame(res.mat)
      
      return(resm)
      
    }
    
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
    
    write.table(length(summaryData$DietaryGuild[summaryData$DietaryGuild=="Fruit"]), 
                paste("Sample_size/checkSampleFruit", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    write.table(length(summaryData$DietaryGuild[summaryData$DietaryGuild=="Leaf"]), 
                paste("Sample_size/checkSampleLeaf", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    write.table(length(summaryData$geographicCode[!is.na(summaryData$geographicCode)]), 
                paste("Sample_size/checkSampleRange", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    write.table(nrow(summaryData[!is.na(summaryData$Brain)&summaryData$DietaryGuild=="Fruit",]), 
                paste("Sample_size/checkSampleBrain", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    write.table(nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Bodymass)&summaryData$DietaryGuild=="Fruit",]), 
                paste("Sample_size/checkSampleEQ", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    write.table(nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Neocortex)&summaryData$DietaryGuild=="Fruit",]), 
                paste("Sample_size/checkSampleNeocortex", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    write.table(nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Hippocampus)&summaryData$DietaryGuild=="Fruit",]), 
                paste("Sample_size/checkSampleHippocampus", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    write.table(nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Cerebellum)&summaryData$DietaryGuild=="Fruit",]), 
                paste("Sample_size/checkSampleCerebellum", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    write.table(nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Striatum)&summaryData$DietaryGuild=="Fruit",]), 
                paste("Sample_size/checkSampleStriatum", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    write.table(nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$MOB)&summaryData$DietaryGuild=="Fruit",]), 
                paste("Sample_size/checkSampleMOB", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    
    ##--------
    # Evolutionary history of diet
    ##--------
    
    #library(RPANDA)
    #write.table(summaryData, paste("OutputEvolModel/Data", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    
    vectorDiet <-  summaryData$DietaryGuild
    names(vectorDiet) <-  summaryData$SpeciesForPhylogeny
    vectorDiet <- vectorDiet[vectorDiet!="Other"&!is.na(summaryData$geographicCode)]
    #Load and save tree corresponding to species with diet
    #    options(warn=1)
    
    phylo <- read.tree("Tree/Tree_biogeobears.nex")
    
    phylo <- drop.tip(phylo,
                      phylo$tip.label[
                        which(phylo$tip.label
                              %nin%names(vectorDiet))]) 
    simmapdiet1 <- make.simmap(tree=phylo, vectorDiet, model="ARD", pi="estimated", nsim=numberSimulations)#inequal and not symmetrical rate of transition from folivory to frugivory etc...
    #simmapdiet2 <- make.simmap(tree=phylo, vectorDiet, model="ARD", pi="equal", nsim=numberSimulations)#inequal and not symmetrical rate of transition from folivory to frugivory etc...
    #simmapdiet3 <- make.simmap(tree=phylo, vectorDiet, model="ARD", pi=c(0,1), nsim=numberSimulations)#inequal and not symmetrical rate of transition from folivory to frugivory etc...
    
    #QtransitionRateEst[which(!is.na(QtransitionRateEst[,1]))[1],] <- as.vector(simmapdiet1[[1]]$Q[,1])
    write.table(as.vector(simmapdiet1[[1]]$Q[,1]), paste("OutputEvolModel/Output_simmap_transition", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
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
    phylo <- read.tree("Tree/Tree_biogeobears.nex")
    
    colnames(summaryData)[colnames(summaryData)=="DietaryGuild"] <- "Guild"
    
    print("Brain")
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
    write.table(resultBrainFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Brain", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    
    print("EQ")
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
    write.table(resultEQFrugivory, paste("OutputEvolModel/Output_evolutionary_history_EQ", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    
    
    #Neocortex
    print("Neocortex")
    summaryData$ratioNeocortex <- summaryData$Neocortex/ summaryData$Brain
    hist(summaryData$ratioNeocortex )
    
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
    write.table(resultNeocortexFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Neocortex", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    
    #Hippocampus
    print("ratioHippocampus")
    summaryData$ratioHippocampus <- summaryData$Hippocampus/ summaryData$Brain
    hist(summaryData$ratioHippocampus )
    
    resultHippocampusFrugivory <- runComparisonModelsCompetition(
      simmap=simmapdiet1,
      data=summaryData[!is.na(summaryData$ratioHippocampus)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
      subgroup="Fruit",
      numberMaps=numberSimulations,
      trait="ratioHippocampus",
      tree=phylo,
      ana_events_tables=BSM_output$RES_ana_events_tables,
      clado_events_tables=BSM_output$RES_clado_events_tables
    )
    write.table(resultHippocampusFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Hippocampus", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    
    #Cerebellum
    print("ratioCerebellum")
    summaryData$ratioCerebellum <- summaryData$Cerebellum/ summaryData$Brain
    hist(summaryData$ratioCerebellum )
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
    write.table(resultCerebellumFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Cerebellum", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    
    
    #Striatum
    print("ratioStriatum")
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
    write.table(resultStriatumFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Striatum", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    
    
    
    #MOB
    print("MOB")
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
    write.table(resultMOBFrugivory, paste("OutputEvolModel/Output_evolutionary_history_MOB", a, "_", b, "_", c, "_", d, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
    
  }
  
  # Run parallel
  mclapply(1:randomSampling, parrallel_run, mc.cores=10, mc.preschedule = F)
  
  # for (test_rep in 1:randomSampling){
  #   print(test_rep)
  #   parrallel_run(test_rep)
  # }
  
}

## END ANALYSIS
##----------------------

##--
## Saving environment
save.image("geography_trait_models2.RData")