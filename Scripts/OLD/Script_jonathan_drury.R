########################################################################################################################################
#                                                  SCRIPT FOR JONATHAN DRURY
########################################################################################################################################

#This scripts will lead to the error as encountered. Further details below.

##---------------
##NOTE

#To make the script work, load the result of the biogeobears (BSM_output_file.RData)
#Run the simmap
#Run the function for evolutionary models. The error that will be printed is mine: it will indicate which model failed. I provide another file "function_RPANDA", where I copied
#functions from RPANDA to track down the error. These functions include the corrections I have mentioned in my first email.

##--------------

##--------------
## Libraries
library(caper)
library(phytools)
library(BioGeoBEARS)
library(RPANDA)
##--------------

##--------------
##Functions

#To: run biogeobears, clean trees, run evolutionary models
source("Functions.R")

##--------------

##--------------
## Environments

#Biogeobearsoutput
load("BSM_output_file2.Rdata")
##--------------

##--------------
#Load data table (reduced for including only variables of interest, so as to be the clearest possible. Please email me for further details.

data.df <- read.delim("table_data_rdc_jonathan_drury.txt")


##------------------------------------------------------------------------------------------------

################
## ANALYSIS
################

##--------------
#Run simmap

#data vector
vectorDiet <-  data.df$Guild
names(vectorDiet) <-  data.df$SpeciesForPhylogeny
vectorDiet <- vectorDiet[vectorDiet!="Other"&!is.na(data.df$geographicCode)]

# #Load and save tree corresponding to species with diet
# phylo_all <-read.nexus("consensusTree_10kTrees_Primates_Version3.nex")
# phylo_init <- phylo_all
# 
# #ultrametric?
# phylo <- force.ultrametric(tree=phylo_init, method="extend")#method="nnls")
# is.ultrametric(phylo)
# #check if this differs (was done outside the loop and was ok)
# plot(phylo_init$edge.length, phylo$edge.length, xlab="Original tree", ylab="Post-chronos tree", pch=20, bty="n")
# abline(a=0, b=1, col="gray", lwd=0.5)
# #binary?
# phylo <- multi2di(phylo)
# is.binary(phylo)
# #remove useless species
# '%nin%' <- Negate('%in%')
# phylo <- drop.tip(phylo,
#                   phylo$tip.label[
#                     which(phylo$tip.label
#                           %nin%data.df$SpeciesForPhylogeny[data.df$Guild!="Other"])]) 
# #reload for branch numerotation
# write.tree(phylo, np("Tree_diet.nex"))
# phylo <- read.tree("Tree_diet.nex")
# #clean it (function from Matzke)
# phylo <- cleanTree(tr=phylo, trfn="Tree_diet.nex")
# phylo <- force.ultrametric(tree=phylo, method="extend")#method="nnls")
# #last check
# is.ultrametric(phylo)
# is.binary(phylo)
# #reload
# write.tree(phylo, np("Tree_diet.nex"))
phylo <- read.tree("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/Tree_biogeobears.nex")
#run simmap
numberSimulations=15
'%nin%' <- Negate('%in%')
phylo <- drop.tip(phylo,
                  phylo$tip.label[
                    which(phylo$tip.label
                          %nin%names(vectorDiet))]) 
simmapdiet1 <- make.simmap(tree=phylo, vectorDiet, model="ARD", pi="estimated", nsim=numberSimulations)#inequal and not symmetrical rate of transition from folivory to frugivory etc...


length(simmapdiet1[[1]]$tip.label)
length(phylo$tip.label)

which(simmapdiet1[[1]]$tip.label %nin% phylo$tip.label)
which(phylo$tip.label %nin% simmapdiet1[[1]]$tip.label)

##---------------------
# Run (non)competitive models

#Reload tree to have same than used for biogeobears
#Fit evolutionary history
source("~/PhD/Meta_analysis/Cognition_metaanalysis/Functions.R")
resultFrugivory <- runComparisonModelsCompetition(
  simmap=simmapdiet1,
  data=data.df[!is.na(data.df$Trait)&!is.na(data.df$geographicCode)&data.df$SpeciesForPhylogeny%in%phylo$tip.label,],
  subgroup="Fruit",
  numberMaps=numberSimulations,
  trait="Trait",
  tree=phylo,
  ana_events_tables=BSM_output$RES_ana_events_tables,
  clado_events_tables=BSM_output$RES_clado_events_tables
)

#End
##-------------