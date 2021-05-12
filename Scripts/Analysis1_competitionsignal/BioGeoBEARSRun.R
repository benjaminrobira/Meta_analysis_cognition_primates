####---------------------------------
#### BioGeoBearsRun
####---------------------------------

# This script allows to fit the biogeobears algorithm to primate ranging data given different threshold to consider that a primate species initially belong to a biogeography area.

###Set working directory
setwd("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates")

#Import environment
rm(list=ls())
load("REnvironments/geography_traits.RData")

##--------
#Home made functions
#To source all phylogenetics functions (biogeobears + models of evolution)
source("Scripts/Functions.R")

#My toolkit
source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
##--------

##--------
#Libraries

#Phylogeny
library(caper)
library(MCMCglmm)
library(RPANDA)
library(BioGeoBEARS)
library(phytools)
library(ape)
library(geiger)
library(optimix)

##--------

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

##----
##Run the biogeobears for speeding up after

##--
# THE TREE

#Phylogenetic tree: force it to be ultrametric, problem was minor because doesn't change branch length actually
# phylo_all <-read.nexus("Raw_data/Tree/TreeBlock_10kTrees_Primates_Version3.nex")
# phylo_init <- phylo_all[[e]]

phylo_all <-read.nexus("Raw_data/Tree/consensusTree_10kTrees_Primates_Version3.nex")
phylo_init <- phylo_all

phylo <- force.ultrametric(tree=phylo_init, method="extend")#method="nnls")
is.ultrametric(phylo)

#check if this differ (was done outside the loop and was ok)
plot(phylo_init$edge.length, phylo$edge.length, xlab="Original tree", ylab="Post-chronos tree", pch=20, bty="n")
abline(a=0, b=1, col="gray", lwd=0.5)

phylo <- multi2di(phylo)
is.binary(phylo)

phylo <- drop.tip(phylo,
                  phylo$tip.label[
                    which(phylo$tip.label
                          %nin%summaryData$SpeciesForPhylogeny[!is.na(summaryData$geographicCode)])]) 


write.tree(phylo, np("Raw_data/Tree/Tree_biogeobears.nex"))
phylo <- read.tree("Raw_data/Tree/Tree_biogeobears.nex")

phylo <- cleanTree(tr=phylo, trfn="Tree_biogeobears.nex")
phylo <- force.ultrametric(tree=phylo, method="extend")#method="nnls")

is.ultrametric(phylo)
is.binary(phylo)

write.tree(phylo, np("Raw_data/Tree/Tree_biogeobears.nex"))
phylo <- read.tree("Raw_data/Tree/Tree_biogeobears.nex")

#numstates_from_numareas(numareas=12, maxareas=3, include_null_range=TRUE)#check how big is the transition matrix to deal with

frugivoryThresholdVector <- seq(from=20, to=40, by=20)
folivoryThresholdVector <- seq(from=40, to=60, by=20)
geographicThresholdVector <- c(10,30)/100
randomSampling=10
numberSimulations=10
numberTrees=1
for(c in 1:length(geographicThresholdVector)){
  
  ###------------------
  ##Adding the co-occurence
  ###------------------
  
  summaryData$geographicCode <- matrixRangingSensitivity[match(summaryData$SpeciesForPhylogeny,matrixRangingSensitivity$SpeciesForPhylogeny),which(thresholdPresenceRange==geographicThresholdVector[c])]   
  
  ##--------
  # Evolutionary history of range
  ##--------
  
  ##----------
  ##Create geography table for phylogenetic fit
  ##----------
  
  geographyDataTableTest <- c(
    paste(
      length(summaryData$geographicCode[!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label]),
      length(areaName), 
      paste("(", paste(LETTERS[1:length(areaName)], collapse=" "), ")", sep=""), sep="\t"),
    paste(summaryData$SpeciesForPhylogeny[!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label], 
          summaryData$geographicCode[!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label], sep="\t")
  )
  
  #Save table with a tabulation as necessary
  write.table(x=geographyDataTableTest,
              file="test_geo.data", 
              sep="\t",col.names=FALSE, row.names=FALSE, na="", quote = FALSE)
  
  
  BSM_output_from_function <- runBioGeoBearsandBSM_DEC(
    pathToTree="Raw_data/Tree/Tree_biogeobears.nex",
    pathToGeo="test_geo.data",
    maxCooccurrenceGeo=3,
    mapNumber=numberSimulations,
    coreNumber=6,
    pathToSaveRdata=paste("BioGeoBEARS/BioGeoBEARS_", c, ".Rdata", sep=""),
    pathToSaveFigureGeo=paste("Brain_DEC_v1", c, ".pdf", sep=""),
    pathToBSMMapList=paste("BioGeoBEARS/BSM_output_file", c, ".Rdata", sep=""),
    pathToBSMoutput=paste("BioGeoBEARS/BSM_output_file", c, ".Rdata", sep=""),
    pathToFileStockingPdfForBSMAnalysis="BioGeoBEARS/"
  )
  dev.off()
} 

save.image("BioGeoBEARS/geography_traits_biogeobears.RData", version=2)
