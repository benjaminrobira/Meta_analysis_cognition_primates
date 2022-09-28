##########
## Testing for interdepence of evolution of the different areas
##########

rm(list=ls())
load("REnvironments/PGLSdirectionSelection.RData")

library(mvMORPH)
library(tidyr)
library(ape)
library(phylolm)
library(phytools)
library(caper)
library(geiger)
library(BioGeoBEARS)

# Preparing the data ------------------------------------------------------

#Select the continuous trait of brain sizes
summaryData_formvBM <- summaryData %>% 
  dplyr::select(SpeciesForPhylogeny, EQ.log, ratioNeocortex.log, ratioHippocampus.log, ratioCerebellum.log, ratioStriatum.log) %>% #, ratioMOB.log) %>% 
  tidyr::drop_na()

#Subset tree to get only species included in dataset
phyloConsensus_rdc <- phylo

phyloConsensus_rdc <- drop.tip(phyloConsensus_rdc,
                               phyloConsensus_rdc$tip.label[
                    which(phyloConsensus_rdc$tip.label
                          %nin%summaryData_formvBM$SpeciesForPhylogeny)]) 

write.tree(phyloConsensus_rdc, np("Raw_data/Tree/Tree_mvBM.nex"))
phyloConsensus_rdc <- read.tree("Raw_data/Tree/Tree_mvBM.nex")

phyloConsensus_rdc <- cleanTree(tr=phyloConsensus_rdc, trfn="Tree_mvBM.nex")
phyloConsensus_rdc <- force.ultrametric(tree=phyloConsensus_rdc, method="extend")#method="nnls")

is.ultrametric(phyloConsensus_rdc)
is.binary(phyloConsensus_rdc)

write.tree(phyloConsensus_rdc, np("Raw_data/Tree/Tree_mvBM.nex"))
phyloConsensus_rdc <- read.tree("Raw_data/Tree/Tree_mvBM.nex")

# Fitting the multivariate brownian motion --------------------------------

rownames(summaryData_formvBM) <- summaryData_formvBM$SpeciesForPhylogeny
#Delete species not in tree
summaryData_formvBM <- summaryData_formvBM[summaryData_formvBM$SpeciesForPhylogeny %in% phyloConsensus_rdc$tip.label,]
#Reorder df following tree
summaryData_formvBM <- summaryData_formvBM %>%
  arrange(factor(SpeciesForPhylogeny, levels = phyloConsensus_rdc$tip.label))

#Sample size
nrow(summaryData_formvBM)

#Fit multivariate model
multivariateBrownianModel <- mvBM(phyloConsensus_rdc, summaryData_formvBM[,-c(1,2)], model="BM1", method="inverse")
multivariateBrownianModel$sigma

multivariateOUModel <- mvOU(phyloConsensus_rdc, summaryData_formvBM[,-c(1,2)], model="OU1", method="inverse")
multivariateOUModel$sigma


