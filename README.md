# Meta_analysis_cognition_primates
This repository stores all material linked to "Primate Cognition and the Red Queen Intelligence Hypothesis". It is divided into several folders:

1) Raw_data:  
It includes a) brain and associated date (each within a folder named after their study of origin) b) phylogenetic trees (from the 10kTrees project, "Tree" file).

2) Scripts:   
It includes: 
* The script for the analysis of the evolutionary history of brain in primates (data cleaning and merging; "DataCleaningMerging")
* The script(s) for running the analysis on the cluster ("RunClusterX") 
* The script with all necessary functions to make RunClusterX work ("Function").
* The scripts for extracting and plotting results from cluster run ("EvolutionaryHistoryResults") 
* The script for conducting the PGLS analysis to deepen results from contrasting competitive vs non-competitive scenarios ("PGLSdirectionSelection") 
* The script for running the diversification analysis ("DiversificationAnalysis")

3) Text:  
It includes the R markdown document that is actually the written version of the article ("Article.Rmd") and the associated bibtex file for article reference only (intern compilation for adding *R* packages; "biliograpphyarticle.bib"). Supplementary files (the .lua and the .csl) are added so as to allow correct compilation. They must be stored together with the .Rmd file.
