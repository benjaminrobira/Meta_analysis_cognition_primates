# Meta analysis of brain size in primates :monkey:

This repository stores all material linked to "**Species sympatry shapes brain size evolution in Primates**" by *Benjamin Robira* and *Benoît Perez-Lamarque*.

It is divided into several folders:

## [:file_folder: **Raw data**](Raw_data)

It includes a) brain and associated data (each within a folder named after their study of origin) b) phylogenetic trees (from the 10kTrees project/consensus tree, larger than the 10kTrees one, used for diversification analysis; "Tree" file). IUCN data should be directly requested to the IUCN.

## [:file_folder: **Scripts**](Scripts)

It includes: 

  - The script with all necessary functions to make the EvolutionaryHistoryClusterRun_X work (`Function`).
* [:file_folder: **Analysis1_competitionsignal:**](Scripts/Analysis1_competitionsignal) Analysis to see whether brain size evolution in frugivorous primates is better predicted by models accounting or not for sympatry.
  - The script for cleaning, merging, and processing ranging data and attribute/trait data (`DataCleaning`).
  - The script with prior run of the biogeobear algorithm to retrace biogeography history (`BioGeoBEARSRun`).
  - The script(s) for running the analysis on the cluster (e.g. `EvolutionaryHistoryClusterRun_a1b1`). 
  - The script for extracting and plotting results from cluster run (`EvolutionaryHistoryAnalysis`). 
  - The script for checking the evolutionary dependence of brain areas (`Checking_interdependence_evolution_areas`).
* [:file_folder: **Analysis2_competitiondirection:**](Scripts/Analysis2_competitiondirection) Analysis based on PGLS to see the direction of the effect of sympatry on brain size
  - The script for extracting the necessary spatial covariates: number of sympatric species, and average overlap (`SpatialProcessingSympatry`). 
  - The script for conducting the PGLS analysis on selection direction (`PGLSdirectionSelection`).

* [:file_folder: **Analysis3_diversification:**](Scripts/Analysis3_diversification) Analysis to see whether how the diversification rate covaries with brain size or sympatry
  - The script for running the diversification algorithm ClaDS on the cluster (`script_diversification_primates`). 
  - The script for conducting the PGLS analysis on diversification and brain size (`PGLSdiversification_bayesian` - non-Bayesian version also available, but suffered from bias).
  - The script for conducting the PGLS analysis on diversification and sympatry (`PGLSdiversificationAndSympatry`).
  - The file with all inputs/outputs from the diversification analysis.

:question: Please if you encounter any issue or see any error in the code, have questions or would like to discuss about this topic, feel free to email me at [:e-mail:](mailto:benjamin.robira@normalesup.org) benjamin.robira@normalesup.org  

My apology for English errors I have made,  
Benjamin.