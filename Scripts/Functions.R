####---------------------------------
#### Function script
####---------------------------------


# This script resumes the necessary functions for running the meta-analyses to investigate whether between species competition (direct or indirect), exploitation or cooperation better describe evolutionary scenarios of brain size evolution.
# It includes:

##---
#runBioGeoBearsandBSM_DEC

# This function allows to run the BioGeoBears algorithm (version DEC) so as to compute the history of primate biogeography

##---
#runComparisonModelsCompetition

# Runs the different evolutionary models excluding for competition (BM, OU, EB), or including competition (MC, DDlin, DDexp)

##---
#cleanTree

# Allows cleaning tree (looking for abnormalities with branch length, etc...)


##------------------------------------------------------------------------------------------------------------------------------------------------------------------

runBioGeoBearsandBSM_DEC <- 
  function(
    ###################
    ## Run the BioGeoBears analysis followed by stocchastic mapping as described in http://phylo.wikidot.com/biogeobears#script. It uses the DEC algotrithm, non-stratified. This can be freely adapted to DEC+J etc...
    ###################
    
    pathToTree,
    pathToGeo,
    maxCooccurrenceGeo,
    mapNumber,
    coreNumber,
    pathToSaveRdata,
    pathToSaveFigureGeo,
    pathToBSMMapList,
    pathToBSMoutput,
    pathToFileStockingPdfForBSMAnalysis
  ){

#######################################################
# Installing BioGeoBEARS
#######################################################

# Install optimx
#install.packages("optimx", dependencies=TRUE, repos="http://cran.rstudio.com")
library(optimx)

# Also get snow (for parallel processing)
#install.packages("snow")
library(snow)

# Install phylobase
#install.packages("phylobase", dependencies=TRUE, repos="http://cran.rstudio.com")

library(devtools)
#Intasll BioGeoBEARS:
#devtools::install_github(repo="nmatzke/BioGeoBEARS")

#Load necessary libraries
library(rexpokit)
library(cladoRcpp)
library(BioGeoBEARS)
# Load the package (after installation, see above).
library(GenSA)    # GenSA is better than optimx (although somewhat slower)
library(FD)       # for FD::maxent() (make sure this is up-to-date)
library(snow)     # (if you want to use multicore functionality; some systems/R versions prefer library(parallel), try either)
library(parallel)

#######################################################
# Phylogeny file
# Notes: 
# 1. Must be binary/bifurcating: no polytomies
# 2. No negative branchlengths (e.g. BEAST MCC consensus trees sometimes have negative branchlengths)
# 3. Be careful of very short branches, as BioGeoBEARS will interpret ultrashort branches as direct ancestors
# 4. You can use non-ultrametric trees, but BioGeoBEARS will interpret any tips significantly below the 
#    top of the tree as fossils!  This is only a good idea if you actually do have fossils in your tree,
#    as in e.g. Wood, Matzke et al. (2013), Systematic Biology.
# 5. The default settings of BioGeoBEARS make sense for trees where the branchlengths are in units of 
#    millions of years, and the tree is 1-1000 units tall. If you have a tree with a total height of
#    e.g. 0.00001, you will need to adjust e.g. the max values of d and e, or (simpler) multiply all
#    your branchlengths to get them into reasonable units.
# 6. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#######################################################

tree <- read.tree(pathToTree)

print("TREE FEATURES ARE AS FOLLOWS:")
print("Plot")
plot(tree)
print("Max node depth edge length")
print(max(node.depth.edgelength(tree)))
print("Min edge length")
print(min(tree$edge.length))
print("Is binary?")
print(is.binary(tree))
print("Is ultrametric?")
print(is.ultrametric(tree))
print("Is rooted?")
print(is.rooted(tree))

###
# TREE PATH
trfn = np(pathToTree)
moref(trfn)
###

#######################################################
# Geography file
# Notes:
# 1. This is a PHYLIP-formatted file. This means that in the 
#    first line, 
#    - the 1st number equals the number of rows (species)
#    - the 2nd number equals the number of columns (number of areas)
#    - after a tab, put the areas in parentheses, with spaces: (A B C D)
#
# 1.5. Example first line:
#    10    4    (A B C D)
# 
# 2. The second line, and subsequent lines:
#    speciesA    0110
#    speciesB    0111
#    speciesC    0001
#         ...
# 
# 2.5a. This means a TAB between the species name and the area 0/1s
# 2.5b. This also means NO SPACE AND NO TAB between the area 0/1s.
# 
# 3. See example files at:
#    http://phylo.wikidot.com/biogeobears#files
# 
# 4. Make you understand what a PLAIN-TEXT EDITOR is:
#    http://phylo.wikidot.com/biogeobears#texteditors
#
# 3. The PHYLIP format is the same format used for C++ LAGRANGE geography files.
#
# 4. All names in the geography file must match names in the phylogeny file.
#
# 5. DON'T USE SPACES IN SPECIES NAMES, USE E.G. "_"
#
# 6. Operational taxonomic units (OTUs) should ideally be phylogenetic lineages, 
#    i.e. genetically isolated populations.  These may or may not be identical 
#    with species.  You would NOT want to just use specimens, as each specimen 
#    automatically can only live in 1 area, which will typically favor DEC+J 
#    models.  This is fine if the species/lineages really do live in single areas,
#    but you wouldn't want to assume this without thinking about it at least. 
#    In summary, you should collapse multiple specimens into species/lineages if 
#    data indicates they are the same genetic population.
######################################################


#################
# PARAMETERIZATION
#################

# This is the example geography file for Hawaiian Psychotria
# (from Ree & Smith 2008)
geogfn = np(pathToGeo)

# Look at the raw geography text file:
moref(geogfn)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn=pathToGeo)
tipranges

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = maxCooccurrenceGeo

#######################################################
#######################################################
# DEC AND DEC+J ANALYSIS
#######################################################
#######################################################
# NOTE: The BioGeoBEARS "DEC" model is identical with 
# the Lagrange DEC model, and should return identical
# ML estimates of parameters, and the same 
# log-likelihoods, for the same datasets.
#
# Ancestral state probabilities at nodes will be slightly 
# different, since BioGeoBEARS is reporting the 
# ancestral state probabilities under the global ML
# model, and Lagrange is reporting ancestral state
# probabilities after re-optimizing the likelihood
# after fixing the state at each node. These will 
# be similar, but not identical. See Matzke (2014),
# Systematic Biology, for discussion.
#
# Also see Matzke (2014) for presentation of the 
# DEC+J model.
#######################################################
#######################################################

#######################################################
#######################################################

#######################################################
# Run DEC
#######################################################

# Intitialize a default model (DEC model)
BioGeoBEARS_run_object = define_BioGeoBEARS_run()

# Give BioGeoBEARS the location of the phylogeny Newick file
BioGeoBEARS_run_object$trfn = trfn

# Give BioGeoBEARS the location of the geography text file
BioGeoBEARS_run_object$geogfn = geogfn

# Input the maximum range size
BioGeoBEARS_run_object$max_range_size = max_range_size

BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
# 1. Here, un-comment ONLY the files you want to use.
# 2. Also un-comment "BioGeoBEARS_run_object = section_the_tree(...", below.
# 3. For example files see (a) extdata_dir, 
#  or (b) http://phylo.wikidot.com/biogeobears#files
#  and BioGeoBEARS Google Group posts for further hints)
#
# Uncomment files you wish to use in time-stratified analyses:
#BioGeoBEARS_run_object$timesfn = "timeperiods.txt"
#BioGeoBEARS_run_object$dispersal_multipliers_fn = "manual_dispersal_multipliers.txt"
#BioGeoBEARS_run_object$areas_allowed_fn = "areas_allowed.txt"
#BioGeoBEARS_run_object$areas_adjacency_fn = "areas_adjacency.txt"
#BioGeoBEARS_run_object$distsfn = "distances_matrix.txt"
# See notes on the distances model on PhyloWiki's BioGeoBEARS updates page.

# Speed options and multicore processing if desired
BioGeoBEARS_run_object$on_NaN_error = -1e50    # returns very low lnL if parameters produce NaN error (underflow check)
BioGeoBEARS_run_object$speedup = TRUE          # shorcuts to speed ML search; use FALSE if worried (e.g. >3 params)
BioGeoBEARS_run_object$use_optimx = "GenSA"    # if FALSE, use optim() instead of optimx()
BioGeoBEARS_run_object$num_cores_to_use = coreNumber
# (use more cores to speed it up; this requires
# library(parallel) and/or library(snow). The package "parallel" 
# is now default on Macs in R 3.0+, but apparently still 
# has to be typed on some Windows machines. Note: apparently 
# parallel works on Mac command-line R, but not R.app.
# BioGeoBEARS checks for this and resets to 1
# core with R.app)

# Sparse matrix exponentiation is an option for huge numbers of ranges/states (600+)
# I have experimented with sparse matrix exponentiation in EXPOKIT/rexpokit,
# but the results are imprecise and so I haven't explored it further.
# In a Bayesian analysis, it might work OK, but the ML point estimates are
# not identical.
# Also, I have not implemented all functions to work with force_sparse=TRUE.
# Volunteers are welcome to work on it!!
BioGeoBEARS_run_object$force_sparse = FALSE    # force_sparse=TRUE causes pathology & isn't much faster at this scale

# This function loads the dispersal multiplier matrix etc. from the text files into the model object. Required for these to work!
# (It also runs some checks on these inputs for certain errors.)
BioGeoBEARS_run_object = readfiles_BioGeoBEARS_run(BioGeoBEARS_run_object)

# Divide the tree up by timeperiods/strata (uncomment this for stratified analysis)
#BioGeoBEARS_run_object = section_the_tree(inputs=BioGeoBEARS_run_object, make_master_table=TRUE, plot_pieces=FALSE)
# The stratified tree is described in this table:
#BioGeoBEARS_run_object$master_table

# Good default settings to get ancestral states
BioGeoBEARS_run_object$return_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_TTL_loglike_from_condlikes_table = TRUE
BioGeoBEARS_run_object$calc_ancprobs = TRUE    # get ancestral states from optim run

# Set up DEC model
# (nothing to do; defaults)

# Look at the BioGeoBEARS_run_object; it's just a list of settings etc.
BioGeoBEARS_run_object

# This contains the model object
BioGeoBEARS_run_object$BioGeoBEARS_model_object

# This table contains the parameters of the model 
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table

# Run this to check inputs. Read the error messages if you get them!
check_BioGeoBEARS_run(BioGeoBEARS_run_object)

# For a slow analysis, run once, then set runslow=FALSE to just 
# load the saved result.
runslow = TRUE
resfn = pathToSaveRdata
if (runslow)
{
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  resDEC = res
} else {
  # Loads to "res"
  load(resfn)
  resDEC = res
}

#############
#Run stochastic mapping
############

model_name = "DEC"
res = resDEC
tr = tree

pdffn = pathToSaveFigureGeo##paste0("Psychotria_", model_name, "_v1.pdf")
pdf(pdffn, width=6, height=6)

analysis_titletxt = "Geographic area"

# Setup
results_object = res
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.5, statecex=0.5, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)
# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
# cmdstr = paste("open ", pdffn, sep="")
# system(cmdstr) # Plot it


#######################################################
# Stochastic mapping on DEC M3b stratified with islands coming up
#######################################################
clado_events_tables = NULL
ana_events_tables = NULL
lnum = 0

#######################################################
# Get the inputs for Biogeographical Stochastic Mapping
# Note: this can be slow for large state spaces and trees, since 
# the independent likelihoods for each branch are being pre-calculated
# E.g., for 10 areas, this requires calculation of a 1024x1024 matrix
# for each branch.  On a tree with ~800 tips and thus ~1600 branches, this was about 1.6 gigs
# for storage of "BSM_inputs_file.Rdata".
# Update: 2015-09-23 -- now, if you used multicore functionality for the ML analysis,
# the same settings will be used for get_inputs_for_stochastic_mapping().
#######################################################
BSM_inputs_fn = pathToBSMMapList
runInputsSlow = TRUE
if (runInputsSlow)
{
  stochastic_mapping_inputs_list = get_inputs_for_stochastic_mapping(res=res)
  save(stochastic_mapping_inputs_list, file=BSM_inputs_fn)
} else {
  # Loads to "stochastic_mapping_inputs_list"
  load(BSM_inputs_fn)
} # END if (runInputsSlow)

# Check inputs (doesn't work the same on unconstr)
names(stochastic_mapping_inputs_list)
stochastic_mapping_inputs_list$phy2
stochastic_mapping_inputs_list$COO_weights_columnar
stochastic_mapping_inputs_list$unconstr
set.seed(seed=as.numeric(Sys.time()))

runBSMslow = TRUE
if (runBSMslow == TRUE)
{
  # Saves to: RES_clado_events_tables.Rdata
  # Saves to: RES_ana_events_tables.Rdata
  BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=2*mapNumber, nummaps_goal=mapNumber, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
  
  RES_clado_events_tables = BSM_output$RES_clado_events_tables
  RES_ana_events_tables = BSM_output$RES_ana_events_tables
} else {#won't be used here. left in case
  # Load previously saved...
  
  # Loads to: RES_clado_events_tables
  load(file="RES_clado_events_tables.Rdata")
  # Loads to: RES_ana_events_tables
  load(file="RES_ana_events_tables.Rdata")
  BSM_output = NULL
  BSM_output$RES_clado_events_tables = RES_clado_events_tables
  BSM_output$RES_ana_events_tables = RES_ana_events_tables
} # END if (runBSMslow == TRUE)

#Save BSM output
BSM_output_fn = pathToBSMoutput
save(BSM_output, file=BSM_output_fn)

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables
head(clado_events_tables[[1]])
head(ana_events_tables[[1]])
length(clado_events_tables)
length(ana_events_tables)

include_null_range = TRUE
areanames = names(tipranges@df)
areas = areanames
max_range_size = maxCooccurrenceGeo

# Note: If you did something to change the states_list from the default given the number of areas, you would
# have to manually make that change here as well! (e.g., areas_allowed matrix, or manual reduction of the states_list)
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)

colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)

############################################
# Setup for painting a single stochastic map
############################################
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = FALSE
clado_events_table = clado_events_tables[[1]]
ana_events_table = ana_events_tables[[1]]

# cols_to_get = names(clado_events_table[,-ncol(clado_events_table)])
# colnums = match(cols_to_get, names(ana_events_table))
# ana_events_table_cols_to_add = ana_events_table[,colnums]
# anagenetic_events_txt_below_node = rep("none", nrow(ana_events_table_cols_to_add))
# ana_events_table_cols_to_add = cbind(ana_events_table_cols_to_add, anagenetic_events_txt_below_node)
# rows_to_get_TF = ana_events_table_cols_to_add$node <= length(tr$tip.label)
# master_table_cladogenetic_events = rbind(ana_events_table_cols_to_add[rows_to_get_TF,], clado_events_table)

############################################
# Open a PDF
############################################
pdffn = paste0(pathToFileStockingPdfForBSMAnalysis, model_name, "_single_stochastic_map_n1.pdf")
pdf(file=pdffn, width=6, height=6)

# Convert the BSM into a modified res object
master_table_cladogenetic_events = clado_events_tables[[1]]
resmod = stochastic_map_states_into_res(res=res, master_table_cladogenetic_events=master_table_cladogenetic_events, stratified=stratified)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=FALSE, show.tip.label=TRUE)

# Paint on the branch states
paint_stochastic_map_branches(res=resmod, master_table_cladogenetic_events=master_table_cladogenetic_events, colors_list_for_states=colors_list_for_states, lwd=5, lty=par("lty"), root.edge=TRUE, stratified=stratified)

plot_BioGeoBEARS_results(results_object=resmod, analysis_titletxt="Stochastic map", addl_params=list("j"), plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, skiptree=TRUE, show.tip.label=TRUE)

############################################
# Close PDF
############################################
dev.off()
# cmdstr = paste("open ", pdffn, sep="")
# system(cmdstr)

#######################################################
# Plot all mapNumber stochastic maps to PDF
#######################################################
# Setup
include_null_range = include_null_range
areanames = areanames
areas = areanames
max_range_size = max_range_size
states_list_0based = rcpp_areas_list_to_states_list(areas=areas, maxareas=max_range_size, include_null_range=include_null_range)
colors_list_for_states = get_colors_for_states_list_0based(areanames=areanames, states_list_0based=states_list_0based, max_range_size=max_range_size, plot_null_range=TRUE)
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))
stratified = stratified

# Loop through the maps and plot to PDF
pdffn = paste0(pathToFileStockingPdfForBSMAnalysis, model_name, "_", length(clado_events_tables), "BSMs_v1.pdf")
pdf(file=pdffn, width=6, height=6)

nummaps_goal = mapNumber
for (i in 1:nummaps_goal){
  clado_events_table = clado_events_tables[[i]]
  analysis_titletxt = paste0(model_name, " - Stochastic Map #", i, "/", nummaps_goal)
  plot_BSM(results_object=res, clado_events_table=clado_events_table, stratified=stratified, analysis_titletxt=analysis_titletxt, addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, include_null_range=include_null_range)
} # END for (i in 1:nummaps_goal)

dev.off()
# cmdstr = paste("open ", pdffn, sep="")
# system(cmdstr)

#######################################################
# Summarize stochastic map tables
#######################################################
length(clado_events_tables)
length(ana_events_tables)

head(clado_events_tables[[1]][,-20])
tail(clado_events_tables[[1]][,-20])

head(ana_events_tables[[1]])
tail(ana_events_tables[[1]])

areanames = names(tipranges@df)
actual_names = areanames
actual_names

# Get the dmat and times (if any)
dmat_times = get_dmat_times_from_res(res=res, numstates=NULL)
dmat_times

# Extract BSM output
clado_events_tables = BSM_output$RES_clado_events_tables
ana_events_tables = BSM_output$RES_ana_events_tables

# Simulate the source areas
BSMs_w_sourceAreas = simulate_source_areas_ana_clado(res, clado_events_tables, ana_events_tables, areanames)
clado_events_tables = BSMs_w_sourceAreas$clado_events_tables
ana_events_tables = BSMs_w_sourceAreas$ana_events_tables

# Count all anagenetic and cladogenetic events
counts_list = count_ana_clado_events(clado_events_tables, ana_events_tables, areanames, actual_names)

summary_counts_BSMs = counts_list$summary_counts_BSMs
print(conditional_format_table(summary_counts_BSMs))

# Histogram of event counts
hist_event_counts(counts_list, pdffn=paste0(pathToFileStockingPdfForBSMAnalysis, model_name, "_histograms_of_event_counts.pdf"))

#######################################################
# Print counts to files
#######################################################
tmpnames = names(counts_list)
cat("\n\nWriting tables* of counts to tab-delimited text files:\n(* = Tables have dimension=2 (rows and columns). Cubes (dimension 3) and lists (dimension 1) will not be printed to text files.) \n\n")
for (i in 1:length(tmpnames))
{
  cmdtxt = paste0("item = counts_list$", tmpnames[i])
  eval(parse(text=cmdtxt))
  
  # Skip cubes
  if (length(dim(item)) != 2)
  {
    next()
  }
  
  outfn = paste0(tmpnames[i], ".txt")
  if (length(item) == 0)
  {
    cat(outfn, " -- NOT written, *NO* events recorded of this type", sep="")
    cat("\n")
  } else {
    cat(outfn)
    cat("\n")
    write.table(conditional_format_table(item), file=outfn, quote=FALSE, sep="\t", col.names=TRUE, row.names=TRUE)
  } # END if (length(item) == 0)
} # END for (i in 1:length(tmpnames))
cat("...done.\n")

#######################################################
# Check that ML ancestral state/range probabilities and
# the mean of the BSMs approximately line up
#######################################################
library(MultinomialCI)    # For 95% CIs on BSM counts
check_ML_vs_BSM(res, clado_events_tables, model_name, tr=NULL, plot_each_node=FALSE, linreg_plot=TRUE, MultinomialCI=TRUE)

}

##---------------------------------------------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------------------------------------------


runComparisonModelsCompetition <- function(
  
  #######################
  ## Run the phylogenetic comparison between random evolution (BM, OU) or competitive scenario with MC= exclusion from competitive
  #Taxa = mutual exclusion, or DDlin and DDexp, which are positive or negative density dependance to taxa of same "group" (here feeding group)
  #of the evolutionary rate, either linearily (lin), or exponentially (exp)
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
  )
  { 
    # 
    # simmap=simmapdiet1
    # data=summaryData[!is.na(summaryData$Brain.log)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,]
    # subgroup="Fruit"
    # numberMaps=numberSimulations
    # trait="Brain.log"
    # tree=phylo
    # ana_events_tables=BSM_output$RES_ana_events_tables
    # clado_events_tables=BSM_output$RES_clado_events_tables

    #source("~/PhD/Meta_analysis/Cognition_metaanalysis/function_RPANDA.R")
    #extracted from: https://journals.plos.org/plosbiology/article?rev=1&id=10.1371/journal.pbio.2003563#pbio.2003563.ref028
    
    #load required R packages
    require(RPANDA)
    require(geiger)
    
    library(svMisc)
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
    
    for(i in 1:N){
      #print(i)
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
      {o2<-fitContinuous(subtree,M,model="BM", ncores=1)
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
      print("BM")
      tryCatch(
      {
      o3<-fitContinuous(subtree,M,model="OU", ncores=1)
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
      print("OU")
      
      tryCatch(
        {
          o32<-fitContinuous(subtree,M,model="EB", ncores=1)
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
      print("EB")
      
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
      print("MC")
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
      print("DDexp")
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
      print("DDlin")
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
      progress(i/N*100)
    }
    
    res.mat <- as.data.frame(res.mat)
    resm <- data.frame(res.mat)
   
    return(resm)

}

##---------------------------------------------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------------------------------------------
##---------------------------------------------------------------------------------------------------------------------------------------------------------

cleanTree <- function(tr, trfn){
  
#From: http://phylo.wikidot.com/testing-for-fixing-common-tree-file-problems-in-biogeobears
  
  # Is the tree dichotomous?  And rooted?
  if (is.binary.tree(tr) == FALSE)
  {
    stop("Stopping, because tree is not binary. (In other words, you have polytomies. BioGeoBEARS cannot handle polytomies.)")
  }
  if (is.rooted(tr) == FALSE)
  {
    stop("Stopping, because tree is not rooted. (BioGeoBEARS requires rooted trees.)")
  }
  
  #######################################################
  # Check for negative and 0-length branches
  # Note: in BioGeoBEARS, tips with branches very close to length 0 (but not exactly 0)
  #       are treated as direct ancestors.  But, negative branchlengths, and internal
  #       branchlengths=0, all have to be fixed.
  #######################################################
  # Minimum branchlength for any new internal branches
  min_brlen = 0.01
  
  # Are there 0-length branches?
  min(tr$edge.length)
  
  # Remove branchlengths less than 0
  if (min(tr$edge.length) < 0)
  {
    cat("\n\n")
    cat("Negative branchlengths detected at these nodes:\n\n")
    trtable = prt(tr, printflag=FALSE)
    negative_brlens_TF = trtable$edge.length < 0
    negative_brlens_TF[is.na(negative_brlens_TF)] = FALSE
    negative_BL_rows = trtable[negative_brlens_TF, ]
    print(negative_BL_rows)
    cat("\n")
    
    cat("Correcting, using impose_min_brlen()")
    tr = impose_min_brlen(phy=tr, min_brlen=min_brlen, leave_BL0_terminals=TRUE)
    
    if (grepl(pattern="\\.newick", x=trfn) == TRUE)
    {
      new_trfn = gsub(pattern="\\.newick", replacement="_noNeg.newick", x=trfn)
    } else {
      new_trfn = paste0(trfn, "_noNeg.newick")
    }
    trfn = new_trfn
    
    write.tree(tr, file=new_trfn)
    cat("\n\nFixed negative branchlengths, and saved to ", new_trfn, "\n", sep="")
  } # END if (min(tr$edge.length) < 0)
  min(tr$edge.length)
  sum(tr$edge.length == min(tr$edge.length))
  
  #######################################################
  # *Internal* branches of 0 length?
  #######################################################
  trtable = prt(tr, printflag=FALSE)
  internal_TF = trtable$node.type == "internal"
  edges_BL0_TF = trtable$edge.length == 0
  sum_TFs = (internal_TF + edges_BL0_TF)
  sum_TFs[is.na(sum_TFs)] = 0
  internal_BL0_TF = (sum_TFs == 2)
  sum(internal_BL0_TF)
  
  # Edit the branchlengths, if needed
  if (sum(internal_BL0_TF) > 0)
  {
    internal_BL0_rows = trtable[internal_BL0_TF, ]
    internal_BL0_rows
    nodes_to_change = internal_BL0_rows$node
    edges_to_change = internal_BL0_rows$parent_br
    edges_to_change
    
    cat("\n\n")
    cat("Internal branches of length 0 detected at these nodes:\n\n")
    print(internal_BL0_rows)
    cat("\n")
    cat("Changing to length min_brlen=", min_brlen, "...")
    tr$edge.length[edges_to_change] = min_brlen
    
    # New filename
    if (grepl(pattern="\\.newick", x=trfn) == TRUE)
    {
      new_trfn = gsub(pattern="\\.newick", replacement="_minBL.newick", x=trfn)
    } else {
      new_trfn = paste0(trfn, "_minBL.newick")
    }
    trfn = new_trfn
    
    write.tree(tr, file=new_trfn)
    cat("...saved to ", new_trfn, "\n", sep="")
  } # END if (sum(internal_BL0_TF) > 0)
  
  ##################################################
  # *External* branches of 0 length?
  ##################################################
  external_TF = trtable$node.type == "tip"
  edges_BL0_TF = trtable$edge.length == 0
  sum_TFs = (external_TF + edges_BL0_TF)
  sum_TFs[is.na(sum_TFs)] = 0
  external_BL0_TF = (sum_TFs == 2)
  sum(external_BL0_TF)
  
  tip_brlen = 0.0000001
  if (sum(external_BL0_TF) > 0)
  {
    cat("\n\n")
    cat("Tip branches of length zero detected. Changing to length ", tip_brlen, sep="")
    cat("\n\n")
    edgenums_to_change = trtable$parent_br[external_BL0_TF]
    edgenums_to_change
    
    # Change them
    tr$edge.length[edgenums_to_change] = tip_brlen
    
    # New filename
    if (grepl(pattern="\\.newick", x=trfn) == TRUE)
    {
      new_trfn = gsub(pattern="\\.newick", replacement="_tipsNo0.newick", x=trfn)
    } else {
      new_trfn = paste0(trfn, "_tipsNo0.newick")
    }
    trfn = new_trfn
    
    write.tree(tr, file=new_trfn)
    cat("...saved to ", new_trfn, "\n", sep="")
  } # END if (sum(external_BL0_TF) > 0)
  
return(tr)
}

















