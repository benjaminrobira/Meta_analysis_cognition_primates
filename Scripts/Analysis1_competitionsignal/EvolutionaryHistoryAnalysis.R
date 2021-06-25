####---------------------------------
#### EvolutionaryHistoryAnalysis
####---------------------------------

# This script allows processing the result of the estimation with different evolutionary scenario of primate brain size evolution. It is still in construction and is therefore a huge mess.

# It is divided into several parts (to be completed):

#1) Creation of the map and phylogenetic interaction tree
#2) Creation of the brain + raw data figure
#3) Creation of the model fit figure

###Set working directory
setwd("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates")

#Import environment
rm(list=ls())

##--------
#Home made functions
#To source all phylogenetics functions (biogeobears + models of evolution)
source("Functions.R")

#My toolkit
source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
##--------

##--------
#Libraries

#Plot
library(RColorBrewer)
library(tidyr)
library(stringr)
library(svMisc)
library(plotrix)
library(circlize)

#Spatial
library(rworldmap) # World map
library(cleangeo) #to clean it otherwise issues with intersection
library(maps)
library(rgeos) #for readOGR; gArea/gCentroid...
library(sf) #for intersection
library(rgdal)
library(geosphere)

#Phylogenetics
library(caper)
library(MCMCglmm)
library(RPANDA)
library(BioGeoBEARS)
library(phytools)
library(ape)
library(geiger)
library(optimix)

#Parallelizing
library(snow)
library(foreach)
library(doParallel)

##--------

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

###----------------------
## Plotting results of phylogenetic history
###----------------------

tree <-read.nexus("Raw_data/Tree/TreeBlock_10kTrees_Primates_Version3.nex")
tree <- tree[[1]]
tree <- force.ultrametric(tree=tree, method="extend")#method="nnls")
is.ultrametric(tree)
is.rooted(tree)

'%nin%' <- Negate('%in%')
drop.tip(tree,
         tree$tip.label[
           which(tree$tip.label
                 %nin%unique(summaryData$SpeciesForPhylo[summaryData$DietaryGuild=="Fruit"]))]) 

#Write and reload tree so branches are correctly numbered
write.tree(phylo, np("Raw_data/Tree/Tree_diet.nex"))
tree <- read.tree("Raw_data/Tree/Tree_diet.nex")

numberOrder=6
hc = as.hclust(tree)
labels = hc$labels

labels.tordc <- as.data.frame(labels)
colnames(labels.tordc) <- "Name"
labels.tordc <- separate(labels.tordc, col="Name", into=c("Name1", "Name2", "Name3"), sep="_")

labels.rdc <- apply(labels.tordc, 1, function(x){
  if(!is.na(x[3])){
    paste(toupper(substr(x[1], 1, 1)), ". ", substr(x[2], 1, 3), ". ", substr(x[3], 1, 1), ".", sep="")
  } else{paste(toupper(substr(x[1], 1, 1)), ". ", substr(x[2], 1, 3), sep="")
  }
}
)

#Link between species/groups:
speciesLabels <- hc$labels#Should be in the tree order
#groupLabels <- rep(seq(from=1, to=23, by=2), each=2)[-1]

########
## Fig 1: map and interaction
########

#Reimport areas with cropping
centroid <- matrix(NA, ncol=2, nrow=length(areaName))
for(i in 1:length(areaName)){
  areaTransitory <- readOGR(dsn=paste("T:/IUCN_data_primate/Geographic_areas/Shapefiles/",areaName[i],".shp",sep=""))
  areaTransitory = clgeo_Clean(areaTransitory)
  areaTransitory <- spTransform(areaTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  #Have mercator for intersection
  areaTransitory <- spTransform(areaTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  areaTransitory <- gIntersection(areaTransitory, worldMap_mercator, byid=FALSE)
  areaTransitory = clgeo_Clean(areaTransitory)
  #Reunite polygon in case
  areaTransitory <- gBuffer(areaTransitory, byid=F, width=0)
  
  #back transform to long/lat
  areaTransitory <- spTransform(areaTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  assign(paste("area", i, sep="_"), areaTransitory)
  if(i==1){
    centroid[i,] <- c(summary(areaTransitory)$bbox[1,2] + 5, summary(areaTransitory)$bbox[2,1]) 
  }
  else if (i==2){
    centroid[i,] <- c(summary(areaTransitory)$bbox[1,2] - 5, summary(areaTransitory)$bbox[2,1] - 5) 
  }
  else{
    centroid[i,] <- geosphere::centroid(areaTransitory)
  }
}
warnings()

#Just to see the difference as pinpointed by the warning
# geosphere::areaPolygon(areaTransitory)
# gArea(areaTransitory)
# really false
# geosphere::centroid(areaTransitory)
# gCentroid(areaTransitory)
#ok

layout(mat=t(c(1,2)), widths=c(5,10), heights=c(10))
par(mar=c(0, 0, 0, 0), mgp=c(2, 0.5, 0), xpd=TRUE)


#Create the map of the geographic area
#Have background
map("world", fill=TRUE, col="lightgray", bg="white", border=NA)#, ylim=c(-60, 50))

for(i in 1:length(areaName)){
  plot(get(paste("area", i, sep="_")), col=colourArea[i], border=colourArea[i], add=TRUE) #border="black",
}
#Have borders
#plot(worldMap, col=NA, border="white",bg="white", lwd=0.1, add=TRUE)
# for(i in 1:length(areaName)){
#   plot(get(paste("area", i, sep="_")), col=NA, border="black",  add=TRUE)
# }
points(x=centroid[,1], y=centroid[,2], pch=19, col=c(colourArea[1], colourArea[2], rep("white", times=10)), cex=1.15)
points(x=centroid[,1], y=centroid[,2], cex=1.15, col=c("white", "white", rep("black", times=10)))
text(x=centroid[,1], y=centroid[,2], labels=1:length(areaName), cex=0.4, col=c("white", "white", rep("black", times=10)), adj=c(0.5,0.5))

####
## Fig interaction
####

#Import data
summaryDataForPlot <- read_delim("C:/Users/robira/Documents/PhD/Meta_analysis/Cognition_metaanalysis/OutputEvolModel/Dataplot.txt","\t", escape_double = FALSE, trim_ws = TRUE)

#Tree
tree <- read.tree("C:/Users/robira/Documents/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/Tree_diet.nex")

#Create geography df
geoBinary <- as.data.frame(summaryDataForPlot[,c(3,which(colnames(summaryDataForPlot)=="geographicCode"))])
colnames(geoBinary) <- c("SpeciesPhylo", "Loc")

#Create species ID
hc = as.hclust(tree)#bird.orders)
labels = hc$labels  # name of birds

labels.tordc <- as.data.frame(labels)
colnames(labels.tordc) <- "Name"
labels.tordc <- separate(labels.tordc, col="Name", into=c("Name1", "Name2", "Name3"), sep="_")

labels.rdc <- apply(labels.tordc, 1, function(x){
  if(!is.na(x[3])){
    paste(toupper(substr(x[1], 1, 1)), ". ", substr(x[2], 1, 3), ". ", substr(x[3], 1, 1), ".", sep="")
  } else{paste(toupper(substr(x[1], 1, 1)), ". ", substr(x[2], 1, 3), sep="")
  }
}
)

#Match to have right order for geography
locationSpecies <- geoBinary$Loc[match(labels, geoBinary$SpeciesPhylo)]
colLoc <- colourArea

#Match to have diet
dietSpecies <- summaryDataForPlot$DietaryGuild[match(labels, summaryDataForPlot$SpeciesForPhylogeny)]

#Getting species labels abbreviated
speciesLabels <- hc$labels#Should be in the tree order

#Create the circos plot linking species based on their diet and geography
circos.clear()
circos.par(gap.degree=0, gap.after=0, cell.padding=c(0,0,0,0), track.margin = c(0, 0), "canvas.xlim" = c(-1.1, 1.1), "canvas.ylim" = c(-1.1, 1.1))
circos.initialize(speciesLabels, xlim = c(0, 1))

#Species name + area it belongs to
circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1, track.margin=c(0.01, 0.01),
             panel.fun = function(x, y) {
               i=CELL_META$sector.numeric.index
               circos.text(CELL_META$xcenter, 1, labels.rdc[i], adj = c(0, 0), 
                           facing = "clockwise", niceFacing = TRUE,
                           col = "black", cex = 0.35, font=3)
               
               geo <- as.numcharac(unlist(strsplit(locationSpecies[i], "")))
               for(g in 1:length(geo)){
                 if(geo[g]==1){
                   #circos.points(CELL_META$xcenter,0.75/length(geo)*(2*g-1)/2, col=as.character(colLoc[g]), pch=19, cex=0.2)
                   circos.rect(CELL_META$cell.xlim[1],0.75/length(geo)*(g-1), CELL_META$cell.xlim[2], 0.75/length(geo)*(g), col=as.character(colLoc[g]), border=NA)
                 }
               }
             })

#Plot the geographic links
for(i in 1:length(speciesLabels)){
  #locI <- which(strsplit(locationSpecies[i], "")==1)
  for(j in i:length(speciesLabels)){
    #locJ <- which(strsplit(locationSpecies[j], "")==1)
    product <- as.numcharac(unlist(strsplit(locationSpecies[j], "")))*as.numcharac(unlist(strsplit(locationSpecies[i], "")))
    if(i==j|(length(unique(product))==1&unique(product)[1]==0)){
      #Do nothing
    }
    else{
      if(dietSpecies[i]=="Fruit"&dietSpecies[i]==dietSpecies[j]){
        # colour <- as.data.frame(table(colLoc[which(product==1)]))
        # if(is.finite(max(colour$Freq))){
        # }else{
        #   print(c(i,j))
        # }
        # colour <- colour[colour$Freq==max(colour$Freq),1][1]
        circos.link(speciesLabels[i], runif(1, 0, 1), speciesLabels[j], runif(1, 0, 1), lwd=0.2, col="lightgray")#adjustcolor(as.character(colour), alpha.f=0.9))
      }
      else{
        #circos.link(speciesLabels[i], runif(1, 0, 1), speciesLabels[j], runif(1, 0, 1), lwd=1, col="lightgray")
      }
    }
  }
}
circos.clear()

#Plot the phylogenetic tree in a new circular plot
n = length(labels)  # number of species
dend = as.dendrogram(hc)

par(new = TRUE) # <- magic
circos.par("canvas.xlim" = c(-1.05, 1.05), "canvas.ylim" = c(-1.25, 1.25))
circos.initialize("a", xlim = c(0, n)) # only one sector
# circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.3, 
#              panel.fun = function(x, y) {
#                for(i in seq_len(n)) {
#                  circos.text(i-0.5, 0, labels.rdc[i], adj = c(0, 0.5), 
#                              facing = "clockwise", niceFacing = TRUE,
#                              col = "black", cex = 0.2, font=3)
#                }
#              })


#suppressPackageStartupMessages(library(dendextend))
#dend = color_branches(dend, k = 6, col = 1:6)
dend_height = attr(dend, "height")
circos.track(ylim = c(0, dend_height), bg.border = NA,
             track.height = 0.95, panel.fun = function(x, y) {
               circos.dendrogram(dend)
             })
circos.clear()

########
## Fig 2: brain and measures
########

library(RColorBrewer)
colourVector <- c("darkgrey", brewer.pal(n = 5, name = "Set1")[1:5])
colourVectorbis <- c("lightgray", brewer.pal(n = 5, name = "Pastel1")[1:5])
colour.circle.points <- c("black", "darkred", "darkblue", "darkgreen", "purple4", "orange4")

###
## Fig brain values / circular
###

summaryDataForPlot$EQ <- summaryDataForPlot$Brain*1.036*(10**-3)/(0.085*summaryDataForPlot$Bodymass**0.775)
summaryDataForPlot$ratioNeocortex <- summaryDataForPlot$Neocortex/summaryDataForPlot$Brain
summaryDataForPlot$ratioHippocampus <- summaryDataForPlot$Hippocampus/summaryDataForPlot$Brain
summaryDataForPlot$ratioCerebellum <- summaryDataForPlot$Cerebellum/summaryDataForPlot$Brain
summaryDataForPlot$ratioStriatum <- summaryDataForPlot$Striatum/summaryDataForPlot$Brain
summaryDataForPlot$ratioMOB <- summaryDataForPlot$MOB/summaryDataForPlot$Brain

#Brain data
relativeValueEQ <- summaryDataForPlot$EQ[match(speciesLabels, summaryDataForPlot$SpeciesForPhylogeny)] - 1#runif(length(speciesLabels), -1, 1)
relativeValueNeocortex <- scale(summaryDataForPlot$ratioNeocortex[match(speciesLabels, summaryDataForPlot$SpeciesForPhylogeny)])
relativeValueHippocampus <- scale(summaryDataForPlot$ratioHippocampus[match(speciesLabels, summaryDataForPlot$SpeciesForPhylogeny)])
relativeValueCerebellum <- scale(summaryDataForPlot$ratioCerebellum[match(speciesLabels, summaryDataForPlot$SpeciesForPhylogeny)])
relativeValueStriatum <- scale(summaryDataForPlot$ratioStriatum[match(speciesLabels, summaryDataForPlot$SpeciesForPhylogeny)])
relativeValueMOB <- scale(summaryDataForPlot$ratioMOB[match(speciesLabels, summaryDataForPlot$SpeciesForPhylogeny)])


library(circlize)
circos.clear()
circos.par(gap.degree=0, gap.after=0, cell.padding=c(0,0,0,0), track.margin=c(0, 0))
circos.initialize(speciesLabels, xlim = c(0, 1))

#Species name
circos.track(ylim = c(0, 20), bg.border = NA, track.height = 0.05, track.margin=c(0.01, 0.1),
             panel.fun = function(x, y) {
               i=CELL_META$sector.numeric.index
               circos.text(CELL_META$xcenter, 0, labels.rdc[i], adj = c(0, 0),
                           facing = "clockwise", niceFacing = TRUE,
                           col = "black", cex = 0.35, font=3)
             })

#Background
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  circos.rect(0, 0, 1, 1, col=colourVector[1], border=colourVector[1])
}, track.height = 1/15)

circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  circos.rect(0, 0, 1, 1, col=colourVectorbis[1], border=colourVectorbis[1])
}, track.height = 1/15)

#Background
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  circos.rect(0, 0, 1, 1, col=colourVector[2], border=colourVector[2])
}, track.height = 1/15)

circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  circos.rect(0, 0, 1, 1, col=colourVectorbis[2], border=colourVectorbis[2])
}, track.height = 1/15)

#Background
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  circos.rect(0, 0, 1, 1, col=colourVector[3], border=colourVector[3])
}, track.height = 1/15)

circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  circos.rect(0, 0, 1, 1, col=colourVectorbis[3], border=colourVectorbis[3])
}, track.height = 1/15)

#Background
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  circos.rect(0, 0, 1, 1, col=colourVector[4], border=colourVector[4])
}, track.height = 1/15)

circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  circos.rect(0, 0, 1, 1, col=colourVectorbis[4], border=colourVectorbis[4])
}, track.height = 1/15)

#Background
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  circos.rect(0, 0, 1, 1, col=colourVector[5], border=colourVector[5])
}, track.height = 1/15)

circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  circos.rect(0, 0, 1, 1, col=colourVectorbis[5], border=colourVectorbis[5])
}, track.height = 1/15)

#Background
circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  circos.rect(0, 0, 1, 1, col=colourVector[6], border=colourVector[6])
}, track.height = 1/15)

circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
  circos.rect(0, 0, 1, 1, col=colourVectorbis[6], border=colourVectorbis[6])
}, track.height = 1/15)


library(plotrix)
#Main circle
for(i in 1:15){
  draw.circle(x=0,y=0,1-1/15-(i-1)*1/15, col=NA, border="white")
}

#increment of 0.5
for(i in 1:30){
  draw.circle(x=0,y=0,1-1/15/2-(i-1)*1/15/2, col=NA, border="white", lty=2)
}

#Value

#EQ
absMax <- max(abs(relativeValueEQ), na.rm=TRUE)
circos.track(ylim = c(0, 1), bg.border = NA, track.index=2, panel.fun = function(x, y) {
  i=CELL_META$sector.numeric.index
  #circos.rect(0, 0, 1, 1, col=colourPositive, border=colourPositive)
  if(is.na(relativeValueEQ[i])){}  else{
    if(relativeValueEQ[i] > 0 & dietSpecies[i]=="Fruit"){
      circos.points(CELL_META$xcenter, relativeValueEQ[i]/absMax, pch=19, col=colour.circle.points[1], cex=0.7)
      circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueEQ[i]/absMax, lty=3)
    }
    else if(relativeValueEQ[i] > 0 & dietSpecies[i]=="Leaf"){
      circos.points(CELL_META$xcenter, relativeValueEQ[i]/absMax, pch=21, col=colour.circle.points[1], bg="white", cex=0.7)
      circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueEQ[i]/absMax, lty=3)
    }
    else{}
  }
}, track.height = 0.1)


circos.track(ylim = c(0, 1), bg.border = NA, track.index=3,  panel.fun = function(x, y) {
  i=CELL_META$sector.numeric.index
  if(is.na(relativeValueEQ[i])){}  else{
    #circos.rect(0, 0, 1, 1, col=colourNegative, border=colourNegative)
    if(relativeValueEQ[i] <= 0 & dietSpecies[i]=="Fruit"){
      circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueEQ[i]/absMax, lty=3)
      circos.points(CELL_META$xcenter, 1 + relativeValueEQ[i]/absMax, pch=19, col=colour.circle.points[1], cex=0.7)
    }
    else if(relativeValueEQ[i] <= 0 & dietSpecies[i]=="Leaf"){
      circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueEQ[i]/absMax, lty=3)
      circos.points(CELL_META$xcenter, 1 + relativeValueEQ[i]/absMax, pch=21, col=colour.circle.points[1], bg="white", cex=0.7)
    }
    else{}
  }
}, track.height = 0.1)


#Striatum
absMax <- max(abs(relativeValueStriatum), na.rm=TRUE)
circos.track(ylim = c(0, 1), bg.border = NA, track.index=4, panel.fun = function(x, y) {
  i=CELL_META$sector.numeric.index
  #circos.rect(0, 0, 1, 1, col=colourPositive, border=colourPositive)
  if(is.na(relativeValueStriatum[i])){}  else{
    if(relativeValueStriatum[i] > 0 & dietSpecies[i]=="Fruit"){
      circos.points(CELL_META$xcenter, relativeValueStriatum[i]/absMax, pch=19, col=colour.circle.points[2], cex=0.65)
      circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueStriatum[i]/absMax, lty=3)
    }
    else if(relativeValueStriatum[i] > 0 & dietSpecies[i]=="Leaf"){
      circos.points(CELL_META$xcenter, relativeValueStriatum[i]/absMax, pch=21, col=colour.circle.points[2], bg="white", cex=0.65)
      circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueStriatum[i]/absMax, lty=3)
    }
    else{}
  }
}, track.height = 0.1)


circos.track(ylim = c(0, 1), bg.border = NA, track.index=5,  panel.fun = function(x, y) {
  i=CELL_META$sector.numeric.index
  if(is.na(relativeValueStriatum[i])){}  else{
    #circos.rect(0, 0, 1, 1, col=colourNegative, border=colourNegative)
    if(relativeValueStriatum[i] <= 0 & dietSpecies[i]=="Fruit"){
      circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueStriatum[i]/absMax, lty=3)
      circos.points(CELL_META$xcenter, 1 + relativeValueStriatum[i]/absMax, pch=19, col=colour.circle.points[2], cex=0.65)
    }
    else if(relativeValueStriatum[i] <= 0 & dietSpecies[i]=="Leaf"){
      circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueStriatum[i]/absMax, lty=3)
      circos.points(CELL_META$xcenter, 1 + relativeValueStriatum[i]/absMax, pch=21, col=colour.circle.points[2], bg="white", cex=0.65)
    }
    else{}
  }
}, track.height = 0.1)

#MOB
absMax <- max(abs(relativeValueMOB), na.rm=TRUE)
circos.track(ylim = c(0, 1), bg.border = NA, track.index=6, panel.fun = function(x, y) {
  i=CELL_META$sector.numeric.index
  #circos.rect(0, 0, 1, 1, col=colourPositive, border=colourPositive)
  if(is.na(relativeValueMOB[i])){}  else{
    if(relativeValueMOB[i] > 0 & dietSpecies[i]=="Fruit"){
      circos.points(CELL_META$xcenter, relativeValueMOB[i]/absMax, pch=19, col=colour.circle.points[3], cex=0.65)
      circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueMOB[i]/absMax, lty=3)
    }
    else if(relativeValueMOB[i] > 0 & dietSpecies[i]=="Leaf"){
      circos.points(CELL_META$xcenter, relativeValueMOB[i]/absMax, pch=21, col=colour.circle.points[3], bg="white", cex=0.65)
      circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueMOB[i]/absMax, lty=3)
    }
    else{}
  }
}, track.height = 0.1)


circos.track(ylim = c(0, 1), bg.border = NA, track.index=7,  panel.fun = function(x, y) {
  i=CELL_META$sector.numeric.index
  if(is.na(relativeValueMOB[i])){}  else{
    #circos.rect(0, 0, 1, 1, col=colourNegative, border=colourNegative)
    if(relativeValueMOB[i] <= 0 & dietSpecies[i]=="Fruit"){
      circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueMOB[i]/absMax, lty=3)
      circos.points(CELL_META$xcenter, 1 + relativeValueMOB[i]/absMax, pch=19, col=colour.circle.points[3], cex=0.65)
    }
    else if(relativeValueMOB[i] <= 0 & dietSpecies[i]=="Leaf"){
      circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueMOB[i]/absMax, lty=3)
      circos.points(CELL_META$xcenter, 1 + relativeValueMOB[i]/absMax, pch=21, col=colour.circle.points[3], bg="white", cex=0.65)
    }
    else{}
  }
}, track.height = 0.1)

#Neocortex
absMax <- max(abs(relativeValueNeocortex), na.rm=TRUE)
circos.track(ylim = c(0, 1), bg.border = NA, track.index=8, panel.fun = function(x, y) {
  i=CELL_META$sector.numeric.index
  #circos.rect(0, 0, 1, 1, col=colourPositive, border=colourPositive)
  if(is.na(relativeValueNeocortex[i])){}  else{
    if(relativeValueNeocortex[i] > 0 & dietSpecies[i]=="Fruit"){
      circos.points(CELL_META$xcenter, relativeValueNeocortex[i]/absMax, pch=19, col=colour.circle.points[4], cex=0.6)
      circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueNeocortex[i]/absMax, lty=3)
    }
    else if(relativeValueNeocortex[i] > 0 & dietSpecies[i]=="Leaf"){
      circos.points(CELL_META$xcenter, relativeValueNeocortex[i]/absMax, pch=21, col=colour.circle.points[4], bg="white", cex=0.6)
      circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueNeocortex[i]/absMax, lty=3)
    }
    else{}
  }
}, track.height = 0.1)


circos.track(ylim = c(0, 1), bg.border = NA, track.index=9,  panel.fun = function(x, y) {
  i=CELL_META$sector.numeric.index
  if(is.na(relativeValueNeocortex[i])){}  else{
    #circos.rect(0, 0, 1, 1, col=colourNegative, border=colourNegative)
    if(relativeValueNeocortex[i] <= 0 & dietSpecies[i]=="Fruit"){
      circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueNeocortex[i]/absMax, lty=3)
      circos.points(CELL_META$xcenter, 1 + relativeValueNeocortex[i]/absMax, pch=19, col=colour.circle.points[4], cex=0.6)
    }
    else if(relativeValueNeocortex[i] <= 0 & dietSpecies[i]=="Leaf"){
      circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueNeocortex[i]/absMax, lty=3)
      circos.points(CELL_META$xcenter, 1 + relativeValueNeocortex[i]/absMax, pch=21, col=colour.circle.points[4], bg="white", cex=0.6)
    }
    else{}
  }
}, track.height = 0.1)

#Hippocampus
absMax <- max(abs(relativeValueHippocampus), na.rm=TRUE)
circos.track(ylim = c(0, 1), bg.border = NA, track.index=10, panel.fun = function(x, y) {
  i=CELL_META$sector.numeric.index
  #circos.rect(0, 0, 1, 1, col=colourPositive, border=colourPositive)
  if(is.na(relativeValueHippocampus[i])){}  else{
    if(relativeValueHippocampus[i] > 0 & dietSpecies[i]=="Fruit"){
      circos.points(CELL_META$xcenter, relativeValueHippocampus[i]/absMax, pch=19, col=colour.circle.points[5], cex=0.55)
      circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueHippocampus[i]/absMax, lty=3)
    }
    else if(relativeValueHippocampus[i] > 0 & dietSpecies[i]=="Leaf"){
      circos.points(CELL_META$xcenter, relativeValueHippocampus[i]/absMax, pch=21, col=colour.circle.points[5], bg="white", cex=0.55)
      circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueHippocampus[i]/absMax, lty=3)
    }
    else{}
  }
}, track.height = 0.1)


circos.track(ylim = c(0, 1), bg.border = NA, track.index=11,  panel.fun = function(x, y) {
  i=CELL_META$sector.numeric.index
  if(is.na(relativeValueHippocampus[i])){}  else{
    #circos.rect(0, 0, 1, 1, col=colourNegative, border=colourNegative)
    if(relativeValueHippocampus[i] <= 0 & dietSpecies[i]=="Fruit"){
      circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueHippocampus[i]/absMax, lty=3)
      circos.points(CELL_META$xcenter, 1 + relativeValueHippocampus[i]/absMax, pch=19, col=colour.circle.points[5], cex=0.55)
    }
    else if(relativeValueHippocampus[i] <= 0 & dietSpecies[i]=="Leaf"){
      circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueHippocampus[i]/absMax, lty=3)
      circos.points(CELL_META$xcenter, 1 + relativeValueHippocampus[i]/absMax, pch=21, col=colour.circle.points[5], bg="white", cex=0.55)
    }
    else{}
  }
}, track.height = 0.1)

#Cerebellum
absMax <- max(abs(relativeValueCerebellum), na.rm=TRUE)
circos.track(ylim = c(0, 1), bg.border = NA, track.index=12, panel.fun = function(x, y) {
  i=CELL_META$sector.numeric.index
  #circos.rect(0, 0, 1, 1, col=colourPositive, border=colourPositive)
  if(is.na(relativeValueCerebellum[i])){}  else{
    if(relativeValueCerebellum[i] > 0 & dietSpecies[i]=="Fruit"){
      circos.points(CELL_META$xcenter, relativeValueCerebellum[i]/absMax, pch=19, col=colour.circle.points[6], cex=0.5)
      circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueCerebellum[i]/absMax, lty=3)
    }
    else if(relativeValueCerebellum[i] > 0 & dietSpecies[i]=="Leaf"){
      circos.points(CELL_META$xcenter, relativeValueCerebellum[i]/absMax, pch=21, col=colour.circle.points[6], bg="white", cex=0.5)
      circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueCerebellum[i]/absMax, lty=3)
    }
    else{}
  }
}, track.height = 0.1)


circos.track(ylim = c(0, 1), bg.border = NA, track.index=13,  panel.fun = function(x, y) {
  i=CELL_META$sector.numeric.index
  if(is.na(relativeValueCerebellum[i])){}  else{
    #circos.rect(0, 0, 1, 1, col=colourNegative, border=colourNegative)
    if(relativeValueCerebellum[i] <= 0 & dietSpecies[i]=="Fruit"){
      circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueCerebellum[i]/absMax, lty=3)
      circos.points(CELL_META$xcenter, 1 + relativeValueCerebellum[i]/absMax, pch=19, col=colour.circle.points[6], cex=0.5)
    }
    else if(relativeValueCerebellum[i] <= 0 & dietSpecies[i]=="Leaf"){
      circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueCerebellum[i]/absMax, lty=3)
      circos.points(CELL_META$xcenter, 1 + relativeValueCerebellum[i]/absMax, pch=21, col=colour.circle.points[6], bg="white", cex=0.5)
    }
    else{}
  }
}, track.height = 0.1)

###
## Fig brain
###

#from: https://www.r-bloggers.com/2018/10/how-to-highlight-3d-brain-regions/

library(rgl)
library(misc3d)
library(neurobase)
if (!requireNamespace("aal")) {
  devtools::install_github("muschellij2/aal")
} else {
  library(aal)
}
if (!requireNamespace("MNITemplate")) {
  devtools::install_github("jfortin1/MNITemplate")
} else {
  library(MNITemplate)
}
library(aal)
library(MNITemplate)
img = aal_image()
template = readMNI(res = "2mm")
cut <- 4500
dtemp <- dim(template)
# All of the sections you can label
labs = aal_get_labels()

# Pick the region of the brain you would like to highlight - in this case the hippocamus_L
hippocampus = labs$index[grep("Hippocampus", labs$name)]
cerebellum = labs$index[grep("Cerebellum", labs$name)]
striatum = c(labs$index[grep("Caudate", labs$name)], labs$index[grep("Putamen", labs$name)], labs$index[grep("Pallidum", labs$name)])
MOB = labs$index[grep("Olfactory_L", labs$name)]
#neocortex = c(labs$index[grep("Frontal", labs$name)], labs$index[grep("Parietal", labs$name)], labs$index[grep("Temporal", labs$name)],
#labs$index[grep("Occipital", labs$name)], labs$index[grep("Cingulate", labs$name)], labs$index[grep("Insula", labs$name)]) #incomplete


colourHip=colourVector[4]
colourCereb=colourVector[6]
colourOlf=colourVector[3]
colourStri=colourVector[2]

### this would be the ``activation'' or surface you want to render 

#Set up the view manually the first time to extract coords
#um <- par3d()$userMatrix
#I obtained those as a correct view
um <- matrix(
  c(-0.4987500,-0.80934894,-0.3101658,0,-0.4182398,-0.08870265,0.9039953,0,-0.7591600,0.58059108,-0.2942616,0,0.0000000,0.00000000,0.0000000,1), 
  ncol=4, nrow=4, byrow=TRUE)
view3d(userMatrix = um)
#c(-0.6391721,-0.76426947,-0.08573965,0,-0.2281326,0.08195215,0.97017491,0,-0.7344485,0.63966882,-0.22673632,0,0.0000000,0.00000000,0.00000000,1)
#c(-0.7031180,0.4854550,-0.5195753,0.0000000,-0.5022486,0.1781951,0.8461637,0.0000000,0.5033601,0.8559089,0.1185269,0.0000000,0.0000000,0.0000000,0.0000000,1.0000000)


contour3d(template, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], level = cut, alpha = 0.03, draw = TRUE)

mask = remake_img(vec = img %in% hippocampus, img = img)
contour3d(mask, level = c(0.5), alpha = c(0.4), add = TRUE, color=colourHip)
mask = remake_img(vec = img %in% cerebellum, img = img)
contour3d(mask, level = c(0.5), alpha = c(0.4), add = TRUE, color=colourCereb)
mask = remake_img(vec = img %in% MOB, img = img)
contour3d(mask, level = c(0.5), alpha = c(0.4), add = TRUE, color=colourOlf)
#mask = remake_img(vec = img %in% neocortex, img = img)
#contour3d(mask, level = c(0.5), alpha = c(0.4), add = TRUE, color=c("orange"))
mask = remake_img(vec = img %in% striatum, img = img)
contour3d(mask, level = c(0.5), alpha = c(0.4), add = TRUE, color=colourStri)


### add text
#text3d(x=dtemp[1]/2, y=dtemp[2]/2, z = dtemp[3]*1.1, text="Top")
#text3d(x=-0.98, y=dtemp[2]/2, z = dtemp[3]/2, text="Right")
#text3d(x=dtemp[1]/2, y=dtemp[2], z = dtemp[3]/2, text="Front")
#text3d(x=dtemp[1]/2, y=-1, z = dtemp[3]/2, text="Back")

#segmentBottomTop
arrow3d(p0=c(dtemp[1]/2,dtemp[2]/2,dtemp[3]), p1=c(dtemp[1]/2,dtemp[2]/2,1.17*dtemp[3]), type="extrusion")
#segmentBackFront
arrow3d(p0=c(dtemp[1]/2,dtemp[2],dtemp[3]/2), p1=c(dtemp[1]/2,1.1*dtemp[2],dtemp[3]/2), type="extrusion")

#legend3d(x=0.1, y=0.2, legend = c("Cerebellum", "Hippocampus", "MOB", "Striatum"), fill = c(colourCereb, colourHip, colourOlf, colourStri), cex=1, bty="n", ncol=2)
#rglwidget()


#from: https://r-graphics.org/recipe-miscgraph-3d-save
#rgl.snapshot('3dplot.png', fmt = 'png', top = TRUE, width = 300)
#snapshot3d('Plots/3dplot.png', fmt = 'png', top = TRUE, webshot = FALSE, width=32000, height=32000)
#rgl.postscript('3dplot.pdf', fmt = 'pdf')

#remotes::install_github("rstudio/chromote")
#remotes::install_github("rstudio/webshot2")

library(webshot2)
snapshot3d('Plots/3dplot.png', fmt = 'png')

writeASY()
rgl.postscript('3dplot', fmt="svg")

plot(x=0, y=0,
     xlab="", ylab="", 
     xlim=c(0,1), ylim=c(0,1),
     las=1,
     type="n", 
     tcl=-0.25, frame.plot=FALSE, 
     xaxt="n",xaxs="i",yaxs="i", yaxt="n"
)
library(png)
brainIMG <- readPNG("C:/Users/robira/Documents/PhD/Meta_analysis/Cognition_metaanalysis/3dplot.png")
addImg(brainIMG, x = 0.5, y = 0.5, width = 1)


########
## Fig 3:
########

##---------
## Check sample size
repetition=2*2*2*10#length(frugivoryThresholdVector)*length(folivoryThresholdVector)*length(geographicThresholdVector)*randomSampling
  
transitionMatrix <- matrix(NA, nrow=repetition, ncol=2)

for(a in 1:2){
  for(b in 1:2){
    for(c in 1:2){
      for(d in 1:10){
        start=which(is.na(transitionMatrix[,1]))[1]
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_simmap_transition",a,"_",b,"_",c,"_",d, ".txt", sep=""))
          transitionMatrix[start,1] <- toAdd[1,1]
          transitionMatrix[start,2] <- toAdd[2,1]
          }, error=function(e){
            #Do nothing
          }
        )
      }
    }
  }
}

minProba.v <- apply(transitionMatrix, 2, min)
maxProba.v <- apply(transitionMatrix, 2, max)


##---------
## Check transition matrix


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

counter=0
for(a in 1:2){
  for(b in 1:2){
    for(c in 1:2){
      for(d in 1:10){
        counter=counter+1
        tryCatch(
        {toAdd <- read.delim(paste("Processed_data/Sample_size/checkSampleFruit",a,"_",b,"_",c,"_",d, ".txt", sep=""))
         checkSampleFruit[counter] <- toAdd[1]
        }, error=function(e){
          #Do nothing
        }
        )
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/Sample_size/checkSampleLeaf",a,"_",b,"_",c,"_",d, ".txt", sep=""))
           checkSampleLeaf[counter] <- toAdd[1]
          }, error=function(e){
            #Do nothing
          }
        )
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/Sample_size/checkSampleRange",a,"_",b,"_",c,"_",d, ".txt", sep=""))
           checkSampleRange[counter] <- toAdd[1]
          }, error=function(e){
            #Do nothing
          }
        )
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/Sample_size/checkSampleBrain",a,"_",b,"_",c,"_",d, ".txt", sep=""))
           checkSampleBrain[counter] <- toAdd[1]
          }, error=function(e){
            #Do nothing
          }
        )
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/Sample_size/checkSampleEQ",a,"_",b,"_",c,"_",d, ".txt", sep=""))
           checkSampleEQ[counter] <- toAdd[1]
          }, error=function(e){
            #Do nothing
          }
        )
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/Sample_size/checkSampleNeocortex",a,"_",b,"_",c,"_",d, ".txt", sep=""))
           checkSampleNeocortex[counter] <- toAdd[1]
          }, error=function(e){
            #Do nothing
          }
        )
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/Sample_size/checkSampleHippocampus",a,"_",b,"_",c,"_",d, ".txt", sep=""))
           checkSampleHippocampus[counter] <- toAdd[1]
          }, error=function(e){
            #Do nothing
          }
        )
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/Sample_size/checkSampleCerebellum",a,"_",b,"_",c,"_",d, ".txt", sep=""))
           checkSampleCerebellum[counter] <- toAdd[1]
          }, error=function(e){
            #Do nothing
          }
        )
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/Sample_size/checkSampleStriatum",a,"_",b,"_",c,"_",d, ".txt", sep=""))
           checkSampleStriatum[counter] <- toAdd[1]
          }, error=function(e){
            #Do nothing
          }
        )
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/Sample_size/checkSampleMOB",a,"_",b,"_",c,"_",d, ".txt", sep=""))
           checkSampleMOB[counter] <- toAdd[1]
          }, error=function(e){
            #Do nothing
          }
        )
      }
    }
  }
}

checkSampleFruit <- unlist(checkSampleFruit)
checkSampleLeaf <- unlist(checkSampleLeaf)
checkSampleRange <- unlist(checkSampleRange)
checkSampleBrain <-  unlist(checkSampleBrain)
checkSampleEQ <-  unlist(checkSampleEQ)
checkSampleNeocortex <-  unlist(checkSampleNeocortex)
checkSampleHippocampus <- unlist(checkSampleHippocampus)
checkSampleCerebellum <- unlist(checkSampleCerebellum)
checkSampleStriatum <- unlist(checkSampleStriatum)
checkSampleMOB <- unlist(checkSampleMOB)


#Min-max sample
paste(c(min(checkSampleFruit), max(checkSampleFruit)), collapse="-")
paste(c(min(checkSampleLeaf), max(checkSampleLeaf)), collapse="-")
paste(c(min(checkSampleRange), max(checkSampleRange)), collapse="-")
paste(c(min(checkSampleBrain), max(checkSampleBrain)), collapse="-")
paste(c(min(checkSampleEQ), max(checkSampleEQ)), collapse="-")
paste(c(min(checkSampleNeocortex), max(checkSampleNeocortex)), collapse="-")
paste(c(min(checkSampleHippocampus), max(checkSampleHippocampus)), collapse="-")
paste(c(min(checkSampleCerebellum), max(checkSampleCerebellum)), collapse="-")
paste(c(min(checkSampleStriatum), max(checkSampleStriatum)), collapse="-")
paste(c(min(checkSampleMOB), max(checkSampleMOB)), collapse="-")

##----------
## Checking the effect of scaling: Neocortex as an example
summaryNeocortexFrugivory <- as.data.frame(matrix(NA, nrow=10*(repetition+1), ncol=53))
summaryNeocortexRawFrugivory <- as.data.frame(matrix(NA, nrow=10*(repetition+1), ncol=53))
summaryNeocortexBodymassFrugivory <- as.data.frame(matrix(NA, nrow=10*(repetition+1), ncol=53))
summaryNeocortexBodymassRawFrugivory <- as.data.frame(matrix(NA, nrow=10*(repetition+1), ncol=53))

numberSimulations=10

for(i in 1:((10))){
  
  start=which(is.na(summaryNeocortexFrugivory[,1]))[1]
  end=which(is.na(summaryNeocortexFrugivory[,1]))[1] + numberSimulations - 1
  
  tryCatch(
    {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_Neocortex1_1_1_",i, ".txt", sep=""))
    summaryNeocortexFrugivory[start:end,] <- toAdd
    }, error=function(e){
      #Do nothing
    }
  )
  
  tryCatch(
    {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_Neocortex1_1_1_",i, ".txt", sep=""))
    summaryNeocortexFrugivory[start:end,] <- toAdd
    }, error=function(e){
      #Do nothing
    }
  )
  
  tryCatch(
    {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_NeocortexRaw1_1_1_",i, ".txt", sep=""))
    summaryNeocortexRawFrugivory[start:end,] <- toAdd
    }, error=function(e){
      #Do nothing
    }
  )
  
  tryCatch(
    {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_NeocortexBodymass1_1_1_",i, ".txt", sep=""))
    summaryNeocortexBodymassFrugivory[start:end,] <- toAdd
    }, error=function(e){
      #Do nothing
    }
  )
  
  tryCatch(
    {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_NeocortexBodymassRaw1_1_1_",i, ".txt", sep=""))
    summaryNeocortexBodymassRawFrugivory[start:end,] <- toAdd
    }, error=function(e){
      #Do nothing
    }
  )
}

colnames(summaryNeocortexFrugivory) <- colnames(toAdd)
colnames(summaryNeocortexRawFrugivory) <- colnames(toAdd)
colnames(summaryNeocortexBodymassFrugivory) <- colnames(toAdd)
colnames(summaryNeocortexBodymassRawFrugivory) <- colnames(toAdd)

colNum <- c("darkgray", brewer.pal(n = 5, name = "Set1")[5], brewer.pal(n = 4, name = "Set1"))

models <- c("BM", "OU", "EB", "MC", "DDlin", "DDexp")
colourModels <- brewer.pal(n = 6, name = "Set1")

layout(mat=rbind(c(1,2), c(3,4)), widths=c(5,5), heights=c(5,5))
par(mar=c(4, 4, 2, 1), mgp=c(2, 0.5, 0), xpd=TRUE)

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6, main="Ratio brain",
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
for(i in 1:6){
  
  meanPt <- mean(as.numcharac(summaryNeocortexFrugivory[, ncol(summaryNeocortexFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryNeocortexFrugivory[, ncol(summaryNeocortexFrugivory)-6+i]))
  #sd <- sd/nrow(summaryNeocortexFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
draw.circle(x=0.3,y=1.1,0.2, col=colNum[2], border=NA)
text(x=0.3, y=1.1, labels="5", xpd=TRUE, col="white", font=2)

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(summaryNeocortexFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryNeocortexFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)


plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6, main="Raw",
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
for(i in 1:6){
  
  meanPt <- mean(as.numcharac(summaryNeocortexRawFrugivory[, ncol(summaryNeocortexRawFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryNeocortexRawFrugivory[, ncol(summaryNeocortexRawFrugivory)-6+i]))
  #sd <- sd/nrow(summaryNeocortexRawFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
draw.circle(x=0.3,y=1.1,0.2, col=colNum[2], border=NA)
text(x=0.3, y=1.1, labels="5", xpd=TRUE, col="white", font=2)

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(summaryNeocortexRawFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryNeocortexRawFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6, main="Bodymass EQ",
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
for(i in 1:6){
  
  meanPt <- mean(as.numcharac(summaryNeocortexBodymassFrugivory[, ncol(summaryNeocortexBodymassFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryNeocortexBodymassFrugivory[, ncol(summaryNeocortexBodymassFrugivory)-6+i]))
  #sd <- sd/nrow(summaryNeocortexBodymassFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
draw.circle(x=0.3,y=1.1,0.2, col=colNum[2], border=NA)
text(x=0.3, y=1.1, labels="5", xpd=TRUE, col="white", font=2)

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(summaryNeocortexBodymassFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryNeocortexBodymassFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)


plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6, main="Bodymass raw",
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
for(i in 1:6){
  
  meanPt <- mean(as.numcharac(summaryNeocortexBodymassRawFrugivory[, ncol(summaryNeocortexBodymassRawFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryNeocortexBodymassRawFrugivory[, ncol(summaryNeocortexBodymassRawFrugivory)-6+i]))
  #sd <- sd/nrow(summaryNeocortexBodymassRawFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
draw.circle(x=0.3,y=1.1,0.2, col=colNum[2], border=NA)
text(x=0.3, y=1.1, labels="5", xpd=TRUE, col="white", font=2)

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(summaryNeocortexBodymassRawFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryNeocortexBodymassRawFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)
##----------


##---------
## Load the data

summaryBrainFrugivory <- as.data.frame(matrix(NA, nrow=10*(repetition+1), ncol=53))
summaryEQFrugivory <- as.data.frame(matrix(NA, nrow=10*(repetition+1), ncol=53))
summaryNeocortexFrugivory <- as.data.frame(matrix(NA, nrow=10*(repetition+1), ncol=53))
summaryHippocampusFrugivory <- as.data.frame(matrix(NA, nrow=10*(repetition+1), ncol=53))
summaryCerebellumFrugivory <- as.data.frame(matrix(NA, nrow=10*(repetition+1), ncol=53))
summaryStriatumFrugivory <- as.data.frame(matrix(NA, nrow=10*(repetition+1), ncol=53))
summaryMOBFrugivory <- as.data.frame(matrix(NA, nrow=10*(repetition+1), ncol=53))

counter=0
start=counter
end=counter
numberSimulations=10
for(a in 1:2){
  for(b in 1:2){
    for(c in 1:2){
      for(d in 1:10){
        counter=end+1
        start=counter
        end=counter + numberSimulations - 1
              
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_BrainBodymassRaw",a,"_",b,"_",c,"_",d,".txt", sep=""))
        summaryBrainFrugivory[start:end,] <- as.data.frame(toAdd)
        }, error=function(e){
          #Do nothing
          }
        )
        
        tryCatch(
        {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_EQ",a,"_",b,"_",c,"_",d,".txt", sep=""))
        summaryEQFrugivory[start:end,] <- toAdd
        }, error=function(e){
          #Do nothing
        }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_NeocortexBodymassRaw",a,"_",b,"_",c,"_",d,".txt", sep=""))
        summaryNeocortexFrugivory[start:end,] <- toAdd
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_HippocampusBodymassRaw",a,"_",b,"_",c,"_",d,".txt", sep=""))
        summaryHippocampusFrugivory[start:end,] <- toAdd
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_CerebellumBodymassRaw",a,"_",b,"_",c,"_",d,".txt", sep=""))
        summaryCerebellumFrugivory[start:end,] <- toAdd
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_StriatumBodymassRaw",a,"_",b,"_",c,"_",d,".txt", sep=""))
        summaryStriatumFrugivory[start:end,] <- toAdd
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_MOBBodymassRaw",a,"_",b,"_",c,"_",d,".txt", sep=""))
        summaryMOBFrugivory[start:end,] <- toAdd
          }, error=function(e){
            #Do nothing
          }
        )
      }
    }
  }
}

summaryBrainFrugivory <- summaryBrainFrugivory[!is.na(summaryBrainFrugivory[,1]),]
summaryEQFrugivory <- summaryEQFrugivory[!is.na(summaryEQFrugivory[,1]),]
summaryNeocortexFrugivory <- summaryNeocortexFrugivory[!is.na(summaryNeocortexFrugivory[,1]),]
summaryHippocampusFrugivory <- summaryHippocampusFrugivory[!is.na(summaryHippocampusFrugivory[,1]),]
summaryCerebellumFrugivory <- summaryCerebellumFrugivory[!is.na(summaryCerebellumFrugivory[,1]),]
summaryStriatumFrugivory <- summaryStriatumFrugivory[!is.na(summaryStriatumFrugivory[,1]),]
summaryMOBFrugivory <- summaryMOBFrugivory[!is.na(summaryMOBFrugivory[,1]),]

colnames(summaryBrainFrugivory) <- colnames(toAdd)
colnames(summaryEQFrugivory) <- colnames(toAdd)
colnames(summaryNeocortexFrugivory) <- colnames(toAdd)
colnames(summaryHippocampusFrugivory) <- colnames(toAdd)
colnames(summaryCerebellumFrugivory) <- colnames(toAdd)
colnames(summaryStriatumFrugivory) <- colnames(toAdd)
colnames(summaryMOBFrugivory) <- colnames(toAdd)

##----

colNum <-c("darkgrey", brewer.pal(n = 5, name = "Set1")[1:5])

models <- c("BM", "OU", "EB", "MC", "DDlin", "DDexp")
colourModels <- brewer.pal(n = 6, name = "Set1")
# 
# ## Brain
# 
# plot(
#   x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6, 
#   xlim=c(0,7), ylim=c(0,1.2),
#   las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
#   xaxt="n",xaxs="i",yaxs="i", yaxt="n")
# 
# addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
# axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
# text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
# for(i in 1:6){
#   meanPt <- mean(as.numcharac(summaryBrainFrugivory[, ncol(summaryBrainFrugivory)-6+i]))
#   sd <- sd(as.numcharac(summaryBrainFrugivory[, ncol(summaryBrainFrugivory)-6+i]))
#   #sd <- sd/nrow(summaryBrainFrugivory) #error not sd
#   errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
#   points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
# }
# 
# ## EQ
# 
# plot(
#   x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
#   xlim=c(0,7), ylim=c(0,1.2),
#   las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
#   xaxt="n",xaxs="i",yaxs="i", yaxt="n")
# 
# addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
# axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
# text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
# for(i in 1:6){
#   meanPt <- mean(as.numcharac(summaryEQFrugivory[, ncol(summaryEQFrugivory)-6+i]))
#   sd <- sd(as.numcharac(summaryEQFrugivory[, ncol(summaryEQFrugivory)-6+i]))
#   #sd <- sd/nrow(summaryEQFrugivory) #error not sd
#   errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
#   points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
#   
# }
# 
# #b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
# #add their values below:
# 
# text(x=c(5,6), y=c(-0.2, -0.2),
#      labels=c(
#        paste("r~", round(mean(as.numcharac(summaryEQFrugivory$DDlingeo.b)), digit=3), sep=""),
#        paste("r~", round(mean(as.numcharac(summaryEQFrugivory$DDexpgeo.r)), digit=3), sep="")
#      ), xpd=TRUE)
# 
# draw.circle(x=0.3,y=1.1,0.2, col=colNum[1], border=NA)
# text(x=0.3, y=1.1, labels="1", xpd=TRUE, col="white", font=2, cex=1.5)
# 
# ##----------------------


layout(mat=rbind(c(1,2,3), c(4,5,6)), widths=c(5,5,5), heights=c(5,5))
par(mar=c(4, 4, 2, 1), mgp=c(2, 0.5, 0), xpd=TRUE)
#note: 1= second run for frugivory 20%
#note: _2= first run for frugivory 20%

## EQ

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
for(i in 1:6){
  meanPt <- mean(as.numcharac(summaryEQFrugivory[, ncol(summaryEQFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryEQFrugivory[, ncol(summaryEQFrugivory)-6+i]))
  #sd <- sd/nrow(summaryEQFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}

#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("r~", round(mean(as.numcharac(summaryEQFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryEQFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

draw.circle(x=0.3,y=1.1,0.2, col=colNum[1], border=NA)
text(x=0.3, y=1.1, labels="1", xpd=TRUE, col="white", font=2)
text(x=3.5, y=1.1, labels="EQ", xpd=TRUE, col="black", font=2, cex=1.5)

##-------------

##------------
#Striatum

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
for(i in 1:6){
  
  meanPt <- mean(as.numcharac(summaryStriatumFrugivory[, ncol(summaryStriatumFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryStriatumFrugivory[, ncol(summaryStriatumFrugivory)-6+i]))
  #sd <- sd/nrow(summaryStriatumFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
draw.circle(x=0.3,y=1.1,0.2, col=colNum[2], border=NA)
text(x=0.3, y=1.1, labels="2", xpd=TRUE, col="white", font=2)
text(x=3.5, y=1.1, labels="Striatum", xpd=TRUE, col="black", font=2, cex=1.5)


#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("r~", round(mean(as.numcharac(summaryStriatumFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryStriatumFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)


##------------

##-------------
#MOB

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
for(i in 1:6){
  
  meanPt <- mean(as.numcharac(summaryMOBFrugivory[, ncol(summaryMOBFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryMOBFrugivory[, ncol(summaryMOBFrugivory)-6+i]))
  #sd <- sd/nrow(summaryMOBFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
draw.circle(x=0.3,y=1.1,0.2, col=colNum[3], border=NA)
text(x=0.3, y=1.1, labels="3", xpd=TRUE, col="white", font=2)
text(x=3.5, y=1.1, labels="MOB", xpd=TRUE, col="black", font=2, cex=1.5)

#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("r~", round(mean(as.numcharac(summaryMOBFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryMOBFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

##------------

##------------

#Hippocampus

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
for(i in 1:6){
  
  meanPt <- mean(as.numcharac(summaryHippocampusFrugivory[, ncol(summaryHippocampusFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryHippocampusFrugivory[, ncol(summaryHippocampusFrugivory)-6+i]))
  #sd <- sd/nrow(summaryHippocampusFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
draw.circle(x=0.3,y=1.1,0.2, col=colNum[4], border=NA)
text(x=0.3, y=1.1, labels="6", xpd=TRUE, col="white", font=2)
text(x=3.5, y=1.1, labels="Hippocampus", xpd=TRUE, col="black", font=2, cex=1.5)


#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("r~", round(mean(as.numcharac(summaryHippocampusFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryHippocampusFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

##------------

##-------------
#Neocortex

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
for(i in 1:6){
  
  meanPt <- mean(as.numcharac(summaryNeocortexFrugivory[, ncol(summaryNeocortexFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryNeocortexFrugivory[, ncol(summaryNeocortexFrugivory)-6+i]))
  #sd <- sd/nrow(summaryNeocortexFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
draw.circle(x=0.3,y=1.1,0.2, col=colNum[5], border=NA)
text(x=0.3, y=1.1, labels="5", xpd=TRUE, col="white", font=2)
text(x=3.5, y=1.1, labels="Neocortex", xpd=TRUE, col="black", font=2, cex=1.5)


#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("r~", round(mean(as.numcharac(summaryNeocortexFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryNeocortexFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)


##-------------

##-------------
#Cerebellum 

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
for(i in 1:6){
  
  meanPt <- mean(as.numcharac(summaryCerebellumFrugivory[, ncol(summaryCerebellumFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryCerebellumFrugivory[, ncol(summaryCerebellumFrugivory)-6+i]))
  #sd <- sd/nrow(summaryCerebellumFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
draw.circle(x=0.3,y=1.1,0.2, col=colNum[6], border=NA)
text(x=0.3, y=1.1, labels="7", xpd=TRUE, col="white", font=2)
text(x=3.5, y=1.1, labels="Cerebellum", xpd=TRUE, col="black", font=2, cex=1.5)


#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("r~", round(mean(as.numcharac(summaryCerebellumFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryCerebellumFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

###----------------------


#############
#############
## GARBAGE
# 
# 
# #--------------------------------------------------------------------------------
# 
# ##########
# ## STATISTICAL ANALYSES: IS COGNITION PROMOTED BY INTERSPECIES COMPETITION?
# ##########
# 
# 
# nameForPhylo <- as.data.frame(phylo$tip.label)
# colnames(nameForPhylo) <- "Species"
# nameForPhylo$Species <- as.character(nameForPhylo$Species)
# 
# options(warn=1)
# nameForPhylo <- separate(nameForPhylo, col=Species, into=c("Name1", "Name2"), sep="_", remove=FALSE)#some have three "_", bot a big deal for now
# nameForPhylo$Species_abbrv <- paste(nameForPhylo$Name1, substr(nameForPhylo$Name2,1,4), sep="_")
# summaryKamilar$SpeciesPhylo <- nameForPhylo$Species[match(summaryKamilar$Species_abbrv,nameForPhylo$Species_abbrv)]
# 
# #Phylogenetical tree
# 
# library(ape)
# library(phytools)
# #Consensus tree from 10k Trees
# phylo_all <-read.nexus("Raw_data/Tree/consensusTree_10kTrees_Primates_Version3.nex")
# phylo_init <- phylo_all
# 
# is.ultrametric(phylo_init)
# #Tree is not ultrametric while it should
# phylo <- force.ultrametric(tree=phylo_init, method="extend")#method="nnls")
# is.ultrametric(phylo)
# 
# #Check changes
# plot(phylo_init$edge.length, phylo$edge.length, xlab="Original tree", ylab="Post-chronos tree", pch=20, bty="n")#no much changes; good
# abline(a=0, b=1, col="gray", lwd=0.5)
# 
# #Matching name of dataset to phylogenetics name of tree
# nameForPhylo <- as.data.frame(phylo$tip.label)
# colnames(nameForPhylo) <- "Species"
# nameForPhylo$Species <- as.character(nameForPhylo$Species)
# 
# library(tidyr)
# options(warn=1)
# nameForPhylo <- separate(nameForPhylo, col=Species, into=c("Name1", "Name2"), sep="_", remove=FALSE)#some have three "_", bot a big deal for now
# nameForPhylo$Species_abbrv <- paste(nameForPhylo$Name1, substr(nameForPhylo$Name2,1,4), sep="_")
# summaryKamilar$SpeciesPhylo <- nameForPhylo$Species[match(summaryKamilar$Species_abbrv,nameForPhylo$Species_abbrv)]
# 
# dataForBrain <- summaryKamilar[, c(1,2,3,5,6,7,10,13,14,15,16,17,18,ncol(summaryKamilar))]
# dataForBrain <- dataForBrain[!is.na(dataForBrain$OverlapHR)&
#                                !is.na(dataForBrain$Brain)&
#                                !is.na(dataForBrain$Bodymass),]
# dataForBrain <- dataForBrain[!is.na(dataForBrain$SpeciesPhylo),]
# nrow(dataForBrain)#lost one
# 
# '%nin%' <- Negate('%in%')
# 
# dataForBrain_rdc <- dataForBrain[
#   dataForBrain$DietGuild=="Fruit"|
#     dataForBrain$DietGuild=="Leaves",]
# 
# phylo <- drop.tip(phylo,
#                   phylo$tip.label[
#                     which(phylo$tip.label
#                           %nin%unique(dataForBrain_rdc$SpeciesPhylo))]) 
# 
# 
# 
# 
# ####
# ## Fig interaction
# ####
# 
# #TO MODIFY WITH OWN DATA, CONSTRUCTED ON A TEST DATASET
# library(ape)
# #data(bird.orders)
# 
# numberOrder=6
# hc = as.hclust(tree)#bird.orders)
# labels = hc$labels  # name of birds
# 
# labels.tordc <- as.data.frame(labels)
# colnames(labels.tordc) <- "Name"
# labels.tordc <- separate(labels.tordc, col="Name", into=c("Name1", "Name2", "Name3"), sep="_")
# 
# labels.rdc <- apply(labels.tordc, 1, function(x){
#   if(!is.na(x[3])){
#     paste(toupper(substr(x[1], 1, 1)), ". ", substr(x[2], 1, 3), ". ", substr(x[3], 1, 1), ".", sep="")
#   } else{paste(toupper(substr(x[1], 1, 1)), ". ", substr(x[2], 1, 3), sep="")
#   }
# }
# )
# 
# #Link between species/groups:
# speciesLabels <- hc$labels#Should be in the tree order
# #groupLabels <- rep(seq(from=1, to=23, by=2), each=2)[-1]
# 
# 
# #Geo + diet
# geoBinary <- as.data.frame(summaryDataForPlot)
# colnames(geoBinary) <- c("SpeciesPhylo", "Loc")
# locationSpecies <- geoBinary$Loc[match(speciesLabels, geoBinary$SpeciesPhylo)]
# dietSpecies <- summaryKamilar$Guild[match(speciesLabels, summaryKamilar$SpeciesPhylo)]
# colLoc <- c(brewer.pal(n = 9, name = "Pastel1"), "gray93")#1:length(unique(tableSiteGeo$Locality))#brewer.pal(n = length(unique(tableSiteGeo$Locality)), name = "Set1")
# loc.v <- tableSiteGeo$Locality
# 
# #Brain data
# relativeValueBrain <- summaryKamilar$EQ[match(speciesLabels, summaryKamilar$SpeciesPhylo)] - 1#runif(length(speciesLabels), -1, 1)
# #relativeValueBrain <- scale(log(relativeValueBrain))
# positiveCol="magenta4"#pastellize(x="purple", p=0.2)
# negativeCol="darkgreen"#pastellize(x="green", p=0.2)
# 
# 
# library(circlize)
# circos.clear()
# circos.par(gap.degree=0, gap.after=0, cell.padding=c(0,0,0,0), track.margin = c(0, 0))
# circos.initialize(speciesLabels, xlim = c(0, 1))
# 
# 
# ###TRIAL 1
# #Brain size + diet
# 
# absMax <- max(abs(relativeValueBrain))
# 
# library(RColorBrewer)
# colourPositive=brewer.pal(n = 5, name = "Pastel1")[1]
# colourNegative=brewer.pal(n = 5, name = "Pastel1")[2]
# 
# #Plot relative brain size
# #Background
# circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
#   circos.rect(0, 0, 1, 1, col=colourPositive, border=colourPositive)
# }, track.height = 0.1)
# 
# circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
#   circos.rect(0, 0, 1, 1, col=colourNegative, border=colourNegative)
# }, track.height = 0.1)
# 
# library(plotrix)
# #Main circle
# draw.circle(x=0,y=0,0.9, col=NA, border="white")
# 
# #increment of 0.5
# for(i in 1:5){
#   draw.circle(x=0,y=0,1-(i-1)*0.05, col=NA, border="white", lty=2)
# }
# 
# 
# #Value
# #for(i in 1:length(speciesLabels)){
# circos.track(ylim = c(0, 1), bg.border = NA, track.index=1, panel.fun = function(x, y) {
#   i=CELL_META$sector.numeric.index
#   #circos.rect(0, 0, 1, 1, col=colourPositive, border=colourPositive)
#   if(relativeValueBrain[i] > 0 & dietSpecies[i]=="Fruit"){
#     circos.points(CELL_META$xcenter, relativeValueBrain[i]/absMax, pch=19, col="black")
#     circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueBrain[i]/absMax, lty=3)
#   }
#   else if(relativeValueBrain[i] > 0 & dietSpecies[i]=="Leaves"){
#     circos.points(CELL_META$xcenter, relativeValueBrain[i]/absMax, pch=21, col="black", bg="white")
#     circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueBrain[i]/absMax, lty=3)
#   }
#   else{}
# }, track.height = 0.1)
# 
# circos.track(ylim = c(0, 1), bg.border = NA, track.index=2,  panel.fun = function(x, y) {
#   i=CELL_META$sector.numeric.index
#   #circos.rect(0, 0, 1, 1, col=colourNegative, border=colourNegative)
#   if(relativeValueBrain[i] <= 0 & dietSpecies[i]=="Fruit"){
#     circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueBrain[i]/absMax, lty=3)
#     circos.points(CELL_META$xcenter, 1 + relativeValueBrain[i]/absMax, pch=19, col="black")
#   }
#   else if(relativeValueBrain[i] <= 0 & dietSpecies[i]=="Leaves"){
#     circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueBrain[i]/absMax, lty=3)
#     circos.points(CELL_META$xcenter, 1 + relativeValueBrain[i]/absMax, pch=21, col="black", bg="white")
#   }
#   else{}
# }, track.height = 0.1)
# 
# #Species name
# circos.track(ylim = c(0, 20), bg.border = NA, track.height = 0.1, track.margin=c(0.01, 0.1),
#              panel.fun = function(x, y) {
#                i=CELL_META$sector.numeric.index
#                circos.text(CELL_META$xcenter, 0, labels.rdc[i], adj = c(0, 0), 
#                            facing = "clockwise", niceFacing = TRUE,
#                            col = "black", cex = 0.5, font=3)
#              })
# 
# 
# #Plot the geographic links
# for(i in 1:length(speciesLabels)){
#   #locI <- which(strsplit(locationSpecies[i], "")==1)
#   for(j in i:length(speciesLabels)){
#     #locJ <- which(strsplit(locationSpecies[j], "")==1)
#     product <- as.numcharac(unlist(strsplit(locationSpecies[j], "")))*as.numcharac(unlist(strsplit(locationSpecies[i], "")))
#     if(i==j|(length(unique(product))==1&unique(product)[1]==0)){
#       #Do nothing
#     }
#     else{
#       if(dietSpecies[i]==dietSpecies[j]){
#         colour <- as.data.frame(table(colLoc[which(product==1)]))
#         colour <- colour[colour$Freq==max(colour$Freq),1][1]
#         circos.link(speciesLabels[i], runif(1, 0, 1), speciesLabels[j], runif(1, 0, 1), lwd=1, col=as.character(colour))
#       }
#       else{
#         #circos.link(speciesLabels[i], runif(1, 0, 1), speciesLabels[j], runif(1, 0, 1), lwd=1, col="lightgray")
#       }
#     }
#   }
# }
# circos.clear()
# 
# 
# #Plot the phylogenetic tree in a new circular plot
# ct = cutree(hc, numberOrder)  # cut tree into order
# n = length(labels)  # number of bird species
# dend = as.dendrogram(hc)
# 
# library(circlize)
# par(new = TRUE) # <- magic
# circos.par("canvas.xlim" = c(-1.75, 1.75), "canvas.ylim" = c(-1.75, 1.75))
# circos.initialize("a", xlim = c(0, n)) # only one sector
# # circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.3, 
# #              panel.fun = function(x, y) {
# #                for(i in seq_len(n)) {
# #                  circos.text(i-0.5, 0, labels.rdc[i], adj = c(0, 0.5), 
# #                              facing = "clockwise", niceFacing = TRUE,
# #                              col = "black", cex = 0.2, font=3)
# #                }
# #              })
# 
# 
# #suppressPackageStartupMessages(library(dendextend))
# #dend = color_branches(dend, k = 6, col = 1:6)
# dend_height = attr(dend, "height")
# circos.track(ylim = c(0, dend_height), bg.border = NA, 
#              track.height = 0.95, panel.fun = function(x, y) {
#                circos.dendrogram(dend)
#              })
# circos.clear()
# 
# 
# 
# 
# ###TRIAL 2
# library(circlize)
# circos.clear()
# circos.par(gap.degree=0, gap.after=0, cell.padding=c(0,0,0,0))
# circos.initialize(speciesLabels, xlim = c(0, 1))
# 
# # #Plot the family
# # circos.track(ylim = c(0, 1), panel.fun = function(x,y) {
# #     i=CELL_META$sector.numeric.index
# #     #print(i)
# #     if(i!=1&i!=length(groupLabels)){
# #       before=groupLabels[i-1]
# #       now=groupLabels[i]
# #       after=groupLabels[i+1]
# #     } else {
# #       if(i==1){
# #         before=groupLabels[length(groupLabels)]
# #         now=groupLabels[i]
# #         after=groupLabels[i+1]
# #       } else{
# #         before=groupLabels[i-1]
# #         now=groupLabels[i]
# #         after=groupLabels[1]
# #       }
# #     }
# #     if(before==now & after==now){
# #       circos.rect(0, 0, 1, 1, col="lightgray", border="lightgray")
# #     }
# #     else{
# #       if(before!=now){
# #         circos.rect(0.1, 0, 1, 1, col="lightgray", border="lightgray")
# #       } else{
# #         circos.rect(0, 0, 0.9, 1, col="lightgray", border="lightgray")
# #       }
# #     }
# # }, bg.border = NA, track.height = 0.1)
# # 
# 
# #Plot the pictogram for diet
# 
# library(png)
# fruit <- readPNG("Metaanalysis/fruit3.png")
# leaf <- readPNG("Metaanalysis/leaf.png")
# dim(fruit)
# dim(leaf)
# 
# fruit.raster <- as.raster(fruit)
# leaf.raster <- as.raster(leaf)
# 
# circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) { 
#   print(CELL_META$sector.numeric.index)
#   i=CELL_META$sector.numeric.index
#   if(dietSpecies[i]=="Leaves"){
#     circos.raster(leaf.raster, CELL_META$xcenter, CELL_META$ycenter, 
#                   width = "0.25cm",#CELL_META$xrange,
#                   #height = dim(leaf)[1]/dim(leaf)[2]*CELL_META$xrange,
#                   facing = "inside")
#     # circos.raster(leaf.raster, CELL_META$xcenter, CELL_META$ycenter, 
#     #               width = CELL_META$xrange,
#     #               height = dim(leaf)[1]/dim(leaf)[2]*CELL_META$xrange,
#     #               facing = "bending.inside")
#     #addImg(leaf,  CELL_META$xcenter, CELL_META$ycenter, width = dim(leaf)[2]/dim(leaf)[1]*CELL_META$xrange)
#     #plot(leaf.raster)
#   } else{
#     #print("fruit")
#     circos.raster(fruit.raster, CELL_META$xcenter, CELL_META$ycenter, 
#                   width = "0.25cm",#CELL_META$xrange,
#                   #height = dim(leaf)[1]/dim(leaf)[2]*CELL_META$xrange,
#                   facing = "inside")
#     #circos.raster(fruit.raster, CELL_META$xcenter, CELL_META$ycenter, 
#     # width = CELL_META$xrange, 
#     # height = dim(fruit)[1]/dim(fruit)[2]*CELL_META$xrange, 
#     # facing = "bending.inside")
#   }
# }, track.height = 0.1)
# 
# 
# #Species name
# circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1, 
#              panel.fun = function(x, y) {
#                i=CELL_META$sector.numeric.index
#                circos.text(i-0.5, 0, labels.rdc[i], adj = c(0, 0), 
#                            facing = "clockwise", niceFacing = TRUE,
#                            col = "black", cex = 0.5, font=3)
#              })
# 
# #Plot the relative brain size
# circos.track(ylim = c(-4, 4), bg.border = NA, panel.fun = function(x, y) {
#   i=CELL_META$sector.numeric.index
#   circos.barplot(relativeValueBrain[i], 0, col=ifelse(relativeValueBrain[i] > 0, positiveCol, negativeCol), 
#                  border=ifelse(relativeValueBrain[i] > 0, positiveCol, negativeCol))
# }, track.height = 0.3)
# 
# library(plotrix)
# draw.circle(x=0,y=0,0.725, lty=2)
# 
# 
# #Plot the geographic links
# for(i in 1:length(speciesLabels)){
#   #locI <- which(strsplit(locationSpecies[i], "")==1)
#   for(j in 1:length(speciesLabels)){
#     #locJ <- which(strsplit(locationSpecies[j], "")==1)
#     product <- as.numcharac(unlist(strsplit(locationSpecies[j], "")))*as.numcharac(unlist(strsplit(locationSpecies[i], "")))
#     if(i==j|(length(unique(product))==1&unique(product)[1]==0)){
#       #Do nothing
#     }
#     else{
#       if(dietSpecies[i]==dietSpecies[j]){
#         colour <- as.data.frame(table(colLoc[which(product==1)]))
#         colour <- colour[colour$Freq==max(colour$Freq),1][1]
#         circos.link(speciesLabels[i], runif(1, 0, 1), speciesLabels[j], runif(1, 0, 1), lwd=1, col=colour)
#       }
#       else{
#         circos.link(speciesLabels[i], runif(1, 0, 1), speciesLabels[j], runif(1, 0, 1), lwd=1, col="lightgray")
#       }
#     }
#   }
# }
# circos.clear()
# 
# 
# #Plot the phylogenetic tree in a new circular plot
# ct = cutree(hc, numberOrder)  # cut tree into order
# n = length(labels)  # number of bird species
# dend = as.dendrogram(hc)
# 
# library(circlize)
# par(new = TRUE) # <- magic
# circos.par("canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-2, 2))
# circos.initialize("a", xlim = c(0, n)) # only one sector
# # circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.3, 
# #              panel.fun = function(x, y) {
# #                for(i in seq_len(n)) {
# #                  circos.text(i-0.5, 0, labels.rdc[i], adj = c(0, 0.5), 
# #                              facing = "clockwise", niceFacing = TRUE,
# #                              col = "black", cex = 0.2, font=3)
# #                }
# #              })
# 
# 
# #suppressPackageStartupMessages(library(dendextend))
# #dend = color_branches(dend, k = 6, col = 1:6)
# dend_height = attr(dend, "height")
# circos.track(ylim = c(0, dend_height), bg.border = NA, 
#              track.height = 0.5, panel.fun = function(x, y) {
#                circos.dendrogram(dend, col="lightgray")
#              })
# circos.clear()
# 
# ##----------------------
# 
# 
# 
# 
# ###
# ## Fig brain
# ###
# 
# #from: https://www.r-bloggers.com/2018/10/how-to-highlight-3d-brain-regions/
# 
# library(rgl)
# library(misc3d)
# library(neurobase)
# if (!requireNamespace("aal")) {
#   devtools::install_github("muschellij2/aal")
# } else {
#   library(aal)
# }
# if (!requireNamespace("MNITemplate")) {
#   devtools::install_github("jfortin1/MNITemplate")
# } else {
#   library(MNITemplate)
# }
# img = aal_image()
# template = readMNI(res = "2mm")
# cut <- 4500
# dtemp <- dim(template)
# # All of the sections you can label
# labs = aal_get_labels()
# 
# 
# # Pick the region of the brain you would like to highlight - in this case the hippocamus_L
# hippocampus = labs$index[grep("Hippocampus", labs$name)]
# cerebellum = labs$index[grep("Cerebellum", labs$name)]
# striatum = c(labs$index[grep("Caudate", labs$name)], labs$index[grep("Putamen", labs$name)], labs$index[grep("Pallidum", labs$name)])
# MOB = labs$index[grep("Olfactory_L", labs$name)]
# #neocortex = c(labs$index[grep("Frontal", labs$name)], labs$index[grep("Parietal", labs$name)], labs$index[grep("Temporal", labs$name)],
# #labs$index[grep("Occipital", labs$name)], labs$index[grep("Cingulate", labs$name)], labs$index[grep("Insula", labs$name)]) #incomplete
# 
# 
# library(RColorBrewer)
# colourVector <- brewer.pal(n = 5, name = "Set1")
# 
# colourHip=colourVector[1]
# colourCereb=colourVector[2]
# colourOlf=colourVector[3]
# colourStri=colourVector[4]
# 
# ### this would be the ``activation'' or surface you want to render 
# 
# #Set up the view manually the first time to extract coords
# #um <- par3d()$userMatrix
# #I obtained those as a correct view
# um <- matrix(c(-0.6391721,-0.76426947,-0.08573965,0,-0.2281326,0.08195215,0.97017491,0,-0.7344485,0.63966882,-0.22673632,0,0.0000000,0.00000000,0.00000000,1), ncol=4, nrow=4, byrow=TRUE)
# view3d(userMatrix = um)
# 
# contour3d(template, x=1:dtemp[1], y=1:dtemp[2], z=1:dtemp[3], level = cut, alpha = 0.03, draw = TRUE)
# 
# mask = remake_img(vec = img %in% hippocampus, img = img)
# contour3d(mask, level = c(0.5), alpha = c(0.4), add = TRUE, color=colourHip)
# mask = remake_img(vec = img %in% cerebellum, img = img)
# contour3d(mask, level = c(0.5), alpha = c(0.4), add = TRUE, color=colourCereb)
# mask = remake_img(vec = img %in% MOB, img = img)
# contour3d(mask, level = c(0.5), alpha = c(0.4), add = TRUE, color=colourOlf)
# #mask = remake_img(vec = img %in% neocortex, img = img)
# #contour3d(mask, level = c(0.5), alpha = c(0.4), add = TRUE, color=c("orange"))
# mask = remake_img(vec = img %in% striatum, img = img)
# contour3d(mask, level = c(0.5), alpha = c(0.4), add = TRUE, color=colourStri)
# 
# 
# ### add text
# #text3d(x=dtemp[1]/2, y=dtemp[2]/2, z = dtemp[3]*1.1, text="Top")
# #text3d(x=-0.98, y=dtemp[2]/2, z = dtemp[3]/2, text="Right")
# #text3d(x=dtemp[1]/2, y=dtemp[2], z = dtemp[3]/2, text="Front")
# #text3d(x=dtemp[1]/2, y=-1, z = dtemp[3]/2, text="Back")
# 
# #segmentBottomTop
# arrow3d(p0=c(dtemp[1]/2,dtemp[2]/2,dtemp[3]), p1=c(dtemp[1]/2,dtemp[2]/2,1.17*dtemp[3]), type="extrusion")
# #segmentBackFront
# arrow3d(p0=c(dtemp[1]/2,dtemp[2],dtemp[3]/2), p1=c(dtemp[1]/2,1.1*dtemp[2],dtemp[3]/2), type="extrusion")
# 
# legend3d(x=0.1, y=0.2, legend = c("Cerebellum", "Hippocampus", "Olfactory bulb", "Striatum"), fill = c(colourCereb, colourHip, colourOlf, colourStri), cex=1, bty="n", ncol=2)
# rglwidget()
# 
# #from: https://r-graphics.org/recipe-miscgraph-3d-save
# #rgl.snapshot('3dplot.png', fmt = 'png')
# rgl.postscript('3dplot.pdf', fmt = 'pdf')
# 
# 
# 
# ##----------------------
# 
