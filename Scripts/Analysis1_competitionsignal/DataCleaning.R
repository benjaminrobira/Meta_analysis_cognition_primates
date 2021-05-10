####---------------------------------
#### EvolutionaryHistoryAnalysis
####---------------------------------

# This script allows pre-processing necessary data for conducting the meta analysis on primate brain size.
# It will merge, clean, and adjust data relative to 1) primate biogeography (from the IUCN) and 2) data relative to brain size and other attributes (e.g. diet, body mass)

# It is divided into several parts (to be completed):

#1) Ranging data
#2) Brain data

###Set working directory
setwd("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates")

#Import environment
rm(list=ls())

##--------
#Home made functions
#To source all phylogenetics functions (biogeobears + models of evolution)
source("~/PhD/Meta_analysis/Cognition_metaanalysis/Functions.R")

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


####------------------
## ANALYSIS
####------------------

##////
# Data preparation
##////

#####
## Primate ranging
####
areaName <- c("Madagascar_east", "Madagascar_west", "Africa_west", "Africa_central", "Africa_south_east", 
              "America_central", "America_south_north", "America_south_south", 
              "Asia_west", "Asia_central_east",  "Asia_south", "Asia_island")
 

colourArea <- brewer.pal(n = 12, name="Paired")

#Select only the terrestrial area, because was roughly drawn with google earth
#+Note the shapefile were transformed from kml using: https://products.aspose.app/gis/conversion/kml-to-shapefile
 
worldMap <- getMap()
worldMap = clgeo_Clean(worldMap)

#Remove antartica because it causes problems
for(i in 1:nrow(worldMap)){
  if(!is.na(worldMap[i,]@data$continent[1])&worldMap[i,]@data$continent[1]=="Antarctica"){
    print(i)
  }
}
worldMap <- worldMap[-7,]

#Add the long/lat correct proj
worldMap <- spTransform(worldMap, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
#See if mercator proj working this time
worldMap_mercator <- spTransform(worldMap,CRS=CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))


##Import geographic areas: uncroped for now (because it seems that some islands are not well represented of the global map and are absent? I observed all absence for species on somes islands: near the cameroon, tanzania, or southeast indonesia. 
#I guess it is due to bad polygons from the world map
for(i in 1:length(areaName)){
  areaTransitory <- readOGR(dsn=paste("T:/IUCN_data_primate/Geographic_areas/Shapefiles/",areaName[i],".shp",sep=""))
  areaTransitory = clgeo_Clean(areaTransitory)
  areaTransitory <- spTransform(areaTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
  assign(paste("area", i, sep="_"), areaTransitory)
}

###
## Get geographic localisation for each species
primateSpeciesRange <- readOGR(dsn="T:/IUCN_data_primate/data_0.shp")
primateSpeciesRange <- spTransform(primateSpeciesRange, CRS(proj4string(worldMap)))

#Get the species
speciesGeo <- rep(NA, times=length(primateSpeciesRange))

for(i in 1:length(primateSpeciesRange)){
  if(!is.na(primateSpeciesRange[i,]@data$SUBSPECIES[1])){
    speciesGeo[i] <- paste(primateSpeciesRange[i,]@data$BINOMIAL[1], primateSpeciesRange[i,]@data$SUBSPECIES[1], sep="_")
  } else{
    speciesGeo[i] <- primateSpeciesRange[i,]@data$BINOMIAL[1]
  }
}
speciesGeo <- gsub(" ", "_",  speciesGeo)

#Extract name of phylogeny

treeForName <-read.nexus("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3.nex")
speciesLabels <-  as.data.frame(treeForName$tip.label)

#Determine the belonging geographic area based on overlap
thresholdPresenceRange <- seq(from=5, to=30, by=5)/100
matrixRangingSensitivity <- matrix(NA, nrow=nrow(speciesLabels), ncol=length(thresholdPresenceRange))
for(r in 1:length(thresholdPresenceRange)){#for all thresholds
  print(r)
  thresholdPresence=thresholdPresenceRange[r]
  matrixPrimateRange <- matrix(NA, nrow=nrow(speciesLabels), ncol=length(areaName)+1)
  for(i in 1:nrow(matrixPrimateRange)){#for all species for which phylogeny is known
    matrixPrimateRange[i,1] <- speciesLabels[i,1]
    lengthString <- str_length(speciesLabels[i,1])
    substrSpeciesGeo <- sapply(speciesGeo, function (x){substr(x, 1, lengthString)})#I do it in case there are subsubspecies etc while the phylogenetic tree does not got into such details
    polygonSpeciesId <- which(substrSpeciesGeo==speciesLabels[i,1])
    
    if(!is.na(polygonSpeciesId[1])){#if data on ranging are available
      speciesRangeTransitory <- primateSpeciesRange[polygonSpeciesId,]
      speciesRangeTransitory <- clgeo_Clean(speciesRangeTransitory)
      speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
      speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
    
      #Unit polygons: you can use:
      # unionSpatialPolygons from maptools: unionSpatialPolygons(pols,rep(1, length(pols))))
      # gUnaryUnion from rgeos:  gUnaryUnion(pols, id = pols@data$NAME_0))
      # gBuffer from rgeos: gBuffer(pols, byid=F, width=0)) -> the fastest https://gis.stackexchange.com/questions/278818/fastest-way-to-union-a-set-of-polygons-in-r
      #st_nuion from st:
      #pols_sf <- st_as_sf(pols)
      #st_union(pols_sf))
      

      speciesRangeTransitory <- gBuffer(speciesRangeTransitory, byid=F, width=0)
      #speciesRangeTransitory <- spTransform( test , CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj

      for(j in 1:length(areaName)){
        #transform biogeographic area to mercator
        areaTransitory <- get(paste("area", j, sep="_"))
        areaTransitory <- clgeo_Clean(areaTransitory)
        areaTransitory <- spTransform(areaTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
        areaTransitory <- spTransform(areaTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        areaTransitory <- clgeo_Clean(areaTransitory)
        
        #transform species area to mercator
        #speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        
        #compute intersection
        speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        intersectingArea <- gIntersection(speciesRangeTransitory, areaTransitory, byid=FALSE)

        #back transform to long/lat
        areaTransitory <- spTransform(areaTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
        speciesRangeTransitory <- spTransform(speciesRangeTransitory,   CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

        if(!is.null(intersectingArea)){
          intersectingArea <- gBuffer(intersectingArea, byid=F, width=0)
          intersectingArea  <- clgeo_Clean(intersectingArea)
          intersectingArea <- spTransform(intersectingArea, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
          overlap <- geosphere::areaPolygon(intersectingArea)/geosphere::areaPolygon(speciesRangeTransitory)
        } else{
          overlap <- 0
        }
        if(length(overlap) > 1){print(c(i,j))}
        if(overlap >= thresholdPresence){
          matrixPrimateRange[i, j+1] <- 1 
        } else{
          matrixPrimateRange[i, j+1] <- 0
        }
      }
      matrixRangingSensitivity[i,r] <- paste(matrixPrimateRange[i, 2:ncol(matrixPrimateRange)], collapse="")
    } else{
      matrixRangingSensitivity[i,r] <- NA
      }
      progress(i/nrow(matrixPrimateRange)*100)
  }
}

check <- apply(matrixPrimateRange[,2:ncol(matrixPrimateRange)], 1, function(x) sort(unique(x))[length(unique(x))])
which(check==0)
matrixPrimateRange[which(check==0),]
#remove values for Macaca sylvanus (Macaque de Barbarie)
matrixRangingSensitivity[matrixPrimateRange[,1]=="Macaca_sylvanus",] <- NA#Because still no species name for the summary matrix


ref <- matrixRangingSensitivity[,3]
howManyDifferent <- rep(NA, times=ncol(matrixRangingSensitivity))
for(i in 1:ncol(matrixRangingSensitivity)){
  howManyDifferent[i] <- length(which(matrixRangingSensitivity[,i]!=ref&!is.na(matrixRangingSensitivity[,i])))#matrixRangingSensitivity[,i]!=paste(rep(NA, times=12), collapse="")))
}
howManyDifferent <- howManyDifferent/length(!is.na(ref))#ref[ref!=paste(rep(NA, times=12), collapse="")])

nrow(matrixRangingSensitivity[!is.na(matrixRangingSensitivity[,1]),])

#Add species name
matrixRangingSensitivity <- as.data.frame(matrixRangingSensitivity)
matrixRangingSensitivity$SpeciesForPhylogeny <- speciesLabels[,1]

########
## Check which species distrib is NA and correct it using the IUCN red list (due to the fact that the the tree has too much precision

sort(matrixRangingSensitivity$SpeciesForPhylogeny[is.na(matrixRangingSensitivity[,1])])
#Those four are problem of too detailed info

#Eulemur_fulvus#Madagascar west + north/east "11000000000" 

speciesRangeTransitory <- primateSpeciesRange[which(sapply(as.vector(speciesGeo), function(x) substr(x, 1, 14))=="Eulemur_fulvus"),]
speciesRangeTransitory <- clgeo_Clean(speciesRangeTransitory)
speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
speciesRangeTransitory <- gBuffer(speciesRangeTransitory, byid=F, width=0)

#east
areaTransitory <- get(paste("area", 1, sep="_"))
areaTransitory <- clgeo_Clean(areaTransitory)
areaTransitory <- spTransform(areaTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
areaTransitory <- spTransform(areaTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
areaTransitory <- clgeo_Clean(areaTransitory)

#compute intersection
speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
intersectingArea <- gIntersection(speciesRangeTransitory, areaTransitory, byid=FALSE)

#back transform to long/lat
areaTransitory <- spTransform(areaTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
speciesRangeTransitory <- spTransform(speciesRangeTransitory,   CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
intersectingArea <- gBuffer(intersectingArea, byid=F, width=0) 
intersectingArea  <- clgeo_Clean(intersectingArea)
intersectingArea <- spTransform(intersectingArea, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
overlap <- geosphere::areaPolygon(intersectingArea)/geosphere::areaPolygon(speciesRangeTransitory)
overlap

#Eulemur_macaco#Madagascar north/east "10000000000"

speciesRangeTransitory <- primateSpeciesRange[which(sapply(as.vector(speciesGeo), function(x) substr(x, 1, 14))=="Eulemur_macaco"),]
speciesRangeTransitory <- clgeo_Clean(speciesRangeTransitory)
speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
speciesRangeTransitory <- gBuffer(speciesRangeTransitory, byid=F, width=0)

#east
areaTransitory <- get(paste("area", 1, sep="_"))
areaTransitory <- clgeo_Clean(areaTransitory)
areaTransitory <- spTransform(areaTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
areaTransitory <- spTransform(areaTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
areaTransitory <- clgeo_Clean(areaTransitory)

#compute intersection
speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
intersectingArea <- gIntersection(speciesRangeTransitory, areaTransitory, byid=FALSE)

#back transform to long/lat
areaTransitory <- spTransform(areaTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
speciesRangeTransitory <- spTransform(speciesRangeTransitory,   CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
intersectingArea <- gBuffer(intersectingArea, byid=F, width=0) 
intersectingArea  <- clgeo_Clean(intersectingArea) 
intersectingArea <- spTransform(intersectingArea, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
overlap <- geosphere::areaPolygon(intersectingArea)/geosphere::areaPolygon(speciesRangeTransitory)
overlap

Hapalemur_griseus#Madagascar west&north/east "11000000000"
##Check for it whether the two are always overlaped moren than 30%

which(apply(speciesLabels, 1, function(x) substr(x, 1, 17))=="Hapalemur_griseus")
# 54 55 56 57 58

speciesRangeTransitory <- primateSpeciesRange[which(sapply(as.vector(speciesGeo), function(x) substr(x, 1, 17))=="Hapalemur_griseus"),]
speciesRangeTransitory <- clgeo_Clean(speciesRangeTransitory)
speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
speciesRangeTransitory <- gBuffer(speciesRangeTransitory, byid=F, width=0)

#east
areaTransitory <- get(paste("area", 1, sep="_"))
areaTransitory <- clgeo_Clean(areaTransitory)
areaTransitory <- spTransform(areaTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
areaTransitory <- spTransform(areaTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
areaTransitory <- clgeo_Clean(areaTransitory)

#compute intersection
speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
intersectingArea <- gIntersection(speciesRangeTransitory, areaTransitory, byid=FALSE)

#back transform to long/lat
areaTransitory <- spTransform(areaTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
speciesRangeTransitory <- spTransform(speciesRangeTransitory,   CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
intersectingArea <- gBuffer(intersectingArea, byid=F, width=0) 
intersectingArea  <- clgeo_Clean(intersectingArea) 
intersectingArea <- spTransform(intersectingArea, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
overlap <- geosphere::areaPolygon(intersectingArea)/geosphere::areaPolygon(speciesRangeTransitory)
overlap

#Macaca_nemestrina#Asia island "000000000001"

#Then manually modify
vectorSpeciesToCorrect <- c("Eulemur_fulvus", "Eulemur_macaco", "Hapalemur_griseus", "Macaca_nemestrina")
vectorRangeToCorrect <- c("110000000000","100000000000","110000000000","0000000000001")

for(i in 1:length(vectorSpeciesToCorrect)){
  toChange <- which(substr(matrixRangingSensitivity$SpeciesForPhylogeny, 1, 14)==substr(vectorSpeciesToCorrect[i], 1, 14))
  matrixRangingSensitivity[toChange,1:(ncol(matrixRangingSensitivity)-1)] <- vectorRangeToCorrect[1]
}

###----
## Plotting the variability of biogeography classification in function of the spatial overlap threshold that is chosen

plot(
  x=0, y=0, xlab="Overlap threshold", ylab="Variation percent (relative to 15%)", 
  xlim=c(thresholdPresenceRange[1],thresholdPresenceRange[length(thresholdPresenceRange)]), ylim=c(0,0.4),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=5/100, xmax=30/100, xintsmall=2.5/100, xintbig=5/100, ymin=0, ymax=0.4, yintsmall=0.025, yintbig=0.1, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=1, tcl=-0.25)
axis(side=1, at=thresholdPresenceRange, labels=thresholdPresenceRange, las=1, tcl=-0.25)

points(x=thresholdPresenceRange, y=howManyDifferent, pch=19)
lines(x=thresholdPresenceRange, y=howManyDifferent)

## End of processing ranging data
##---

##----------------------------------------------------------------------------------------------------------------------------------------------------------------
##----------------------------------------------------------------------------------------------------------------------------------------------------------------

################################
## Working on traits + necessary attribute (diet, body mass etc...)
################################

##----------------
#Loading data
##----------------

options(warn=1)
#Powell
Data_powell <- read.delim("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Powell/Data_powell.txt")
Data_powell_mosaic <- read.delim("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Powell/Data_powell_mosaic.txt")
Data_powell2 <- read.delim("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Powell/Data_powell2.txt")
Data_powell2<- Data_powell2[Data_powell2$MSW05_Order=="Primates",]

#Todorov
Data_todorov <-  read.delim("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Todorov/data.txt")

#Willems
Data_willems <- read.csv("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Willems/willems_data.csv", sep=";")

#Decasien
Data_decasien_body <- read.delim("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Decasien/data_decasien1_body.txt")
Data_decasien_brain <- read.delim("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Decasien/data_decasien1_brain.txt")
Data_decasien_diet <- read.delim("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Decasien/data_decasien1_diet.txt")
Data_decasien_groupsize <- read.delim("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Decasien/data_decasien1_groupsize.txt")
Data_decasien_social <- read.delim2("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Decasien/data_decasien1_social.txt")
Data_decasien_mosaic <- read.delim2("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Decasien/Decasien_mosaic.txt")

#Kamilar
Data_kamilar1 <- read.csv("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Kamilar/kamilar_data1.csv", sep=";")
Data_kamilar2 <- read.csv("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Kamilar/kamilar_data2.csv", sep=";")
Data_kamilar3 <- read.csv("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Kamilar/kamilar_data3.csv", sep=";")
Data_kamilar4 <- read.csv("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Kamilar/kamilar_data4.csv", sep=";")
Data_kamilar5 <- read.csv("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Kamilar/kamilar_data5.csv", sep=";")
Data_kamilar6 <- read.csv("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Kamilar/kamilar_data6.csv", sep=";")
Data_kamilar7 <- read.csv("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Kamilar/kamilar_data7.csv", sep=";")
Data_kamilar8 <- read.csv("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Kamilar/kamilar_data8.csv", sep=";")
Data_kamilar9 <- read.csv("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Kamilar/kamilar_data9.csv", sep=";")
Data_kamilar10 <- read.csv("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Kamilar/kamilar_data10.csv", sep=";")

#Pearce
Data_pearce <- read.csv("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Pearce_data.csv", sep=";")

#Wrangham
Data_wrangham <- read.csv("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/wranghametal_data.csv", sep=";")

#Grueter
Data_grueter <- read.delim("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Grueter/grueter.txt")

#Navarrete
Data_navarrete <- read.delim("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Navarrete/Data_navarrete_copied_fromTable1_main.txt")
Data_navarrete <- aggregate(Data_navarrete[, c(4:ncol(Data_navarrete))], by=list(Data_navarrete$Species), FUN=mean)

#Error seen 
Data_grueter$Species <- as.character(Data_grueter$Species)
Data_grueter$Species[Data_grueter$Species=="paella"] <- "apella"

Data_grueter$CompleteName <- paste(Data_grueter[,1], Data_grueter[,2], sep="_")
Data_grueter$CompleteName <- gsub(" ", "_", Data_grueter$CompleteName)

#End loading
#........................................................................................

##----------------------------------------------------
# Creating a unique "species" ID between datasets
##----------------------------------------------------

#Uniformize name for everyone
dataset <-
c("Data_willems",
"Data_pearce",
"Data_wrangham",
"Data_grueter",
"Data_powell",
"Data_powell_mosaic",
"Data_powell2",
"Data_todorov",
"Data_decasien_mosaic",
"Data_decasien_diet",
"Data_decasien_body",
"Data_navarrete",
"Data_kamilar1",
"Data_kamilar2",
"Data_kamilar3",
"Data_kamilar4",
"Data_kamilar5",
"Data_kamilar6",
"Data_kamilar7",
"Data_kamilar8",
"Data_kamilar9",
"Data_kamilar10")

whereColSpecies <-
c(
1,
2,
2,
ncol(Data_grueter),
1,
1,
5,
1,
2,
1,
1,
1,
rep(1, times=10)
)

vectorSpecies <- c()

for(a in 1:length(dataset)){
  test <- get(dataset[a])
  colnames(test)[colnames(test)=="Species"] <- "Species_init"
  colnames(test)[whereColSpecies[a]] <- "Species"
  test$Species <- as.character(test$Species)
  test$Species <- gsub(" ", "_", test$Species)
  test$Species <- tolower(test$Species)
  test$Species <- firstup(test$Species)
  test <- separate(test, col="Species", into=c("Name1", "Name2"), remove=FALSE)
  test$Species_abbrv <- paste(test$Name1, substr(test$Name2, 1,4), sep="_")
  vectorSpecies <- c(vectorSpecies, as.character(test$Species_abbrv))
  
  #First run without the row below to check if duplicates still present
  assign(dataset[a], test)
}#warnings ok

vectorSpecies <- sort(unique(vectorSpecies))
length(vectorSpecies)
vectorSpecies

#End naming
#........................................................................................

##-------------------------------
# Creating the summary dataset
##-------------------------------

#Summary dataset: based on all species from all datasets
summaryData <- as.data.frame(vectorSpecies)
colnames(summaryData) <- "Species_abbrv"
summaryData$Family <- Data_powell2$MSW05_Family[match(summaryData$Species_abbrv,Data_powell2$Species_abbrv)]
summaryData <- summaryData[!is.na(summaryData[,1])&summaryData[,1]!="_NA",]

summaryData <- separate(summaryData, col=Species_abbrv, into=c("Part1", "Part2"), sep="_", remove=FALSE)
dfTransitory <- unique(summaryData[,c(2,4)])
dfTransitory <- dfTransitory[!is.na(dfTransitory[,2]),]
summaryData$Family <- dfTransitory$Family[match(summaryData$Part1,dfTransitory$Part1)]

summaryData <- summaryData[,c(1,4)]
summaryData[is.na(summaryData[,2]),]
#Don't care for missing ones now because might be deleted if missing data afterwards

#Matching correct name for phylogeny
#summaryData$SpeciesForPhylogeny <- Data_powell$Species.name.adjusted.to.10kTrees[match(summaryData$Species_abbrv, Data_powell$Species_abbrv)]


treeForName <-read.nexus("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3.nex")

speciesLabels <-  as.data.frame(treeForName$tip.label)
colnames(speciesLabels) <- "SpeciesTree"
speciesLabels <- speciesLabels [order(speciesLabels$SpeciesTree),]
speciesLabels <- as.data.frame(speciesLabels)
colnames(speciesLabels) <- "SpeciesTree"
speciesLabels <- separate(speciesLabels, col=SpeciesTree, into=c("Part1", "Part2"), sep="_", remove=FALSE)
speciesLabels$Species_abbrv <- paste(speciesLabels$Part1, substr(speciesLabels$Part2, 1, 4), sep="_")
summaryData$SpeciesForPhylogeny <- speciesLabels$SpeciesTree[match(summaryData$Species_abbrv, speciesLabels$Species_abbrv)]
#summaryData$SpeciesForPhylogeny <- Data_powell$Species.name.adjusted.to.10kTrees[match(summaryData$Species_abbrv, Data_powell$Species_abbrv)]


###
# This is transitory and is only used for filling covariates then species for which phylogeny is unknown will be discarded later
###

#if not available, add initial date. it will then be modified once you selected species with available data
toCorrectSpecies <- as.data.frame(summaryData$Species_abbrv[which(is.na(summaryData$SpeciesForPhylogeny))])
colnames(toCorrectSpecies) <- "SpeciesInit"
toCorrectSpecies <- separate(toCorrectSpecies, col=SpeciesInit, into=c("Item1", "Item2"), sep="_", remove=FALSE) 
toCorrectSpecies$newName <- paste(toCorrectSpecies$Item1, toCorrectSpecies$Item2, sep="_")

summaryData[,1] <- as.character(summaryData[,1])
summaryData[,2] <- as.character(summaryData[,2])
summaryData[,3] <- as.character(summaryData[,3])

summaryData$SpeciesForPhylogeny[which(is.na(summaryData$SpeciesForPhylogeny))] <- toCorrectSpecies$newName 
summaryData <- summaryData[,c(2,1,3)]
head(summaryData)
tail(summaryData)
#I will need in the end to recheck the name inserted here with those of the phylogenetic tree

#Identification of three pb: delete
summaryData[grep("fascicularis", summaryData$SpeciesForPhylogeny),]
"Macacafascicularis_NA"

summaryData[grep("fuscicollis", summaryData$SpeciesForPhylogeny),]
"Saguinusfuscicollis_NA"

summaryData[grep("fulvus", summaryData$SpeciesForPhylogeny),]
"Lemurfulvus_NA"

summaryData <- summaryData[summaryData$SpeciesForPhylogeny %nin% c("Macacafascicularis_NA", "Saguinusfuscicollis_NA", "Lemurfulvus_NA"),]

##------------------
##Adding the covariates
##------------------

###Brain data: Overall/Hyppocampus/Neocortex/Cerebellum/Insula/Hypothalamus
##Brain data are from Decasien_mosaic or powell2 or todorov

#Decasien
head(Data_decasien_mosaic)
str(Data_decasien_mosaic)
#First, need to aggregate

Subdata_decasien_mosaic <- Data_decasien_mosaic[, c(7,9,13,17,20,21,29,ncol(Data_decasien_mosaic))]
Subdata_decasien_mosaic$BV <- gsub(" " , "",as.character(Subdata_decasien_mosaic$BV))
meanData <- aggregate(as.numeric(as.character(Subdata_decasien_mosaic$BV)), by=list(Subdata_decasien_mosaic$Species_abbrv), FUN=mean, na.rm=TRUE)
summaryData$Brain_volume_mm3_decasien <- meanData$x[match(summaryData$Species_abbrv, meanData$Group.1)]
  
Subdata_decasien_mosaic$Hippocampus <- gsub(" " , "",as.character(Subdata_decasien_mosaic$Hippocampus))
meanData <- aggregate(as.numeric(as.character(Subdata_decasien_mosaic$Hippocampus)), by=list(Subdata_decasien_mosaic$Species_abbrv), FUN=mean, na.rm=TRUE)
summaryData$Hippocampus_volume_mm3_decasien <- meanData$x[match(summaryData$Species_abbrv, meanData$Group.1)]

Subdata_decasien_mosaic$Neocortex..GM.WM <- gsub(" " , "",as.character(Subdata_decasien_mosaic$Neocortex..GM.WM))
meanData <- aggregate(as.numeric(as.character(Subdata_decasien_mosaic$Neocortex..GM.WM)), by=list(Subdata_decasien_mosaic$Species_abbrv), FUN=mean, na.rm=TRUE)
summaryData$Neocortex_volume_mm3_decasien <- meanData$x[match(summaryData$Species_abbrv, meanData$Group.1)]

Subdata_decasien_mosaic$Cerebellum <- gsub(" " , "",as.character(Subdata_decasien_mosaic$Cerebellum))
meanData <- aggregate(as.numeric(as.character(Subdata_decasien_mosaic$Cerebellum)), by=list(Subdata_decasien_mosaic$Species_abbrv), FUN=mean, na.rm=TRUE)
summaryData$Cerebellum_volume_mm3_decasien <- meanData$x[match(summaryData$Species_abbrv, meanData$Group.1)]

Subdata_decasien_mosaic$Striatum <- gsub(" " , "",as.character(Subdata_decasien_mosaic$Striatum))
meanData <- aggregate(as.numeric(as.character(Subdata_decasien_mosaic$Striatum)), by=list(Subdata_decasien_mosaic$Species_abbrv), FUN=mean, na.rm=TRUE)
summaryData$Striatum_volume_mm3_decasien <- meanData$x[match(summaryData$Species_abbrv, meanData$Group.1)]

Subdata_decasien_mosaic$MOB <- gsub(" " , "",as.character(Subdata_decasien_mosaic$MOB))
meanData <- aggregate(as.numeric(as.character(Subdata_decasien_mosaic$MOB)), by=list(Subdata_decasien_mosaic$Species_abbrv), FUN=mean, na.rm=TRUE)
summaryData$MOB_volume_mm3_decasien <- meanData$x[match(summaryData$Species_abbrv, meanData$Group.1)]
#MOB=main olfactive bulb

Subdata_decasien_mosaic$Optic.Tract <- gsub(" " , "",as.character(Subdata_decasien_mosaic$Optic.Tract))
meanData <- aggregate(as.numeric(as.character(Subdata_decasien_mosaic$Optic.Tract)), by=list(Subdata_decasien_mosaic$Species_abbrv), FUN=mean, na.rm=TRUE)
summaryData$Optic.Tract_volume_mm3_decasien <- meanData$x[match(summaryData$Species_abbrv, meanData$Group.1)]

#Check how many data per area
length(unique(Data_decasien_mosaic$Species_abbrv[!is.na(Data_decasien_mosaic$Insula..GM.)&Data_decasien_mosaic$Insula..GM.!=""]))
length(unique(Data_decasien_mosaic$Species_abbrv[!is.na(Data_decasien_mosaic$Thalamus)&Data_decasien_mosaic$Thalamus!=""]))
length(unique(Data_decasien_mosaic$Species_abbrv[!is.na(Data_decasien_mosaic$Striatum)&Data_decasien_mosaic$Striatum!=""]))
length(unique(Data_decasien_mosaic$Species_abbrv[!is.na(Data_decasien_mosaic$Hypothalamus)&Data_decasien_mosaic$Hypothalamus!=""]))
length(unique(Data_decasien_mosaic$Species_abbrv[!is.na(Data_decasien_mosaic$Amygdala)&Data_decasien_mosaic$Amygdala!=""]))
length(unique(Data_decasien_mosaic$Species_abbrv[!is.na(Data_decasien_mosaic$Hippocampus)&Data_decasien_mosaic$Hippocampus!=""]))
length(unique(Data_decasien_mosaic$Species_abbrv[!is.na(Data_decasien_mosaic$Neocortex..GM.WM)&Data_decasien_mosaic$Neocortex..GM.WM!=""]))
length(unique(Data_decasien_mosaic$Species_abbrv[!is.na(Data_decasien_mosaic$Piriform.Lobe)&Data_decasien_mosaic$Piriform.Lobe!=""]))
length(unique(Data_decasien_mosaic$Species_abbrv[!is.na(Data_decasien_mosaic$Cerebellum)&Data_decasien_mosaic$Cerebellum!=""]))
length(unique(Data_decasien_mosaic$Species_abbrv[!is.na(Data_decasien_mosaic$Septum)&Data_decasien_mosaic$Septum!=""]))
length(unique(Data_decasien_mosaic$Species_abbrv[!is.na(Data_decasien_mosaic$MOB)&Data_decasien_mosaic$MOB!=""]))
length(unique(Data_decasien_mosaic$V1..GM.[!is.na(Data_decasien_mosaic$V1..GM.)&Data_decasien_mosaic$V1..GM.!=""]))
length(unique(Data_decasien_mosaic$LGN[!is.na(Data_decasien_mosaic$LGN)&Data_decasien_mosaic$LGN!=""]))
length(unique(Data_decasien_mosaic$Mesencephalon[!is.na(Data_decasien_mosaic$Mesencephalon)&Data_decasien_mosaic$Mesencephalon!=""]))
length(unique(Data_decasien_mosaic$Optic.Tract[!is.na(Data_decasien_mosaic$Optic.Tract)&Data_decasien_mosaic$Optic.Tract!=""]))

#Powell
summaryData$Brain_volume_mm3_powell_mosaic <- Data_powell_mosaic$Brain[match(summaryData$Species_abbrv,Data_powell_mosaic$Species_abbrv)]
summaryData$Neocortex_volume_mm3_powell_mosaic <- Data_powell_mosaic$Neocortex[match(summaryData$Species_abbrv,Data_powell_mosaic$Species_abbrv)]
summaryData$Cerebellum_volume_mm3_powell_mosaic <- Data_powell_mosaic$Cerebellum[match(summaryData$Species_abbrv,Data_powell_mosaic$Species_abbrv)]

#Todorov
summaryData$Brain_volume_mm3_todorov <- Data_todorov$Brainvol[match(summaryData$Species_abbrv,Data_todorov$Species_abbrv)]
summaryData$Hippocampus_volume_mm3_todorov <- Data_todorov$Hippocampus.total[match(summaryData$Species_abbrv,Data_todorov$Species_abbrv)]
summaryData$Neocortex_volume_mm3_todorov <- Data_todorov$Neo[match(summaryData$Species_abbrv,Data_todorov$Species_abbrv)]

#Navarrete
summaryData$Brain_volume_mm3_navarrete <- Data_navarrete$Totalbrain[match(summaryData$Species_abbrv,Data_navarrete$Species_abbrv)]
summaryData$Hippocampus_volume_mm3_navarrete <- Data_navarrete$Hippocampus[match(summaryData$Species_abbrv,Data_navarrete$Species_abbrv)]
summaryData$Neocortex_volume_mm3_navarrete <- Data_navarrete$Neocortex[match(summaryData$Species_abbrv,Data_navarrete$Species_abbrv)]
summaryData$Cerebellum_volume_mm3_navarrete <- Data_navarrete$Cerebellum[match(summaryData$Species_abbrv,Data_navarrete$Species_abbrv)]
summaryData$Striatum_volume_mm3_navarrete <- Data_navarrete$Striatum[match(summaryData$Species_abbrv,Data_navarrete$Species_abbrv)]

#Compare BV with grueter and powell, to see difference
summaryData$Brain_volume_mm3_grueter <- Data_grueter$Endocranial_volume_cm3[match(summaryData$Species_abbrv,Data_grueter$Species_abbrv)]*(10**3)
summaryData$Brain_volume_mm3_powell <- Data_powell$ECV[match(summaryData$Species_abbrv,Data_powell$Species_abbrv)]*(10**3)

#Add body mass
Data_pearce_rdcmasse <- Data_pearce[,c(19,15)]
Data_pearce_rdcmasse$Species_abbrv <- gsub(" ","_",Data_pearce_rdcmasse$Species_abbrv)
Data_pearce_rdcmasse <- unique(Data_pearce_rdcmasse)
Data_pearce_rdcmasse <- aggregate(Data_pearce_rdcmasse$mass, by=list(Data_pearce_rdcmasse[,1]), FUN=mean)
colnames(Data_pearce_rdcmasse) <- c("Species_abbrv", "masse")

summaryData$Body_mass_g_pearce <- Data_pearce_rdcmasse$masse[match(summaryData$Species_abbrv,Data_pearce_rdcmasse$Species_abbrv)]
summaryData$Body_mass_g_powell <- Data_powell$Body.mass[match(summaryData$Species_abbrv,Data_powell$Species_abbrv)]
summaryData$Body_mass_g_grueter <- Data_grueter$Body_mass_g[match(summaryData$Species_abbrv,Data_grueter$Species_abbrv)]
summaryData$Body_mass_g_decasien <- Data_decasien_body$Final.Body.Weight..g.[match(summaryData$Species_abbrv,Data_decasien_body$Species_abbrv)]
summaryData$Body_mass_g_decasien <- as.numeric(as.character(gsub(" " , "",as.character(summaryData$Body_mass_g_decasien))))

#Add diet
#Breadth:
unique(Data_powell2$X6.1_DietBreadth)
summaryData$Dietbreadth_powell2 <- Data_powell2$X6.1_DietBreadth[match(summaryData$Species_abbrv,Data_powell2$Species_abbrv)]
summaryData$Dietbreadth_powell2[summaryData$Dietbreadth_powell2==-999] <- NA

#%
summaryData$Diet_frug_powell <- Data_powell$X.fruit[match(summaryData$Species_abbrv,Data_powell$Species_abbrv)]
summaryData$Diet_leaves_powell <- Data_powell$X.leaves[match(summaryData$Species_abbrv,Data_powell$Species_abbrv)]
summaryData$Diet_youngleaves_powell <- Data_powell$X.young.leaves[match(summaryData$Species_abbrv,Data_powell$Species_abbrv)]
summaryData$Diet_matureleaves_powell <- Data_powell$X.mature.leaves[match(summaryData$Species_abbrv,Data_powell$Species_abbrv)]

summaryData$Diet_frug_decasien <- Data_decasien_diet$X..Fruit[match(summaryData$Species_abbrv,Data_decasien_diet$Species_abbrv)]
summaryData$Diet_frug_decasien[summaryData$Diet_frug_decasien=="n/a"] <- NA
summaryData$Diet_frug_decasien <- as.numeric(as.character(summaryData$Diet_frug_decasien))

summaryData$Diet_leaves_willems <- Data_willems$PropLeaves[match(summaryData$Species_abbrv,Data_willems$Species_abbrv)]*100

#Having an idea of what are the variables that can be considered
summarySampleVar <- apply(summaryData, 2, function(x){length(x[!is.na(x)])})
summarySampleVar

#End dataset creation
#........................................................................................
# -> use Hyppocampus (spatio-temporal) and Neocortex (information processing) and Cerebellum (spatio-temporal)
# -> use Striatum (motor action) and MOB (olfactory) as "null" area?


####
#Correlation analysis of covariates from different dataset: eyeing to data robustness
####

# Brain: 4 11 14 17 22 23 #do not use 11 bc powell x2
# Hippo: 5 15 18
# Neo 6 12 16 19
# Cereb 7 13 20
# Striatum 8 21
# Body 24 25 26 27
# Frug 29 33
# Fol 30 34
# grep("frug", names(summarySampleVar))
# names(summarySampleVar)[grep("Brain", names(summarySampleVar))]
# 
#I removed one column since then, so substract one

#Col to correlate #####◙A RENUMEROTER
colNumTest <- c(4, 4, 4, 4, 14, 14, 14, 17, 17, 22,
                5, 5, 15,
                6, 6, 6, 12, 12, 16,
                7, 7, 13,
                8, 
                24, 24, 24, 25, 25, 26,
                29,
                30
                ) 

colNumToCompare <- c(14, 17, 22, 23, 17, 22, 23, 22, 23, 23,
                     15, 18, 18,
                     12, 16, 19, 16, 19, 19,
                     13, 20, 20,
                     21,
                     25, 26, 27, 26, 27, 27,
                     33,
                     34)

cbind(colnames(summaryData[colNumTest]), colnames(summaryData[colNumToCompare]))
#Vectors to save results
barLower <- rep(NA, times=length(colNumTest))
barUpper <- rep(NA, times=length(colNumTest))
meanCoeff <- rep(NA, times=length(colNumTest))
N <- rep(NA, times=length(colNumTest))

for (i in 1:length(colNumTest)){
  test <- #abs(as.numeric(as.character(summaryData[,colNumTest[i]])) - as.numeric(as.character(summaryData[,colNumToCompare[i]])))
    cor.test(as.numeric(as.character(summaryData[,colNumTest[i]])), 
                   as.numeric(as.character( summaryData[,colNumToCompare[i]])), method="pearson")
  barLower[i] <- test$conf.int[1]
  barUpper[i] <- test$conf.int[2]
  meanCoeff[i] <- test$estimate[1]
  N[i] <- nrow(summaryData[!is.na(summaryData[,colNumTest[i]])&!is.na(summaryData[,colNumToCompare[i]]),])
}

#Plot

plot(
  x=0, y=0, xlab="", ylab="Coefficient of correlation", 
  xlim=c(0,length(meanCoeff)+1), ylim=c(0.6,1),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=0, xmax=length(meanCoeff), xintsmall=0.25, xintbig=1, ymin=0.6, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=2, by=0.2), labels=seq(from=0, to=2, by=0.2), las=2, tcl=-0.25)

#Comparison
whatCompared <- c(
  rep("Brain", times=10),
  "Hippocampus",
  "Hippocampus",
  "Hippocampus",
  "Neocortex",
  "Neocortex",
  "Neocortex",
  "Neocortex",
  "Neocortex",
  "Neocortex",
  "Cerebellum",
  "Cerebellum",
  "Cerebellum",
  "Striatum",
  rep("Body", times=6),
  "Frug.",
  "Fol."
)

#Plot legend of what is compared in coloured rectangles
whereToPlot <- as.data.frame(table(whatCompared))
whereToPlot$loc <- whereToPlot$Freq/2

colourWhatCompared <- c(brewer.pal(n = length(unique(whatCompared)) - 1, name = "Pastel1"), "darkgrey")

#Colour rectangle to indicate what is compared
refLoc=0
for (i in 1:length(whatCompared)){
rect(xleft=i-1,
     xright=i,
     ybottom=0.55,#-0.05,
     ytop=0.6,#0,
     border=NA,
     col=colourWhatCompared[which(unique(whatCompared)==whatCompared[i])],
     xpd=TRUE
)
errorBars(location=i-0.5, meanPt=meanCoeff[i], refUnit=1, col="black", minValue=0, maxValue=1, upperBarValue=barUpper[i], lowerBarValue=barLower[i], horiz=FALSE, symmetrical=FALSE)
points(x=i-0.5, y=meanCoeff[i], pch=19, col=colourWhatCompared[which(unique(whatCompared)==whatCompared[i])],
        xpd=TRUE)
text(x=i-0.5, y=barUpper[i]+0.015, labels=N[i], pch=19, col=colourWhatCompared[which(unique(whatCompared)==whatCompared[i])], cex=0.6,xpd=TRUE)

if(i==length(whatCompared)|whatCompared[i]!=whatCompared[i+1]){
  refLoc=refLoc+whereToPlot$loc[whereToPlot[,1]==whatCompared[i]]
  if(whereToPlot$loc[whereToPlot[,1]==whatCompared[i]]<=1){
    segments(x0=refLoc, x1=refLoc, y0=0.575,#-0.025,
             y1=0.525,#-0.075,
             col=colourWhatCompared[which(unique(whatCompared)==whatCompared[i])], xpd=TRUE)
    text(x=refLoc, y=0.5,#-0.1, 
         labels=whereToPlot[whereToPlot[,1]==whatCompared[i],1], col=colourWhatCompared[which(unique(whatCompared)==whatCompared[i])], cex=0.5, xpd=TRUE)
  }
  else{
    text(x=refLoc, y=0.575,#-0.025, 
         labels=whereToPlot[whereToPlot[,1]==whatCompared[i],1], col="black", cex=0.5, xpd=TRUE)
  }
  refLoc=refLoc+whereToPlot$loc[whereToPlot[,1]==whatCompared[i]]#add second time for having complete rectangle
}
}

cbind(colnames(summaryData[colNumTest]), colnames(summaryData[colNumToCompare]))
#Vectors to save results
barLower <- rep(NA, times=length(colNumTest))
barUpper <- rep(NA, times=length(colNumTest))
meanCoeff <- rep(NA, times=length(colNumTest))
N <- rep(NA, times=length(colNumTest))

for (i in 1:length(colNumTest)){
  transitoryinit <- as.data.frame(cbind(as.numeric(as.character(summaryData[,colNumTest[i]])),as.numeric(as.character(summaryData[,colNumToCompare[i]]))))
  transitoryinit <- transitoryinit[!is.na(transitoryinit[,1])&!is.na(transitoryinit[,2]),]
  transitory <- apply(transitoryinit, 1, function(v){abs(v[1]**2 - v[2]**2)/2/v[1]/v[2]}) #note= (abs((v1-v2))/v1 + abs((v2-v1))/v2)/2
  
  #When rate is 0 for both, gives NA, so to transform to 0
  transitory[is.na(transitory)] <- 0
  
  transitory[!is.finite(transitory)] <- apply(transitoryinit[!is.finite(transitory),], 1, max)
  
  barLower[i] <- min(transitory)
  barUpper[i] <- max(transitory)
  meanCoeff[i] <- mean(transitory)
  N[i] <- length(transitory)
}

ymax <- round((barUpper)/10)*10


plot(
  x=0, y=0, xlab="", ylab="Variability", 
  xlim=c(0,length(meanCoeff)+1), ylim=c(0,1),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=0, xmax=length(meanCoeff), xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)

#Comparison
whatCompared <- c(
  rep("Brain", times=10),
  "Hippocampus",
  "Hippocampus",
  "Hippocampus",
  "Neocortex",
  "Neocortex",
  "Neocortex",
  "Neocortex",
  "Neocortex",
  "Neocortex",
  "Cerebellum",
  "Cerebellum",
  "Cerebellum",
  "Striatum",
  rep("Body", times=6),
  "Frug.",
  "Fol."
)

#Plot legend of what is compared in coloured rectangles
whereToPlot <- as.data.frame(table(whatCompared))
whereToPlot$loc <- whereToPlot$Freq/2

colourWhatCompared <- c(brewer.pal(n = length(unique(whatCompared)) - 1, name = "Pastel1"), "darkgrey")

#Colour rectangle to indicate what is compared
refLoc=0
for (i in 1:length(whatCompared)){
  rect(xleft=i-1,
       xright=i,
       ybottom=-0.05,
       ytop=0,
       border=NA,
       col=colourWhatCompared[which(unique(whatCompared)==whatCompared[i])],
       xpd=TRUE
  )
  errorBars(location=i-0.5, meanPt=meanCoeff[i], refUnit=1, col="black", minValue=0, maxValue=1, upperBarValue=barUpper[i], lowerBarValue=barLower[i], horiz=FALSE, symmetrical=FALSE)
  points(x=i-0.5, y=meanCoeff[i], pch=19, col=colourWhatCompared[which(unique(whatCompared)==whatCompared[i])],
         xpd=TRUE)
  text(x=i-0.5, y=barUpper[i]+0.05, labels=N[i], pch=19, col=colourWhatCompared[which(unique(whatCompared)==whatCompared[i])], cex=0.6,xpd=TRUE)
  
  if(i==length(whatCompared)|whatCompared[i]!=whatCompared[i+1]){
    refLoc=refLoc+whereToPlot$loc[whereToPlot[,1]==whatCompared[i]]
    if(whereToPlot$loc[whereToPlot[,1]==whatCompared[i]]<=1){
      segments(x0=refLoc, x1=refLoc, y0=-0.025,
               y1=-0.075,
               col=colourWhatCompared[which(unique(whatCompared)==whatCompared[i])], xpd=TRUE)
      text(x=refLoc, y=-0.1, 
           labels=whereToPlot[whereToPlot[,1]==whatCompared[i],1], col=colourWhatCompared[which(unique(whatCompared)==whatCompared[i])], cex=0.5, xpd=TRUE)
    }
    else{
      text(x=refLoc, y=-0.025, 
           labels=whereToPlot[whereToPlot[,1]==whatCompared[i],1], col="black", cex=0.5, xpd=TRUE)
    }
    refLoc=refLoc+whereToPlot$loc[whereToPlot[,1]==whatCompared[i]]#add second time for having complete rectangle
  }
}

#End 
#........................................................................................

##-----
#Visual representation of the sample

#Adding geographical
summaryData$geographicCode <- matrixRangingSensitivity[match(summaryData$SpeciesForPhylogeny,matrixRangingSensitivity$SpeciesForPhylogeny),which(thresholdPresenceRange==min(thresholdPresenceRange))]   

#Adding diet as binary
summaryData$Guild_init_decasien <- Data_decasien_diet$Diet.Category[match(summaryData$Species_abbrv,Data_decasien_diet$Species_abbrv)]
summaryData$Guild_init_decasien[summaryData$Guild_init_decasien=="Fol"] <- "Leaf"
summaryData$Guild_init_decasien[summaryData$Guild_init_decasien=="Frug"] <- "Fruit"

dataForSample <- summaryData
labels.tordc <- as.data.frame(dataForSample$SpeciesForPhylogeny)
colnames(labels.tordc) <- "Name"
labels.tordc <- separate(labels.tordc, col="Name", into=c("Name1", "Name2", "Name3"), sep="_")

labels.rdc <- apply(labels.tordc, 1, function(x){
  if(!is.na(x[3])){
    paste(toupper(substr(x[1], 1, 1)), ". ", substr(x[2], 1, 3), ". ", substr(x[3], 1, 1), ".", sep="")
  } else{paste(toupper(substr(x[1], 1, 1)), ". ", substr(x[2], 1, 3), sep="")
  }
}
)

dataForSample$Species <- labels.rdc 
  
dataForSample$Bodymass <- NA
dataForSample$Brain <- NA
dataForSample$Neocortex <- NA
dataForSample$Hippocampus <- NA
dataForSample$Cerebellum <- NA
dataForSample$Striatum <- NA
dataForSample$MOB <- NA
#dataForSample$Optic_tract <- NA
dataForSample$Range <- NA
dataForSample$Frug <- NA
dataForSample$Fol <- NA

dataForSample <- dataForSample[order(dataForSample$SpeciesForPhylo, decreasing=TRUE),]

for(i in 1:nrow(dataForSample)){
v.transitory <- c(dataForSample$Body_mass_g_decasien[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
         dataForSample$Body_mass_g_powell[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
         dataForSample$Body_mass_g_pearce[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
         dataForSample$Body_mass_g_grueter[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]])

dataForSample$Bodymass[i] <- ifelse(length(unique(v.transitory))==1&is.na(unique(v.transitory)[1]), 0, 1)

v.transitory <- c(dataForSample$Brain_volume_mm3_decasien[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
                  dataForSample$Brain_volume_mm3_powell[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
                  dataForSample$Brain_volume_mm3_todorov[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
                  dataForSample$Brain_volume_mm3_grueter[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
                  dataForSample$Brain_volume_mm3_navarrete[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]])

dataForSample$Brain[i] <- ifelse(length(unique(v.transitory))==1&is.na(unique(v.transitory)[1]), 0, 1)

v.transitory <- c(dataForSample$Neocortex_volume_mm3_decasien[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
                  dataForSample$Neocortex_volume_mm3_powell_mosaic[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
                  dataForSample$Neocortex_volume_mm3_todorov[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
                  dataForSample$Neocortex_volume_mm3_navarrete[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]])

dataForSample$Neocortex[i] <- ifelse(length(unique(v.transitory))==1&is.na(unique(v.transitory)[1]), 0, 1)

v.transitory <- c(dataForSample$Hippocampus_volume_mm3_decasien[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
                  dataForSample$Hippocampus_volume_mm3_todorov[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
                  dataForSample$Hippocampus_volume_mm3_navarrete[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]])

dataForSample$Hippocampus[i] <- ifelse(length(unique(v.transitory))==1&is.na(unique(v.transitory)[1]), 0, 1)

v.transitory <- c(dataForSample$Cerebellum_volume_mm3_decasien[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
                  dataForSample$Cerebellum_volume_mm3_powell_mosaic[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
                  dataForSample$Cerebellum_volume_mm3_navarette[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]])

dataForSample$Cerebellum[i] <- ifelse(length(unique(v.transitory))==1&is.na(unique(v.transitory)[1]), 0, 1)

v.transitory <- c(dataForSample$Striatum_volume_mm3_decasien[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
                  dataForSample$Striatum_volume_mm3_navarrete[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]])

dataForSample$Striatum[i] <- ifelse(length(unique(v.transitory))==1&is.na(unique(v.transitory)[1]), 0, 1)

v.transitory <- c(dataForSample$MOB_volume_mm3_decasien[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]])

dataForSample$MOB[i] <- ifelse(length(unique(v.transitory))==1&is.na(unique(v.transitory)[1]), 0, 1)

#v.transitory <- c(dataForSample$Optic.Tract_volume_mm3_decasien[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]])

#dataForSample$Optic_tract[i] <- ifelse(length(unique(v.transitory))==1&is.na(unique(v.transitory)[1]), 0, 1)

v.transitory <- c(dataForSample$geographicCode)

dataForSample$Range[i] <- ifelse(length(unique(v.transitory))==1&is.na(unique(v.transitory)[1]), 0, 1)
if(dataForSample$SpeciesForPhylogeny[i] =="Macaca_sylvanus"){#Range available only not used here because isolated area
  dataForSample$Range[i] <- 1
}

v.transitory <- c(dataForSample$Diet_frug_powell[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
                  dataForSample$Diet_frug_decasien[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]])

dataForSample$Frug[i] <- ifelse(length(unique(v.transitory))==1&is.na(unique(v.transitory)[1]), 0, 1)

v.transitory <- c(dataForSample$Diet_leaves_powell[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]],
                  dataForSample$Diet_leaves_willems[dataForSample$Species_abbrv==dataForSample$Species_abbrv[i]])

dataForSample$Fol[i] <- ifelse(length(unique(v.transitory))==1&is.na(unique(v.transitory)[1]), 0, 1)

dataForSample$Diet_guild[i] <- ifelse(length(unique(dataForSample$Guild_init_decasien))==1&is.na(unique(dataForSample$Guild_init_decasien)[1]), 0, 1)

}

dataForSample <- dataForSample[,c(3, (ncol(dataForSample)-11):ncol(dataForSample))]
dataForSample <- dataForSample[order(dataForSample[,1], decreasing=TRUE),]

plot(
  x=0, y=0, xlab="", ylab="", cex.sub=1.6,
  xlim=c(-10,ncol(dataForSample)-1), ylim=c(0,nrow(dataForSample)),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

text(x=rep(-2, times=nrow(dataForSample)), y=1:nrow(dataForSample)-0.5, labels=dataForSample$SpeciesForPhylo, xpd=TRUE, cex=0.4)
text(x=rep(-0, times=nrow(dataForSample)), y=1:nrow(dataForSample)-0.5, labels=dataForSample$Species, xpd=TRUE, cex=0.4)
text(x=3:ncol(dataForSample)-1.5, y=rep(nrow(dataForSample)+3, times=length(3:ncol(dataForSample))), labels=colnames(dataForSample)[3:ncol(dataForSample)], xpd=TRUE, cex=0.4, srt=45)

for(i in 1:nrow(dataForSample)){
  for(j in 3:ncol(dataForSample)){
    if(!is.na(dataForSample[i,j])&dataForSample[i,j]==1){
      rect(
        xleft=j-2,
        xright=j-1,
        ybottom=i-1,
        ytop=i,
        border=NA,
        col="black"
      )
    } else{
      rect(
        xleft=j-2,
        xright=j-1,
        ybottom=i-1,
        ytop=i,
        border=NA,
        col="lightgrey"
      )
    }
  }
}

addGrid(xmin=1, xmax=ncol(dataForSample), xintsmall=1, xintbig=1, ymin=0, ymax=nrow(dataForSample), yintsmall=1, yintbig=1, colsmall="white", colbig="white", axisPlot=FALSE)

# segments(x0=c(2,8), x1=c(2,8), y0=c(0,0), y1=c(nrow(dataForSample),nrow(dataForSample)), col=brewer.pal(n = 5, name = "Set1")[2], xpd=TRUE)
# segments(x0=c(9,12), x1=c(9,12), y0=c(0,0), y1=c(nrow(dataForSample),nrow(dataForSample)), col=brewer.pal(n = 5, name = "Set1")[3], xpd=TRUE)

## END ANALYSIS
##----------------------

##--
## Saving environment
save.image("REnvironments/geography_traits.RData")