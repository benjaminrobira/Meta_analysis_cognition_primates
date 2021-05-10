# install.packages(c("RColorBrewer",
# "rworldmap",
# "cleangeo",
# "rgdal",
# "ape",
# "stringr",
# "svMisc",
# "rgeos",
# "tidyr",
# "RPANDA",
# "ape",
# "caper",
# "MCMCglmm",
# "phytools",
# "BioGeoBEARS",
# "geosphere"), dependencies=TRUE)


###
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
###


###
#Functions
#To source all phylogenetics functions (biogeobears + models of evolution)
source("~/PhD/Meta_analysis/Cognition_metaanalysis/Functions.R")

#My toolkit
source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")

#Capitalize first letter
firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}
###


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

Eulemur_fulvus#Madagascar west + north/east "11000000000"


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


Eulemur_macaco#Madagascar north/east "10000000000"


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

Macaca_nemestrina#Asia island "000000000001"

#Then manually modify
vectorSpeciesToCorrect <- c("Eulemur_fulvus", "Eulemur_macaco", "Hapalemur_griseus", "Macaca_nemestrina")
vectorRangeToCorrect <- c("110000000000","100000000000","110000000000","0000000000001")

for(i in 1:length(vectorSpeciesToCorrect)){
  toChange <- which(substr(matrixRangingSensitivity$SpeciesForPhylogeny, 1, 14)==substr(vectorSpeciesToCorrect[i], 1, 14))
  matrixRangingSensitivity[toChange,1:(ncol(matrixRangingSensitivity)-1)] <- vectorRangeToCorrect[1]
}

pdf(file="Figure/Sensitivity_geography.pdf")

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

dev.off()

save.image("geography.RData")

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
pdf(file="Figure/Sensitivity_data.pdf")

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

dev.off()

#End 
#........................................................................................

##
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

pdf(file="sample.pdf", height=30, width=10)

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

dev.off()
##----------------

save.image("geography_traits.RData")

##////
# Data analysis
##////


###
## Phylogeny analysis of the evolutionary history: looping to take into account all uncertainty
##

frugivoryThresholdVector <- seq(from=20, to=40, by=20)
folivoryThresholdVector <- seq(from=40, to=60, by=20)
geographicThresholdVector <- c(10,30)/100
randomSampling=10
numberSimulations=10
numberTrees=1

totModels=randomSampling*numberSimulations*numberTrees*length(frugivoryThresholdVector)*length(folivoryThresholdVector)*length(geographicThresholdVector)

QtransitionRateEst <- as.data.frame(matrix(NA, ncol=2, nrow=totModels))
#QtransitionRateEqual <- as.data.frame(matrix(NA, ncol=2, nrow=totModels))
#QtransitionRatePrior <- as.data.frame(matrix(NA, ncol=2, nrow=totModels))

progression=0

##----
##Run the biogeobears for speeding up after

##--
# THE TREE

#Phylogenetic tree: force it to be ultrametric, problem was minor because doesn't change branch length actually
# phylo_all <-read.nexus("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/TreeBlock_10kTrees_Primates_Version3.nex")
# phylo_init <- phylo_all[[e]]

phylo_all <-read.nexus("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3.nex")
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


write.tree(phylo, np("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/Tree_biogeobears.nex"))
phylo <- read.tree("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/Tree_biogeobears.nex")

phylo <- cleanTree(tr=phylo, trfn="Tree_biogeobears.nex")
phylo <- force.ultrametric(tree=phylo, method="extend")#method="nnls")

is.ultrametric(phylo)
is.binary(phylo)

write.tree(phylo, np("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/Tree_biogeobears.nex"))
phylo <- read.tree("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/Tree_biogeobears.nex")

#numstates_from_numareas(numareas=12, maxareas=3, include_null_range=TRUE)#check how big is the transition matrix to deal with

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
              file="~/PhD/Meta_analysis/Cognition_metaanalysis/test_geo.data", 
              sep="\t",col.names=FALSE, row.names=FALSE, na="", quote = FALSE)
  
  
  BSM_output_from_function <- runBioGeoBearsandBSM_DEC(
    pathToTree="~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/Tree_biogeobears.nex",
    pathToGeo="test_geo.data",
    maxCooccurrenceGeo=3,
    mapNumber=numberSimulations,
    coreNumber=9,
    pathToSaveRdata=paste("BioGeoBEARS/BioGeoBEARS_", c, ".Rdata", sep=""),
    pathToSaveFigureGeo=paste("~/PhD/Meta_analysis/Cognition_metaanalysis/Brain_DEC_v1", c, ".pdf", sep=""),
    pathToBSMMapList=paste("BioGeoBEARS/BSM_output_file", c, ".Rdata", sep=""),
    pathToBSMoutput=paste("BioGeoBEARS/BSM_output_file", c, ".Rdata", sep=""),
    pathToFileStockingPdfForBSMAnalysis="~/PhD/Meta_analysis/Cognition_metaanalysis/BioGeoBEARS/"
  )
  dev.off()
} 

save.image("geography_traits_biogeobears.RData")

load("geography_traits_biogeobears.RData")

# summaryBrainFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))
# summaryEQFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))
# summaryNeocortexFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))
# summaryHippocampusFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))
# summaryCerebellumFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))
# summaryStriatumFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))
# summaryMOBFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))
# summaryOpticTractFrugivory <- as.data.frame(matrix(NA, ncol=6, nrow=totModels))

# repetition=length(frugivoryThresholdVector)*length(folivoryThresholdVector)*length(geographicThresholdVector)*randomSampling
# checkSampleFruit <- rep(NA, times=repetition)
# checkSampleLeaf <- rep(NA, times=repetition)
# checkSampleRange <- rep(NA, times=repetition)
# checkSampleBrain <-  rep(NA, times=repetition)
# checkSampleEQ <-  rep(NA, times=repetition)
# checkSampleNeocortex <-  rep(NA, times=repetition)
# checkSampleHippocampus <- rep(NA, times=repetition)
# checkSampleCerebellum <- rep(NA, times=repetition)
# checkSampleStriatum <- rep(NA, times=repetition)
# checkSampleMOB <- rep(NA, times=repetition)
# #checkSampleOptic_tract <- rep(NA, times=)
# checkSampleRange <- rep(NA, times=repetition)

counter=0

##-----
## Run the different evolutionary models // RUN ON THE CLUSTER

#' for(a in 1:length(frugivoryThresholdVector)){
#'   
#'   frugivoryThreshold=frugivoryThresholdVector[a]
#'   
#'   for(b in 1:length(folivoryThresholdVector)){
#'     
#'     folivoryThreshold=folivoryThresholdVector[b] 
#'     
#'     for(c in 1:length(geographicThresholdVector)){
#'       
#'       ###------------------
#'       ##Adding the co-occurence
#'       ###------------------
#'       
#'       summaryData$geographicCode <- matrixRangingSensitivity[match(summaryData$SpeciesForPhylogeny,matrixRangingSensitivity$SpeciesForPhylogeny),which(thresholdPresenceRange==geographicThresholdVector[c])]   
#'       
#'       load(paste("BioGeoBEARS/BSM_output_file", c, ".Rdata", sep=""))
#'       
#'       ##--------
#'       #for(e in 1:numberTrees){
#'         
#'           #Parallelizing loop
#'       
#'           #for(d in 1:randomSampling){
#'           cores=detectCores()
#'           cl <- makeCluster(cores[1]-2) #not to overload your computer
#'           registerDoParallel(cl)
#'           foreach(d=1:randomSampling#, .packages=c('caper',
#'                   #'RPANDA',
#'                   #'BioGeoBEARS',
#'                   #'phytools',
#'                   #'ape',
#'                   #'phytools')
#'                   ) %dopar% {
#' 
#'             setwd("C:/Users/robira/Documents/PhD/Meta_analysis/Cognition_metaanalysis")
#'             number= (a-1)*randomSampling + (b-1)*randomSampling + (c-1)*randomSampling + d-1 #C'est faux ce comptage...
#'             
#'             #Source functions & packages
#'             
#'             ###
#'             #Functions
#'             #To source all phylogenetics functions (biogeobears + models of evolution)
#'             source("~/PhD/Meta_analysis/Cognition_metaanalysis/Functions.R")
#'             
#'             #My toolkit
#'             source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
#'             
#'             #Capitalize first letter
#'             firstup <- function(x) {
#'               substr(x, 1, 1) <- toupper(substr(x, 1, 1))
#'               x
#'             }
#'             ###
#'             
#'             # test <- matrix(runif(10, 0, 10), ncol=2, nrow=5)
#'             # write.table(test, paste("OutputEvolModel/Data", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
#'             # 
#'             # # #Phylogenetics
#'             # tryCatch(library(caper), error=function(e) {
#'               # install.packages("caper", dependencies=TRUE, repos="https://ftp.igh.cnrs.fr/pub/CRAN/")
#'               # library(caper)
#'               # })
#'             # tryCatch(library(MCMCglmm), error=function(e) {
#'               # install.packages("MCMCglmm", dependencies=TRUE, repos="https://ftp.igh.cnrs.fr/pub/CRAN/")
#'               # library(MCMCglmm)
#'             # })
#'             # tryCatch(library(RPANDA), error=function(e) {
#'               # install.packages("devtools", dependencies=TRUE, repos="https://ftp.igh.cnrs.fr/pub/CRAN/")
#'               # library(devtools)
#'               
#'               # install_github("hmorlon/PANDA", dependencies = TRUE)
#'               # library(RPANDA)
#'             # })
#'             # tryCatch(library(BioGeoBEARS), error=function(e) {
#'               
#'               # install.packages("optimx", dependencies=TRUE, repos="http://cran.rstudio.com")
#'               # library(optimx)
#'               
#'               # install.packages("snow", dependencies=TRUE, repos="https://ftp.igh.cnrs.fr/pub/CRAN/")
#'               # library(snow)
#'               
#'               
#'               # library(devtools)
#'                               
#'               # devtools::install_github(repo="nmatzke/BioGeoBEARS")
#'               # library(BioGeoBEARS)
#'             # })
#'             # tryCatch(library(phytools), error=function(e) {
#'               # install.packages("phytools", dependencies=TRUE, repos="https://ftp.igh.cnrs.fr/pub/CRAN/")
#'               # library(phytools)
#'             # })
#'             # tryCatch(library(ape), error=function(e) {
#'               # install.packages("ape", dependencies=TRUE, repos="https://ftp.igh.cnrs.fr/pub/CRAN/")
#'               # library(ape)
#'             # })
#'             # tryCatch(library(geiger), error=function(e) {
#'               # install.packages("geiger", dependencies=TRUE, repos="https://ftp.igh.cnrs.fr/pub/CRAN/")
#'               # library(geiger)
#'             # })  
#'       			#Phylogenetics
#'       			library(caper)
#'       			library(MCMCglmm)
#'       			library(RPANDA)
#'       			library(BioGeoBEARS)
#'       			library(phytools)
#'       			library(ape)
#'       			library(phytools)
#'       			library(geiger)
#'       			library(optimx)
#'             
#'             #test <- matrix(runif(10, 0, 10), ncol=2, nrow=5)
#'             #write.table(test, paste("OutputEvolModel/Data_packages", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
#'             
#'             #write.table(summaryData, "SummaryData.txt", row.names=FALSE, col.names=TRUE, sep="\t")
#'             
#'             #####---------------
#'             #Random sampling for covariate in case of multiple sources
#'             #####---------------
#'             summaryData$Family<- NA
#'             summaryData$DietaryGuild <- NA
#'             summaryData$FrugivoryPercent <- NA
#'             summaryData$FolivoryPercent <- NA
#'             summaryData$Brain <- NA
#'             summaryData$Neocortex <- NA 
#'             summaryData$Cerebellum <- NA
#'             summaryData$Hippocampus <- NA
#'             summaryData$Striatum <- NA
#'             summaryData$MOB <- NA 
#'             summaryData$Optic.Tract <- NA 
#'             summaryData$Bodymass <- NA
#'             
#'             for(i in 1:nrow(summaryData)){
#'               #Frugivory 
#'               value <- c(summaryData$Diet_frug_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                          summaryData$Diet_frug_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
#'               
#'               value <- value[!is.na(value)]
#'               #Because there are issues if of length 1
#'               if(length(value)==1){
#'                 value <- c(value, value)
#'               }
#'               
#'               if(length(value)>0){
#'                 summaryData$FrugivoryPercent[i] <- sample(value, 1)
#'               }
#'               else{
#'                 summaryData$FrugivoryPercent[i] <- NA
#'               }
#'               
#'               #Folivory 
#'               value <- c(summaryData$Diet_leaves_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                          summaryData$Diet_leaves_willems[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
#'               
#'               value <- value[!is.na(value)]
#'               #Because there are issues if of length 1
#'               if(length(value)==1){
#'                 value <- c(value, value)
#'               }
#'               
#'               if(length(value)>0){
#'                 summaryData$FolivoryPercent[i] <- sample(value, 1)
#'               }
#'               else{
#'                 summaryData$FolivoryPercent[i] <- NA
#'               }
#'               
#'               ##--------------
#'               #Associate the dietary guild
#'               if(!is.na(as.numeric(as.character(summaryData$FrugivoryPercent[i])))&as.numeric(as.character(summaryData$FrugivoryPercent[i]))>=frugivoryThreshold){
#'                 summaryData$DietaryGuild[i] <- "Fruit"
#'               } else if(!is.na(as.numeric(as.character(summaryData$FolivoryPercent[i])))&as.numeric(as.character(summaryData$FolivoryPercent[i]))>=folivoryThreshold){
#'                 summaryData$DietaryGuild[i] <- "Leaf"
#'               }
#'               else if(!is.na(summaryData$Guild_init_decasien[i])){
#'                 summaryData$DietaryGuild[i] <- summaryData$Guild_init_decasien[i]
#'                 if(summaryData$DietaryGuild[i] =="Om"|summaryData$DietaryGuild[i]=="Frug/Fol"){
#'                   summaryData$DietaryGuild[i] <- "Fruit"
#'                 }
#'                 }
#'               else {
#'                 summaryData$DietaryGuild[i] <- "Other"
#'               }
#'               ##--------------
#'             
#'               #Brain volume
#'               value <- c(
#'                 summaryData$Brain_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                 summaryData$Brain_volume_mm3_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                 summaryData$Brain_volume_mm3_todorov[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                 summaryData$Brain_volume_mm3_grueter[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                 summaryData$Brain_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
#'               
#'               value <- value[!is.na(value)]
#'               #Because there are issues if of length 1
#'               if(length(value)==1){
#'                 value <- c(value, value)
#'               }
#'               
#'               if(length(value)>0){
#'                 summaryData$Brain[i] <- sample(value, 1)
#'               }
#'               else{
#'                 summaryData$Brain[i] <- NA
#'               }
#'               
#'               #Striatum
#'               value <- c(
#'                 summaryData$Striatum_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                 summaryData$Striatum_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
#'               
#'               value <- value[!is.na(value)]
#'               #Because there are issues if of length 1
#'               if(length(value)==1){
#'                 value <- c(value, value)
#'               }
#'               
#'               if(length(value)>0){
#'                 summaryData$Striatum[i] <- sample(value, 1)
#'               }
#'               else{
#'                 summaryData$Striatum[i] <- NA
#'               }
#'               
#'               #Neocortex
#'               value <- c(summaryData$Neocortex_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                          summaryData$Neocortex_volume_mm3_powell_mosaic[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                          summaryData$Neocortex_volume_mm3_todorov[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                          summaryData$Neocortex_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
#'               
#'               value <- value[!is.na(value)]
#'               #Because there are issues if of length 1
#'               if(length(value)==1){
#'                 value <- c(value, value)
#'               }
#'               
#'               if(length(value)>0){
#'                 summaryData$Neocortex[i] <- sample(value, 1)
#'               }
#'               else{
#'                 summaryData$Neocortex[i] <- NA
#'               }
#'               
#'               #Cerebellum
#'               value <- c(summaryData$Cerebellum_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                          summaryData$Cerebellum_volume_mm3_powell_mosaic[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                          summaryData$Cerebellum_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
#'               
#'               value <- value[!is.na(value)]
#'               #Because there are issues if of length 1
#'               if(length(value)==1){
#'                 value <- c(value, value)
#'               }
#'               
#'               if(length(value)>0){
#'                 summaryData$Cerebellum[i] <- sample(value, 1)
#'               }
#'               else{
#'                 summaryData$Cerebellum[i] <- NA
#'               }
#'               
#'               #Hippocampus
#'               value <- c(summaryData$Hippocampus_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                          summaryData$Hippocampus_volume_mm3_todorov[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                          summaryData$Hippocampus_volume_mm3_navarrete[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
#'               
#'               value <- value[!is.na(value)]
#'               #Because there are issues if of length 1
#'               if(length(value)==1){
#'                 value <- c(value, value)
#'               }
#'               
#'               if(length(value)>0){
#'                 summaryData$Hippocampus[i] <- sample(value, 1)
#'               }
#'               else{
#'                 summaryData$Hippocampus[i] <- NA
#'               }
#'               
#'               #MOB
#'               value <- c(
#'                 summaryData$MOB_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
#'               
#'               value <- value[!is.na(value)]
#'               #Because there are issues if of length 1
#'               if(length(value)==1){
#'                 value <- c(value, value)
#'               }
#'               
#'               if(length(value)>0){
#'                 summaryData$MOB[i] <- sample(value, 1)
#'               }
#'               else{
#'                 summaryData$MOB[i] <- NA
#'               }
#'               
#'               # #Optic.Tract #Too low sample
#'               # value <- c(
#'               #   summaryData$Optic.Tract_volume_mm3_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
#'               # 
#'               # value <- value[!is.na(value)]
#'               # #Because there are issues if of length 1
#'               # if(length(value)==1){
#'               #   value <- c(value, value)
#'               # }
#'               # 
#'               # if(length(value)>0){
#'               #   summaryData$Optic.Tract[i] <- sample(value, 1)
#'               # }
#'               # else{
#'               #   summaryData$Optic.Tract[i] <- NA
#'               # }
#'               
#'               #Bodymass
#'               value <- c(summaryData$Body_mass_g_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                          summaryData$Body_mass_g_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                          summaryData$Body_mass_g_pearce[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
#'                          summaryData$Body_mass_g_grueter[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
#'               
#'               value <- value[!is.na(value)]
#'               #Because there are issues if of length 1
#'               if(length(value)==1){
#'                 value <- c(value, value)
#'               }
#'               
#'               if(length(value)>0){
#'                 summaryData$Bodymass[i] <- sample(value, 1)
#'               }
#'               else{
#'                 summaryData$Bodymass[i] <- NA
#'               }
#'               
#'             }
#'             
#'             summaryData$Family<- summaryData$Family[match(summaryData$Species_abbrv,summaryData$Species_abbrv)]
#'     
#'             #write.table(summaryData, paste("OutputEvolModel/Data", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
#'             #save data for plotting
#'             
#'             if(a==1&b==1&c==1&d==1){
#'               summaryDataForPlot <- summaryData
#'               write.table(summaryData, "OutputEvolModel/Dataplot.txt", row.names=FALSE, col.names=TRUE, sep="\t")
#'             }
#'             #Sample size
#'             counter=counter+1
#'             checkSampleFruit[counter] <- length(summaryData$DietaryGuild[summaryData$DietaryGuild=="Fruit"])
#'             checkSampleLeaf[counter]  <- length(summaryData$DietaryGuild[summaryData$DietaryGuild=="Leaf"])
#'             checkSampleRange[counter]  <- length(summaryData$geographicCode[!is.na(summaryData$geographicCode)])
#'             checkSampleBrain[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&summaryData$DietaryGuild=="Fruit",])
#'             checkSampleEQ[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Bodymass)&summaryData$DietaryGuild=="Fruit",])
#'             checkSampleNeocortex[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Neocortex)&summaryData$DietaryGuild=="Fruit",]) 
#'             checkSampleHippocampus[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Hippocampus)&summaryData$DietaryGuild=="Fruit",])
#'             checkSampleCerebellum[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Cerebellum)&summaryData$DietaryGuild=="Fruit",])
#'             checkSampleStriatum[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Striatum)&summaryData$DietaryGuild=="Fruit",])
#'             checkSampleMOB[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$MOB)&summaryData$DietaryGuild=="Fruit",])
#'             #checkSampleOptic_tract[counter]  <- nrow(summaryData[!is.na(summaryData$Brain)&!is.na(summaryData$Optic.Tract)&summaryData$DietaryGuild=="Fruit",])#too low sample
#' 
#'             ##--------
#'             # Evolutionary history of diet
#'             ##--------
#'             
#'             #library(RPANDA)
#'             write.table(summaryData, paste("OutputEvolModel/Data", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
#'             
#'             vectorDiet <-  summaryData$DietaryGuild
#'             names(vectorDiet) <-  summaryData$SpeciesForPhylogeny
#'             vectorDiet <- vectorDiet[vectorDiet!="Other"&!is.na(summaryData$geographicCode)]
#'             #Load and save tree corresponding to species with diet
#'             options(warn=1)
#' 
#'             #Phylogenetic tree: force it to be ultrametric, problem was minor because doesn't change branch length actually
#'             # phylo_all <-read.nexus("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/TreeBlock_10kTrees_Primates_Version3.nex")
#'             # phylo_init <- phylo_all[[e]]
#'             # 
#'             
#'             #run simmap
#'             # phylo <- drop.tip(phylo,
#'             #                   phylo$tip.label[
#'             #                     which(phylo$tip.label
#'             #                           %nin%names(vectorDiet))]) 
#'             
#'             
#'             #write.tree(phylo, np("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/Tree_diet.nex"))
#'             phylo <- read.tree("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/Tree_biogeobears.nex")
#'             
#'             phylo <- drop.tip(phylo,
#'                               phylo$tip.label[
#'                                 which(phylo$tip.label
#'                                       %nin%names(vectorDiet))]) 
#'             simmapdiet1 <- make.simmap(tree=phylo, vectorDiet, model="ARD", pi="estimated", nsim=numberSimulations)#inequal and not symmetrical rate of transition from folivory to frugivory etc...
#'             #simmapdiet2 <- make.simmap(tree=phylo, vectorDiet, model="ARD", pi="equal", nsim=numberSimulations)#inequal and not symmetrical rate of transition from folivory to frugivory etc...
#'             #simmapdiet3 <- make.simmap(tree=phylo, vectorDiet, model="ARD", pi=c(0,1), nsim=numberSimulations)#inequal and not symmetrical rate of transition from folivory to frugivory etc...
#'             
#'             #QtransitionRateEst[which(!is.na(QtransitionRateEst[,1]))[1],] <- as.vector(simmapdiet1[[1]]$Q[,1])
#'             write.table(as.vector(simmapdiet1[[1]]$Q[,1]), paste("OutputEvolModel/Output_simmap_transition", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
#'             #QtransitionRateEqual[which(!is.na(QtransitionRateEqual[,1]))[1],] <- as.vector(simmapdiet2[[1]]$Q[,1])
#'             #QtransitionRatePrior[which(!is.na(QtransitionRatePrior[,1]))[1],] <- as.vector(simmapdiet3[[1]]$Q[,1])
#' 
#'             # Evolutionary history of traits (~brain size) with and without competition
#'             #following Drury et al.'s approach to see whether one trait (here one brain area size) phylogenetical history is better described if considering competition
#'             #https://academic.oup.com/sysbio/article/65/4/700/1753588
#'             #https://journals.plos.org/plosbiology/article?rev=1&id=10.1371/journal.pbio.2003563#pbio.2003563.ref028
#'             ##--------
#'             
#'             #Create variable of rinterest
#'             summaryData$ratioBrain <- summaryData$Brain*1.036*(10**-3)/summaryData$Bodymass #Following decasien for multiplication by 1.036
#'             summaryData$EQ <- summaryData$Brain*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775) #Following decasien, according to #Jerison, H. J. Evolution of the Brain and Intelligence (Academic, 1973).
#'             
#'             #Make it having symmetrical (and if possible gaussian distribution, since it seems to be a prerequisite of the analysis)
#'             hist(summaryData$ratioBrain)
#'             hist(log(summaryData$ratioBrain))#good enough normal
#'             #summaryData$ratioBrain.log <- log(summaryData$ratioBrain)
#'             
#'             hist(summaryData$EQ)
#'             #summaryData$EQ.log <- log(summaryData$EQ)
#'             #hist(summaryData$EQ.log)
#' 
#'             hist(summaryData$Brain)
#'             summaryData$Brain.log <- log(summaryData$Brain)
#'             hist(summaryData$Brain.log)
#'             
#'             # cor.test(summaryData$EQ, summaryData$Brain)
#'             # cor.test(summaryData$EQ, summaryData$ratioBrain)
#'             # cor.test(summaryData$Brain, summaryData$ratioBrain)
#'             
#'             #Reload tree to have same than used for biogeobears
#'             phylo <- read.tree("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/Tree_biogeobears.nex")
#'             
#'             colnames(summaryData)[colnames(summaryData)=="DietaryGuild"] <- "Guild"
#'             # tableDrury <- cbind(summaryData$SpeciesForPhylogeny,summaryData$EQ,summaryData$Guild,summaryData$geographicCode) 
#'             # colnames(tableDrury) <- c("SpeciesForPhylogeny", "Trait", "Guild", "geographicCode")
#'             # 
#'             # write.table(tableDrury , "table_data_rdc_jonathan_drury.txt", sep="\t")
#'             # 
#'             resultBrainFrugivory <- runComparisonModelsCompetition(
#'               simmap=simmapdiet1,
#'               data=summaryData[!is.na(summaryData$Brain.log)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
#'               subgroup="Fruit",
#'               numberMaps=numberSimulations,
#'               trait="Brain.log",
#'               tree=phylo,
#'               ana_events_tables=BSM_output$RES_ana_events_tables,
#'               clado_events_tables=BSM_output$RES_clado_events_tables
#'             )
#'             
#'             resultEQFrugivory <- runComparisonModelsCompetition(
#'               simmap=simmapdiet1,
#'               data=summaryData[!is.na(summaryData$EQ)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
#'               subgroup="Fruit",
#'               numberMaps=numberSimulations,
#'               trait="EQ",
#'               tree=phylo,
#'               ana_events_tables=BSM_output$RES_ana_events_tables,
#'               clado_events_tables=BSM_output$RES_clado_events_tables
#'             )
#'             
#'             #Neocortex
#'             summaryData$ratioNeocortex <- summaryData$Neocortex/ summaryData$Brain
#'             hist(summaryData$ratioNeocortex )
#'             #summaryData$ratioNeocortex.log <- log(summaryData$ratioNeocortex)
#'             #hist(summaryData$ratioNeocortex.log)
#'             
#'             resultNeocortexFrugivory <- runComparisonModelsCompetition(
#'               simmap=simmapdiet1,
#'               data=summaryData[!is.na(summaryData$ratioNeocortex)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
#'               subgroup="Fruit",
#'               numberMaps=numberSimulations,
#'               trait="ratioNeocortex",
#'               tree=phylo,
#'               ana_events_tables=BSM_output$RES_ana_events_tables,
#'               clado_events_tables=BSM_output$RES_clado_events_tables
#'             )
#'             
#'             #Hippocampus
#'             summaryData$ratioHippocampus <- summaryData$Hippocampus/ summaryData$Brain
#'             hist(summaryData$ratioHippocampus )
#'             
#'             resultHippocampusFrugivory <- runComparisonModelsCompetition(
#'               simmap=simmapdiet1,
#'               data=summaryData[!is.na(summaryData$ratioHippocampus)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
#'               subgroup="Fruit",
#'               numberMaps=numberSimulations,
#'               trait="ratioHippocampus",
#'               tree=phylo,
#'               ana_events_tables=BSM_output$RES_ana_events_tables,
#'               clado_events_tables=BSM_output$RES_clado_events_tables
#'             )
#'             
#'             #Cerebellum
#'             summaryData$ratioCerebellum <- summaryData$Cerebellum/ summaryData$Brain
#'             hist(summaryData$ratioCerebellum )
#'             #summaryData$ratioCerebellum.log <- log(summaryData$ratioCerebellum)
#'             #hist(summaryData$ratioCerebellum.log)
#'             
#'             resultCerebellumFrugivory <- runComparisonModelsCompetition(
#'               simmap=simmapdiet1,
#'               data=summaryData[!is.na(summaryData$ratioCerebellum)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
#'               subgroup="Fruit",
#'               numberMaps=numberSimulations,
#'               trait="ratioCerebellum",
#'               tree=phylo,
#'               ana_events_tables=BSM_output$RES_ana_events_tables,
#'               clado_events_tables=BSM_output$RES_clado_events_tables
#'             )
#'             
#'             #Striatum
#'             summaryData$ratioStriatum <- summaryData$Striatum/ summaryData$Brain
#'             hist(summaryData$ratioStriatum)
#'             
#'             resultStriatumFrugivory <- runComparisonModelsCompetition(
#'               simmap=simmapdiet1,
#'               data=summaryData[!is.na(summaryData$ratioStriatum)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
#'               subgroup="Fruit",
#'               numberMaps=numberSimulations,
#'               trait="ratioStriatum",
#'               tree=phylo,
#'               ana_events_tables=BSM_output$RES_ana_events_tables,
#'               clado_events_tables=BSM_output$RES_clado_events_tables
#'             )
#'             
#'             #MOB
#'             summaryData$ratioMOB <- summaryData$MOB/ summaryData$Brain
#'             hist(summaryData$ratioMOB)
#'             summaryData$ratioMOB.log <- log(summaryData$ratioMOB)
#'             hist(summaryData$ratioMOB.log)
#'             
#'             resultMOBFrugivory <- runComparisonModelsCompetition(
#'               simmap=simmapdiet1,
#'               data=summaryData[!is.na(summaryData$ratioMOB)&!is.na(summaryData$geographicCode)&summaryData$SpeciesForPhylogeny%in%phylo$tip.label,],
#'               subgroup="Fruit",
#'               numberMaps=numberSimulations,
#'               trait="ratioMOB.log",
#'               tree=phylo,
#'               ana_events_tables=BSM_output$RES_ana_events_tables,
#'               clado_events_tables=BSM_output$RES_clado_events_tables
#'             )
#'             
#'             # #Optic.Tract
#'             # runBioGeoBearsandBSM_DEC(
#'             #   pathToTree="~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3_rdcbrain.nex",
#'             #   pathToGeo="test_geo.data",
#'             #   maxCooccurrenceGeo=4,
#'             #   coreNumber=8,
#'             #   pathToSaveRdata="BioGeoBEARS_brain.Rdata",
#'             #   pathToSaveFigureGeo="~/PhD/Meta_analysis/Cognition_metaanalysis/Brain_DEC_v1.pdf",
#'             #   pathToBSMMapList="BSM_inputs_file.Rdata",
#'             #   pathToBSMoutput="BSM_output_file.Rdata",
#'             #   pathToFileStockingPdfForBSMAnalysis="~/PhD/Meta_analysis/Cognition_metaanalysis/"
#'             # )
#'             # 
#'             # summaryData$ratioOptic.Tract <- summaryData$Optic.Tract/ summaryData$Brain
#'             # hist(summaryData$ratioOptic.Tract)
#'             # resultOptic.TractFrugivory <- runComparisonModelsCompetition(
#'             #   simmap=simmapdiet3,
#'             #   data=summaryData,
#'             #   subgroup="Fruit",
#'             #   numberMaps=numberSimulations,
#'             #   trait="ratioOptic.Tract",
#'             #   tree=tree,
#'             #   ana_events_tables=BSM_output$RES_ana_events_tables,
#'             #   clado_events_tables=BSM_output$RES_clado_events_tables
#'             # )
#'             
#'             start=which(is.na(summaryBrainFrugivory[,1]))[1]
#'             end=which(is.na(summaryBrainFrugivory[,1]))[1] + numberSimulations - 1
#'               
#'             
#'             summaryBrainFrugivory[start:end,1:6] <-
#'               resultBrainFrugivory[,(ncol(resultBrainFrugivory)-5):ncol(resultBrainFrugivory)]
#'             summaryEQFrugivory[start:end,1:6] <-
#'               resultEQFrugivory[,(ncol(resultEQFrugivory)-5):ncol(resultEQFrugivory)]
#'             summaryNeocortexFrugivory[start:end,1:6] <-
#'               resultNeocortexFrugivory[, (ncol(resultNeocortexFrugivory)-5):ncol(resultNeocortexFrugivory)]
#'             summaryHippocampusFrugivory[start:end,1:6] <-
#'               resultHippocampusFrugivory[, (ncol(resultHippocampusFrugivory)-5):ncol(resultHippocampusFrugivory)]
#'             summaryCerebellumFrugivory[start:end,1:6] <-
#'               resultCerebellumFrugivory[, (ncol(resultCerebellumFrugivory)-5):ncol(resultCerebellumFrugivory)]
#'             summaryStriatumFrugivory[start:end,1:6] <-
#'               resultStriatumFrugivory[, (ncol(resultStriatumFrugivory)-5):ncol(resultStriatumFrugivory)]
#'             summaryMOBFrugivory[start:end,1:6] <-
#'               resultMOBFrugivory[, (ncol(resultMOBFrugivory)-5):ncol(resultMOBFrugivory)]
#'             summaryOpticTractFrugivory[start:end,1:6] <-
#'               resultOptic.TractFrugivory[, (ncol(resultOptic.TractFrugivory)-5):ncol(resultOptic.TractFrugivory)]
#'             
#'             
#'             #Part added because now parallelizing code // to avoid to many changes
#'             summaryBrainFrugivory <- summaryBrainFrugivory[start:end,]
#'             summaryEQFrugivory <- summaryEQFrugivory[start:end,]
#'             summaryNeocortexFrugivory <- summaryNeocortexFrugivory[start:end,]
#'             summaryHippocampusFrugivory <- summaryHippocampusFrugivory[start:end,]
#'             summaryCerebellumFrugivory <- summaryCerebellumFrugivory[start:end,]
#'             summaryStriatumFrugivory <- summaryStriatumFrugivory[start:end,]
#'             summaryMOBFrugivory <- summaryMOBFrugivory[start:end,]
#'             
#'             write.table(summaryBrainFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Brain", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
#'             write.table(summaryEQFrugivory, paste("OutputEvolModel/Output_evolutionary_history_EQ", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
#'             write.table(summaryNeocortexFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Neocortex", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
#'             write.table(summaryHippocampusFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Hippocampus", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
#'             write.table(summaryCerebellumFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Cerebellum", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
#'             write.table(summaryStriatumFrugivory, paste("OutputEvolModel/Output_evolutionary_history_Striatum", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
#'             write.table(summaryMOBFrugivory, paste("OutputEvolModel/Output_evolutionary_history_MOB", number, ".txt", sep=""), row.names=FALSE, col.names=TRUE, sep="\t")
#'             #write.table(summaryOpticTractFrugivory, "Output_evolutionary_history_OpticTract.txt", row.names=FALSE, col.names=TRUE, sep="\t")
#'             
#'             
#'             #progression=progression+numberSimulations
#'             #progress(progression*100/totModels)
#'           }
#'           #stop cluster
#'           stopCluster(cl)
#'       #}
#'     }
#'   }
#' }
#' 
#' save.image("geography_trait_models.RData")


###----------------------
## Plotting results of phylogenetic history
###----------------------

tree <-read.nexus("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/TreeBlock_10kTrees_Primates_Version3.nex")
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
write.tree(phylo, np("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/Tree_diet.nex"))
tree <- read.tree("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/Tree_diet.nex")


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


colourHip=colourVector[5]
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

legend3d(x=0.1, y=0.2, legend = c("Cerebellum", "Hippocampus", "MOB", "Striatum"), fill = c(colourCereb, colourHip, colourOlf, colourStri), cex=1, bty="n", ncol=2)
#rglwidget()


#from: https://r-graphics.org/recipe-miscgraph-3d-save
#rgl.snapshot('3dplot.png', fmt = 'png', top = TRUE, width = 300)
snapshot3d('3dplot.png', fmt = 'png', top = TRUE, webshot = FALSE, width=32000, height=32000)
#rgl.postscript('3dplot.pdf', fmt = 'pdf')

#remotes::install_github("rstudio/chromote")
#remotes::install_github("rstudio/webshot2")

library(webshot2)
snapshot3d('3dplot.png', fmt = 'png')


writeASY()
rgl.postscript('3dplot', fmt="svg")



, width = 300, height = 300


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
## Fig 3: results of the models -------->>>>>>>>>>>> ADD THE EB RESULTS !!!!!!!!
########

##---------
## Check sample size

transitionMatrix <- matrix(NA, nrow=(9), ncol=2)

for(i in 0:(9)){
  start=which(is.na(transitionMatrix[,1]))[1]
  
  toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/OutputEvolModel_test_ENS/Output_simmap_transition",i, ".txt", sep=""))
  transitionMatrix[start,1] <- toAdd[1,1]
  transitionMatrix[start,2] <- toAdd[2,1]
}

minProba.v <- apply(transitionMatrix, 2, min)
maxProba.v <- apply(transitionMatrix, 2, max)


##---------
## Check transition matrix

repetition=(9)#length(frugivoryThresholdVector)*length(folivoryThresholdVector)*length(geographicThresholdVector)*randomSampling
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


for(i in 0:(9)){

  toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/Sample_size/checkSampleFruit",i, ".txt", sep=""))
  checkSampleFruit[i] <- toAdd[1]
  toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/Sample_size/checkSampleLeaf",i, ".txt", sep=""))
  checkSampleLeaf[i] <- toAdd[1]
  toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/Sample_size/checkSampleRange",i, ".txt", sep=""))
  checkSampleRange[i] <- toAdd[1]
  toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/Sample_size/checkSampleBrain",i, ".txt", sep=""))
  checkSampleBrain[i] <- toAdd[1]
  toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/Sample_size/checkSampleEQ",i, ".txt", sep=""))
  checkSampleEQ[i] <- toAdd[1]
  toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/Sample_size/checkSampleNeocortex",i, ".txt", sep=""))
  checkSampleNeocortex[i] <- toAdd[1]
  toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/Sample_size/checkSampleHippocampus",i, ".txt", sep=""))
  checkSampleHippocampus[i] <- toAdd[1]
  toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/Sample_size/checkSampleCerebellum",i, ".txt", sep=""))
  checkSampleCerebellum[i] <- toAdd[1]
  toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/Sample_size/checkSampleStriatum",i, ".txt", sep=""))
  checkSampleStriatum[i] <- toAdd[1]
  toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/Sample_size/checkSampleMOB",i, ".txt", sep=""))
  checkSampleMOB[i] <- toAdd[1]

}

#Min-max sample
sapply(checkSampleFruit, function(x) paste(c(min(x), max(x)), collapse="-"))
sapply(checkSampleLeaf, function(x) paste(c(min(x), max(x)), collapse="-"))
sapply(checkSampleRange, function(x) paste(c(min(x), max(x)), collapse="-"))
sapply(checkSampleBrain, function(x) paste(c(min(x), max(x)), collapse="-"))
sapply(checkSampleEQ, function(x) paste(c(min(x), max(x)), collapse="-"))
sapply(checkSampleNeocortex, function(x) paste(c(min(x), max(x)), collapse="-"))
sapply(checkSampleHippocampus, function(x) paste(c(min(x), max(x)), collapse="-"))
sapply(checkSampleCerebellum, function(x) paste(c(min(x), max(x)), collapse="-"))
sapply(checkSampleStriatum, function(x) paste(c(min(x), max(x)), collapse="-"))
sapply(checkSampleMOB, function(x) paste(c(min(x), max(x)), collapse="-"))
sapply(checkSampleRange, function(x) paste(c(min(x), max(x)), collapse="-"))


##---------
## Load the data

summaryBrainFrugivory <- as.data.frame(matrix(NA, nrow=10*((9)+1), ncol=53))
summaryEQFrugivory <- as.data.frame(matrix(NA, nrow=10*((9)+1), ncol=53))
summaryNeocortexFrugivory <- as.data.frame(matrix(NA, nrow=10*((9)+1), ncol=53))
summaryHippocampusFrugivory <- as.data.frame(matrix(NA, nrow=10*((9)+1), ncol=53))
summaryCerebellumFrugivory <- as.data.frame(matrix(NA, nrow=10*((9)+1), ncol=53))
summaryStriatumFrugivory <- as.data.frame(matrix(NA, nrow=10*((9)+1), ncol=53))
summaryMOBFrugivory <- as.data.frame(matrix(NA, nrow=10*((9)+1), ncol=53))

for(i in 0:((9))){
  start=which(is.na(summaryBrainFrugivory[,1]))[1]
  end=which(is.na(summaryBrainFrugivory[,1]))[1] + numberSimulations - 1
  
  tryCatch(
    {toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/OutputEvolModel_test_ENS/Output_evolutionary_history_Brain",i, ".txt", sep=""))
  summaryBrainFrugivory[start:end,] <- as.data.frame(toAdd)
  }, error=function(e){
    #Do nothing
    }
  )
  
  tryCatch(
  {toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/OutputEvolModel_test_ENS/Output_evolutionary_history_EQ",i, ".txt", sep=""))
  summaryEQFrugivory[start:end,] <- toAdd
  }, error=function(e){
    #Do nothing
  }
  )
  
  tryCatch(
    {toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/OutputEvolModel_test_ENS/Output_evolutionary_history_Neocortex",i, ".txt", sep=""))
  summaryNeocortexFrugivory[start:end,] <- toAdd
    }, error=function(e){
      #Do nothing
    }
  )
  
  tryCatch(
    {toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/OutputEvolModel_test_ENS/Output_evolutionary_history_Hippocampus",i, ".txt", sep=""))
  summaryHippocampusFrugivory[start:end,] <- toAdd
    }, error=function(e){
      #Do nothing
    }
  )
  
  tryCatch(
    {toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/OutputEvolModel_test_ENS/Output_evolutionary_history_Cerebellum",i, ".txt", sep=""))
  summaryCerebellumFrugivory[start:end,] <- toAdd
    }, error=function(e){
      #Do nothing
    }
  )
  
  tryCatch(
    {toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/OutputEvolModel_test_ENS/Output_evolutionary_history_Striatum",i, ".txt", sep=""))
  summaryStriatumFrugivory[start:end,] <- toAdd
    }, error=function(e){
      #Do nothing
    }
  )
  
  tryCatch(
    {toAdd <- read.delim(paste("~/PhD/Meta_analysis/Cognition_metaanalysis/OutputEvolModel_test_ENS/Output_evolutionary_history_MOB",i, ".txt", sep=""))
  summaryMOBFrugivory[start:end,] <- toAdd
    }, error=function(e){
      #Do nothing
    }
  )
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

colNum <- c("darkgray", brewer.pal(n = 5, name = "Set1")[5], brewer.pal(n = 4, name = "Set1"))

models <- c("BM", "OU", "EB", "MC", "DDlin", "DDexp")
colourModels <- brewer.pal(n = 6, name = "Set1")

## Brain

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
for(i in 1:6){
  meanPt <- mean(as.numcharac(summaryBrainFrugivory[, ncol(summaryBrainFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryBrainFrugivory[, ncol(summaryBrainFrugivory)-6+i]))
  #sd <- sd/nrow(summaryBrainFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
}

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
       paste("b~", round(mean(as.numcharac(summaryEQFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryEQFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

draw.circle(x=0.3,y=1.1,0.2, col=colNum[1], border=NA)
text(x=0.3, y=1.1, labels="1", xpd=TRUE, col="white", font=2, cex=1.5)

##----------------------


layout(mat=rbind(c(1,2,3), c(4,5,6)), widths=c(5,5,5), heights=c(5,5))
par(mar=c(4, 4, 2, 1), mgp=c(2, 0.5, 0), xpd=TRUE)
#note: 1= second run for frugivory 20%
#note: _2= first run for frugivory 20%

# 
# plot(
#   x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
#   xlim=c(0,7), ylim=c(0,1.2),
#   las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
#   xaxt="n",xaxs="i",yaxs="i", yaxt="n")
# 
# source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
# addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
# axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
# text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
# for(i in 1:6){
#   
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
#        paste("b~", round(mean(as.numcharac(summaryEQFrugivory$DDlingeo.b)), digit=3), sep=""),
#        paste("r~", round(mean(as.numcharac(summaryEQFrugivory$DDexpgeo.r)), digit=3), sep="")
#        ), xpd=TRUE)
# 
# library(plotrix)
# draw.circle(x=0.3,y=1.1,0.2, col=colNum[1], border=NA)
# text(x=0.3, y=1.1, labels="1", xpd=TRUE, col="white", font=2)
# 
# ##----------------------

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
draw.circle(x=0.3,y=1.1,0.2, col=colNum[6], border=NA)
text(x=0.3, y=1.1, labels="2", xpd=TRUE, col="white", font=2)


#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(summaryStriatumFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryStriatumFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)


##------------

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
draw.circle(x=0.3,y=1.1,0.2, col=colNum[5], border=NA)
text(x=0.3, y=1.1, labels="3", xpd=TRUE, col="white", font=2)

#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(summaryMOBFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryMOBFrugivory$DDexpgeo.r)), digit=3), sep="")
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
draw.circle(x=0.3,y=1.1,0.2, col=colNum[2], border=NA)
text(x=0.3, y=1.1, labels="5", xpd=TRUE, col="white", font=2)


#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(summaryNeocortexFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryNeocortexFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)


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
draw.circle(x=0.3,y=1.1,0.2, col=colNum[3], border=NA)
text(x=0.3, y=1.1, labels="6", xpd=TRUE, col="white", font=2)


#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(summaryHippocampusFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryHippocampusFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

##------------

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
draw.circle(x=0.3,y=1.1,0.2, col=colNum[4], border=NA)
text(x=0.3, y=1.1, labels="7", xpd=TRUE, col="white", font=2)

#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(summaryCerebellumFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryCerebellumFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

###----------------------


###
## Phylogeny analysis of the current observed pattern: phylogenetic regression
##






#############
#############
## GARBAGE


    #--------------------------------------------------------------------------------
    
    ##########
    ## STATISTICAL ANALYSES: IS COGNITION PROMOTED BY INTERSPECIES COMPETITION?
    ##########
    
    
    nameForPhylo <- as.data.frame(phylo$tip.label)
    colnames(nameForPhylo) <- "Species"
    nameForPhylo$Species <- as.character(nameForPhylo$Species)
    
    options(warn=1)
    nameForPhylo <- separate(nameForPhylo, col=Species, into=c("Name1", "Name2"), sep="_", remove=FALSE)#some have three "_", bot a big deal for now
    nameForPhylo$Species_abbrv <- paste(nameForPhylo$Name1, substr(nameForPhylo$Name2,1,4), sep="_")
    summaryKamilar$SpeciesPhylo <- nameForPhylo$Species[match(summaryKamilar$Species_abbrv,nameForPhylo$Species_abbrv)]
  
    #Phylogenetical tree
    
    library(ape)
    library(phytools)
    #Consensus tree from 10k Trees
    phylo_all <-read.nexus("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3.nex")
    phylo_init <- phylo_all
    
    is.ultrametric(phylo_init)
    #Tree is not ultrametric while it should
    phylo <- force.ultrametric(tree=phylo_init, method="extend")#method="nnls")
    is.ultrametric(phylo)
    
    #Check changes
    plot(phylo_init$edge.length, phylo$edge.length, xlab="Original tree", ylab="Post-chronos tree", pch=20, bty="n")#no much changes; good
    abline(a=0, b=1, col="gray", lwd=0.5)
    
    #Matching name of dataset to phylogenetics name of tree
    nameForPhylo <- as.data.frame(phylo$tip.label)
    colnames(nameForPhylo) <- "Species"
    nameForPhylo$Species <- as.character(nameForPhylo$Species)
    
    library(tidyr)
    options(warn=1)
    nameForPhylo <- separate(nameForPhylo, col=Species, into=c("Name1", "Name2"), sep="_", remove=FALSE)#some have three "_", bot a big deal for now
    nameForPhylo$Species_abbrv <- paste(nameForPhylo$Name1, substr(nameForPhylo$Name2,1,4), sep="_")
    summaryKamilar$SpeciesPhylo <- nameForPhylo$Species[match(summaryKamilar$Species_abbrv,nameForPhylo$Species_abbrv)]
    
    dataForBrain <- summaryKamilar[, c(1,2,3,5,6,7,10,13,14,15,16,17,18,ncol(summaryKamilar))]
    dataForBrain <- dataForBrain[!is.na(dataForBrain$OverlapHR)&
                                   !is.na(dataForBrain$Brain)&
                                   !is.na(dataForBrain$Bodymass),]
    dataForBrain <- dataForBrain[!is.na(dataForBrain$SpeciesPhylo),]
    nrow(dataForBrain)#lost one
    
    '%nin%' <- Negate('%in%')
    
    dataForBrain_rdc <- dataForBrain[
      dataForBrain$DietGuild=="Fruit"|
        dataForBrain$DietGuild=="Leaves",]
    
    phylo <- drop.tip(phylo,
                      phylo$tip.label[
                        which(phylo$tip.label
                              %nin%unique(dataForBrain_rdc$SpeciesPhylo))]) 
    


    
    ####
    ## Fig interaction
    ####
    
    #TO MODIFY WITH OWN DATA, CONSTRUCTED ON A TEST DATASET
    library(ape)
    #data(bird.orders)
    
    numberOrder=6
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
    
    #Link between species/groups:
    speciesLabels <- hc$labels#Should be in the tree order
    #groupLabels <- rep(seq(from=1, to=23, by=2), each=2)[-1]
    
    
    #Geo + diet
    geoBinary <- as.data.frame(summaryDataForPlot)
    colnames(geoBinary) <- c("SpeciesPhylo", "Loc")
    locationSpecies <- geoBinary$Loc[match(speciesLabels, geoBinary$SpeciesPhylo)]
    dietSpecies <- summaryKamilar$Guild[match(speciesLabels, summaryKamilar$SpeciesPhylo)]
    colLoc <- c(brewer.pal(n = 9, name = "Pastel1"), "gray93")#1:length(unique(tableSiteGeo$Locality))#brewer.pal(n = length(unique(tableSiteGeo$Locality)), name = "Set1")
    loc.v <- tableSiteGeo$Locality
    
    #Brain data
    relativeValueBrain <- summaryKamilar$EQ[match(speciesLabels, summaryKamilar$SpeciesPhylo)] - 1#runif(length(speciesLabels), -1, 1)
    #relativeValueBrain <- scale(log(relativeValueBrain))
    positiveCol="magenta4"#pastellize(x="purple", p=0.2)
    negativeCol="darkgreen"#pastellize(x="green", p=0.2)
    
    
    library(circlize)
    circos.clear()
    circos.par(gap.degree=0, gap.after=0, cell.padding=c(0,0,0,0), track.margin = c(0, 0))
    circos.initialize(speciesLabels, xlim = c(0, 1))
    
    
    ###TRIAL 1
    #Brain size + diet
    
    absMax <- max(abs(relativeValueBrain))
    
    library(RColorBrewer)
    colourPositive=brewer.pal(n = 5, name = "Pastel1")[1]
    colourNegative=brewer.pal(n = 5, name = "Pastel1")[2]
    
    #Plot relative brain size
    #Background
    circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
      circos.rect(0, 0, 1, 1, col=colourPositive, border=colourPositive)
    }, track.height = 0.1)
    
    circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
      circos.rect(0, 0, 1, 1, col=colourNegative, border=colourNegative)
    }, track.height = 0.1)
    
    library(plotrix)
    #Main circle
    draw.circle(x=0,y=0,0.9, col=NA, border="white")
    
    #increment of 0.5
    for(i in 1:5){
      draw.circle(x=0,y=0,1-(i-1)*0.05, col=NA, border="white", lty=2)
    }
    
    
    #Value
    #for(i in 1:length(speciesLabels)){
    circos.track(ylim = c(0, 1), bg.border = NA, track.index=1, panel.fun = function(x, y) {
      i=CELL_META$sector.numeric.index
      #circos.rect(0, 0, 1, 1, col=colourPositive, border=colourPositive)
      if(relativeValueBrain[i] > 0 & dietSpecies[i]=="Fruit"){
        circos.points(CELL_META$xcenter, relativeValueBrain[i]/absMax, pch=19, col="black")
        circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueBrain[i]/absMax, lty=3)
      }
      else if(relativeValueBrain[i] > 0 & dietSpecies[i]=="Leaves"){
        circos.points(CELL_META$xcenter, relativeValueBrain[i]/absMax, pch=21, col="black", bg="white")
        circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueBrain[i]/absMax, lty=3)
      }
      else{}
    }, track.height = 0.1)
    
    circos.track(ylim = c(0, 1), bg.border = NA, track.index=2,  panel.fun = function(x, y) {
      i=CELL_META$sector.numeric.index
      #circos.rect(0, 0, 1, 1, col=colourNegative, border=colourNegative)
      if(relativeValueBrain[i] <= 0 & dietSpecies[i]=="Fruit"){
        circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueBrain[i]/absMax, lty=3)
        circos.points(CELL_META$xcenter, 1 + relativeValueBrain[i]/absMax, pch=19, col="black")
      }
      else if(relativeValueBrain[i] <= 0 & dietSpecies[i]=="Leaves"){
        circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueBrain[i]/absMax, lty=3)
        circos.points(CELL_META$xcenter, 1 + relativeValueBrain[i]/absMax, pch=21, col="black", bg="white")
      }
      else{}
    }, track.height = 0.1)
    
    #Species name
    circos.track(ylim = c(0, 20), bg.border = NA, track.height = 0.1, track.margin=c(0.01, 0.1),
                 panel.fun = function(x, y) {
                   i=CELL_META$sector.numeric.index
                   circos.text(CELL_META$xcenter, 0, labels.rdc[i], adj = c(0, 0), 
                               facing = "clockwise", niceFacing = TRUE,
                               col = "black", cex = 0.5, font=3)
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
          if(dietSpecies[i]==dietSpecies[j]){
            colour <- as.data.frame(table(colLoc[which(product==1)]))
            colour <- colour[colour$Freq==max(colour$Freq),1][1]
            circos.link(speciesLabels[i], runif(1, 0, 1), speciesLabels[j], runif(1, 0, 1), lwd=1, col=as.character(colour))
          }
          else{
            #circos.link(speciesLabels[i], runif(1, 0, 1), speciesLabels[j], runif(1, 0, 1), lwd=1, col="lightgray")
          }
        }
      }
    }
    circos.clear()
    
    
    #Plot the phylogenetic tree in a new circular plot
    ct = cutree(hc, numberOrder)  # cut tree into order
    n = length(labels)  # number of bird species
    dend = as.dendrogram(hc)
    
    library(circlize)
    par(new = TRUE) # <- magic
    circos.par("canvas.xlim" = c(-1.75, 1.75), "canvas.ylim" = c(-1.75, 1.75))
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
    
    
    
    
    ###TRIAL 2
    library(circlize)
    circos.clear()
    circos.par(gap.degree=0, gap.after=0, cell.padding=c(0,0,0,0))
    circos.initialize(speciesLabels, xlim = c(0, 1))
    
    # #Plot the family
    # circos.track(ylim = c(0, 1), panel.fun = function(x,y) {
    #     i=CELL_META$sector.numeric.index
    #     #print(i)
    #     if(i!=1&i!=length(groupLabels)){
    #       before=groupLabels[i-1]
    #       now=groupLabels[i]
    #       after=groupLabels[i+1]
    #     } else {
    #       if(i==1){
    #         before=groupLabels[length(groupLabels)]
    #         now=groupLabels[i]
    #         after=groupLabels[i+1]
    #       } else{
    #         before=groupLabels[i-1]
    #         now=groupLabels[i]
    #         after=groupLabels[1]
    #       }
    #     }
    #     if(before==now & after==now){
    #       circos.rect(0, 0, 1, 1, col="lightgray", border="lightgray")
    #     }
    #     else{
    #       if(before!=now){
    #         circos.rect(0.1, 0, 1, 1, col="lightgray", border="lightgray")
    #       } else{
    #         circos.rect(0, 0, 0.9, 1, col="lightgray", border="lightgray")
    #       }
    #     }
    # }, bg.border = NA, track.height = 0.1)
    # 
    
    #Plot the pictogram for diet
    
    library(png)
    fruit <- readPNG("~/PhD/Meta_analysis/Cognition_metaanalysis/Metaanalysis/fruit3.png")
    leaf <- readPNG("~/PhD/Meta_analysis/Cognition_metaanalysis/Metaanalysis/leaf.png")
    dim(fruit)
    dim(leaf)
    
    fruit.raster <- as.raster(fruit)
    leaf.raster <- as.raster(leaf)
    
    circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) { 
      print(CELL_META$sector.numeric.index)
      i=CELL_META$sector.numeric.index
      if(dietSpecies[i]=="Leaves"){
        circos.raster(leaf.raster, CELL_META$xcenter, CELL_META$ycenter, 
                      width = "0.25cm",#CELL_META$xrange,
                      #height = dim(leaf)[1]/dim(leaf)[2]*CELL_META$xrange,
                      facing = "inside")
        # circos.raster(leaf.raster, CELL_META$xcenter, CELL_META$ycenter, 
        #               width = CELL_META$xrange,
        #               height = dim(leaf)[1]/dim(leaf)[2]*CELL_META$xrange,
        #               facing = "bending.inside")
        #addImg(leaf,  CELL_META$xcenter, CELL_META$ycenter, width = dim(leaf)[2]/dim(leaf)[1]*CELL_META$xrange)
        #plot(leaf.raster)
      } else{
        #print("fruit")
        circos.raster(fruit.raster, CELL_META$xcenter, CELL_META$ycenter, 
                      width = "0.25cm",#CELL_META$xrange,
                      #height = dim(leaf)[1]/dim(leaf)[2]*CELL_META$xrange,
                      facing = "inside")
        #circos.raster(fruit.raster, CELL_META$xcenter, CELL_META$ycenter, 
        # width = CELL_META$xrange, 
        # height = dim(fruit)[1]/dim(fruit)[2]*CELL_META$xrange, 
        # facing = "bending.inside")
      }
    }, track.height = 0.1)
    
    
    #Species name
    circos.track(ylim = c(0, 1), bg.border = NA, track.height = 0.1, 
                 panel.fun = function(x, y) {
                   i=CELL_META$sector.numeric.index
                   circos.text(i-0.5, 0, labels.rdc[i], adj = c(0, 0), 
                               facing = "clockwise", niceFacing = TRUE,
                               col = "black", cex = 0.5, font=3)
                 })
    
    #Plot the relative brain size
    circos.track(ylim = c(-4, 4), bg.border = NA, panel.fun = function(x, y) {
      i=CELL_META$sector.numeric.index
      circos.barplot(relativeValueBrain[i], 0, col=ifelse(relativeValueBrain[i] > 0, positiveCol, negativeCol), 
                     border=ifelse(relativeValueBrain[i] > 0, positiveCol, negativeCol))
    }, track.height = 0.3)
    
    library(plotrix)
    draw.circle(x=0,y=0,0.725, lty=2)
    
    
    #Plot the geographic links
    for(i in 1:length(speciesLabels)){
      #locI <- which(strsplit(locationSpecies[i], "")==1)
      for(j in 1:length(speciesLabels)){
        #locJ <- which(strsplit(locationSpecies[j], "")==1)
        product <- as.numcharac(unlist(strsplit(locationSpecies[j], "")))*as.numcharac(unlist(strsplit(locationSpecies[i], "")))
        if(i==j|(length(unique(product))==1&unique(product)[1]==0)){
          #Do nothing
        }
        else{
          if(dietSpecies[i]==dietSpecies[j]){
            colour <- as.data.frame(table(colLoc[which(product==1)]))
            colour <- colour[colour$Freq==max(colour$Freq),1][1]
            circos.link(speciesLabels[i], runif(1, 0, 1), speciesLabels[j], runif(1, 0, 1), lwd=1, col=colour)
          }
          else{
            circos.link(speciesLabels[i], runif(1, 0, 1), speciesLabels[j], runif(1, 0, 1), lwd=1, col="lightgray")
          }
        }
      }
    }
    circos.clear()
    
    
    #Plot the phylogenetic tree in a new circular plot
    ct = cutree(hc, numberOrder)  # cut tree into order
    n = length(labels)  # number of bird species
    dend = as.dendrogram(hc)
    
    library(circlize)
    par(new = TRUE) # <- magic
    circos.par("canvas.xlim" = c(-2, 2), "canvas.ylim" = c(-2, 2))
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
                 track.height = 0.5, panel.fun = function(x, y) {
                   circos.dendrogram(dend, col="lightgray")
                 })
    circos.clear()
    
    ##----------------------
 



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
    
    
    library(RColorBrewer)
    colourVector <- brewer.pal(n = 5, name = "Set1")
    
    colourHip=colourVector[1]
    colourCereb=colourVector[2]
    colourOlf=colourVector[3]
    colourStri=colourVector[4]
    
    ### this would be the ``activation'' or surface you want to render 
    
    #Set up the view manually the first time to extract coords
    #um <- par3d()$userMatrix
    #I obtained those as a correct view
    um <- matrix(c(-0.6391721,-0.76426947,-0.08573965,0,-0.2281326,0.08195215,0.97017491,0,-0.7344485,0.63966882,-0.22673632,0,0.0000000,0.00000000,0.00000000,1), ncol=4, nrow=4, byrow=TRUE)
    view3d(userMatrix = um)
    
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
    
    legend3d(x=0.1, y=0.2, legend = c("Cerebellum", "Hippocampus", "Olfactory bulb", "Striatum"), fill = c(colourCereb, colourHip, colourOlf, colourStri), cex=1, bty="n", ncol=2)
    rglwidget()
    
    #from: https://r-graphics.org/recipe-miscgraph-3d-save
    #rgl.snapshot('3dplot.png', fmt = 'png')
    rgl.postscript('3dplot.pdf', fmt = 'pdf')
    
    
    
    ##----------------------
    
    
    
    
    