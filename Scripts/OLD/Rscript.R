
areaName <- c("Madagascar_east", "Madagascar_west", "Africa_west", "Africa_central", "Africa_south_east", 
              "America_central", "America_south_north", "America_south_south", 
              "Asia_west", "Asia_central_east",  "Asia_south", "Asia_island")
 
library(RColorBrewer)
colourArea <- brewer.pal(n = 12, name="Paired")

#Select only the terrestrial area, because was roughly drawn with google earth

#+Note the shapefile were transformed from kml using: https://products.aspose.app/gis/conversion/kml-to-shapefile

library(rworldmap) # World map 
worldMap <- getMap()
library(cleangeo) #to clean it otherwise issues with intersection
worldMap = clgeo_Clean(worldMap)

library(rgeos) #for readOGR; gArea/gCentroid...
library(sf) #for intersection
library(rgdal)
centroid <- matrix(NA, ncol=2, nrow=length(areaName))
for(i in 1:length(areaName)){
  areaTransitory <- readOGR(dsn=paste("T:/IUCN_data_primate/Geographic_areas/Shapefiles/",areaName[i],".shp",sep=""))
  areaTransitory = clgeo_Clean(areaTransitory)
  areaTransitory <- spTransform(areaTransitory, CRS(proj4string(worldMap)))
  areaTransitory <- gIntersection(areaTransitory, worldMap)
  assign(paste("area", i, sep="_"), areaTransitory)
  if(i==1){
    centroid[i,] <- c(summary(areaTransitory)$bbox[1,2] + 5, summary(areaTransitory)$bbox[2,1]) 
  }
  else if (i==2){
    centroid[i,] <- c(summary(areaTransitory)$bbox[1,2] - 5, summary(areaTransitory)$bbox[2,1] - 5) 
  }
  else{
    centroid[i,] <- gCentroid(areaTransitory)@coords
  }
}
warnings()


#Create the map of the geographic area

#Have background
library(maps)
map("world", fill=TRUE, col="lightgray", bg="white", border=NA, ylim=c(-60, 50))

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
library(ape)
treeForName <-read.nexus("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3.nex")
speciesLabels <-  as.data.frame(treeForName$tip.label)

thresholdPresenceRange <- seq(from=5, to=40, by=2.5)/100
matrixRangingSensitivity <- matrix(NA, nrow=nrow(speciesLabels), ncol=length(thresholdPresenceRange)+1)
for(r in 1:length(thresholdPresenceRange)){
  print(r)
  thresholdPresence=thresholdPresenceRange[r]
  matrixPrimateRange <- matrix(NA, nrow=nrow(speciesLabels), ncol=length(areaName)+1)
  library(stringr)
  library(svMisc)
  for(i in 1:nrow(matrixPrimateRange)){
    matrixPrimateRange[i,1] <- speciesLabels[i,1]
    lengthString <- str_length(speciesLabels[i,1])
    substrSpeciesGeo <- sapply(speciesGeo, function (x){substr(x, 1, lengthString)})
    polygonSpeciesId <- which(substrSpeciesGeo==speciesLabels[i,1])
    if(!is.na(polygonSpeciesId[1])){
      speciesRangeTransitory <- primateSpeciesRange[polygonSpeciesId,]
      speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS(proj4string(worldMap)))
      
      for(j in 1:length(areaName)){
        if(!is.null(gIntersection(speciesRangeTransitory, get(paste("area", j, sep="_"))))){
          overlap <- gArea(gIntersection(speciesRangeTransitory, get(paste("area", j, sep="_"))))/gArea(speciesRangeTransitory)
        } else{
          overlap <- 0
        }
        if(overlap >= thresholdPresence){
          matrixPrimateRange[i, j+1] <- 1 
        } else{
          matrixPrimateRange[i, j+1] <- 0
        }
      }
      matrixRangingSensitivity[i,r] <- paste(matrixPrimateRange[i, 2:ncol(matrixPrimateRange)], collapse="")
    }
    progress(i/nrow(matrixPrimateRange)*100)
  }
}


ref <- matrixRangingSensitivity[,5]
howManyDifferent <- rep(NA, times=ncol(matrixRangingSensitivity))
for(i in 1:ncol(matrixRangingSensitivity)){
  howManyDifferent[i] <- length(which(matrixRangingSensitivity[,i]!=ref&!is.na(matrixRangingSensitivity[,i])))#matrixRangingSensitivity[,i]!=paste(rep(NA, times=12), collapse="")))
}
howManyDifferent <- howManyDifferent/length(!is.na(ref))#ref[ref!=paste(rep(NA, times=12), collapse="")])

nrow(matrixRangingSensitivity[!is.na(matrixRangingSensitivity[,5]),])


################################
## First: exploring dataset: having correct format and checking whether information is consensual or not
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
rep(1, times=10)
)

vectorSpecies <- c()

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}


library(tidyr)
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

library(ape)
treeForName <-read.nexus("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3.nex")

speciesLabels <-  as.data.frame(treeForName$tip.label)
colnames(speciesLabels) <- "SpeciesTree"
speciesLabels <- speciesLabels [order(speciesLabels$SpeciesTree),]
speciesLabels <- as.data.frame( speciesLabels)
colnames(speciesLabels) <- "SpeciesTree"
speciesLabels <- separate(speciesLabels, col=SpeciesTree, into=c("Part1", "Part2"), sep="_", remove=FALSE)
speciesLabels$Species_abbrv <- paste(speciesLabels$Part1, substr(speciesLabels$Part2, 1, 4), sep="_")
summaryData$SpeciesForPhylogeny <- speciesLabels$SpeciesTree[match(summaryData$Species_abbrv, speciesLabels$Species_abbrv)]
#summaryData$SpeciesForPhylogeny <- Data_powell$Species.name.adjusted.to.10kTrees[match(summaryData$Species_abbrv, Data_powell$Species_abbrv)]

  
#if not available, add initial date. it will then be modified once you selected species with available data
toCorrectSpecies <- as.data.frame(summaryData$Species_abbrv[which(is.na(summaryData$SpeciesForPhylogeny))])
colnames(toCorrectSpecies) <- "SpeciesInit"
library(tidyr)
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

Subdata_decasien_mosaic$Optic.Tract <- gsub(" " , "",as.character(Subdata_decasien_mosaic$Optic.Tract))
meanData <- aggregate(as.numeric(as.character(Subdata_decasien_mosaic$Optic.Tract)), by=list(Subdata_decasien_mosaic$Species_abbrv), FUN=mean, na.rm=TRUE)
summaryData$Optic.Tract_volume_mm3_decasien <- meanData$x[match(summaryData$Species_abbrv, meanData$Group.1)]

#MOB=main olfactive bulb

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

#Add the HR overlap from: Pearce and Willems

Data_pearce_rdcoverlap <- Data_pearce[,c(19,10)]
Data_pearce_rdcoverlap$Species_abbrv <- gsub(" ","_",Data_pearce_rdcoverlap$Species_abbrv)
Data_pearce_rdcoverlap <- unique(Data_pearce_rdcoverlap)
Data_pearce_rdcoverlap <- aggregate(Data_pearce_rdcoverlap$ov, by=list(Data_pearce_rdcoverlap[,1]), FUN=mean)
colnames(Data_pearce_rdcoverlap) <- c("Species_abbrv", "Overlap")

summaryData$HR_overlap_pearce <- Data_pearce_rdcoverlap$Overlap[match(summaryData$Species_abbrv,Data_pearce_rdcoverlap$Species_abbrv)]
summaryData$HR_overlap_willems <- Data_willems$HRoverlap[match(summaryData$Species_abbrv,Data_willems$Species_abbrv)]

#Add pop density: pearce, powell2, wrangham
Data_pearce_rdcpop <- Data_pearce[,c(19,7)]
Data_pearce_rdcpop$Species_abbrv <- gsub(" ","_",Data_pearce_rdcpop$Species_abbrv)
Data_pearce_rdcpop <- unique(Data_pearce_rdcpop)
Data_pearce_rdcpop <- aggregate(Data_pearce_rdcpop$PD, by=list(Data_pearce_rdcpop[,1]), FUN=mean)
colnames(Data_pearce_rdcpop) <- c("Species_abbrv", "PD")

summaryData$Pop_density_n_km2_pearce <- Data_pearce_rdcpop$PD[match(summaryData$Species_abbrv,Data_pearce_rdcpop$Species_abbrv)]
summaryData$Pop_density_n_km2_powell2 <- Data_powell2$X21.1_PopulationDensity_n.km2[match(summaryData$Species_abbrv,Data_pearce_rdcoverlap$Species_abbrv)]
summaryData$Pop_density_n_km2_powell2[summaryData$Pop_density_n_km2_powell2=="-999"] <- NA
summaryData$Pop_density_n_km2_wrangham <- Data_wrangham$PD...km.z.[match(summaryData$Species_abbrv,Data_pearce_rdcoverlap$Species_abbrv)]
 
#Calculate the M-index and add the one from Willems
summaryData$Mindex_willems <- Data_willems$Mindex[match(summaryData$Species_abbrv,Data_willems$Species_abbrv)]

#Add HRsize and Daily path length
#HR
summaryData$HR_km2_powell2 <- Data_powell2$X22.1_HomeRange_km2[match(summaryData$Species_abbrv,Data_powell2$Species_abbrv)]
summaryData$HR_km2_powell2[summaryData$HR_km2_powell2==-999] <- NA
summaryData$HR_km2_powell <- Data_powell$HR.size.average[match(summaryData$Species_abbrv,Data_powell$Species_abbrv)]

Data_pearce_rdcHR <- Data_pearce[,c(19,9)]
Data_pearce_rdcHR$Species_abbrv <- gsub(" ","_",Data_pearce_rdcHR$Species_abbrv)
Data_pearce_rdcHR <- unique(Data_pearce_rdcHR)
Data_pearce_rdcHR <- aggregate(Data_pearce_rdcHR$HR, by=list(Data_pearce_rdcHR[,1]), FUN=mean)
colnames(Data_pearce_rdcHR) <- c("Species_abbrv", "HR")
summaryData$HR_km2_pearce <- Data_pearce_rdcHR$HR[match(summaryData$Species_abbrv,Data_pearce_rdcHR$Species_abbrv)]

#DPL
summaryData$DPL_m_powell2 <- Data_powell2$Day.range..m.[match(summaryData$Species_abbrv,Data_willems$Species_abbrv)]
Data_pearce_rdcDPL <- Data_pearce[,c(19,8)]
Data_pearce_rdcDPL$Species_abbrv <- gsub(" ","_",Data_pearce_rdcDPL$Species_abbrv)
Data_pearce_rdcDPL <- unique(Data_pearce_rdcDPL)
Data_pearce_rdcDPL <- aggregate(Data_pearce_rdcDPL$DR, by=list(Data_pearce_rdcDPL[,1]), FUN=mean)
colnames(Data_pearce_rdcDPL) <- c("Species_abbrv", "DPL")
summaryData$DPL_m_pearce <- Data_pearce_rdcDPL$DPL[match(summaryData$Species_abbrv,Data_pearce_rdcDPL$Species_abbrv)]*1000

#powell2: human pop HuPopDen_Mean_n.km2
summaryData$Human_density_n_km2_powell2 <-Data_powell2$HuPopDen_Mean_n.km2[match(summaryData$Species_abbrv,Data_powell2$Species_abbrv)]

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
summaryData$Diet_leaves_willems <- Data_willems$PropLeaves[match(summaryData$Species_abbrv,Data_willems$Species_abbrv)]

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

#Col to correlate
colNumTest <- c(4,4,14,14,17,5,6,6,12,7,  
                20,20,20,19,19,21,
                25,25,26,23,
                29,29,30,34,35)
colNumToCompare <- c(14,18,17,18,18,15,12,16,16,13,
                     19,21,22,21,22,22,
                     26,27,27,24,
                     30,31,31,38,39)

cbind(colnames(summaryData[colNumTest]), colnames(summaryData[colNumToCompare]))
#Vectors to save results
CIlower <- rep(NA, times=length(colNumTest))
CIupper <- rep(NA, times=length(colNumTest))
meanCoeff <- rep(NA, times=length(colNumTest))
N <- rep(NA, times=length(colNumTest))


for (i in 1:length(colNumTest)){
  test <- cor.test(as.numeric(as.character(summaryData[,colNumTest[i]])), 
                   as.numeric(as.character( summaryData[,colNumToCompare[i]])), method="pearson")
  CIlower[i] <- test$conf.int[1]
  CIupper[i] <- test$conf.int[2]
  meanCoeff[i] <- test$estimate[1]
  N[i] <- nrow(summaryData[!is.na(summaryData[,colNumTest[i]])&!is.na(summaryData[,colNumToCompare[i]]),])
}


#Plot

plot(
  x=0, y=0, xlab="", ylab="Coefficient of correlation", 
  xlim=c(0,length(meanCoeff)+1), ylim=c(0,1),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
addGrid(xmin=0, xmax=length(meanCoeff), xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=2, by=0.2), labels=seq(from=0, to=2, by=0.2), las=2, tcl=-0.25)

#Comparison
whatCompared <- c(
  rep("Brain", times=4),
  "Hippocampus",
  "Neocortex",
  "Neocortex",
  "Neocortex",
  "Cerebellum",
  rep("Body", times=7),
  rep("Pop. density", times=3),
  "HR overlap",
  rep("HR size", times=3),
  "Frug.",
  "Fol."
)

#Plot legend of what is compared in coloured rectangles
whereToPlot <- as.data.frame(table(whatCompared))
whereToPlot$loc <- whereToPlot$Freq/2

library(RColorBrewer)
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
errorBars(location=i-0.5, meanPt=meanCoeff[i], refUnit=1, col="black", minValue=0, maxValue=1, upperBarValue=CIupper[i], lowerBarValue=CIlower[i], horiz=FALSE, symmetrical=FALSE)
points(x=i-0.5, y=meanCoeff[i], pch=19, col=colourWhatCompared[which(unique(whatCompared)==whatCompared[i])],
        xpd=TRUE)
text(x=i-0.5, y=CIupper[i]+0.05, labels=N[i], pch=19, col=colourWhatCompared[which(unique(whatCompared)==whatCompared[i])], cex=0.6,xpd=TRUE)

if(i==length(whatCompared)|whatCompared[i]!=whatCompared[i+1]){
  refLoc=refLoc+whereToPlot$loc[whereToPlot[,1]==whatCompared[i]]
  if(whereToPlot$loc[whereToPlot[,1]==whatCompared[i]]<=1){
    segments(x0=refLoc, x1=refLoc, y0=-0.025, y1=-0.075, col=colourWhatCompared[which(unique(whatCompared)==whatCompared[i])], xpd=TRUE)
    text(x=refLoc, y=-0.1, labels=whereToPlot[whereToPlot[,1]==whatCompared[i],1], col=colourWhatCompared[which(unique(whatCompared)==whatCompared[i])], cex=0.5, xpd=TRUE)
  }
  else{
    text(x=refLoc, y=-0.025, labels=whereToPlot[whereToPlot[,1]==whatCompared[i],1], col="black", cex=0.5, xpd=TRUE)
  }
  refLoc=refLoc+whereToPlot$loc[whereToPlot[,1]==whatCompared[i]]#add second time for having complete rectangle
}
}

#Astonishing low robustness: check
# i=18
# summaryData[!is.na(summaryData[,colNumTest[i]])&!is.na(summaryData[,colNumToCompare[i]]),c(colNumTest[i],colNumToCompare[i])]
# 
# i=19
# summaryData[!is.na(summaryData[,colNumTest[i]])&!is.na(summaryData[,colNumToCompare[i]]),c(colNumTest[i],colNumToCompare[i])]
# 
# i=20
# summaryData[!is.na(summaryData[,colNumTest[i]])&!is.na(summaryData[,colNumToCompare[i]]),c(colNumTest[i],colNumToCompare[i])]

#End correlation
#........................................................................................

##############
## Creating co-occurence dataset and merge to covariate dataset
##############

summaryData$Diet_frug_decasien[summaryData$Diet_frug_decasien=="n/a"] <- NA
summaryData$Diet_frug_decasien <- as.numeric(as.character(summaryData$Diet_frug_decasien))

#Summarising Kamilar data: counting species in co-occurrence (from all or same dietary guild)
dataset[12:21]
summaryKamilar <- as.data.frame(matrix(NA, nrow=300, ncol=9))
options(warn=2)

thresholdFrugivory=30
noDataPercentageFrug=0
for(i in 12:length(dataset)){
  datasetTransitory <- get(dataset[i])
  datasetTransitory$Guild <- as.character(datasetTransitory$Guild)
  colWithPresence <- (which(colnames(datasetTransitory)=="Guild")+1):(ncol(datasetTransitory)-1)
  
  for(k in 1:nrow(datasetTransitory)){
    if(datasetTransitory$Species_abbrv[k] %in% summaryData$Species_abbrv[!is.na(summaryData$Diet_frug_powell)]){
      whichGuild <- summaryData$Diet_frug_powell[summaryData$Species_abbrv==datasetTransitory$Species_abbrv[k]]
      if(!is.na(whichGuild)&whichGuild>thresholdFrugivory){
        whichGuild <- "Fruit"
      }
      else{
        whichGuild <- "Leaves"
      }
    }
    else if(datasetTransitory$Species_abbrv[k] %in% summaryData$Species_abbrv[!is.na(summaryData$Diet_frug_decasien)]){
      whichGuild <- summaryData$Diet_frug_decasien[summaryData$Species_abbrv==datasetTransitory$Species_abbrv[k]]
      if(!is.na(whichGuild)&whichGuild>thresholdFrugivory){
        whichGuild <- "Fruit"
      }
      else{
        whichGuild <- "Leaves"
      }
    }
    else if(datasetTransitory$Species_abbrv[k] %in% summaryData$Species_abbrv[!is.na(summaryData$Diet_leaves_willems)]){
      whichGuild <- summaryData$Diet_leaves_willems[summaryData$Species_abbrv==datasetTransitory$Species_abbrv[k]]
      if(!is.na(whichGuild)&whichGuild<100-thresholdFrugivory){
        whichGuild <- "Fruit"
      }
      else{
        whichGuild <- "Leaves"
      }
    }
    else{
      #print("Problem")
      noDataPercentageFrug=noDataPercentageFrug+1
      whichGuild <- as.character(datasetTransitory$Guild[k])
      if(whichGuild=="Fruit-Insects"){
        whichGuild <- "Fruit"
      }
      else{
        
      }
    }
    if(as.character(get(dataset[i])$Guild[k])!="Insects" & as.character(get(dataset[i])$Guild[k])!="Exudates"){
      datasetTransitory$Guild[k] <- as.character(whichGuild)
      
    }
  }
  
  for(j in 1:nrow(datasetTransitory)){
    rowNb=which(is.na(summaryKamilar[,1]))[1]
    summaryKamilar[rowNb,1] <- datasetTransitory[j,ncol(datasetTransitory)]
    summaryKamilar[rowNb,2] <- paste("Site", i, sep="_")
    summaryKamilar[rowNb,3] <- nrow(datasetTransitory)
    colOfInterest <- colWithPresence[which(datasetTransitory[j, colWithPresence]==1)]
    summaryKamilar[rowNb,4] <- length(colOfInterest)
    dfReduced <- as.data.frame(datasetTransitory[-j,colOfInterest])
    summaryKamilar[rowNb,6] <- median(apply(dfReduced, 2, sum))#take the median to be less biased by extreme rare values

    whichGuild <- as.character(datasetTransitory$Guild[j])
    if(whichGuild!=as.character(get(dataset[i])$Guild[j])){
      summaryKamilar[rowNb,8] <- "differentGuild"
      summaryKamilar[rowNb,9] <- as.character(get(dataset[i])$Guild[j])
    }
    else{
      summaryKamilar[rowNb,8] <- "sameGuild"
    }
    
    summaryKamilar[rowNb,5] <- as.character(whichGuild)
    dfReduced <- datasetTransitory[-j,]
    dfReduced <- as.data.frame(dfReduced[as.character(dfReduced$Guild)==as.character(whichGuild),colOfInterest])
    summaryKamilar[rowNb,7] <- median(apply(dfReduced, 2, sum))#take the median to be less biased by extreme rare values
    }
}
summaryKamilar <- summaryKamilar[!is.na(summaryKamilar[,1]),]
colnames(summaryKamilar) <- c("Species_abbrv", "SiteID", "Nspeciessite", "Nlocalities", "DietGuild", "NspeciesCooccurrence", "NspeciesCooccurrence_sameGuild",
"sameGuildAsInitial", "InitialGuild")
str(summaryKamilar)
plot(summaryKamilar[,6], summaryKamilar[,7])

#Some species are present multiple times: pondered average for all values

nbTimePresentInData <- as.data.frame(table(summaryKamilar$Species_abbrv))
nbTimePresentInData <- nbTimePresentInData[nbTimePresentInData$Freq > 1,]
for(a in 1:nrow(nbTimePresentInData)){
  dfTransitory <- summaryKamilar[summaryKamilar$Species_abbrv==nbTimePresentInData$Var1[a],]
  toAdd <-
    unlist(c(
      dfTransitory[1,1],
      paste("Multiple",paste(dfTransitory$SiteID, collapse=";"), sep="-"),
      round((as.numeric(as.character(dfTransitory[1,4]))*as.numeric(as.character(dfTransitory[1,3])) + 
           as.numeric(as.character(dfTransitory[2,4]))*as.numeric(as.character(dfTransitory[2,3])))/(sum(as.numeric(as.character(dfTransitory[,4])))), digit=2),
      sum(as.numeric(as.character(dfTransitory[,4]))),
      dfTransitory[1,5],
      round((as.numeric(as.character(dfTransitory[1,4]))*as.numeric(as.character(dfTransitory[1,6])) + 
         as.numeric(as.character(dfTransitory[2,4]))*as.numeric(as.character(dfTransitory[2,6])))/(sum(as.numeric(as.character(dfTransitory[,4])))), digit=2),
         round((as.numeric(as.character(dfTransitory[1,4]))*as.numeric(as.character(dfTransitory[1,7])) + 
         as.numeric(as.character(dfTransitory[2,4]))*as.numeric(as.character(dfTransitory[2,7])))/(sum(as.numeric(as.character(dfTransitory[,4])))), digit=2),
      dfTransitory[1,8:ncol(summaryKamilar)]
    ))
  names(toAdd) <- colnames(dfTransitory)
  
  #Remove species from dataset
  summaryKamilar <- summaryKamilar[summaryKamilar$Species_abbrv!=nbTimePresentInData$Var1[a],]
  
  #Readd pondered mean
  summaryKamilar <- rbind(summaryKamilar, toAdd)
}

## Adding covariates

#Add HR overlap from grueter. Normally based on Willems 2013 but was complemented for a couple of species (hence not to incude in correlation analysis)
summaryData[is.na(summaryData$HR_overlap_willems)&!is.na(summaryData$Brain_volume_mm3_grueter),]
summaryData$HR_overlap_grueter <- Data_grueter$Overlap[match(summaryData$Species_abbrv, Data_grueter$Species_abbr)]
summaryData$HR_overlap_grueter[!is.na(summaryData$HR_overlap_willems)] <- NA

summaryBrainFrugivory <- matrix(NA, ncol=5, nrow=repetitionNumber)
summaryNeocortexFrugivory <- matrix(NA, ncol=5, nrow=repetitionNumber)
summaryHippocampusFrugivory <- matrix(NA, ncol=5, nrow=repetitionNumber)
summaryCerebellumFrugivory <- matrix(NA, ncol=5, nrow=repetitionNumber)
summaryStriatumFrugivory <- matrix(NA, ncol=5, nrow=repetitionNumber)
summaryMOBFrugivory <- matrix(NA, ncol=5, nrow=repetitionNumber)
summaryOpticTractFrugivory <- matrix(NA, ncol=5, nrow=repetitionNumber)

    ##
    #Random sampling for covariate in case of multiple sources
    ##
    summaryKamilar$Hippocampus <- NA
    summaryKamilar$Cerebellum <- NA
    summaryKamilar$Neocortex <- NA 
    summaryKamilar$Brain <- NA
    summaryKamilar$Family<- NA
    summaryKamilar$OverlapHR <- NA 
    summaryKamilar$Bodymass <- NA 
    
    summaryKamilar$Striatum <- NA
    summaryKamilar$MOB <- NA 
    summaryKamilar$Optic.Tract <- NA 
    
    for(i in 1:nrow(summaryKamilar)){
      
      #Hippocampus
      value <- c(summaryData$Hippocampus_volume_mm3_decasien[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]],
      summaryData$Hippocampus_volume_mm3_todorov[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value==1)){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryKamilar$Hippocampus[i] <- sample(value, 1)
      }
      else{
        summaryKamilar$Hippocampus[i] <- NA
      }
      #Cerebellum
      value <- c(summaryData$Cerebellum_volume_mm3_decasien[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]],
      summaryData$Cerebellum_volume_mm3_powell_mosaic[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]],
      summaryData$Cerebellum_volume_mm3_todorov[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value==1)){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryKamilar$Cerebellum[i] <- sample(value, 1)
      }
      else{
        summaryKamilar$Cerebellum[i] <- NA
      }
      
      #Neocortex
      value <- c(summaryData$Neocortex_volume_mm3_decasien[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]],
      summaryData$Neocortex_volume_mm3_powell_mosaic[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]],
      summaryData$Neocortex_volume_mm3_todorov[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value==1)){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryKamilar$Neocortex[i] <- sample(value, 1)
      }
      else{
        summaryKamilar$Neocortex[i] <- NA
      }
      
      #Brain volume
      value <- c(
      summaryData$Brain_volume_mm3_decasien[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]],
      summaryData$Brain_volume_mm3_powell[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]],
      summaryData$Brain_volume_mm3_todorov[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]],
      summaryData$Brain_volume_mm3_grueter[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value==1)){
        value <- c(value, value)
      }
      if(length(value)>0){
        summaryKamilar$Brain[i] <- sample(value, 1)
      }
      else{
        summaryKamilar$Brain[i] <- NA
      }
      
      #Striatum
      value <- c(
        summaryData$Striatum_volume_mm3_decasien[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value==1)){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryKamilar$Striatum[i] <- sample(value, 1)
      }
      else{
        summaryKamilar$Striatum[i] <- NA
      }
      
      
      #MOB
      value <- c(
        summaryData$MOB_volume_mm3_decasien[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value==1)){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryKamilar$MOB[i] <- sample(value, 1)
      }
      else{
        summaryKamilar$MOB[i] <- NA
      }
      
      #Optic.Tract
      value <- c(
        summaryData$Optic.Tract_volume_mm3_decasien[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value==1)){
        value <- c(value, value)
      }
      
      if(length(value)>0){
        summaryKamilar$Optic.Tract[i] <- sample(value, 1)
      }
      else{
        summaryKamilar$Optic.Tract[i] <- NA
      }
      
      
      #Overlap
      value <- c(summaryData$HR_overlap_willems[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]],
      summaryData$HR_overlap_pearce[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]]/100,
      summaryData$HR_overlap_grueter[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value==1)){
        value <- c(value, value)
      }
      if(length(value)>0){
        summaryKamilar$OverlapHR[i] <- sample(value, 1)
      }
      else{
        summaryKamilar$OverlapHR[i] <- NA
      }
      
      #Bodymass
      value <- c(summaryData$Body_mass_g_decasien[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]],
      summaryData$Body_mass_g_powell[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]],
      summaryData$Body_mass_g_pearce[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]],
      summaryData$Body_mass_g_grueter[summaryData$Species_abbrv==summaryKamilar$Species_abbrv[i]])
      
      value <- value[!is.na(value)]
      #Because there are issues if of length 1
      if(length(value==1)){
        value <- c(value, value)
      }
       
      if(length(value)>0){
        summaryKamilar$Bodymass[i] <- sample(value, 1)
      }
      else{
        summaryKamilar$Bodymass[i] <- NA
      }
      
    }
    
    # summaryKamilar$Hippocampus <- summaryData$Hippocampus_volume_mm3_decasien[match(summaryKamilar$Species_abbrv,summaryData$Species_abbrv)]
    # summaryKamilar$Cerebellum <- summaryData$Cerebellum_volume_mm3_decasien[match(summaryKamilar$Species_abbrv,summaryData$Species_abbrv)]
    # summaryKamilar$Neocortex <- summaryData$Neocortex_volume_mm3_decasien[match(summaryKamilar$Species_abbrv,summaryData$Species_abbrv)]
    # summaryKamilar$Brain <- summaryData$Brain_volume_mm3_powell_decasien[match(summaryKamilar$Species_abbrv,summaryData$Species_abbrv)]
    summaryKamilar$Family<- summaryData$Family[match(summaryKamilar$Species_abbrv,summaryData$Species_abbrv)]
    # summaryKamilar$OverlapHR <- summaryData$HR_overlap_willems[match(summaryKamilar$Species_abbrv,summaryData$Species_abbrv)]
    # summaryKamilar$OverlapHR[is.na(summaryKamilar$OverlapHR)] <- 
    # summaryData$HR_overlap_pearce[match(summaryKamilar$Species_abbrv[is.na(summaryKamilar$OverlapHR)],summaryData$Species_abbrv)]/100
    # summaryKamilar$Bodymass <- summaryData$Body_mass_g_powell[match(summaryKamilar$Species_abbrv,summaryData$Species_abbrv)]
    
    #Transform to numeric variables when adequate
    
    colToNum <- c(3,4,6,7,10,11,12,13,15,16,17,18,19)
    for(a in 1:length(colToNum)){
      summaryKamilar[,colToNum[a]] <- as.numeric(as.character(summaryKamilar[,colToNum[a]]))  
    }
    
    #Quick check sample
    length(summaryKamilar$Hippocampus[!is.na(summaryKamilar$Hippocampus)])
    length(summaryKamilar$Brain[!is.na(summaryKamilar$Brain )])
    
    #End co-occurrence data
    #........................................................................................
    
    #Quick visualization of effect:
    
    # #checking without phylogeny: to delete after just to see if worth it
    # plot(summaryKamilar$NspeciesCooccurrence, summaryKamilar$Hippocampus/summaryKamilar$Brain)
    # 
    # 
    # hist(summaryKamilar$NspeciesCooccurrence_sameGuild/summaryKamilar$NspeciesCooccurrence)
    # 
    # colourVector <- summaryKamilar$DietGuild
    # colourVector[colourVector=="Fruit"] <- "red"
    # colourVector[colourVector=="Leaves"] <- "green"
    # colourVector[colourVector=="Insects"] <- "black"
    # colourVector[colourVector=="Exudates"] <- "brown"
    # 
    # #Same guild co-occurrence and brain area sizes
    # 
    # #If rate
    # plot(summaryKamilar$NspeciesCooccurrence_sameGuild/summaryKamilar$NspeciesCooccurrence,
    #      summaryKamilar$Hippocampus/summaryKamilar$Brain, pch=19, col=colourVector)
    # plot(summaryKamilar$NspeciesCooccurrence_sameGuild/summaryKamilar$NspeciesCooccurrence,
    #      summaryKamilar$Neocortex/summaryKamilar$Brain, pch=19, col=colourVector)
    # plot(summaryKamilar$NspeciesCooccurrence_sameGuild/summaryKamilar$NspeciesCooccurrence,
    #      summaryKamilar$Cerebellum/summaryKamilar$Brain, pch=19, col=colourVector)
    # plot(summaryKamilar$NspeciesCooccurrence_sameGuild/summaryKamilar$NspeciesCooccurrence,
    #      log(summaryKamilar$Brain/summaryKamilar$Bodymass), pch=19, col=colourVector)
    # 
    # #If count
    # plot(summaryKamilar$NspeciesCooccurrence_sameGuild,
    #      summaryKamilar$Hippocampus/summaryKamilar$Brain, pch=19, col=colourVector)
    # plot(summaryKamilar$NspeciesCooccurrence_sameGuild,
    #      summaryKamilar$Neocortex/summaryKamilar$Brain, pch=19, col=colourVector)
    # plot(summaryKamilar$NspeciesCooccurrence_sameGuild,
    #      summaryKamilar$Cerebellum/summaryKamilar$Brain, pch=19, col=colourVector)
    # plot(summaryKamilar$NspeciesCooccurrence_sameGuild,
    #      log(summaryKamilar$Brain/summaryKamilar$Bodymass), pch=19, col=colourVector)
    # 
    # 
    # #Global co-occurence and brain area sizes
    # plot(summaryKamilar$NspeciesCooccurrence,
    #      summaryKamilar$Hippocampus/summaryKamilar$Brain, pch=19, col=colourVector)
    # plot(summaryKamilar$NspeciesCooccurrence,
    #      summaryKamilar$Neocortex/summaryKamilar$Brain, pch=19, col=colourVector)
    # plot(summaryKamilar$NspeciesCooccurrence,
    #      summaryKamilar$Cerebellum/summaryKamilar$Brain, pch=19, col=colourVector)
    # plot(summaryKamilar$NspeciesCooccurrence,
    #      log(summaryKamilar$Brain/summaryKamilar$Bodymass), pch=19, col=colourVector)
    # 
    # 
    # #Overlap and brain area sizes
    # hist(summaryKamilar$Overlap)
    # plot(sqrt(summaryKamilar$Overlap),
    #      summaryKamilar$Hippocampus/summaryKamilar$Brain, pch=19, col=colourVector)
    # plot(sqrt(summaryKamilar$Overlap),
    #      summaryKamilar$Neocortex/summaryKamilar$Brain, pch=19, col=colourVector)
    # plot(sqrt(summaryKamilar$Overlap),
    #      summaryKamilar$Cerebellum/summaryKamilar$Brain, pch=19, col=colourVector)
    # plot(sqrt(summaryKamilar$Overlap),
    #      log(summaryKamilar$Brain/summaryKamilar$Bodymass), pch=19, col=colourVector)
    # 
    # plot(summaryKamilar$NspeciesCooccurrence_sameGuild, summaryKamilar$Hippocampus/summaryKamilar$Brain)
    # 
    options(warn=1)
    
    
    library(ape)
    library(caper)
    library(MCMCglmm)
    library(phytools)
    
    #--------------------------------------------------------------------------------
    
    ##########
    ## STATISTICAL ANALYSES: IS COGNITION PROMOTED BY INTERSPECIES COMPETITION?
    ##########
    
    table(summaryKamilar$Species_abbrv)
    
    #Consensus tree
    phylo_all <-read.nexus("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3.nex")
    phylo_init <- phylo_all
    phylo <- force.ultrametric(tree=phylo_init, method="extend")#method="nnls")
    is.ultrametric(phylo)
    
    plot(phylo_init$edge.length, phylo$edge.length, xlab="Original tree", ylab="Post-chronos tree", pch=20, bty="n")
    abline(a=0, b=1, col="gray", lwd=0.5)
    
    
    nameForPhylo <- as.data.frame(phylo$tip.label)
    colnames(nameForPhylo) <- "Species"
    nameForPhylo$Species <- as.character(nameForPhylo$Species)
    
    library(tidyr)
    options(warn=1)
    nameForPhylo <- separate(nameForPhylo, col=Species, into=c("Name1", "Name2"), sep="_", remove=FALSE)#some have three "_", bot a big deal for now
    nameForPhylo$Species_abbrv <- paste(nameForPhylo$Name1, substr(nameForPhylo$Name2,1,4), sep="_")
    summaryKamilar$SpeciesPhylo <- nameForPhylo$Species[match(summaryKamilar$Species_abbrv,nameForPhylo$Species_abbrv)]
      
    #Check were additional sampling could be worth it
    
    dataForSample <- summaryKamilar[, c(1,20,13,12,10,11,17,18,19,16)]
    dataForSample$Frug1 <- summaryData$Diet_frug_powell[match(dataForSample$SpeciesPhylo,summaryData$SpeciesForPhylogeny)]
    dataForSample$Frug2<- summaryData$Diet_frug_decasien[match(dataForSample$SpeciesPhylo,summaryData$SpeciesForPhylogeny)]
    dataForSample$Fol <-summaryData$Diet_leaves_willems[match(dataForSample$SpeciesPhylo,summaryData$SpeciesForPhylogeny)]
    
    
    dataForSample <- dataForSample[order(dataForSample[,1], decreasing=TRUE),]
    
    library(RColorBrewer)
    pdf(file="sample.pdf", height=30, width=10)
    
    plot(
      x=0, y=0, xlab="", ylab="", cex.sub=1.6,
      xlim=c(-10,ncol(dataForSample)), ylim=c(0,nrow(dataForSample)),
      las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
      xaxt="n",xaxs="i",yaxs="i", yaxt="n")
    
    text(x=rep(-5, times=nrow(dataForSample)), y=1:nrow(dataForSample)-0.5, labels=dataForSample$SpeciesPhylo, xpd=TRUE, cex=0.4)
    text(x=rep(-1, times=nrow(dataForSample)), y=1:nrow(dataForSample)-0.5, labels=dataForSample$Species_abbrv, xpd=TRUE, cex=0.4)
    text(x=2:ncol(dataForSample), y=rep(nrow(dataForSample)+2, times=length(2:ncol(dataForSample))), labels=colnames(dataForSample)[2:ncol(dataForSample)], xpd=TRUE, cex=0.4, srt=45)
    
    for(i in 1:nrow(dataForSample)){
      for(j in 2:ncol(dataForSample)){
        if(!is.na(dataForSample[i,j])){
          rect(
            xleft=j-1,
            xright=j,
            ybottom=i-1,
            ytop=i,
            border=NA,
            col="black"
          )
        } else{
          rect(
            xleft=j-1,
            xright=j,
            ybottom=i-1,
            ytop=i,
            border=NA,
            col="lightgrey"
          )
        }
        if(j==2&is.na(dataForSample[i,2])&!is.na(dataForSample[i,3])){
          points(x=j-0.5, y=i-0.5, pch=19, cex=1, col=brewer.pal(n = 5, name = "Set1")[1])
        }
        if(j==3&!is.na(dataForSample[i,3])&is.na(dataForSample[i,9])){
          points(x=j-0.5, y=i-0.5, pch=19, cex=1, col=brewer.pal(n = 5, name = "Set1")[1])
        } 
        if(j==3&!is.na(dataForSample[i,3])&length(which(!is.na(dataForSample[i,10:12])))==0){
          text(x=j-0.5, y=i-0.5, labels="D", pch=19, cex=1, col="white")
        }
        if(!is.na(dataForSample[i,j])&is.na(dataForSample[i,3])){
          points(x=j-0.5, y=i-0.5, pch=19, cex=1, col=brewer.pal(n = 5, name = "Set1")[1])
        }
      }
    }
    
    source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
    addGrid(xmin=2, xmax=ncol(dataForSample), xintsmall=1, xintbig=1, ymin=0, ymax=nrow(dataForSample), yintsmall=1, yintbig=1, colsmall="white", colbig="white", axisPlot=FALSE)
    
    segments(x0=c(2,8), x1=c(2,8), y0=c(0,0), y1=c(nrow(dataForSample),nrow(dataForSample)), col=brewer.pal(n = 5, name = "Set1")[2], xpd=TRUE)
    segments(x0=c(9,12), x1=c(9,12), y0=c(0,0), y1=c(nrow(dataForSample),nrow(dataForSample)), col=brewer.pal(n = 5, name = "Set1")[3], xpd=TRUE)
    
    dev.off()
    ###########------------------------------------------------------------------------------
    ###PGLS RUN
    
    # dataForBrain <- summaryKamilar[, c(1,2,3,5,6,7,10,13,14,15,16,17)]
    # dataForBrain <- dataForBrain[!is.na(dataForBrain$OverlapHR)&
    #                                !is.na(dataForBrain$Brain)&
    #                                !is.na(dataForBrain$Bodymass),]
    # dataForBrain <- dataForBrain[!is.na(dataForBrain$SpeciesPhylo),]
    # nrow(dataForBrain)#lost one
    # 
    # 
    # comp_data <- comparative.data(phy = phylo, data=dataForBrain[,c(3,5,6,8,10,11,12)],
    #                               names.col = SpeciesPhylo, vcv = TRUE)
    # 
    # comp_data$data$OverlapHR
    # hist(comp_data$data$OverlapHR)
    # hist(sqrt(comp_data$data$OverlapHR))
    # 
    # 
    # comp_data$data$OverlapHR.sqrt <- sqrt(comp_data$data$OverlapHR)
    # 
    # modelBrain <- pgls(formula = Brain ~ Bodymass + OverlapHR.sqrt, data = comp_data, lambda = "ML")
    # 
    # dataForHippocampus <- summaryKamilar[, c(1,2,3,5,6,7,10,13,14,17)]
    # dataForHippocampus <- dataForHippocampus[!is.na(dataForHippocampus$NspeciesCooccurrence_sameGuild)&
    #                      !is.na(dataForHippocampus$Hippocampus)&
    #                      !is.na(dataForHippocampus$Brain),]
    # 
    # library(lme4)
    # library(lmerTest)
    # dataForHippocampus$Hippocampus.ratio <- dataForHippocampus$Hippocampus / dataForHippocampus$Brain
    # dataForHippocampus$Cooccurrence.ratio <- dataForHippocampus$NspeciesCooccurrence_sameGuild / dataForHippocampus$NspeciesCooccurrence
    # test <- lmer(Hippocampus.ratio ~ Cooccurrence.ratio + (1|Family), data=dataForHippocampus)
    # summary(test)
    # 
    # test <- lmer(Hippocampus ~ Brain + NspeciesCooccurrence + NspeciesCooccurrence_sameGuild + (1|Family), data=dataForHippocampus)
    # summary(test)
    # 
    # dataForHippocampus
    # nrow(dataForHippocampus)
    # Data_powell[Data_powell$Species_abbrv=="Pithecia_mona",]
    # dataForHippocampus <- dataForHippocampus[!is.na(dataForHippocampus$SpeciesPhylo),]
    # nrow(dataForHippocampus)#lost one
    # 
    # comp_data <- comparative.data(phy = phylo, data=dataForHippocampus,
    #                               names.col = SpeciesPhylo, vcv = TRUE)
    # modelHippocampus <- pgls(formula = Hippocampus ~ Brain + NspeciesCooccurrence_sameGuild, data = comp_data, lambda = "ML")
    # 
    # 
    # dataForNeocortex <- summaryKamilar[, c(1,2,3,5,6,7,12,13,14,17)]
    # dataForNeocortex <- dataForNeocortex[!is.na(dataForNeocortex$NspeciesCooccurrence_sameGuild)&
    #                                            !is.na(dataForNeocortex$Neocortex)&
    #                                            !is.na(dataForNeocortex$Brain),]
    # dataForNeocortex
    # nrow(dataForNeocortex)
    # Data_powell[Data_powell$Species_abbrv=="Pithecia_mona",]
    # dataForNeocortex <- dataForNeocortex[!is.na(dataForNeocortex$SpeciesPhylo),]
    # nrow(dataForNeocortex)#lost one
    # 
    # comp_data <- comparative.data(phy = phylo, data=dataForNeocortex,
    #                               names.col = SpeciesPhylo, vcv = TRUE)
    # modelNeocortex <- pgls(formula = Neocortex ~ Brain + NspeciesCooccurrence_sameGuild, data = comp_data, lambda = "ML")
    # 
    # 
    # 
    # dataForCerebellum <- summaryKamilar[, c(1,2,3,5,6,7,11,13,14,17)]
    # dataForCerebellum <- dataForCerebellum[!is.na(dataForCerebellum$NspeciesCooccurrence_sameGuild)&
    #                                            !is.na(dataForCerebellum$Cerebellum)&
    #                                            !is.na(dataForCerebellum$Brain),]
    # dataForCerebellum
    # nrow(dataForCerebellum)
    # Data_powell[Data_powell$Species_abbrv=="Pithecia_mona",]
    # dataForCerebellum <- dataForCerebellum[!is.na(dataForCerebellum$SpeciesPhylo),]
    # nrow(dataForCerebellum)#lost one
    # comp_data <- comparative.data(phy = phylo, data=dataForCerebellum,
    #                               names.col = SpeciesPhylo, vcv = TRUE)
    # modelCerebellum <- pgls(formula = Cerebellum ~ Brain + NspeciesCooccurrence_sameGuild, data = comp_data, lambda = "ML")
    # 
    # 
    # 
    # 
    # inv.phylo<-inverseA(phylo,nodes="TIPS",scale=TRUE)
    # prior<-list(G=list(G1=list(V=1,nu=0.02)),R=list(V=1,nu=0.02))
    # model_simple<-MCMCglmm(Hippocampus ~ Brain + NspeciesCooccurrence_sameGuild,random=~phylo,
    #                        family="gaussian",ginverse=list(phylo=inv.phylo$Ainv),prior=prior,
    #                        data=data,nitt=5000000,burnin=1000,thin=500)
    # summary(model_simple)
    ###########------------------------------------------------------------------------------
    
    
    ###############
    # You stopped to re-work on home-ranges
    #What remains to do?
    #You still have to finishing the match and check correlation between all.
    #recalculate the Mindex (so you might increase your dataset)
    
    #Then you have can calculate the mean "brain index" of competitors ?
    #You also have to filter guilds ? Remove insects/exhudates ?
    #Control for diet ?
    
    #you will then have to run a model for Hippo + Neo + Cereb
    
    #You then have to subsample to see robustness to phylogenetics data
    #Analysis of consistency values between dataset (and maybe e.g. brain size, resampling among value of different studies to test robustness)
    
    
    
    
    
    #####################
    ## PHYLOGENY WITH COMPETITION; 
    #following Drury et al.'s approach to see whether one trait (here one brain area size) phylogenetical history is better described if considering competition
    #https://academic.oup.com/sysbio/article/65/4/700/1753588
    #https://journals.plos.org/plosbiology/article?rev=1&id=10.1371/journal.pbio.2003563#pbio.2003563.ref028
    #####################
    
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
    
    library(RPANDA)
    vectorDiet <- dataForBrain_rdc$DietGuild
    names(vectorDiet) <- dataForBrain_rdc$SpeciesPhylo
    library(phytools)
    simmapdiet1 <- make.simmap(tree=phylo, vectorDiet, model="ARD", pi="estimated", nsim=50)#inequal and not symmetrical rate of transition from folivory to frugivory etc...
    simmapdiet2 <- make.simmap(tree=phylo, vectorDiet, model="ARD", pi="equal", nsim=50)#inequal and not symmetrical rate of transition from folivory to frugivory etc...
    simmapdiet3 <- make.simmap(tree=phylo, vectorDiet, model="ARD", pi=c(1,0), nsim=50)#inequal and not symmetrical rate of transition from folivory to frugivory etc...
    # 
    # plot(simmapdiet1[[1]])
    # 
    # #Variability intra quick check
    # plotSimmap(simmapdiet1[[1]])
    # plotSimmap(simmapdiet1[[2]])
    # plotSimmap(simmapdiet1[[3]])
    # plotSimmap(simmapdiet1[[4]])
    # plotSimmap(simmapdiet1[[5]])
    # 
    # #Var inter quick check
    # plot(simmapdiet2)
    # plot(simmapdiet3)
    # dev.off()
    # 
    # vectorGeo <- dataForBrain_rdc$GeoArea
    # names(vectorGeo) <- dataForBrain_rdc$SpeciesPhylo
    # library(phytools)
    # simmapGeo <- make.simmap(tree=phylo, vectorGeo, model="ARD")
    # 
    
    ######
    #Table to save geography data
    ######
    
    #Site for geography
    tableSiteGeo <- as.data.frame(
      cbind(
        c(
          "Site_12",
          "Site_13",
          "Site_14",
          "Site_15",
          "Site_16",
          "Site_17",
          "Site_18",
          "Site_19",
          "Site_20",
          "Site_21"
        ),
        c(
          "Africa_central",  
          "Africa_east",  
          "Africa_west", 
          #"Africa",
          #"Africa",
          #"Africa",
          "Amazonia_central",
          "Amazonia_west",
          "Amazonia_guyana",
          #"Amazonia",
          #"Amazonia",
          #"Amazonia",
          "Madagascar_east",
          "Madagascar_west",
          #"Madagascar",
          #"Madagascar",
          "Asia_islands",
          "Asia_mainland"
        )
      )
    )
    colnames(tableSiteGeo) <- c("Site", "Locality")
    
    #Get all site presence
    newDataForMultiple <- gsub("Multiple-", "",dataForBrain_rdc$SiteID)
    newDataForMultiple <- as.data.frame(newDataForMultiple)
    colnames(newDataForMultiple) <- "Site"
    library(tidyr)
    newDataForMultiple  <- separate(newDataForMultiple, col="Site", into=c("Site1", "Site2", "Site3"), sep=";", remove=FALSE)
    
    #If using all the areas
    head(newDataForMultiple)
    vectorOutputLocalityID <- rep("0000000000", times=nrow(newDataForMultiple))
    
    #Create binary sequence for locality
    for(i in 1:nrow(newDataForMultiple)){
      toChangeTo1 <- which(tableSiteGeo$Site %in% newDataForMultiple[i,])
      for(j in 1:length(toChangeTo1)){
        substr(vectorOutputLocalityID[i], toChangeTo1[j], toChangeTo1[j]) <- "1"
      }
    }
    # 
    # head(newDataForMultiple)
    # vectorOutputLocalityID <- rep("00000", times=nrow(newDataForMultiple))
    # 
    # #Create binary sequence for locality
    # for(i in 1:nrow(newDataForMultiple)){
    #   toChangeTo1 <- which(tableSiteGeo$Site %in% newDataForMultiple[i,])
    #   if(1 %in% toChangeTo1| 2%in% toChangeTo1| 3 %in% toChangeTo1){
    #     toChangeTo1 <- 1
    #   } else if(4 %in% toChangeTo1| 5 %in% toChangeTo1| 6 %in% toChangeTo1){
    #     toChangeTo1 <- 2
    #   }else if(7 %in% toChangeTo1| 8 %in% toChangeTo1){
    #     toChangeTo1 <- 3
    #   }else if(9 %in% toChangeTo1){
    #     toChangeTo1 <- 4
    #   } else{
    #     toChangeTo1 <- 5
    #   }
    #   for(j in 1:length(toChangeTo1)){
    #     substr(vectorOutputLocalityID[i], toChangeTo1, toChangeTo1) <- "1"
    #   }
    # }
    
    vectorOutputLocalityID_brain <- vectorOutputLocalityID
    toFollow <-   cbind(
      dataForBrain_rdc$SpeciesPhylo,
      vectorOutputLocalityID_brain
    )
    toFollow <- as.data.frame(toFollow)
    colnames(toFollow) <- c("Species", "Locality")
    toFollow <- toFollow[order(toFollow$Species),]
    colnames(toFollow) <- NULL
    
    toFollow[,1] <- as.character(toFollow[,1])
    toFollow[,2] <- as.character(toFollow[,2])
    
    toFollow.v <- paste(toFollow[,1], toFollow[,2], sep="\t")
    
    geographyDataTable <- c(
      paste(as.character(length(vectorOutputLocalityID_brain)), as.character(length(unique(tableSiteGeo$Locality))), paste("(", paste(LETTERS[1:length(unique(tableSiteGeo$Locality))], collapse=" "),")",sep=""), sep="\t"),#paste(tableSiteGeo$Locality, collapse=" "),                                                                             ")", sep="")),
      toFollow.v
    )
    
    
    #Test file (in case you want to subsample)
    test <- unique(dataForBrain_rdc$SpeciesPhylo)#sample(unique(dataForBrain_rdc$SpeciesPhylo), 20)#
    phylotest <- drop.tip(phylo,
                          phylo$tip.label[
                            which(phylo$tip.label
                                  %nin%test)]) 
    
    library(BioGeoBEARS)
    write.tree(phylotest, np("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3_rdcbrain.nex"))
    tree <- read.tree("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3_rdcbrain.nex")
    
    geographyDataTableTest <- c(
      geographyDataTable[1],
      geographyDataTable[which(toFollow[,1]%in%test)+1]
    )
    #geographyDataTableTest[1] <- "20\t10\t(A B C D E F G H I J)"
    
    #Save table with a tabulation as necessary
    write.table(x=geographyDataTableTest,
                file="C:/Users/robira/Documents/PhD/Meta_analysis/Cognition_metaanalysis/test_geo.data", 
                sep="\t",col.names=FALSE, row.names=FALSE, na="", quote = FALSE)
    
    
    source("~/PhD/Meta_analysis/Cognition_metaanalysis/Functions.R")
    runBioGeoBearsandBSM_DEC(
    pathToTree="~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3_rdcbrain.nex",
    pathToGeo="test_geo.data",
    maxCooccurrenceGeo=4,
    mapNumber=50,
    coreNumber=8,
    pathToSaveRdata="BioGeoBEARS_brain.Rdata",
    pathToSaveFigureGeo="C:/Users/robira/Documents/PhD/Meta_analysis/Cognition_metaanalysis/Brain_DEC_v1.pdf",
    pathToBSMMapList="BSM_inputs_file.Rdata",
    pathToBSMoutput="BSM_output_file.Rdata",
    pathToFileStockingPdfForBSMAnalysis="C:/Users/robira/Documents/PhD/Meta_analysis/Cognition_metaanalysis/"
    )
    
    save.image("Frugivory_30.RData")
    #Create variable of rinterest
    summaryKamilar$ratioBrain <- summaryKamilar$Brain*1.036*(10**-3)/ summaryKamilar$Bodymass #Following decasien for multiplication by 1.036
    summaryKamilar$EQ <- summaryKamilar$Brain*1.036*(10**-3)/ (0.085*summaryKamilar$Bodymass**0.775) #Following decasien, according to #Jerison, H. J. Evolution of the Brain and Intelligence (Academic, 1973).
    
    
    #Make it having symmetrical (and if possible gaussian distribution, since it seems to be a prerequisite of the analysis)
    hist(summaryKamilar$ratioBrain)
    hist(log(summaryKamilar$ratioBrain))#good enough normal
    #summaryKamilar$ratioBrain.log <- log(summaryKamilar$ratioBrain)
    
    hist(summaryKamilar$EQ)
    #summaryKamilar$ratioBrain.log <- log(summaryKamilar$ratioBrain)
    
    
    load("BSM_output_file.Rdata")
    colnames(summaryKamilar)[colnames(summaryKamilar)=="DietGuild"] <- "Guild"
    resultBrainFrugivory <- runComparisonModelsCompetition(
      simmap=simmapdiet3,
      data=summaryKamilar,
      subgroup="Fruit",
      numberMaps=50,
      trait="EQ",
      tree=tree,
      ana_events_tables=BSM_output$RES_ana_events_tables,
      clado_events_tables=BSM_output$RES_clado_events_tables
      )
    
    #Neocortex
    summaryKamilar$ratioNeocortex <- summaryKamilar$Neocortex/ summaryKamilar$Brain
    hist(summaryKamilar$ratioNeocortex )
    summaryKamilar$ratioNeocortex.log <- log(summaryKamilar$ratioNeocortex)
    hist(summaryKamilar$ratioNeocortex.log)
    
    
    resultNeocortexFrugivory <- runComparisonModelsCompetition(
      simmap=simmapdiet3,
      data=summaryKamilar,
      subgroup="Fruit",
      numberMaps=50,
      trait="ratioNeocortex.log",
      tree=tree,
      ana_events_tables=BSM_output$RES_ana_events_tables,
      clado_events_tables=BSM_output$RES_clado_events_tables
    )
    
    #Hippocampus
    summaryKamilar$ratioHippocampus <- summaryKamilar$Hippocampus/ summaryKamilar$Brain
    hist(summaryKamilar$ratioHippocampus )
    
    resultHippocampusFrugivory <- runComparisonModelsCompetition(
      simmap=simmapdiet3,
      data=summaryKamilar,
      subgroup="Fruit",
      numberMaps=50,
      trait="ratioHippocampus",
      tree=tree,
      ana_events_tables=BSM_output$RES_ana_events_tables,
      clado_events_tables=BSM_output$RES_clado_events_tables
    )
    
    #Cerebellum
    summaryKamilar$ratioCerebellum <- summaryKamilar$Cerebellum/ summaryKamilar$Brain
    hist(summaryKamilar$ratioCerebellum )
    summaryKamilar$ratioCerebellum.log <- log(summaryKamilar$ratioCerebellum)
    hist(summaryKamilar$ratioCerebellum.log)
    
    
    resultCerebellumFrugivory <- runComparisonModelsCompetition(
      simmap=simmapdiet3,
      data=summaryKamilar,
      subgroup="Fruit",
      numberMaps=50,
      trait="ratioCerebellum.log",
      tree=tree,
      ana_events_tables=BSM_output$RES_ana_events_tables,
      clado_events_tables=BSM_output$RES_clado_events_tables
    )
    
    #Striatum
    summaryKamilar$ratioStriatum <- summaryKamilar$Striatum/ summaryKamilar$Brain
    hist(summaryKamilar$ratioStriatum)
    
    resultStriatumFrugivory <- runComparisonModelsCompetition(
      simmap=simmapdiet3,
      data=summaryKamilar,
      subgroup="Fruit",
      numberMaps=50,
      trait="ratioStriatum",
      tree=tree,
      ana_events_tables=BSM_output$RES_ana_events_tables,
      clado_events_tables=BSM_output$RES_clado_events_tables
    )
    
    #MOB
    summaryKamilar$ratioMOB <- summaryKamilar$MOB/ summaryKamilar$Brain
    hist(summaryKamilar$ratioMOB)
    summaryKamilar$ratioMOB.log <- log(summaryKamilar$ratioMOB)
    hist(summaryKamilar$ratioMOB.log)
    
    resultMOBFrugivory <- runComparisonModelsCompetition(
      simmap=simmapdiet3,
      data=summaryKamilar,
      subgroup="Fruit",
      numberMaps=50,
      trait="ratioMOB.log",
      tree=tree,
      ana_events_tables=BSM_output$RES_ana_events_tables,
      clado_events_tables=BSM_output$RES_clado_events_tables
    )
    
    
    #Optic.Tract
    source("~/PhD/Meta_analysis/Cognition_metaanalysis/Functions.R")
    runBioGeoBearsandBSM_DEC(
      pathToTree="~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3_rdcbrain.nex",
      pathToGeo="test_geo.data",
      maxCooccurrenceGeo=4,
      coreNumber=8,
      pathToSaveRdata="BioGeoBEARS_brain.Rdata",
      pathToSaveFigureGeo="C:/Users/robira/Documents/PhD/Meta_analysis/Cognition_metaanalysis/Brain_DEC_v1.pdf",
      pathToBSMMapList="BSM_inputs_file.Rdata",
      pathToBSMoutput="BSM_output_file.Rdata",
      pathToFileStockingPdfForBSMAnalysis="C:/Users/robira/Documents/PhD/Meta_analysis/Cognition_metaanalysis/"
    )
    
    summaryKamilar$ratioOptic.Tract <- summaryKamilar$Optic.Tract/ summaryKamilar$Brain
    hist(summaryKamilar$ratioOptic.Tract)
    s
    resultOptic.TractFrugivory <- runComparisonModelsCompetition(
      simmap=simmapdiet3,
      data=summaryKamilar,
      subgroup="Fruit",
      numberMaps=50,
      trait="ratioOptic.Tract",
      tree=tree,
      ana_events_tables=BSM_output$RES_ana_events_tables,
      clado_events_tables=BSM_output$RES_clado_events_tables
    )
    
    summaryBrainFrugivory[r,] <-
      apply(resultBrainFrugivory[, (ncol(resultBrainFrugivory)-4):ncol(resultBrainFrugivory)], 2, function(x) mean(as.numcharac(x), na.rm=TRUE))
    summaryNeocortexFrugivory[r,] <-
      apply(resultNeocortexFrugivory[, (ncol(resultNeocortexFrugivory)-4):ncol(resultNeocortexFrugivory)], 2, function(x) mean(as.numcharac(x), na.rm=TRUE))
    summaryHippocampusFrugivory[r,] <-
      apply(resultHippocampusFrugivory[, (ncol(resultHippocampusFrugivory)-4):ncol(resultHippocampusFrugivory)], 2, function(x) mean(as.numcharac(x), na.rm=TRUE))
    summaryCerebellumFrugivory[r,] <-
      apply(resultCerebellumFrugivory[, (ncol(resultCerebellumFrugivory)-4):ncol(resultCerebellumFrugivory)], 2, function(x) mean(as.numcharac(x), na.rm=TRUE))
    summaryStriatumFrugivory[r,] <-
      apply(resultStriatumFrugivory[, (ncol(resultStriatumFrugivory)-4):ncol(resultStriatumFrugivory)], 2, function(x) mean(as.numcharac(x), na.rm=TRUE))
    summaryMOBFrugivory[r,] <-
      apply(resultMOBFrugivory[, (ncol(resultMOBFrugivory)-4):ncol(resultMOBFrugivory)], 2, function(x) mean(as.numcharac(x)))
    summaryOpticTractFrugivory[r,] <-
      apply(resultOptic.TractFrugivory[, (ncol(resultOptic.TractFrugivory)-4):ncol(resultOptic.TractFrugivory)], 2, function(x) mean(as.numcharac(x), na.rm=TRUE))
    

save.image("Frugivory_30_rep_50.RData")



models <- c("BM", "OU", "MC", "DDlin", "DDexp")
library(RColorBrewer)
colourModels <- brewer.pal(n = 5, name = "Set1")



# simmap=simmapdiet3
# data=summaryKamilar
# subgroup="Fruit"
# numberMaps=50
# trait="ratioBrain.log"
# tree=tree
# ana_events_tables=BSM_output$RES_ana_events_tables
# clado_events_tables=BSM_output$RES_clado_events_tables

####
## Plot results
####

dev.off()

pdf(file="test_frugivory_30_3.pdf", height=10, width=15)

layout(mat=rbind(c(1,2,3), c(4,5,6)), widths=c(5,5,5), heights=c(5,5))
par(mar=c(4, 4, 2, 1), mgp=c(2, 0.5, 0), xpd=TRUE)
#note: 1= second run for frugivory 20%
#note: _2= first run for frugivory 20%

colNum <- c("darkgray", brewer.pal(n = 5, name = "Set1")[5], brewer.pal(n = 4, name = "Set1"))

models <- c("BM", "OU", "MC", "DDlin", "DDexp")
library(RColorBrewer)
colourModels <- brewer.pal(n = 5, name = "Set1")

# 
# plot(
#   x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
#   xlim=c(0,6), ylim=c(0,1.2),
#   las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
#   xaxt="n",xaxs="i",yaxs="i", yaxt="n")
# 
# source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
# addGrid(xmin=0, xmax=6, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
# axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
# text(x=1:5, y=rep(-0.1, times=5), labels=models, cex=1.6,  xpd=TRUE)
# for(i in 1:5){
#   
#   meanPt <- mean(as.numcharac(resultBrainFrugivory[, ncol(resultBrainFrugivory)-5+i]))
#   sd <- sd(as.numcharac(resultBrainFrugivory[, ncol(resultBrainFrugivory)-5+i]))
#   #sd <- sd/nrow(resultBrainFrugivory) #error not sd
#   errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
#   points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
#   
# }
# 
# #b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
# #add their values below:
# 
# text(x=c(4,5), y=c(-0.2, -0.2),
#      labels=c(
#        paste("b~", round(mean(as.numcharac(resultBrainFrugivory$DDlingeo.b)), digit=3), sep=""),
#        paste("r~", round(mean(as.numcharac(resultBrainFrugivory$DDexpgeo.r)), digit=3), sep="")
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
  xlim=c(0,6), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
addGrid(xmin=0, xmax=6, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:5, y=rep(-0.1, times=5), labels=models, cex=1.6,  xpd=TRUE)
for(i in 1:5){
  
  meanPt <- mean(as.numcharac(resultStriatumFrugivory[, ncol(resultStriatumFrugivory)-5+i]))
  sd <- sd(as.numcharac(resultStriatumFrugivory[, ncol(resultStriatumFrugivory)-5+i]))
  #sd <- sd/nrow(resultStriatumFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
library(plotrix)
draw.circle(x=0.3,y=1.1,0.2, col=colNum[6], border=NA)
text(x=0.3, y=1.1, labels="2", xpd=TRUE, col="white", font=2)


#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(4,5), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(resultStriatumFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(resultStriatumFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)


##------------

#MOB

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,6), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
addGrid(xmin=0, xmax=6, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:5, y=rep(-0.1, times=5), labels=models, cex=1.6,  xpd=TRUE)
for(i in 1:5){
  
  meanPt <- mean(as.numcharac(resultMOBFrugivory[, ncol(resultMOBFrugivory)-5+i]))
  sd <- sd(as.numcharac(resultMOBFrugivory[, ncol(resultMOBFrugivory)-5+i]))
  #sd <- sd/nrow(resultMOBFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
library(plotrix)
draw.circle(x=0.3,y=1.1,0.2, col=colNum[5], border=NA)
text(x=0.3, y=1.1, labels="3", xpd=TRUE, col="white", font=2)

#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(4,5), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(resultMOBFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(resultMOBFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

##------------

##-------------
##optic tract
plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,6), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
addGrid(xmin=0, xmax=6, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:5, y=rep(-0.1, times=5), labels=models, cex=1.6,  xpd=TRUE)
for(i in 1:5){
  
  meanPt <- mean(as.numcharac(resultOptic.TractFrugivory[, ncol(resultOptic.TractFrugivory)-5+i]))
  sd <- sd(as.numcharac(resultOptic.TractFrugivory[, ncol(resultOptic.TractFrugivory)-5+i]))
  #sd <- sd/nrow(resultOptic.TractFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
library(plotrix)
draw.circle(x=0.3,y=1.1,0.2, col="black", border=NA)
text(x=0.3, y=1.1, labels="4", xpd=TRUE, col="white", font=2)

#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(4,5), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(resultOptic.TractFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(resultOptic.TractFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

##------------


#Neocortex

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,6), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
addGrid(xmin=0, xmax=6, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:5, y=rep(-0.1, times=5), labels=models, cex=1.6,  xpd=TRUE)
for(i in 1:5){
  
  meanPt <- mean(as.numcharac(resultNeocortexFrugivory[, ncol(resultNeocortexFrugivory)-5+i]))
  sd <- sd(as.numcharac(resultNeocortexFrugivory[, ncol(resultNeocortexFrugivory)-5+i]))
  #sd <- sd/nrow(resultNeocortexFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
library(plotrix)
draw.circle(x=0.3,y=1.1,0.2, col=colNum[2], border=NA)
text(x=0.3, y=1.1, labels="5", xpd=TRUE, col="white", font=2)


#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(4,5), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(resultNeocortexFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(resultNeocortexFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)


##------------

#Hippocampus


plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,6), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
addGrid(xmin=0, xmax=6, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:5, y=rep(-0.1, times=5), labels=models, cex=1.6,  xpd=TRUE)
for(i in 1:5){
  
  meanPt <- mean(as.numcharac(resultHippocampusFrugivory[, ncol(resultHippocampusFrugivory)-5+i]))
  sd <- sd(as.numcharac(resultHippocampusFrugivory[, ncol(resultHippocampusFrugivory)-5+i]))
  #sd <- sd/nrow(resultHippocampusFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
library(plotrix)
draw.circle(x=0.3,y=1.1,0.2, col=colNum[3], border=NA)
text(x=0.3, y=1.1, labels="6", xpd=TRUE, col="white", font=2)


#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(4,5), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(resultHippocampusFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(resultHippocampusFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

##------------

#Cerebellum 

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,6), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
addGrid(xmin=0, xmax=6, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:5, y=rep(-0.1, times=5), labels=models, cex=1.6,  xpd=TRUE)
for(i in 1:5){
  
  meanPt <- mean(as.numcharac(resultCerebellumFrugivory[, ncol(resultCerebellumFrugivory)-5+i]))
  sd <- sd(as.numcharac(resultCerebellumFrugivory[, ncol(resultCerebellumFrugivory)-5+i]))
  #sd <- sd/nrow(resultCerebellumFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
library(plotrix)
draw.circle(x=0.3,y=1.1,0.2, col=colNum[4], border=NA)
text(x=0.3, y=1.1, labels="7", xpd=TRUE, col="white", font=2)

#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(4,5), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(resultCerebellumFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(resultCerebellumFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

dev.off()








dev.off()
pdf(file="test_frugivory_30_2.pdf", height=10, width=15)

layout(mat=rbind(c(1,2,3), c(4,5,6)), widths=c(5,5,5), heights=c(5,5))
par(mar=c(4, 4, 2, 1), mgp=c(2, 0.5, 0), xpd=TRUE)
#note: 1= second run for frugivory 20%
#note: _2= first run for frugivory 20%

colNum <- c("darkgray", brewer.pal(n = 5, name = "Set1")[5], brewer.pal(n = 4, name = "Set1"))

models <- c("BM", "OU", "MC", "DDlin", "DDexp")
library(RColorBrewer)
colourModels <- brewer.pal(n = 5, name = "Set1")
# 
# 
# plot(
#   x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
#   xlim=c(0,6), ylim=c(0,1.2),
#   las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
#   xaxt="n",xaxs="i",yaxs="i", yaxt="n")
# 
# source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
# addGrid(xmin=0, xmax=6, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
# axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
# text(x=1:5, y=rep(-0.1, times=5), labels=models, cex=1.6,  xpd=TRUE)
# for(i in 1:5){
#   
#   meanPt <- mean(as.numcharac(resultBrainFrugivory[, ncol(resultBrainFrugivory)-5+i]))
#   sd <- sd(as.numcharac(resultBrainFrugivory[, ncol(resultBrainFrugivory)-5+i]))
#   #sd <- sd/nrow(resultBrainFrugivory) #error not sd
#   errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
#   points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
#   
# }
# 
# #b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
# #add their values below:
# 
# text(x=c(4,5), y=c(-0.2, -0.2),
#      labels=c(
#        paste("b~", round(mean(as.numcharac(resultBrainFrugivory$DDlingeo.b)), digit=3), sep=""),
#        paste("r~", round(mean(as.numcharac(resultBrainFrugivory$DDexpgeo.r)), digit=3), sep="")
#      ), xpd=TRUE)
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
  xlim=c(0,6), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
addGrid(xmin=0, xmax=6, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:5, y=rep(-0.1, times=5), labels=models, cex=1.6,  xpd=TRUE)
for(i in 1:5){
  
  meanPt <- mean(as.numcharac(resultStriatumFrugivory[, ncol(resultStriatumFrugivory)-5+i]))
  sd <- sd(as.numcharac(resultStriatumFrugivory[, ncol(resultStriatumFrugivory)-5+i]))
  #sd <- sd/nrow(resultStriatumFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
library(plotrix)
draw.circle(x=0.3,y=1.1,0.2, col=colNum[6], border=NA)
text(x=0.3, y=1.1, labels="2", xpd=TRUE, col="white", font=2)


#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(4,5), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(resultStriatumFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(resultStriatumFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)


##------------

#MOB

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,6), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
addGrid(xmin=0, xmax=6, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:5, y=rep(-0.1, times=5), labels=models, cex=1.6,  xpd=TRUE)
for(i in 1:5){
  
  meanPt <- mean(as.numcharac(resultMOBFrugivory[, ncol(resultMOBFrugivory)-5+i]))
  sd <- sd(as.numcharac(resultMOBFrugivory[, ncol(resultMOBFrugivory)-5+i]))
  #sd <- sd/nrow(resultMOBFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
library(plotrix)
draw.circle(x=0.3,y=1.1,0.2, col=colNum[5], border=NA)
text(x=0.3, y=1.1, labels="3", xpd=TRUE, col="white", font=2)

#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(4,5), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(resultMOBFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(resultMOBFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

##------------


plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,6), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
addGrid(xmin=0, xmax=6, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:5, y=rep(-0.1, times=5), labels=models, cex=1.6,  xpd=TRUE)
for(i in 1:5){
  
  meanPt <- mean(as.numcharac(resultOptic.TractFrugivory[, ncol(resultOptic.TractFrugivory)-5+i]))
  sd <- sd(as.numcharac(resultOptic.TractFrugivory[, ncol(resultOptic.TractFrugivory)-5+i]))
  #sd <- sd/nrow(resultOptic.TractFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
library(plotrix)
draw.circle(x=0.3,y=1.1,0.2, col="black", border=NA)
text(x=0.3, y=1.1, labels="3", xpd=TRUE, col="white", font=2)

#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(4,5), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(resultOptic.TractFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(resultOptic.TractFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

##------------


dev.off()


pdf(file="test_frugivory_30_1.pdf", height=10, width=15)
layout(mat=rbind(c(1,2,3), c(4,5,6)), widths=c(5,5,5), heights=c(5,5))
par(mar=c(4, 4, 2, 1), mgp=c(2, 0.5, 0), xpd=TRUE)
#note: 1= second run for frugivory 20%
#note: _2= first run for frugivory 20%

colNum <- c("darkgray", brewer.pal(n = 5, name = "Set1")[5], brewer.pal(n = 4, name = "Set1"))

models <- c("BM", "OU", "MC", "DDlin", "DDexp")
library(RColorBrewer)
colourModels <- brewer.pal(n = 5, name = "Set1")

#Striatum

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,6), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
addGrid(xmin=0, xmax=6, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:5, y=rep(-0.1, times=5), labels=models, cex=1.6,  xpd=TRUE)
for(i in 1:5){
  
  meanPt <- mean(as.numcharac(resultStriatumFrugivory[, ncol(resultStriatumFrugivory)-5+i]))
  sd <- sd(as.numcharac(resultStriatumFrugivory[, ncol(resultStriatumFrugivory)-5+i]))
  #sd <- sd/nrow(resultStriatumFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}
library(plotrix)
draw.circle(x=0.3,y=1.1,0.2, col=colNum[6], border=NA)
text(x=0.3, y=1.1, labels="3", xpd=TRUE, col="white", font=2)


#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(4,5), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(resultStriatumFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(resultStriatumFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

##----------------------
dev.off()

dev.off()
pdf(file="test_frugivory_30_brain.pdf", height=10, width=15)

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,6), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
addGrid(xmin=0, xmax=6, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:5, y=rep(-0.1, times=5), labels=models, cex=1.6,  xpd=TRUE)
for(i in 1:5){
  
  meanPt <- mean(as.numcharac(resultBrainFrugivory[, ncol(resultBrainFrugivory)-5+i]))
  sd <- sd(as.numcharac(resultBrainFrugivory[, ncol(resultBrainFrugivory)-5+i]))
  #sd <- sd/nrow(resultBrainFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}

#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(4,5), y=c(-0.2, -0.2),
     labels=c(
       paste("b~", round(mean(as.numcharac(resultBrainFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(resultBrainFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

library(plotrix)
draw.circle(x=0.3,y=1.1,0.2, col=colNum[1], border=NA)
text(x=0.3, y=1.1, labels="1", xpd=TRUE, col="white", font=2, cex=3)

##----------------------
dev.off()


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
geoBinary <- as.data.frame(toFollow)
colnames(geoBinary) <- c("SpeciesPhylo", "Loc")
locationSpecies <- geoBinary$Loc[match(speciesLabels, geoBinary$SpeciesPhylo)]
dietSpecies <- summaryKamilar$Guild[match(speciesLabels, summaryKamilar$SpeciesPhylo)]
colLoc <- c(brewer.pal(n = 9, name = "Pastel1"), "gray93")#1:length(unique(tableSiteGeo$Locality))#brewer.pal(n = length(unique(tableSiteGeo$Locality)), name = "Set1")
loc.v <- tableSiteGeo$Locality

toFollow


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
## Fig brain values / circular
###
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
# geoBinary <- as.data.frame(toFollow)
# colnames(geoBinary) <- c("SpeciesPhylo", "Loc")
# locationSpecies <- geoBinary$Loc[match(speciesLabels, geoBinary$SpeciesPhylo)]
# dietSpecies <- summaryKamilar$Guild[match(speciesLabels, summaryKamilar$SpeciesPhylo)]
# colLoc <- c(brewer.pal(n = 9, name = "Pastel1"), "gray93")#1:length(unique(tableSiteGeo$Locality))#brewer.pal(n = length(unique(tableSiteGeo$Locality)), name = "Set1")
# loc.v <- tableSiteGeo$Locality
# #Brain data
# relativeValueBrain <- summaryKamilar$EQ[match(speciesLabels, summaryKamilar$SpeciesPhylo)] - 1#runif(length(speciesLabels), -1, 1)
# relativeValueNeocortex <- scale(summaryKamilar$ratioNeocortex[match(speciesLabels, summaryKamilar$SpeciesPhylo)])
# relativeValueHippocampus <- scale(summaryKamilar$ratioHippocampus[match(speciesLabels, summaryKamilar$SpeciesPhylo)])
# relativeValueCerebellum <- scale(summaryKamilar$Cerebellum[match(speciesLabels, summaryKamilar$SpeciesPhylo)])
# relativeValueStriatum <- scale(summaryKamilar$Striatum[match(speciesLabels, summaryKamilar$SpeciesPhylo)])
# relativeValueMOB <- scale(summaryKamilar$MOB[match(speciesLabels, summaryKamilar$SpeciesPhylo)])
# 
# colourVector <- c("darkgrey", brewer.pal(n = 5, name = "Set1")[1:4])
# colourVectorbis <- c("lightgray", brewer.pal(n = 5, name = "Pastel1")[1:4])
# 
# library(circlize)
# circos.clear()
# circos.par(gap.degree=0, gap.after=0, cell.padding=c(0,0,0,0), track.margin=c(0, 0))
# circos.initialize(speciesLabels, xlim = c(0, 1))
# 
# #Background
# circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
#   circos.rect(0, 0, 1, 1, col=colourVector[1], border=colourVector[1])
# }, track.height = 0.05)
# 
# circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
#   circos.rect(0, 0, 1, 1, col=colourVectorbis[1], border=colourVectorbis[1])
# }, track.height = 0.05)
# 
# #Background
# circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
#   circos.rect(0, 0, 1, 1, col=colourVector[2], border=colourVector[2])
# }, track.height = 0.05)
# 
# circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
#   circos.rect(0, 0, 1, 1, col=colourVectorbis[2], border=colourVectorbis[2])
# }, track.height = 0.05)
# 
# #Background
# circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
#   circos.rect(0, 0, 1, 1, col=colourVector[3], border=colourVector[3])
# }, track.height = 0.05)
# 
# circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
#   circos.rect(0, 0, 1, 1, col=colourVectorbis[3], border=colourVectorbis[3])
# }, track.height = 0.05)
# 
# #Background
# circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
#   circos.rect(0, 0, 1, 1, col=colourVector[4], border=colourVector[4])
# }, track.height = 0.05)
# 
# circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
#   circos.rect(0, 0, 1, 1, col=colourVectorbis[4], border=colourVectorbis[4])
# }, track.height = 0.05)
# 
# #Background
# circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
#   circos.rect(0, 0, 1, 1, col=colourVector[5], border=colourVector[5])
# }, track.height = 0.05)
# 
# circos.track(ylim = c(0, 1), bg.border = NA, panel.fun = function(x, y) {
#   circos.rect(0, 0, 1, 1, col=colourVectorbis[5], border=colourVectorbis[5])
# }, track.height = 0.05)
# 
# library(plotrix)
# #Main circle
# for(i in 1:6){
# draw.circle(x=0,y=0,0.9-(i-1)*0.1, col=NA, border="white")
# }
# 
# #increment of 0.5
# for(i in 1:26){
#   draw.circle(x=0,y=0,1-(i-1)*0.05, col=NA, border="white", lty=2)
# }
# 
# 
# #Value
# #for(i in 1:length(speciesLabels)){
# circos.track(ylim = c(0, 1), bg.border = NA, track.index=1, panel.fun = function(x, y) {
#   i=CELL_META$sector.numeric.index
#   #circos.rect(0, 0, 1, 1, col=colourPositive, border=colourPositive)
#   if(is.na(relativeValueNeocortex[i])){}  else{
#   if(relativeValueNeocortex[i] > 0 & dietSpecies[i]=="Fruit"){
#     circos.points(CELL_META$xcenter, relativeValueNeocortex[i]/absMax, pch=19, col="black")
#     circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueNeocortex[i]/absMax, lty=3)
#   }
#   else if(relativeValueNeocortex[i] > 0 & dietSpecies[i]=="Leaves"){
#     circos.points(CELL_META$xcenter, relativeValueNeocortex[i]/absMax, pch=21, col="black", bg="white")
#     circos.segments(CELL_META$xcenter, 0, CELL_META$xcenter, relativeValueNeocortex[i]/absMax, lty=3)
#   }
#   else{}
#   }
# }, track.height = 0.1)
# 
# circos.track(ylim = c(0, 1), bg.border = NA, track.index=2,  panel.fun = function(x, y) {
#   i=CELL_META$sector.numeric.index
#   if(is.na(relativeValueNeocortex[i])){}  else{
#   #circos.rect(0, 0, 1, 1, col=colourNegative, border=colourNegative)
#   if(relativeValueNeocortex[i] <= 0 & dietSpecies[i]=="Fruit"){
#     circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueNeocortex[i]/absMax, lty=3)
#     circos.points(CELL_META$xcenter, 1 + relativeValueNeocortex[i]/absMax, pch=19, col="black")
#   }
#   else if(relativeValueNeocortex[i] <= 0 & dietSpecies[i]=="Leaves"){
#     circos.segments(CELL_META$xcenter, 1, CELL_META$xcenter, 1 + relativeValueNeocortex[i]/absMax, lty=3)
#     circos.points(CELL_META$xcenter, 1 + relativeValueNeocortex[i]/absMax, pch=21, col="black", bg="white")
#   }
#   else{}
#   }
# }, track.height = 0.1)
# 
# 
# 
# 
# 
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


#from: http://phylo.wikidot.com/biogeobears#script


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



tree <- read.tree("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3_rdcbrain.nex")
plot(tree)
max(node.depth.edgelength(tree))
is.binary(tree)
min(tree$edge.length)
is.ultrametric(tree)
is.rooted(tree)

###
# TREE PATH
trfn = np("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3_rdcbrain.nex")
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

moref("test_geo.data")

#################
# PARAMETERIZATION
#################

# This is the example geography file for Hawaiian Psychotria
# (from Ree & Smith 2008)
geogfn = np("test_geo.data")

# Look at the raw geography text file:
moref(geogfn)

# Look at your geographic range data:
tipranges = getranges_from_LagrangePHYLIP(lgdata_fn="test_geo.data")
tipranges

# Maximum range size observed:
max(rowSums(dfnums_to_numeric(tipranges@df)))

# Set the maximum number of areas any species may occupy; this cannot be larger 
# than the number of areas you set up, but it can be smaller.
max_range_size = 3

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
BioGeoBEARS_run_object$num_cores_to_use = 8
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
resfn = "BioGeoBEARS_brain.Rdata"
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

pdffn = "C:/Users/robira/Documents/PhD/Meta_analysis/Cognition_metaanalysis/Brain_DEC_v1.pdf"#paste0("Psychotria_", model_name, "_v1.pdf")
pdf(pdffn, width=6, height=6)

analysis_titletxt = "Brain size"

# Setup
results_object = res
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.5, statecex=0.5, splitcex=0.5, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it



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
BSM_inputs_fn = "BSM_inputs_file.Rdata"
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
  BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=100, nummaps_goal=50, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
  
  RES_clado_events_tables = BSM_output$RES_clado_events_tables
  RES_ana_events_tables = BSM_output$RES_ana_events_tables
} else {
  # Load previously saved...
  
  # Loads to: RES_clado_events_tables
  load(file="RES_clado_events_tables.Rdata")
  # Loads to: RES_ana_events_tables
  load(file="RES_ana_events_tables.Rdata")
  BSM_output = NULL
  BSM_output$RES_clado_events_tables = RES_clado_events_tables
  BSM_output$RES_ana_events_tables = RES_ana_events_tables
} # END if (runBSMslow == TRUE)

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
max_range_size = max_range_size

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
pdffn = paste0(model_name, "_single_stochastic_map_n1.pdf")
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
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Plot all 50 stochastic maps to PDF
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
pdffn = paste0(model_name, "_", length(clado_events_tables), "BSMs_v1.pdf")
pdf(file=pdffn, width=6, height=6)

nummaps_goal = 50
for (i in 1:nummaps_goal)
{
  clado_events_table = clado_events_tables[[i]]
  analysis_titletxt = paste0(model_name, " - Stochastic Map #", i, "/", nummaps_goal)
  plot_BSM(results_object=res, clado_events_table=clado_events_table, stratified=stratified, analysis_titletxt=analysis_titletxt, addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, include_null_range=include_null_range)
} # END for (i in 1:nummaps_goal)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

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
hist_event_counts(counts_list, pdffn=paste0(model_name, "_histograms_of_event_counts.pdf"))

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


  
  
#######################
## Run the phylogenetic comparison between random evolution (BM, OU) or competitive scenario with MC= exclusion from competitive
#Taxa = mutual exclusion, or DDlin and DDexp, which are positive or negative density dependance to taxa of same "group" (here feeding group)
#of the evolutionary rate, either linearily (lin), or exponentially (exp)
######################


#extracted from: https://journals.plos.org/plosbiology/article?rev=1&id=10.1371/journal.pbio.2003563#pbio.2003563.ref028

#load required R packages
require(RPANDA)
require(geiger)

#isolate data from subgroup of interest (here, frugivores)

group.map<-simmap
data.grouped<-subset(data, DietGuild==subgroup)
mass<-data.grouped[,ncol(data.grouped)-1]
names(mass)<-data.grouped[,ncol(data.grouped)-1]
nc<-name.check(tree,mass)
data.grouped.tree<-drop.tip(tree,nc$tree_not_data)
subdata<-data.grouped[which(data.grouped$SpeciesPhylo%in%data.grouped.tree$tip.label),]

#set up analyses

N=numberMaps #number of stochastic maps to analyze


T=trait #which trait to analyze
res.mat<-matrix(nrow=N,ncol=45)
colnames(res.mat)<-c("subgroup","trait","N","BM.lnL","BM.sig2","BM.z0","BM.AICc","BM.conv","OU.lnL","OU.sig2","OU.alpha","OU.z0","OU.AICc","OU.conv","MCgeo.dietmap","MCgeo.lnL","MCgeo.sig2","MCgeo.S","MCgeo.z0","MCgeo.AICc","MCgeo.conv","DDlingeo.dietmap","DDlingeo.lnL","DDlingeo.sig2","DDlingeo.b","DDlingeo.z0","DDlingeo.AICc","DDlingeo.conv","DDexpgeo.dietmap","DDexpgeo.lnL","DDexpgeo.sig2","DDexpgeo.r","DDexpgeo.z0","DDexpgeo.AICc","DDexpgeo.conv","BM.delaic","OU.delaic","MCgeo.delaic","DDlingeo.delaic","DDexpgeo.delaic","BM.wi","OU.wi","MCgeo.wi","DDlingeo.wi","DDexpgeo.wi")

#fit BM, OU, MC, DDexp, and DDlin models to subgroup trait data

for(i in 1:N){
  group.map2 <-drop.tip.simmap(group.map[[i]],
                             group.map[[i]]$tip.label[which(!group.map[[i]]$tip.label%in%tree$tip.label)])
  j=which(colnames(subdata)==T)
  M<-subdata[,j]
  subtree<-data.grouped.tree
  names(M)<-subdata$SpeciesPhylo
  M<-subset(M,M!='NA')
  
  nc<-name.check(subtree,M)
  if(is.list(nc)){
    subtree<-drop.tip(subtree,nc$tree_not_data)
  }
  
  o2<-fitContinuous(subtree,M,model="BM")
  BM.log_lik<-o2$opt$lnL
  BM.sig2<-o2$opt$sigsq
  BM.z0<-o2$opt$z0
  BM.aicc<-o2$opt$aicc
  BM.conv<-as.numeric(tail(o2$res[,length(o2$res[1,])],n=1))
  
  o3<-fitContinuous(subtree,M,model="OU")
  OU.log_lik<-o3$opt$lnL
  OU.sig2<-o3$opt$sigsq
  OU.alpha<-o3$opt$alpha
  OU.z0<-o3$opt$z0
  OU.aicc<-o3$opt$aicc
  OU.conv<-as.numeric(tail(o3$res[,length(o3$res[1,])],n=1))
  
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
  
  BM.delaic<-BM.aicc-min(BM.aicc,OU.aicc,MCgeo.aicc,DDlingeo.aicc,DDexpgeo.aicc)
  OU.delaic<-OU.aicc-min(BM.aicc,OU.aicc,MCgeo.aicc,DDlingeo.aicc,DDexpgeo.aicc)
  MCgeo.delaic<-MCgeo.aicc-min(BM.aicc,OU.aicc,MCgeo.aicc,DDlingeo.aicc,DDexpgeo.aicc)
  DDlingeo.delaic<-DDlingeo.aicc-min(BM.aicc,OU.aicc,MCgeo.aicc,DDlingeo.aicc,DDexpgeo.aicc)
  DDexpgeo.delaic<-DDexpgeo.aicc-min(BM.aicc,OU.aicc,MCgeo.aicc,DDlingeo.aicc,DDexpgeo.aicc)
  all=sum(exp(-0.5*BM.delaic),exp(-0.5*OU.delaic),exp(-0.5*MCgeo.delaic),exp(-0.5*DDlingeo.delaic),exp(-0.5*DDexpgeo.delaic))
  BM.wi<-exp(-0.5*BM.delaic)/all
  OU.wi<-exp(-0.5*OU.delaic)/all
  MCgeo.wi<-exp(-0.5*MCgeo.delaic)/all
  DDlingeo.wi<-exp(-0.5*DDlingeo.delaic)/all
  DDexpgeo.wi<-exp(-0.5*DDexpgeo.delaic)/all
  
  
  MCgeo.iter <- NA
  DDlingeo.iter <- NA
  DDexpgeo.iter <- NA
  int<-c(subgroup,names(subdata)[j],length(subtree$tip.label),BM.log_lik,BM.sig2,BM.z0,BM.aicc,BM.conv,OU.log_lik,OU.sig2,OU.alpha,OU.z0,OU.aicc,OU.conv,MCgeo.iter,MCgeo.lnL,MCgeo.sig2,MCgeo.S,MCgeo.z0,MCgeo.aicc,MCgeo.conv,DDlingeo.iter,DDlingeo.lnL,DDlingeo.sig2,DDlingeo.b,DDlingeo.z0,DDlingeo.aicc,DDlingeo.conv,DDexpgeo.iter,DDexpgeo.lnL,DDexpgeo.sig2,DDexpgeo.r,DDexpgeo.z0,DDexpgeo.aicc,DDexpgeo.conv,BM.delaic,OU.delaic,MCgeo.delaic,DDlingeo.delaic,DDexpgeo.delaic,BM.wi,OU.wi,MCgeo.wi,DDlingeo.wi,DDexpgeo.wi)
  res.mat[i,]<-int
  
}

res.mat <- as.data.frame(res.mat)
resm <- data.frame(res.mat)
return(resm)





###################
## Old drafts:

#######################################################
# Run DEC+J
#######################################################
BioGeoBEARS_run_object = define_BioGeoBEARS_run()
BioGeoBEARS_run_object$trfn = trfn
BioGeoBEARS_run_object$geogfn = geogfn
BioGeoBEARS_run_object$max_range_size = max_range_size
BioGeoBEARS_run_object$min_branchlength = 0.000001    # Min to treat tip as a direct ancestor (no speciation event)
BioGeoBEARS_run_object$include_null_range = TRUE    # set to FALSE for e.g. DEC* model, DEC*+J, etc.
# (For DEC* and other "*" models, please cite: Massana, Kathryn A.; Beaulieu, 
#  Jeremy M.; Matzke, Nicholas J.; O’Meara, Brian C. (2015). Non-null Effects of 
#  the Null Range in Biogeographic Models: Exploring Parameter Estimation in the 
#  DEC Model. bioRxiv,  http://biorxiv.org/content/early/2015/09/16/026914 )
# Also: search script on "include_null_range" for other places to change

# Set up a time-stratified analysis:
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
BioGeoBEARS_run_object$num_cores_to_use = 6
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

# Set up DEC+J model
# Get the ML parameter values from the 2-parameter nested model
# (this will ensure that the 3-parameter model always does at least as good)
dstart = resDEC$outputs@params_table["d","est"]
estart = resDEC$outputs@params_table["e","est"]
jstart = 0.0001

# Input starting values for d, e
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","init"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["d","est"] = dstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","init"] = estart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["e","est"] = estart

# Add j as a free parameter
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","type"] = "free"
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","init"] = jstart
BioGeoBEARS_run_object$BioGeoBEARS_model_object@params_table["j","est"] = jstart

check_BioGeoBEARS_run(BioGeoBEARS_run_object)

resfn = "Psychotria_DEC+J_M0_unconstrained_v1.Rdata"
runslow = TRUE
if (runslow)
{
  #sourceall("/Dropbox/_njm/__packages/BioGeoBEARS_setup/")
  
  res = bears_optim_run(BioGeoBEARS_run_object)
  res    
  
  save(res, file=resfn)
  
  resDECj = res
} else {
  # Loads to "res"
  load(resfn)
  resDECj = res
}


results_object = res
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)




#######################################################
# PDF plots
#######################################################
pdffn = "Psychotria_DEC_vs_DEC+J_M0_unconstrained_v1.pdf"
pdf(pdffn, width=6, height=6)

#######################################################
# Plot ancestral states - DEC
#######################################################
analysis_titletxt ="BioGeoBEARS DEC on Psychotria M0_unconstrained"

# Setup
results_object = resDEC
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res2 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

#######################################################
# Plot ancestral states - DECJ
#######################################################
analysis_titletxt ="BioGeoBEARS DEC+J on Psychotria M0_unconstrained"

# Setup
results_object = resDECj
scriptdir = np(system.file("extdata/a_scripts", package="BioGeoBEARS"))

# States
res1 = plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="text", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

# Pie chart
plot_BioGeoBEARS_results(results_object, analysis_titletxt, addl_params=list("j"), plotwhat="pie", label.offset=0.45, tipcex=0.7, statecex=0.7, splitcex=0.6, titlecex=0.8, plotsplits=TRUE, cornercoords_loc=scriptdir, include_null_range=TRUE, tr=tr, tipranges=tipranges)

dev.off()  # Turn off PDF
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr) # Plot it



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
BSM_inputs_fn = "BSM_inputs_file.Rdata"
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
  BSM_output = runBSM(res, stochastic_mapping_inputs_list=stochastic_mapping_inputs_list, maxnum_maps_to_try=100, nummaps_goal=50, maxtries_per_branch=40000, save_after_every_try=TRUE, savedir=getwd(), seedval=12345, wait_before_save=0.01)
  
  RES_clado_events_tables = BSM_output$RES_clado_events_tables
  RES_ana_events_tables = BSM_output$RES_ana_events_tables
} else {
  # Load previously saved...
  
  # Loads to: RES_clado_events_tables
  load(file="RES_clado_events_tables.Rdata")
  # Loads to: RES_ana_events_tables
  load(file="RES_ana_events_tables.Rdata")
  BSM_output = NULL
  BSM_output$RES_clado_events_tables = RES_clado_events_tables
  BSM_output$RES_ana_events_tables = RES_ana_events_tables
} # END if (runBSMslow == TRUE)

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
max_range_size = 4

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
pdffn = paste0(model_name, "_single_stochastic_map_n1.pdf")
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
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

#######################################################
# Plot all 50 stochastic maps to PDF
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
pdffn = paste0(model_name, "_", length(clado_events_tables), "BSMs_v1.pdf")
pdf(file=pdffn, width=6, height=6)

nummaps_goal = 50
for (i in 1:nummaps_goal)
{
  clado_events_table = clado_events_tables[[i]]
  analysis_titletxt = paste0(model_name, " - Stochastic Map #", i, "/", nummaps_goal)
  plot_BSM(results_object=res, clado_events_table=clado_events_table, stratified=stratified, analysis_titletxt=analysis_titletxt, addl_params=list("j"), label.offset=0.5, plotwhat="text", cornercoords_loc=scriptdir, root.edge=TRUE, colors_list_for_states=colors_list_for_states, show.tip.label=TRUE, include_null_range=include_null_range)
} # END for (i in 1:nummaps_goal)

dev.off()
cmdstr = paste("open ", pdffn, sep="")
system(cmdstr)

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
hist_event_counts(counts_list, pdffn=paste0(model_name, "_histograms_of_event_counts.pdf"))

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

data(BGB.examples)

# NOT RUN {
data(Anolis.data)
#Create a geography.object with a modified edge matrix
#First, specify which region each branch belonged to:
Anolis.regions<-c(rep("cuba",14),rep("hispaniola",17),"puerto_rico")
Anolis.map<-cbind(Anolis.data$phylo$edge,Anolis.regions)
CreateGeoObject(Anolis.data$phylo,map=Anolis.map)

#Create a geography.object with a make.simmap object
#First, specify which region each branch belonged to:
require(phytools)
geo<-c(rep("cuba",7),rep("hispaniola",9),"puerto_rico")
names(geo)<-Anolis.data$phylo$tip.label
stochastic.map<-phytools::make.simmap(Anolis.data$phylo, 
                                      geo, model="ER", nsim=1)
CreateGeoObject(Anolis.data$phylo,map=stochastic.map)

# }


#Then you need to run the model
#You can include co-occurence per se
#Or co-occurence with guild
#Or you can also have somehow a measure of the "average brain size or part brain size" of the other primate in competition !

#Shall you include overlap and M index for Napoleonian intelligence?
#How to account for pop density = food depletion?






phylo_all <-read.nexus("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3.nex")
phylo_init <- phylo_all
phylo <- force.ultrametric(tree=phylo_init, method="extend")#method="nnls")
is.ultrametric(phylo)

plot(phylo_init$edge.length, phylo$edge.length, xlab="Original tree", ylab="Post-chronos tree", pch=20, bty="n")
abline(a=0, b=1, col="gray", lwd=0.5)

nameForPhylo <- as.data.frame(phylo$tip.label)
colnames(nameForPhylo) <- "Species"
nameForPhylo$Species <- as.character(nameForPhylo$Species)


library(tidyr)
options(warn=1)
nameForPhylo <- separate(nameForPhylo, col=Species, into=c("Name1", "Name2"), sep="_", remove=FALSE)#some have three "_", bot a big deal for now
nameForPhylo$Species_abbrv <- paste(nameForPhylo$Name1, substr(nameForPhylo$Name2,1,4), sep="_")

Data_grueter_test <- Data_grueter
Data_grueter_test$SpeciesPhylo <- nameForPhylo$Species[match(Data_grueter_test$Species_abbrv,nameForPhylo$Species_abbrv)]

Data_grueter_test[is.na(Data_grueter_test$SpeciesPhylo),]#


Data_grueter_test <- Data_grueter_test[!is.na(Data_grueter_test$SpeciesPhylo),]#I do quick removing it I just want to see whether I have same result



comp_data <- comparative.data(phy = phylo, data=Data_grueter_test,
                              names.col = SpeciesPhylo, vcv = TRUE,vcv.dim=3)

comp_data$Overlap.sqrt <- sqrt(as.numeric(as.character(comp_data$Overlap)))
comp_data$Body_mass_g.log <- log(as.numeric(as.character(comp_data$Body_mass_g)))
comp_data$Endocranial_volume_cm3.log <- log(as.numeric(as.character(comp_data$Endocranial_volume_cm3)))

modelBrain <- pgls(formula = Endocranial_volume_cm3 ~ Overlap, data = comp_data, lambda = "ML")


#is there any NA in the match, and are they truly NA or from mismatches?
summaryData[is.na(summaryData$TerritorialityMIndex),]
sort(unique(Data_willems$TaxonName))

#Add intra overlap
Data_pearce_rdcoverlap <- Data_pearce[,c(2,8)]
Data_pearce_rdcoverlap$MSW3sp <- gsub(" ","_",Data_pearce_rdcoverlap$MSW3sp)
Data_pearce_rdcoverlap <- unique(Data_pearce_rdcoverlap)
Data_pearce_rdcoverlap <- aggregate(Data_pearce_rdcoverlap$ov, by=list(Data_pearce_rdcoverlap[,1]), FUN=mean)
colnames(Data_pearce_rdcoverlap) <- c("MSW3sp", "Overlap")


overlapData <- sort(unique(c(Data_grueter$CompleteName,Data_pearce_rdcoverlap$MSW3sp,Data_willems$TaxonName)))#normally grueter and willems should be equal bc grueter is based on willems
overlapData <- as.data.frame(overlapData)
colnames(overlapData) <- "Species"



overlapData$overlapGrueter <- Data_grueter$Overlap[match(overlapData$Species,Data_grueter$CompleteName)]
overlapData$overlapPearce <- Data_pearce_rdcoverlap$Overlap[match(overlapData$Species,Data_pearce_rdcoverlap$MSW3sp)]/100
overlapData$overlapWillems <- Data_willems$HRoverlap[match(overlapData$Species,Data_willems$TaxonName)]


plot(overlapData$overlapWillems, overlapData$overlapGrueter)#Ok but different for one
overlapData[which(overlapData$overlapWillems!=overlapData$overlapGrueter),]


plot(overlapData$overlapWillems, overlapData$overlapPearce)
cor.test(overlapData$overlapWillems, overlapData$overlapPearce, method="pearson")
# Pearson's product-moment correlation
# 
# data:  overlapData$overlapWillems and overlapData$overlapPearce
# t = 9.4279, df = 71, p-value = 3.82e-14
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.6222899 0.8328152
# sample estimates:
# cor 
# 0.7456094 


summaryData$OverlapWillems <-
summaryData$OverlapPearce <-
  

  
  summaryData[which(summaryData$Species %nin% overlapData$Species &
                      summaryData$SpeciesForPhylogeny %nin% overlapData$Species),]



#Checking species number
nrow(Data_powell)
nrow(Data_willems)
nrow(Data_decasien_diet)

Data_willems$TaxonName <- gsub(" ", "_", Data_willems$TaxonName)

'%nin%' <- Negate('%in%')

Data_willems$TaxonName %nin% Data_powell$Species.Name
Data_willems$TaxonName[which(Data_willems$TaxonName %nin% Data_powell$Species.Name)]

Data_willems$TaxonName %nin% Data_powell$Species.name.adjusted.to.10kTrees
Data_willems$TaxonName[which(Data_willems$TaxonName %nin% Data_powell$Species.name.adjusted.to.10kTrees)]


Data_decasien_diet$Taxon %nin% Data_powell$Species.Name
Data_decasien_diet$Taxon[which(Data_decasien_diet$Taxon %nin% Data_powell$Species.Name)]

Data_decasien_diet$Taxon %nin% Data_powell$Species.name.adjusted.to.10kTrees
Data_decasien_diet$Taxon[which(Data_decasien_diet$Taxon %nin% Data_powell$Species.name.adjusted.to.10kTrees)]

###
# Checking resemblance for data
###

###Check for diet
#Creating merged data for diet: todorov is supposedly based on decasien. The other seems to be from different initial b (but with potential redundance between them, not checked)
dietDataTest <- Data_powell[, c(2,25,33)]
colnames(dietDataTest)[2:3] <- c("Fruit_powell", "Leaves_powell")
dietDataTest$Fruit_decasien <- Data_decasien_diet$X..Fruit[match(dietDataTest$Species.name.adjusted.to.10kTree, Data_decasien_diet$Taxon)]
dietDataTest$Fruit_decasien[dietDataTest$Fruit_decasien=="n/a"] <- NA
dietDataTest$Fruit_decasien <- as.numeric(as.character(dietDataTest$Fruit_decasien))
dietDataTest$Leaves_willems <- Data_willems$PropLeaves[match(dietDataTest$Species.name.adjusted.to.10kTree, Data_willems$TaxonName)]*100
dietDataTest$Fruit_todorov <- Data_todorov$X.Fruit[match(dietDataTest$Species.name.adjusted.to.10kTree, Data_todorov$Species)] 


plot(dietDataTest$Fruit_powell, dietDataTest$Fruit_decasien)
cor.test(dietDataTest$Fruit_powell, dietDataTest$Fruit_decasien, method="pearson")
# Pearson's product-moment correlation
# 
# data:  dietDataTest$Fruit_powell and dietDataTest$Fruit_decasien
# t = 10.036, df = 81, p-value = 7.227e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.6299360 0.8273164
# sample estimates:
# cor 
# 0.7444734 

plot(dietDataTest$Fruit_powell, dietDataTest$Fruit_todorov)#It's funny because Todorov took its data on frugivory rate from Decasien, which generally correlates pproximatively.
cor.test(dietDataTest$Fruit_powell, dietDataTest$Fruit_todorov, method="pearson")
# Pearson's product-moment correlation
# 
# data:  dietDataTest$Fruit_powell and dietDataTest$Fruit_todorov
# t = 268440000, df = 32, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 1 1
# sample estimates:
# cor 
# 1 

cor.test(dietDataTest$Fruit_decasien, dietDataTest$Fruit_todorov, method="pearson")
# Pearson's product-moment correlation
# 
# data:  dietDataTest$Fruit_decasien and dietDataTest$Fruit_todorov
# t = 5.7035, df = 28, p-value = 4.076e-06
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
#  0.5066271 0.8649248
# sample estimates:
#      cor 
# 0.733091 

#So todorov certainly did not take from decasien but powell...

plot(dietDataTest$Leaves_powell, dietDataTest$Leaves_willems)
cor.test(dietDataTest$Leaves_powell, dietDataTest$Leaves_willems, method="pearson")
# 
# Pearson's product-moment correlation
# 
# data:  dietDataTest$Leaves_powell and dietDataTest$Leaves_willems
# t = 22.801, df = 103, p-value < 2.2e-16
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.8751652 0.9405576
# sample estimates:
# cor 
# 0.9135866 

#Just by curiosity, check why difference with decasien:
check <- dietDataTest[abs(dietDataTest$Fruit_powell- dietDataTest$Fruit_decasien) > 15,]
check <- check[!is.na(check[,1]),]
check

#A very rough check suggests that Powell data were more carefully taken. I hence use them. Note that this discrepency cam stem from
# 1) The number of study included
# 2) The season on which these studies were based
# 3) The way it was collected (%of observation, %feeding time, %unique observation (i.e. at the individual or group level, considering feeding events/daysetc...a s unit?))

##Check for density
#Pearce, wrangham and powell

head(Data_pearce)
head(Data_powell2)
v1 <- Data_powell$Species.Name %in% Data_powell2$Species
v2 <- Data_powell$Species.name.adjusted.to.10kTree %in% Data_powell2$Species
v1==v2

head(Data_wrangham)

nrow(Data_pearce)

Data_pearce_rdcpop <- Data_pearce[,c(2,5)]
Data_pearce_rdcpop$MSW3sp <- gsub(" ","_",Data_pearce_rdcpop$MSW3sp)
Data_pearce_rdcpop <- unique(Data_pearce_rdcpop)
Data_pearce_rdcpop <- aggregate(Data_pearce_rdcpop$PD, by=list(Data_pearce_rdcpop[,1]), FUN=mean)
colnames(Data_pearce_rdcpop) <- c("MSW3sp", "PD")

Data_powell2$Species <- gsub(" ","_",Data_powell2$MSW05_Binomial)

Data_wrangham$X.1 <- gsub(" ","_", Data_wrangham$X.1)

popDensityData <- as.data.frame(unique(c(Data_pearce_rdcpop$MSW3sp,Data_wrangham$X.1,Data_powell2$Species[Data_powell2$MSW05_Order=="Primates"])))
colnames(popDensityData) <- "Species"
popDensityData$popPearce <- Data_pearce_rdcpop$PD[match(popDensityData$Species,Data_pearce_rdcpop$MSW3sp)]
popDensityData$popWrangham <- Data_wrangham$PD...km.z.[match(popDensityData$Species,Data_wrangham$X.1)]
popDensityData$popPowell <- Data_powell2$X21.1_PopulationDensity_n.km2[match(popDensityData$Species,Data_powell2$Species)]
popDensityData$popPowell[popDensityData$popPowell==-999] <- NA


plot(popDensityData$popPowell,popDensityData$popPearce)
cor.test(popDensityData$popPowell,popDensityData$popPearce, method="pearson")
# 
# Pearson's product-moment correlation
# 
# data:  popDensityData$popPowell and popDensityData$popPearce
# t = 8.0051, df = 73, p-value = 1.377e-11
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.5406555 0.7883470
# sample estimates:
# cor 
# 0.6837194 

plot(popDensityData$popPowell,popDensityData$popWrangham)
cor.test(popDensityData$popPowell,popDensityData$popWrangham, method="pearson")
# Pearson's product-moment correlation
# 
# data:  popDensityData$popPowell and popDensityData$popWrangham
# t = 4.9377, df = 24, p-value = 4.868e-05
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# 0.4448578 0.8605946
# sample estimates:
# cor 
# 0.7098849 
plot(popDensityData$popWrangham,popDensityData$popPearce)
cor.test(popDensityData$popWrangham,popDensityData$popPearce, method="pearson")
# Pearson's product-moment correlation
# 
# data:  popDensityData$popWrangham and popDensityData$popPearce
# t = 0.69929, df = 24, p-value = 0.4911
# alternative hypothesis: true correlation is not equal to 0
# 95 percent confidence interval:
# -0.2602897  0.5012261
# sample estimates:
# cor 
# 0.1413096 

#check which one you are going to use
popDensityData$isInTodorov <- rep("YES", times=nrow(Data_todorov))[match(popDensityData$Species,Data_todorov$Species)]
check <-popDensityData[popDensityData$isInTodorov=="YES",]
check[!is.na(check[,1]),]
missingSpe <- Data_todorov$Species[Data_todorov$Species%nin% popDensityData$Species]

#density is sometimes highly different
#The case of alouatta palliata is illustrative:
# https://onlinelibrary.wiley.com/doi/epdf/10.1002/ajp.20420
#if accounting for island, then the density is from a decade to a hundred








#Missing spe density:
# 
# Pygathrix_nemaeus
# https://Swww.cambridge.org/core/journals/oryx/article/conservation-of-the-redshanked-douc-pygathrix-nemaeus-in-lao-peoples-democratic-republic-density-estimates-based-on-distance-sampling-and-habitat-suitability-modelling/3A826166675625B7E721F548471FE84D
# + ~p60 http://sutir.sut.ac.th:8080/jspui/bitstream/123456789/3714/2/Fulltext.pdf

listSpeciesKamilar <- c()
for(b in 1:10){
  test <- get(paste("Data_kamilar", b, sep=""))
  listSpeciesKamilar <- c(listSpeciesKamilar, as.character(test$Species)) 
}
listSpeciesKamilar <- sort(unique(listSpeciesKamilar))


listSpeciesIn <- c()
listSpeciesSympatry <- c()
for(a in 1:length(Data_todorov$Species)){
for(b in 1:10){
  test <- get(paste("Data_kamilar", b, sep=""))
  if(Data_todorov$Species[a] %in% test$Species){
    listSpeciesIn <- c(listSpeciesIn, as.character(Data_todorov$Species[a]))
    rowInterest <- which(as.character(test$Species) == as.character(Data_todorov$Species[a]))
    colInterest <- which(test[rowInterest, 3:ncol(test)]==1) + 2
    testrdc <- test[-rowInterest,c(1,2,colInterest)]
    rowToDelete <- which(apply(testrdc[3:ncol(testrdc)], 1, sum)==0)
    listSpeciesSympatry <- c(listSpeciesSympatry, as.character(testrdc$Species[-rowToDelete]))
  }
 }
}
listSpeciesIn <- unique(listSpeciesIn)
listSpeciesSympatry <- unique(listSpeciesSympatry)

Data_todorov$Species[which(Data_todorov$Species %nin% listSpeciesIn)]
listSpeciesKamilar



speciesBodySizeMissing <- c()
speciesDietMissing <- c()
speciesPopDensityMissing <- c()


