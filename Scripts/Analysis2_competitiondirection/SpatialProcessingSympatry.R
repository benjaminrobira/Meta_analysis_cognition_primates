####---------------------------------
#### SpatialProcessingSympatry
####---------------------------------

# This script allows processing the IUCN data to investigate (1) number of sympatric species (2) average overlap. It will save the environment but also create an output txt file.

# Note to discuss with Benoit: on a MAD team meeting where Sophie Monsarrat presented her work, she spoke about the Phylocine database that allowed to consider reconstructed ranging data so to avoid biases with recent extinction/ range shrinkage.
# Worth checking?

###Set working directory
setwd("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates")

#Import environment
rm(list=ls())
load("REnvironments/geography_traits.RData")

##--------
#Home made functions
#To source all phylogenetics functions (biogeobears + models of evolution)
source("~/PhD/Meta_analysis/Cognition_metaanalysis/Functions.R")

#My toolkit
source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R")
##--------

##--------
#Libraries

#Phylogeny
library(ape)

#Plot
library(stringr)
library(svMisc)

#Spatial
library(rgdal)
library(raster)
library(rgeos)
library(cleangeo) #to clean it otherwise issues with intersection

##--------

##~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


#diet table
resultsDiet <- matrix(NA, ncol=90, nrow=nrow(summaryData))

frugivoryThresholdVector <- seq(from=10, to=30, by=10)
folivoryThresholdVector <- seq(from=50, to=70, by=10)

for(a in 1:length(frugivoryThresholdVector)){
  
  frugivoryThreshold=frugivoryThresholdVector[a]
  
  for(b in 1:length(folivoryThresholdVector)){
    
    folivoryThreshold=folivoryThresholdVector[b] 
    
    for( r in 1:10){
      for(i in 1:nrow(summaryData)){
        #Frugivory 
        value <- c(summaryData$Diet_frug_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                   summaryData$Diet_frug_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
        
        value <- value[!is.na(value)]
        #Because there are issues if of length 1
        if(length(value==1)){
          value <- c(value, value)
        }
        
        if(length(value)>0){
          summaryData$FrugivoryPercent[i] <- sample(value, 1)
        }
        else{
          summaryData$FrugivoryPercent[i] <- NA
        }
        
        #Folivory 
        value <- c(summaryData$Diet_leaves_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
                   summaryData$Diet_leaves_willems[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
        
        value <- value[!is.na(value)]
        #Because there are issues if of length 1
        if(length(value==1)){
          value <- c(value, value)
        }
        
        if(length(value)>0){
          summaryData$FolivoryPercent[i] <- sample(value, 1)
        }
        else{
          summaryData$FolivoryPercent[i] <- NA
        }
        
        ##--------------
        #Associate the dietary guild
        if(!is.na(as.numeric(as.character(summaryData$FrugivoryPercent[i])))&as.numeric(as.character(summaryData$FrugivoryPercent[i]))>=frugivoryThreshold){
          summaryData$DietaryGuild[i] <- "Fruit"
        } else if(!is.na(as.numeric(as.character(summaryData$FolivoryPercent[i])))&as.numeric(as.character(summaryData$FolivoryPercent[i]))>=folivoryThreshold){
          summaryData$DietaryGuild[i] <- "Leaf"
        }
        else if(!is.na(summaryData$Guild_init_decasien[i])){
          summaryData$DietaryGuild[i] <- summaryData$Guild_init_decasien[i]
          if(summaryData$DietaryGuild[i] =="Om"|summaryData$DietaryGuild[i]=="Frug/Fol"){
            summaryData$DietaryGuild[i] <- "Fruit"
          }
        }
        else {
          summaryData$DietaryGuild[i] <- "Other"
        }
      }
      resultsDiet[,which(is.na(resultsDiet[1,]))[1]] <- summaryData$DietaryGuild
    }
  }
}

##--------------

#check how many different
whichDiet <- apply(resultsDiet, 1, function(x) length(unique(x)))
length(whichDiet[whichDiet>1])/length(whichDiet)#10% diff
length(whichDiet[whichDiet>2])/length(whichDiet)#2% with all possible diets

#how many from those with differences are frugivorous
whichDietWithFrugivorous <- apply(resultsDiet, 1, function(x) if("Fruit"%in%unique(x)){"YES"}else{"NO"})

length(whichDiet[whichDiet>1&whichDietWithFrugivorous=="YES"])/length(whichDiet)#~10% diff -> problematic species are frugivorous...
length(whichDiet[whichDiet>2&whichDietWithFrugivorous=="YES"])/length(whichDiet)#2% with all possible diets


#check if variation with least stringent data only
a=1
b=1
frugivoryThreshold=frugivoryThresholdVector[a]
folivoryThreshold=folivoryThresholdVector[b] 
resultsDiet2 <- matrix(NA, ncol=10, nrow=nrow(summaryData))

for( r in 1:10){
  for(i in 1:nrow(summaryData)){
    #Frugivory 
    value <- c(summaryData$Diet_frug_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
               summaryData$Diet_frug_decasien[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
    
    value <- value[!is.na(value)]
    #Because there are issues if of length 1
    if(length(value==1)){
      value <- c(value, value)
    }
    
    if(length(value)>0){
      summaryData$FrugivoryPercent[i] <- sample(value, 1)
    }
    else{
      summaryData$FrugivoryPercent[i] <- NA
    }
    
    #Folivory 
    value <- c(summaryData$Diet_leaves_powell[summaryData$Species_abbrv==summaryData$Species_abbrv[i]],
               summaryData$Diet_leaves_willems[summaryData$Species_abbrv==summaryData$Species_abbrv[i]])
    
    value <- value[!is.na(value)]
    #Because there are issues if of length 1
    if(length(value==1)){
      value <- c(value, value)
    }
    
    if(length(value)>0){
      summaryData$FolivoryPercent[i] <- sample(value, 1)
    }
    else{
      summaryData$FolivoryPercent[i] <- NA
    }
    
    ##--------------
    #Associate the dietary guild
    if(!is.na(as.numeric(as.character(summaryData$FrugivoryPercent[i])))&as.numeric(as.character(summaryData$FrugivoryPercent[i]))>=frugivoryThreshold){
      summaryData$DietaryGuild[i] <- "Fruit"
    } else if(!is.na(as.numeric(as.character(summaryData$FolivoryPercent[i])))&as.numeric(as.character(summaryData$FolivoryPercent[i]))>=folivoryThreshold){
      summaryData$DietaryGuild[i] <- "Leaf"
    }
    else if(!is.na(summaryData$Guild_init_decasien[i])){
      summaryData$DietaryGuild[i] <- summaryData$Guild_init_decasien[i]
      if(summaryData$DietaryGuild[i] =="Om"|summaryData$DietaryGuild[i]=="Frug/Fol"){
        summaryData$DietaryGuild[i] <- "Fruit"
      }
    }
    else {
      summaryData$DietaryGuild[i] <- "Other"
    }
  }
  resultsDiet2[,which(is.na(resultsDiet2[1,]))[1]] <- summaryData$DietaryGuild
} 

whichDiet <- apply(resultsDiet2, 1, function(x) length(unique(x)))
length(whichDiet[whichDiet>1])/length(whichDiet)#10% diff
#ok good
summaryData[which(whichDiet>1),]

whichIsFrugivorous <- apply(resultsDiet2, 1, function(x) "Fruit" %in% unique(x))
speciesLabelsFrugivorous <- summaryData$SpeciesForPhylogeny[whichIsFrugivorous]

###
## Get geographic localisation for each species
primateSpeciesRange <- readOGR(dsn="T:/IUCN_data_primate/data_0.shp")
primateSpeciesRange <- spTransform(primateSpeciesRange, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))

#Get the species
speciesGeo <- rep(NA, times=length(primateSpeciesRange))
idNumber <- rep(NA, times=length(primateSpeciesRange))
for(i in 1:length(primateSpeciesRange)){
  if(!is.na(primateSpeciesRange[i,]@data$SUBSPECIES[1])){
    speciesGeo[i] <- paste(primateSpeciesRange[i,]@data$BINOMIAL[1], primateSpeciesRange[i,]@data$SUBSPECIES[1], sep="_")
  } else{
    speciesGeo[i] <- primateSpeciesRange[i,]@data$BINOMIAL[1]
  }
  idNumber[i] <- primateSpeciesRange@polygons[[i]]@ID[1]
}
speciesGeo <- gsub(" ", "_",  speciesGeo)

#Which are frugivorous

speciesGeoDiet <- rep("Nonfruit", times=length(speciesGeo))
for(i in 1:length(speciesLabelsFrugivorous)){
  if(substr(speciesLabelsFrugivorous[i],1,14)%in%sapply(vectorSpeciesToCorrect, function(x) substr(x, 1, 14))){
    speciesSubstr <- substr(speciesLabelsFrugivorous[i],1,14)
    substrSpeciesGeo <- sapply(speciesGeo, function (x){substr(x, 1, 14)})#I do it in case there are subsubspecies etc while the phylogenetic tree does not got into such details
    if(length(which(substrSpeciesGeo==speciesSubstr))){
      speciesGeoDiet[which(substrSpeciesGeo==speciesSubstr)] <- "Fruit"
    }
  }
  else{
    lengthString <- str_length(speciesLabelsFrugivorous[i])
    substrSpeciesGeo <- sapply(speciesGeo, function (x){substr(x, 1, lengthString)})#I do it in case there are subsubspecies etc while the phylogenetic tree does not got into such details
    if(length(which(substrSpeciesGeo==speciesLabelsFrugivorous[i]))){
      speciesGeoDiet[which(substrSpeciesGeo==speciesLabelsFrugivorous[i])] <- "Fruit"
    }
  }
}

#Extract name of phylogeny
treeForName <-read.nexus("~/PhD/Meta_analysis/Cognition_metaanalysis/Raw_data/Tree/consensusTree_10kTrees_Primates_Version3.nex")
speciesLabels <-  as.data.frame(treeForName$tip.label)

dataRangePrimate <- matrix(NA, nrow=nrow(speciesLabels), ncol=4)
dataRangePrimate <- as.data.frame(dataRangePrimate )
dataRangePrimate[,1] <- speciesLabels[,1]

vectorSpeciesToCorrect <- c("Eulemur_fulvus", "Eulemur_macaco", "Hapalemur_griseus", "Macaca_nemestrina")


for(i in 1:nrow(speciesLabels)){#for all species for which phylogeny is known
  
  if(speciesLabels[i,1]%in%speciesLabelsFrugivorous){#substract only the frugivorous species
    if(substr(speciesLabels[i,1],1,14)%in%sapply(vectorSpeciesToCorrect, function(x) substr(x, 1, 14))){
      polygonSpeciesId <- which(substr(speciesGeo, 1, 14)==substr(speciesLabels[i,1],1,14))
    }
    else{
      lengthString <- str_length(speciesLabels[i,1])
      substrSpeciesGeo <- sapply(speciesGeo, function (x){substr(x, 1, lengthString)})#I do it in case there are subsubspecies etc while the phylogenetic tree does not got into such details
      polygonSpeciesId <- which(substrSpeciesGeo==speciesLabels[i,1])
    }
    
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
      
      
      #Crop all primate species range based on max/min long/lat so to reduce computational time
      toCrop <- sp::Polygon(cbind(
        t(speciesRangeTransitory@bbox)[c(1,1,2,2),1],
        t(speciesRangeTransitory@bbox)[c(1,2,2,1),2])
      )
      toCrop <- sp::Polygons(list(toCrop ), ID = "A")                      
      toCrop <- sp::SpatialPolygons(list(toCrop))
      
      primateWithoutSpecies <- primateSpeciesRange[-polygonSpeciesId,]
      
      #primateWithoutSpecies <- clgeo_Clean(primateWithoutSpecies)
      primateWithoutSpecies <- spTransform(primateWithoutSpecies, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
      primateWithoutSpecies <- spTransform(primateWithoutSpecies, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
      
      tryCatch(cropedPrimateWithoutSpecies <- crop(primateWithoutSpecies, toCrop), error=function(e) {
        primateWithoutSpecies <- gBuffer(primateWithoutSpecies, byid=TRUE, width=0)
        cropedPrimateWithoutSpecies <- crop(primateWithoutSpecies, toCrop)
      })
      
      if(!is.null(cropedPrimateWithoutSpecies)){
        cropedPrimateWithoutSpecies <- clgeo_Clean(cropedPrimateWithoutSpecies)
        if(!is.null(cropedPrimateWithoutSpecies)){
          
          #Determine the number of species based on the ID
          
          whichSpeciesForComparison <- rep(NA, times=length(cropedPrimateWithoutSpecies))
          for(c in 1:length(cropedPrimateWithoutSpecies)){
            whichSpeciesForComparison[c] <- cropedPrimateWithoutSpecies@polygons[[c]]@ID[1]
          }
          
          #substract species that overlap and that are frugivorous
          speciesForComparison <- speciesGeo[idNumber%in%whichSpeciesForComparison&speciesGeoDiet=="Fruit"]
          
          if(length(speciesForComparison)>0){
            #Do the pairwise overlap taking care of the potential multiple polygons by species
            overlapPair <- rep(NA, times=length(unique(speciesForComparison)))
            
            for(j in 1:length(unique(speciesForComparison))){
              #print(j)
              if(substr(unique(speciesForComparison)[j],1 ,14) %in% sapply(vectorSpeciesToCorrect, function(x) substr(x, 1, 14))){#specific cases as above
                whichToSelect <- which(whichSpeciesForComparison%in%idNumber[which(sapply(speciesGeo, function(x) substr(x, 1, 14))==substr(unique(speciesForComparison)[j], 1 ,14))])
              }
              else{
                whichToSelect <- which(whichSpeciesForComparison%in%idNumber[which(speciesGeo==unique(speciesForComparison)[j], 1 ,14)])
              }
              
              oneSpeciesToCompare <- cropedPrimateWithoutSpecies[whichToSelect,]
              oneSpeciesToCompare <- gBuffer(oneSpeciesToCompare, byid=F, width=0)  
              oneSpeciesToCompare <- clgeo_Clean(oneSpeciesToCompare)
              oneSpeciesToCompare <- spTransform(oneSpeciesToCompare, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
              oneSpeciesToCompare <- spTransform(oneSpeciesToCompare, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
              
              speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
              
              intersectingArea <- gIntersection(speciesRangeTransitory, oneSpeciesToCompare)
              if(!is.null(intersectingArea)){
                intersectingArea <- gBuffer(intersectingArea, byid=F, width=0) 
                if(!is.null(intersectingArea)){
                  intersectingArea  <- clgeo_Clean(intersectingArea)
                  if(!is.null(intersectingArea)){
                    intersectingArea <- spTransform(intersectingArea, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
                    speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
                    overlapPair[j] <- geosphere::areaPolygon(intersectingArea)/geosphere::areaPolygon(speciesRangeTransitory)#counSpecies
                  }
                }
              }
            }
            
            dataRangePrimate[i,2] <- length(overlapPair[!is.na(overlapPair) & overlapPair > 0.10])#number of species within range (overlap > 10%)
            dataRangePrimate[i,3] <- mean(overlapPair[!is.na(overlapPair) & overlapPair > 0.10])#overlap mean (not counting 0, obviously, and those of less than 10%)
            
            speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=merc +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
          } else{
            dataRangePrimate[i,2] <- 0
            dataRangePrimate[i,3] <- 0
          }
          # intersectingArea <- gIntersection(speciesRangeTransitory, cropedPrimateWithoutSpecies)
          # intersectingArea <- gBuffer(intersectingArea, byid=F, width=0)  
          # intersectingArea  <- clgeo_Clean(intersectingArea)
          # intersectingArea <- spTransform(intersectingArea, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
          # speciesRangeTransitory <- spTransform(speciesRangeTransitory, CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))#the initially correct proj
          # dataRangePrimate[i,4] <-  geosphere::areaPolygon(intersectingArea)/geosphere::areaPolygon(speciesRangeTransitory)#overlap total
        }
      }
    }
  }
  #print(i)
  progress(i/nrow(speciesLabels)*100)
}
#I have checked the warnings and it seems everything is working fine
dataRangePrimate <- dataRangePrimate[,-4]
colnames(dataRangePrimate) <- c("Species", "Number_species_cooccurrence", "Overlap_average")
write.table(dataRangePrimate, "dataRangePrimate.txt", col.names=TRUE)

nrow(dataRangePrimate[!is.na(dataRangePrimate[,1])&!is.na(dataRangePrimate[,2]),])
save.image("REnvironments/Data_spatial_primate.RData")
