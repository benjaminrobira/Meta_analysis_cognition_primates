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

##~~~~~~~~~~~~~

repetition=2*2*2*10#length(frugivoryThresholdVector)*length(folivoryThresholdVector)*length(geographicThresholdVector)*randomSampling


##---------
## Load the data

summaryBodyMassFrugivory <- as.data.frame(matrix(NA, nrow=10*(repetition+1), ncol=53))
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
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel_BS/Output_evolutionary_history_BodyMass_",a,"_",b,"_",c,"_",d,".txt", sep=""))
          summaryBodyMassFrugivory[start:end,] <- toAdd
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel_BS/Output_evolutionary_history_NeocortexBrainSize",a,"_",b,"_",c,"_",d,".txt", sep=""))
          summaryNeocortexFrugivory[start:end,] <- toAdd
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel_BS/Output_evolutionary_history_HippocampusBrainSize",a,"_",b,"_",c,"_",d,".txt", sep=""))
          summaryHippocampusFrugivory[start:end,] <- toAdd
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel_BS/Output_evolutionary_history_CerebellumBrainSize",a,"_",b,"_",c,"_",d,".txt", sep=""))
          summaryCerebellumFrugivory[start:end,] <- toAdd
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel_BS/Output_evolutionary_history_StriatumBrainSize",a,"_",b,"_",c,"_",d,".txt", sep=""))
          summaryStriatumFrugivory[start:end,] <- toAdd
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel_BS/Output_evolutionary_history_MOBBrainSize",a,"_",b,"_",c,"_",d,".txt", sep=""))
          summaryMOBFrugivory[start:end,] <- toAdd
          }, error=function(e){
            #Do nothing
          }
        )
      }
    }
  }
}

summaryBodyMassFrugivory <- summaryBodyMassFrugivory[!is.na(summaryBodyMassFrugivory[,1]),]
summaryNeocortexFrugivory <- summaryNeocortexFrugivory[!is.na(summaryNeocortexFrugivory[,1]),]
summaryHippocampusFrugivory <- summaryHippocampusFrugivory[!is.na(summaryHippocampusFrugivory[,1]),]
summaryCerebellumFrugivory <- summaryCerebellumFrugivory[!is.na(summaryCerebellumFrugivory[,1]),]
summaryStriatumFrugivory <- summaryStriatumFrugivory[!is.na(summaryStriatumFrugivory[,1]),]
summaryMOBFrugivory <- summaryMOBFrugivory[!is.na(summaryMOBFrugivory[,1]),]

colnames(summaryBodyMassFrugivory) <- colnames(toAdd)
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
# ## BodyMass
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
#   meanPt <- mean(as.numcharac(summaryBodyMassFrugivory[, ncol(summaryBodyMassFrugivory)-6+i]))
#   sd <- sd(as.numcharac(summaryBodyMassFrugivory[, ncol(summaryBodyMassFrugivory)-6+i]))
#   #sd <- sd/nrow(summaryBodyMassFrugivory) #error not sd
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
#        paste("r~", round(mean(as.numcharac(summaryBodyMassFrugivory$DDlingeo.b)), digit=3), sep=""),
#        paste("r~", round(mean(as.numcharac(summaryBodyMassFrugivory$DDexpgeo.r)), digit=3), sep="")
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

## BodyMass

plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=1,  xpd=TRUE)
for(i in 1:6){
  meanPt <- mean(as.numcharac(summaryBodyMassFrugivory[, ncol(summaryBodyMassFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryBodyMassFrugivory[, ncol(summaryBodyMassFrugivory)-6+i]))
  #sd <- sd/nrow(summaryBodyMassFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col=colourModels[i], xpd=TRUE)
  
}

#b and r are the rate for density dependance (DD) of the speciation rate. If positive, positive DD, otherwise, negative.
#add their values below:

text(x=c(5,6), y=c(-0.2, -0.2),
     labels=c(
       paste("r~", round(mean(as.numcharac(summaryBodyMassFrugivory$DDlingeo.b)), digit=3), sep=""),
       paste("r~", round(mean(as.numcharac(summaryBodyMassFrugivory$DDexpgeo.r)), digit=3), sep="")
     ), xpd=TRUE)

draw.circle(x=0.3,y=1.1,0.2, col=colNum[1], border=NA)
text(x=0.3, y=1.1, labels="1", xpd=TRUE, col="white", font=2)
text(x=3.5, y=1.1, labels="BodyMass", xpd=TRUE, col="black", font=2, cex=1.5)

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

