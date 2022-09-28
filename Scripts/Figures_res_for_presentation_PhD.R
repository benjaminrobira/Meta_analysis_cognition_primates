#############
## Plotting for presentation PhD
#############

#Load environments
load("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/REnvironments/Data_spatial_primate.RData")
load("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/REnvironments/geography_traits_biogeobears.RData")

load("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/REnvironments/PGLSdiversification_withautocorr.RData")
load("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/REnvironments/PGLSdirectionSelection.RData")
load("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/REnvironments/PGLSdiversificationAndSympatry.RData")

#Import own function
source("T:/Saved_PhD/Empirical_analysis/Scripts&Functions/Functions/toolbox.R", local = knitr::knit_global())

#Library
library(plotrix)


###---------------
## LOADING DATA

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
totModelsWorked=c(0,0,0,0,0,0,0)
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
          totModelsWorked[1]=totModelsWorked[1]+1
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_EQ",a,"_",b,"_",c,"_",d,".txt", sep=""))
          summaryEQFrugivory[start:end,] <- toAdd
          totModelsWorked[2]=totModelsWorked[2]+1
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_NeocortexBodymassRaw",a,"_",b,"_",c,"_",d,".txt", sep=""))
          summaryNeocortexFrugivory[start:end,] <- toAdd
          totModelsWorked[3]=totModelsWorked[3]+1
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_HippocampusBodymassRaw",a,"_",b,"_",c,"_",d,".txt", sep=""))
          summaryHippocampusFrugivory[start:end,] <- toAdd
          totModelsWorked[4]=totModelsWorked[4]+1
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_CerebellumBodymassRaw",a,"_",b,"_",c,"_",d,".txt", sep=""))
          summaryCerebellumFrugivory[start:end,] <- toAdd
          totModelsWorked[5]=totModelsWorked[5]+1
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_StriatumBodymassRaw",a,"_",b,"_",c,"_",d,".txt", sep=""))
          summaryStriatumFrugivory[start:end,] <- toAdd
          totModelsWorked[6]=totModelsWorked[6]+1
          }, error=function(e){
            #Do nothing
          }
        )
        
        tryCatch(
          {toAdd <- read.delim(paste("Processed_data/OutputEvolModel/Output_evolutionary_history_MOBBodymassRaw",a,"_",b,"_",c,"_",d,".txt", sep=""))
          summaryMOBFrugivory[start:end,] <- toAdd
          totModelsWorked[7]=totModelsWorked[7]+1
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
library(RColorBrewer)

colNum <-c("darkgrey", brewer.pal(n = 5, name = "Set1")[1:5])

models <- c("BM", "OU", "EB", "MC", expression(DD[italic(lin)]), expression(DD[italic(exp)]))
colourModels <- brewer.pal(n = 6, name = "Set1")
###--------------
#layout(mat=rbind(c(1,2,3), c(4,5,6)), widths=c(5,5,5), heights=c(5,5))

#note: 1= second run for frugivory 20%
#note: _2= first run for frugivory 20%

pdf("Plots/presentationEQ.pdf")
## EQ
par(mar=c(3.75, 4.5, 2, 0), mgp=c(2.75, 0.5, 0), xpd=TRUE)
plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6, cex.lab=2,font.lab=2, cex.lab=1.5, font.lab=2,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

rect(xleft=0, xright=3.5, ybottom=0, ytop=1, col=adjustcolor("royalblue", alpha.f=0.2), border=NA)
rect(xleft=3.5, xright=7, ybottom=0, ytop=1, col=adjustcolor("blue", alpha.f=0.2), border=NA)

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25, cex.axis=1.5, font=2)
segments(x0 = 3.5, y0=0, x1=3.5, y1=1, lty=2)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=2, font=2, xpd=TRUE)

for(i in 1:6){meanPt <- mean(as.numcharac(summaryEQFrugivory[, ncol(summaryEQFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryEQFrugivory[, ncol(summaryEQFrugivory)-6+i]))
  #sd <- sd/nrow(summaryEQFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col="black", cex=2, xpd=TRUE)
}

draw.circle(x=0.3,y=1.1,0.25, col=colNum[1], border=NA)
#text(x=0.3, y=1.1, labels="1", xpd=TRUE, col="white", font=2)
text(x=0.75, y=1.1, labels="EQ", xpd=TRUE, col="black", font=2, cex=2, adj=0)

dev.off()

pdf("Plots/presentationEQEmpty.pdf")
## EQ
par(mar=c(3.75, 4.5, 2, 0), mgp=c(2.75, 0.5, 0), xpd=TRUE)
plot(
  x=0, y=0, xlab="", ylab="AICc weight", cex.sub=1.6, cex.lab=2,font.lab=2, cex.lab=1.5, font.lab=2,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

rect(xleft=0, xright=3.5, ybottom=0, ytop=1, col=adjustcolor("royalblue", alpha.f=0.2), border=NA)
rect(xleft=3.5, xright=7, ybottom=0, ytop=1, col=adjustcolor("blue", alpha.f=0.2), border=NA)

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25, cex.axis=1.5, font=2)
segments(x0 = 3.5, y0=0, x1=3.5, y1=1, lty=2)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=2, font=2, xpd=TRUE)


draw.circle(x=0.3,y=1.1,0.25, col=colNum[1], border=NA)
#text(x=0.3, y=1.1, labels="1", xpd=TRUE, col="white", font=2)
text(x=0.75, y=1.1, labels="EQ", xpd=TRUE, col="black", font=2, cex=2, adj=0)

dev.off()
##-------------

##------------
#Striatum
pdf("Plots/presentationStriatum.pdf")
par(mar=c(3.75, 4.5, 2, 0), mgp=c(2.75, 0.5, 0), xpd=TRUE)

plot(
  x=0, y=0, xlab="",  ylab="AICc weight", cex.sub=2.5, font.lab=2, cex.lab=2.5, font.lab=2,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

rect(xleft=0, xright=3.5, ybottom=0, ytop=1, col=adjustcolor("royalblue", alpha.f=0.2), border=NA)
rect(xleft=3.5, xright=7, ybottom=0, ytop=1, col=adjustcolor("blue", alpha.f=0.2), border=NA)

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25, cex.axis=1.5, font=2)
#segments(x0 = -1, x1 = -1, y0 = 0, y1 = 1, lty = 2, col = colNum[2])
segments(x0 = 3.5, y0=0, x1=3.5, y1=1, lty=2)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=2, font=2, xpd=TRUE)
for(i in 1:6){
  
  meanPt <- mean(as.numcharac(summaryStriatumFrugivory[, ncol(summaryStriatumFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryStriatumFrugivory[, ncol(summaryStriatumFrugivory)-6+i]))
  #sd <- sd/nrow(summaryStriatumFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col="black", cex=2, xpd=TRUE)
  
}
draw.circle(x=0.3,y=1.1,0.25, col=colNum[2], border=NA)
#text(x=0.3, y=1.1, labels="2", xpd=TRUE, col="white", font=2)
text(x=0.75, y=1.1, labels="Striatum ~ SOCIAL", xpd=TRUE, col="black", font=2, cex=2, adj=0)

dev.off()

pdf("Plots/presentationStriatumEmpty.pdf")
par(mar=c(3.75, 4.5, 2, 0), mgp=c(2.75, 0.5, 0), xpd=TRUE)

plot(
  x=0, y=0, xlab="",  ylab="AICc weight", cex.sub=2.5, font.lab=2, cex.lab=2.5, font.lab=2,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

rect(xleft=0, xright=3.5, ybottom=0, ytop=1, col=adjustcolor("royalblue", alpha.f=0.2), border=NA)
rect(xleft=3.5, xright=7, ybottom=0, ytop=1, col=adjustcolor("blue", alpha.f=0.2), border=NA)

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25, cex.axis=1.5, font=2)
#segments(x0 = -1, x1 = -1, y0 = 0, y1 = 1, lty = 2, col = colNum[2])
segments(x0 = 3.5, y0=0, x1=3.5, y1=1, lty=2)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=2, font=2, xpd=TRUE)

draw.circle(x=0.3,y=1.1,0.25, col=colNum[2], border=NA)
#text(x=0.3, y=1.1, labels="2", xpd=TRUE, col="white", font=2)
text(x=0.75, y=1.1, labels="Striatum ~ SOCIAL", xpd=TRUE, col="black", font=2, cex=2, adj=0)

dev.off()
##------------

##-------------
#MOB
pdf("Plots/presentationMOB.pdf")
par(mar=c(3.75, 4.5, 2, 0), mgp=c(2.75, 0.5, 0), xpd=TRUE)

plot(
  x=0, y=0, xlab="",  ylab="AICc weight", cex.sub=2.5, font.lab=2, cex.lab=2.5, font.lab=2,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

rect(xleft=0, xright=3.5, ybottom=0, ytop=1, col=adjustcolor("royalblue", alpha.f=0.2), border=NA)
rect(xleft=3.5, xright=7, ybottom=0, ytop=1, col=adjustcolor("blue", alpha.f=0.2), border=NA)

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25, cex.axis=1.5, font=2)
#segments(x0 = -1, x1 = -1, y0 = 0, y1 = 1, lty = 2, col = colNum[3])
segments(x0 = 3.5, y0=0, x1=3.5, y1=1, lty=2)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=2, font=2, xpd=TRUE)
for(i in 1:6){
  
  meanPt <- mean(as.numcharac(summaryMOBFrugivory[, ncol(summaryMOBFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryMOBFrugivory[, ncol(summaryMOBFrugivory)-6+i]))
  #sd <- sd/nrow(summaryMOBFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col="black", cex=2, xpd=TRUE)
  
}
draw.circle(x=0.3,y=1.1,0.25, col=colNum[3], border=NA)
#text(x=0.3, y=1.1, labels="3", xpd=TRUE, col="white", font=2)
text(x=0.75, y=1.1, labels="MOB", xpd=TRUE, col="black", font=2, cex=2, adj=0)

dev.off()

pdf("Plots/presentationMOBEmpty.pdf")
par(mar=c(3.75, 4.5, 2, 0), mgp=c(2.75, 0.5, 0), xpd=TRUE)

plot(
  x=0, y=0, xlab="",  ylab="AICc weight", cex.sub=2.5, font.lab=2, cex.lab=2.5, font.lab=2,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

rect(xleft=0, xright=3.5, ybottom=0, ytop=1, col=adjustcolor("royalblue", alpha.f=0.2), border=NA)
rect(xleft=3.5, xright=7, ybottom=0, ytop=1, col=adjustcolor("blue", alpha.f=0.2), border=NA)

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25, cex.axis=1.5, font=2)
#segments(x0 = -1, x1 = -1, y0 = 0, y1 = 1, lty = 2, col = colNum[3])
segments(x0 = 3.5, y0=0, x1=3.5, y1=1, lty=2)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=2, font=2, xpd=TRUE)

draw.circle(x=0.3,y=1.1,0.25, col=colNum[3], border=NA)
#text(x=0.3, y=1.1, labels="3", xpd=TRUE, col="white", font=2)
text(x=0.75, y=1.1, labels="MOB", xpd=TRUE, col="black", font=2, cex=2, adj=0)

dev.off()
##------------

##------------
#Hippocampus
pdf("Plots/presentationHippocampus.pdf")
par(mar=c(3.75, 4.5, 2, 0), mgp=c(2.75, 0.5, 0), xpd=TRUE)

plot(
  x=0, y=0, xlab="",  ylab="AICc weight", cex.sub=2.5, font.lab=2, cex.lab=2.5, font.lab=2,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

rect(xleft=0, xright=3.5, ybottom=0, ytop=1, col=adjustcolor("royalblue", alpha.f=0.2), border=NA)
rect(xleft=3.5, xright=7, ybottom=0, ytop=1, col=adjustcolor("blue", alpha.f=0.2), border=NA)

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
segments(x0 = 3.5, y0=0, x1=3.5, y1=1, lty=2)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25, cex.axis=1.5, font=2)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=2, font=2, xpd=TRUE)
for(i in 1:6){
  
  meanPt <- mean(as.numcharac(summaryHippocampusFrugivory[, ncol(summaryHippocampusFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryHippocampusFrugivory[, ncol(summaryHippocampusFrugivory)-6+i]))
  #sd <- sd/nrow(summaryHippocampusFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col="black", cex=2, xpd=TRUE)
  
}
draw.circle(x=0.3,y=1.1,0.25, col=colNum[4], border=NA)
#text(x=0.3, y=1.1, labels="6", xpd=TRUE, col="white", font=2)
text(x=0.75, y=1.1, labels="Hippocampus ~ FORAGING", xpd=TRUE, col="black", font=2, cex=2, adj=0)

dev.off()

pdf("Plots/presentationHippocampusEmpty.pdf")
par(mar=c(3.75, 4.5, 2, 0), mgp=c(2.75, 0.5, 0), xpd=TRUE)

plot(
  x=0, y=0, xlab="",  ylab="AICc weight", cex.sub=2.5, font.lab=2, cex.lab=2.5, font.lab=2,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

rect(xleft=0, xright=3.5, ybottom=0, ytop=1, col=adjustcolor("royalblue", alpha.f=0.2), border=NA)
rect(xleft=3.5, xright=7, ybottom=0, ytop=1, col=adjustcolor("blue", alpha.f=0.2), border=NA)

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
segments(x0 = 3.5, y0=0, x1=3.5, y1=1, lty=2)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25, cex.axis=1.5, font=2)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=2, font=2, xpd=TRUE)

draw.circle(x=0.3,y=1.1,0.25, col=colNum[4], border=NA)
#text(x=0.3, y=1.1, labels="6", xpd=TRUE, col="white", font=2)
text(x=0.75, y=1.1, labels="Hippocampus ~ FORAGING", xpd=TRUE, col="black", font=2, cex=2, adj=0)

dev.off()
##------------

##-------------
#Neocortex
pdf("Plots/presentationNeocortex.pdf")
par(mar=c(3.75, 4.5, 2, 0), mgp=c(2.75, 0.5, 0), xpd=TRUE)

plot(
  x=0, y=0, xlab="",  ylab="AICc weight", cex.sub=2.5, font.lab=2, cex.lab=2.5, font.lab=2,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

rect(xleft=0, xright=3.5, ybottom=0, ytop=1, col=adjustcolor("royalblue", alpha.f=0.2), border=NA)
rect(xleft=3.5, xright=7, ybottom=0, ytop=1, col=adjustcolor("blue", alpha.f=0.2), border=NA)

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25, cex.axis=1.5, font=2)
#segments(x0 = -1, x1 = -1, y0 = 0, y1 = 1, lty = 2, col = colNum[5])
segments(x0 = 3.5, y0=0, x1=3.5, y1=1, lty=2)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=2, font=2, xpd=TRUE)
for(i in 1:6){
  
  meanPt <- mean(as.numcharac(summaryNeocortexFrugivory[, ncol(summaryNeocortexFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryNeocortexFrugivory[, ncol(summaryNeocortexFrugivory)-6+i]))
  #sd <- sd/nrow(summaryNeocortexFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col="black", cex=2, xpd=TRUE)
  
}
draw.circle(x=0.3,y=1.1,0.25, col=colNum[5], border=NA)
#text(x=0.3, y=1.1, labels="5", xpd=TRUE, col="white", font=2)
text(x=0.75, y=1.1, labels="Neocortex ~ INFO. PROCESSING", xpd=TRUE, col="black", font=2, cex=2, adj=0)

dev.off()

pdf("Plots/presentationNeocortexEmpty.pdf")
par(mar=c(3.75, 4.5, 2, 0), mgp=c(2.75, 0.5, 0), xpd=TRUE)
plot(
  x=0, y=0, xlab="",  ylab="AICc weight", cex.sub=2.5, font.lab=2, cex.lab=2.5, font.lab=2,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

rect(xleft=0, xright=3.5, ybottom=0, ytop=1, col=adjustcolor("royalblue", alpha.f=0.2), border=NA)
rect(xleft=3.5, xright=7, ybottom=0, ytop=1, col=adjustcolor("blue", alpha.f=0.2), border=NA)

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25, cex.axis=1.5, font=2)
segments(x0 = 3.5, y0=0, x1=3.5, y1=1, lty=2)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=2, font=2, xpd=TRUE)

draw.circle(x=0.3,y=1.1,0.25, col=colNum[5], border=NA)
#text(x=0.3, y=1.1, labels="5", xpd=TRUE, col="white", font=2)
text(x=0.75, y=1.1, labels="Neocortex ~ INFO. PROCESSING", xpd=TRUE, col="black", font=2, cex=2, adj=0)

dev.off()
##-------------

##-------------
#Cerebellum 
pdf("Plots/presentationCerebellum.pdf")
par(mar=c(3.75, 4.5, 2, 0), mgp=c(2.75, 0.5, 0), xpd=TRUE)

plot(
  x=0, y=0, xlab="",  ylab="AICc weight", cex.sub=2.5, font.lab=2, cex.lab=2.5, font.lab=2,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

rect(xleft=0, xright=3.5, ybottom=0, ytop=1, col=adjustcolor("royalblue", alpha.f=0.2), border=NA)
rect(xleft=3.5, xright=7, ybottom=0, ytop=1, col=adjustcolor("blue", alpha.f=0.2), border=NA)

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25, cex.axis=1.5, font=2)
#segments(x0 = -1, x1 = -1, y0 = 0, y1 = 1, lty = 2, col = colNum[6])$
segments(x0 = 3.5, y0=0, x1=3.5, y1=1, lty=2)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=2, font=2, xpd=TRUE)
for(i in 1:6){
  
  meanPt <- mean(as.numcharac(summaryCerebellumFrugivory[, ncol(summaryCerebellumFrugivory)-6+i]))
  sd <- sd(as.numcharac(summaryCerebellumFrugivory[, ncol(summaryCerebellumFrugivory)-6+i]))
  #sd <- sd/nrow(summaryCerebellumFrugivory) #error not sd
  errorBars(location=i, meanPt=meanPt, barValue=sd, refUnit=1, col="black", minValue=0, maxValue=1, horiz=FALSE)
  points(x=i, y=meanPt, pch=19, col="black", cex=2, xpd=TRUE)
  
}
draw.circle(x=0.3,y=1.1,0.25, col=colNum[5], border=NA)#PUT SAME COLOUR THAN NEOCORTEX FOR PRESENTATION
#text(x=0.3, y=1.1, labels="7", xpd=TRUE, col="white", font=2)
text(x=0.75, y=1.1, labels="Cerebellum ~ INFO. PROCESSING", xpd=TRUE, col="black", font=2, cex=2, adj=0)

dev.off()

pdf("Plots/presentationCerebellumEmpty.pdf")
par(mar=c(3.75, 4.5, 2, 0), mgp=c(2.75, 0.5, 0), xpd=TRUE)
plot(
  x=0, y=0, xlab="",  ylab="AICc weight", cex.sub=2.5, font.lab=2, cex.lab=2.5, font.lab=2,
  xlim=c(0,7), ylim=c(0,1.2),
  las=1, type="n", tcl=-0.25, frame.plot=FALSE, 
  xaxt="n",xaxs="i",yaxs="i", yaxt="n")

rect(xleft=0, xright=3.5, ybottom=0, ytop=1, col=adjustcolor("royalblue", alpha.f=0.2), border=NA)
rect(xleft=3.5, xright=7, ybottom=0, ytop=1, col=adjustcolor("blue", alpha.f=0.2), border=NA)

addGrid(xmin=0, xmax=7, xintsmall=0.25, xintbig=1, ymin=0, ymax=1, yintsmall=0.05, yintbig=0.2, axisPlot=FALSE)
axis(side=2, at=seq(from=0, to=1, by=0.2), labels=seq(from=0, to=1, by=0.2), las=2, tcl=-0.25, cex.axis=1.5, font=2)
segments(x0 = 3.5, y0=0, x1=3.5, y1=1, lty=2)
text(x=1:6, y=rep(-0.1, times=6), labels=models, cex=2, font=2, xpd=TRUE)

draw.circle(x=0.3,y=1.1,0.25, col=colNum[5], border=NA)
#text(x=0.3, y=1.1, labels="7", xpd=TRUE, col="white", font=2)#PUT SAME COLOUR THAN NEOCORTEX FOR PRESENTATION
text(x=0.75, y=1.1, labels="Cerebellum ~ INFO. PROCESSING", xpd=TRUE, col="black", font=2, cex=2, adj=0)
dev.off()

###----------------------