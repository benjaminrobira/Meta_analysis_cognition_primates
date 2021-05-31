#########################
## TOOLBOX
#########################

##What is it?

#This toolbox gathers many handmade functions (of my own, or burrowed, and potentially modified from the Internet!) that were created because repeatedly used. 
#Note that I have plenty of others linked to movement analyses specifically (in case you need).
#How each function work is described, please enjoy!


## OH NO! :( AN ERROR!
#Found a mistake or an error? I am sad (really, this is my highest fright), and I truly hope this is no big deal. 
#Please contact me at benjamin.robira@normalesup.org so we can solve it together!

##--------------------------------------------------
#Function list:

##----
###Data processing
##----
#firstup
#as.numcharac
#'%nin%'
#getmode 
#cbind.fill

##----
###P-value related
##----
#pvalueToText
#pvalueRound 
#textTestOneVar
#textEstOneVar
#roundIntelligent

##----
###Linear models
##----
#reportGamSmooth
#mainResults
#stabResults 
#diagnostics.plot.dharma

##----
###Graphic related
##----
#addImg
#pastellize 
#addGrid
#addCorner 
#errorBars
#addLabel
#patternLayer
#patternLayerCoordinates

##----
### R markdown
##----
# RmdWords
# citeR
##--------------------------------------------------

##################
## Data processing
##################

##-------------
# firstup
##-------------

#Capitalize the first letter of a string

firstup <- function(x) {
  substr(x, 1, 1) <- toupper(substr(x, 1, 1))
  x
}

#--

##-------------
# as.numcharac
##-------------

#This function allows to combine as.numeric and as.character

as.numcharac <- function(x){
	return(as.numeric(as.character(x)))
}

#--


##-------------
# %nin%
##-------------

#This function allows to select what is not included in a given variable


'%nin%' <- Negate('%in%')

#--

##-------------
# getmode
##-------------

#This function calculates the mode (the value the most represented)

# Create the function.
getmode <- function(
v #a numerical vector
) {
   uniqv <- unique(v)
   return(uniqv[which.max(tabulate(match(v, uniqv)))])
}

#--

##-------------
# Concatenate vector of different size (fill with NA)
##-------------


cbind.fill <- function(...){
    nm <- list(...) 
    nm <- lapply(nm, as.matrix)
    n <- max(sapply(nm, nrow)) 
    do.call(cbind, lapply(nm, function (x) 
        rbind(x, matrix(NA, n-nrow(x), ncol(x))))) 
}


#--

##-------------
# replaceNAtable
##-------------

replaceNAtable <- function(input, replaceWith){

input <- apply(input, 2, function(x){
  toChange <- x
  toChange[is.na(toChange)] <- replaceWith
  toChange <- gsub(" ", "", toChange)
  return(toChange)
  })
 return(input)
}
#--

##################
## P-value related
##################

##-------------
# pvalueToText
##-------------

#This function will transform the p-value into the classical symbols "*" "**" "***" "n.s.". Default is 0.05, 0.01, 0.001. > 0.05

pvalueToText <- function(
p, #The pvalue as numeric
thresholdVerySign=0.001,
thresholdIntermediateSign=0.01,
thresholdWeaklySign=0.05
){
	if(p<=thresholdVerySign){
		return("***");
	}
	else if(p<=thresholdIntermediateSign){
		return("**");
	}
	else if (p<=thresholdWeaklySign){
		return("*");
	}
	else{
		return("n.s.");
	}
}

#--

##-------------
# pvalueRound
##-------------

#This function will return as character p=value of the pvalue. In case it's below the treshold for very significant, in order to avoid very small pvalue it will be displayed as "p<threshold". 


pvalueRound <- function(
p,
thresholdVerySign=0.001,
text=TRUE
){
	if(p<=thresholdVerySign){
		if(text){
			return(paste("p<", thresholdVerySign, sep=""))
		}else{
			return(paste("<", thresholdVerySign, sep=""))
		}
	}
	else{
		if(text){
			return(paste("p=", round(p, digit=3), sep=""))
		}else{
			return(round(p, digit=3))
		}
	}
}

#--

# test <- matrix(NA, nrow=100, ncol=3)
# test[,1] <- runif(100,0,1)
# test[,2] <- runif(100,0,1)
# test[,3] <- rep(c("yes", "no"), times=50)

# test <- as.data.frame(test)
# test


# test <- as.data.frame(test)
# colnames(test) <- c("a", "b", "c")
# test[,1] <- as.numeric(as.character(test[,1]))
# test[,2] <- as.numeric(as.character(test[,2]))

# library(betareg)
# model <- betareg(a ~b + c, data=test)



 # model <- glm(a ~b + c, data=test, family="poisson")
 
 
 # res <- betareg(persistanceVar.crt ~ onTrail + group + morningAfternoon + Season +#differenceAltitude +
 # travelledDistance, data=dataForLinearModel[!is.na(dataForLinearModel$onTrail),])
# null <- betareg(persistanceVar.crt ~ group + morningAfternoon + Season +#differenceAltitude +
 # travelledDistance, data=dataForLinearModel[!is.na(dataForLinearModel$onTrail),])
# library(lmtest)

# lmtest::lrtest(res, null) 
# coeffmodelPersistence <- lmtest::coeftest(res)
# coeffmodelPersistence 


##-------------
# textModelFullNull
##-------------

# This functions will allow to directly write (e.g. for Rmarkdown documents) the statistics, df and p-value of a linear model full vs null comparison.

# comparison: the dataframe of the lmtest::lrtest(null, full)
# Note: for incompatibility, the chisq term has to be inserted prior the function... I am working on it!

textModelFullNull <- function(comparison, digit=2){
  output <- paste(
    #as.character("$\chi$"),
    "=",
    round(comparison[2,4], digit=digit),
    ", df=",
    comparison[2,3],
    ", ",
    pvalueRound(comparison[2,5]),
    sep=""
  )
  return(output)
}

##-------------
# textTestOneVar
##-------------

# This functions will allow to directly write (e.g. for Rmarkdown documents) the statistics and p-value for a predictor test.
#It necessitates:
# -coeffLine: a vector containing the statistics and the p-value (note: use that of drop1 for multiple levels!), and the df (for df version)
# -statistics: the "statistics used: z, t, chisq, F, deviance etc... that you want to be displayed

textTestOneVar <- function(coeffLine, digit=2, statistics){
  output <- paste(
    paste(statistics, "=", sep=""),
	ifelse(abs(coeffLine[1])<1/(10**digit),format(coeffLine[1], digits=digit+1, scientific = T),round(coeffLine[1], digit=digit)),
    ", ",
    pvalueRound(coeffLine[2]),
    sep=""
  )
  return(output)
}

textTestOneVarWithdf <- function(coeffLine, digit=2, statistics){
  output <- paste(
	paste(statistics, "=", sep=""),
	ifelse(abs(coeffLine[1])<1/(10**digit),format(coeffLine[1], digits=digit+1, scientific = T),round(coeffLine[1], digit=digit)),
    ", ",
	"df=",
	coeffLine[3],
	", ",
    pvalueRound(coeffLine[2]),
    sep=""
  )
  return(output)
}

##-------------
# textEstOneVar 
##-------------

# This functions will allow to directly write (e.g. for Rmarkdown documents) the statistics and p-value for a predictor test.
#It necessitates:
# -confintLine: a vector containing the estimate, the lower border of the confint, the upper border of the confint in that order.

textEstOneVar <- function(confintLine, digit=2){
  output <- paste(
    "est.=",
    ifelse(abs(confintLine[1])<1/(10**digit),format(confintLine[1], digits=digit+1, scientific = T),round(confintLine[1], digit=digit)),
    ", CI95%=[",
    ifelse(abs(confintLine[2])<1/(10**digit),format(confintLine[2], digits=digit+1, scientific = T),round(confintLine[2], digit=digit)),
    ",",
    ifelse(abs(confintLine[3])<1/(10**digit),format(confintLine[3], digits=digit+1, scientific = T),round(confintLine[3], digit=digit)),
    "]",
    sep=""
  )
  return(output)
}  

##-------------
# roundIntelligent
##-------------

#To round in a smart way (with scientific writing or not): the rounding will occur for the value before the power of ten in the scientific writing
roundIntelligent <- function(x, digit=2){
	ifelse(abs(x)<1/(10**digit),format(x, digits=digit+1, scientific = T),round(x, digit=digit))
}

##------------



##############
## Linear models
##############

##-------------
# reportGamSmooth: Report the result of smoothed terms for GAMs 
##-------------

#The input vector should be the line displaying result of the smoothed terms. This can be accessed if fitting the gam with the gam function of the mgcv package
#using summary(model)$s.table

reportGamSmooth <- function(vectorSmoothTermResults, digit=2){
  paste(
    "edf=",
    ifelse(abs(vectorSmoothTermResults[1])<1/(10**digit),format(vectorSmoothTermResults[1], digits=digit+1, scientific = T),round(vectorSmoothTermResults[1], digit=digit)),
    ", ",
    "Ref.df=", 
    vectorSmoothTermResults[2],
    ", ",
    "F=",
    ifelse(abs(vectorSmoothTermResults[3])<1/(10**digit),format(vectorSmoothTermResults[3], digits=digit+1, scientific = T),round(vectorSmoothTermResults[3], digit=digit)),
    ", ",
    pvalueRound(vectorSmoothTermResults[4]),
    sep=""
  )
}


##-------------
# mainResults
##-------------

#Extract the results of a (generalized) linear model

mainResults <- function(
model,
modelName
){

	coeffModel <- as.data.frame(summary(model)$coefficients)
	confintModel <- as.data.frame(confint(model))
	vectorVarNb <- nrow(coeffModel)
	obsNumber <- nrow(model$model)
	family <- family(model)

	if(family=="gaussian"){
	drop1Model <- as.data.frame(drop1(model, test="F"))
	}
	else{
	drop1Model <- as.data.frame(drop1(model, test="Chisq"))
	}

	output <- matrix(NA, ncol=10, nrow=vectorVarNb)
	colnames(output) <- c("Model", "Adj. R²","Variables", "Est.", "Sd", "Lower CI (95%)", "Upper CI (95%)", #"RSS", 
	"Stat.", "Df", "p-value")

	output[,1] <- paste(modelName, " (N=", obsNumber,")", sep="")

	library(rsq)
	output[,2] <- rsq(model)
	output[,3] <- rownames(coeffModel)
	output[,4] <- coeffModel[,1]
	output[,5] <- coeffModel[,2]
	output[,6] <- confintModel[,1]
	output[,7] <- confintModel[,2]

	#Determine which variable is factor

	output[1,8] <- drop1Model[1,5]
	output[1,9] <- drop1Model[1,1]
	output[1,10] <- drop1Model[1,6]

	counterFactor=0
	for(i in 2:ncol(model$model)){
		if(!is.factor(model$model[,i])){
			output[i+counterFactor,8] <- drop1Model[i,5]
			output[i+counterFactor,9] <- drop1Model[i,1]
			output[i+counterFactor,10] <- drop1Model[i,6]
		}
		else{
			for(j in 1:(length(levels(model$model[,i]))-1)){
				output[i+counterFactor+j-1,8] <- drop1Model[i,5]
				output[i+counterFactor+j-1,9] <- drop1Model[i,1]
				output[i+counterFactor+j-1,10] <- drop1Model[i,6]
			}
			counterFactor=counterFactor+length(levels(model$model[,i]))-2
		}
	}
	return(output)
}
##------------


##-------------
# stabResults
##-------------

#Extract the stability of a (generalized) linear model

stabResults <- function(
model,
modelName
){
	vectorVarNb <- nrow(coeffModel)
	obsNumber <- nrow(model$model)
	family <- family(model)

	dfBetasmodel <- round(cbind(coefficients(model), coefficients(model)+
				  t(apply(X=dfbeta(model), MARGIN=2, FUN=range))), 5)

	output <- matrix(NA, ncol=10, nrow=vectorVarNb)
	colnames(output) <- c("Model","N", "VIF", "dffits", "lev", "CD", "Disp.","Variable", "dfbeta (min)", "dfbeta (max)")


	output[,1] <- modelName

	library(rsq)
	output[,2] <- obsNumber
	output[,3] <- max(vif(model))
	output[,4] <- max(dffits(model))
	output[,5] <- max(as.vector(influence(model)$hat))
	output[,6] <- max(cooks.distance(model))
	source("Empirical_analysis/Scripts&Functions/Functions/diagnostics.R")
	if(family=="poisson"){output[,7] <- overdisp.test(model)[4]}
	
	output[,8] <-  rownames(summary(model)$coefficients)
	output[,9] <- dfBetasmodel[,1]
	output[,10] <- dfBetasmodel[,2]
	
	return(output)
}
	
##---------------

##-------------
# Hist residuals, QQplot, and Residuals vs fitted values (homogeneity)
##-------------

diagnostics.plot.dharma <-function(mod.res, col=grey(level=0.25, alpha=0.5), breaks.histo=20){
  old.par = par(no.readonly = TRUE)
  par(mfrow=c(2, 2))
  par(mar=c(3, 3, 3, 0.5))
  hist(residuals(mod.res), probability=T, xlab="", ylab="", main="", breaks=breaks.histo)
  mtext(text="Histogram of residuals", side=3, line=0)
  x=seq(min(residuals(mod.res)), max(residuals(mod.res)), length.out=100)
  lines(x, dnorm(x, mean=0, sd=sd(residuals(mod.res))))
  
  library(DHARMa)
  simulationOutput <- simulateResiduals(fittedModel = mod.res, plot = FALSE)
  
  plotQQunif(simulationOutput) # left plot in plot.DHARMa()
  plotResiduals(simulationOutput) # right plot in plot.DHARMa()
  
  #Old way without dharma, from Roger Mundry
  # qqnorm(residuals(mod.res), main="", pch=19)
  # qqline(residuals(mod.res))
  # mtext(text="qq-plot of residuals", side=3, line=0)
  # plot(fitted(mod.res), residuals(mod.res), pch=19, col=col)
  # abline(h=0, lty=2)
  # mtext(text="residuals against fitted values", side=3, line=0)
  par(old.par)
}

##---------------


##############
## Graphic related
##############


##-----------
# addImg
##-----------

#This function is used for optimising display of jpeg or png images directly into the graphs. it was taken from the internet 
#https://stackoverflow.com/questions/27800307/adding-a-picture-to-plot-in-r


addImg <- function(
  obj, # an image file imported as an array (e.g. png::readPNG, jpeg::readJPEG)
  x = NULL, # mid x coordinate for image
  y = NULL, # mid y coordinate for image
  width = NULL, # width of image (in x coordinate units)
  interpolate = TRUE # (passed to graphics::rasterImage) A logical vector (or scalar) indicating whether to apply linear interpolation to the image when drawing. 
){
  if(is.null(x) | is.null(y) | is.null(width)){stop("Must provide args 'x', 'y', and 'width'")}
  USR <- par()$usr # A vector of the form c(x1, x2, y1, y2) giving the extremes of the user coordinates of the plotting region
  PIN <- par()$pin # The current plot dimensions, (width, height), in inches
  DIM <- dim(obj) # number of x-y pixels for the image
  ARp <- DIM[1]/DIM[2] # pixel aspect ratio (y/x)
  WIDi <- width/(USR[2]-USR[1])*PIN[1] # convert width units to inches
  HEIi <- WIDi * ARp # height in inches
  HEIu <- HEIi/PIN[2]*(USR[4]-USR[3]) # height in units
  rasterImage(image = obj, 
              xleft = x-(width/2), xright = x+(width/2),
              ybottom = y-(HEIu/2), ytop = y+(HEIu/2), 
              interpolate = interpolate, xpd=TRUE)
}

#--

##-----------
# pastellizing colors
##-----------

#from: https://datascienceconfidential.github.io/r/graphics/2017/11/23/soothing-pastel-colours-in-r.html
pastellize <- function(x, p){
  
  library(ColorPalette)
  
  # x is a colour
  # p is a number in [0,1]
  # p = 1 will give no pastellization
  
  # convert hex or letter names to rgb
  if (is.character(x)) x <- col2rgb(x)/255
  
  # convert vector to rgb
  if (is.numeric(x)) x <- matrix(x, nr=3)
  
  col <- rgb2hsv(x, maxColorValue=1)
  col[2,1] <- col[2,1]*p
  col <- hsv2rgb(h=col[1,1], s=col[2,1], v=col[3,1])
  
  # return in convenient format for plots
  rgb(col[1], col[2], col[3])
}

#--

##---------------
# addGrid
##---------------

#This functions add a background grid

addGrid <- function(
xmin, #min of abscisse
xmax, #max of abscisse
xintsmall, #interval on abscisse for the thiner lines
xintbig, #interval on abscisse for the thicker lines
ymin, #min or ordinate
ymax, #max or ordinate
yintsmall, #interval on ordinate for the thiner lines
yintbig, #interval on ordinate for the thicker lines
colsmall="gray97", #colour for the thiner lines
colbig="gray93", #colour for the thicker lines
axisPlot=TRUE, #redraw or not the axes based on the sequence given for the ordinate and abscisse
lty=1,
lwdsmall=1,
lwdbig=1,
round=FALSE,#Rounding labels for axis?
digit=c(0,0)#Rounding number for x then y
){

segments(x0=xmin, x1=xmax, 
         y0=seq(from=ymin, to=ymax, by=yintsmall), 
         y1=seq(from=ymin, to=ymax, by=yintsmall),
         col=colsmall, lty=lty, lwd=lwdsmall)	


segments(x0=seq(from=xmin, to=xmax, by=xintsmall), x1=seq(from=xmin, to=xmax, by=xintsmall), 
         y0=ymin, 
         y1=ymax,
         col=colsmall, lty=lty, lwd=lwdsmall)
		

segments(x0=xmin, x1=xmax, 
         y0=seq(from=ymin,to=ymax, by=yintbig), 
         y1=seq(from=ymin,to=ymax, by=yintbig),
         col=colbig, lty=lty, lwd=lwdbig)


segments(x0=seq(from=xmin, to=xmax, by=xintbig), x1=seq(from=xmin, to=xmax, by=xintbig), 
         y0=ymin, 
         y1=ymax,
         col=colbig, lty=lty, lwd=lwdbig)
		 
if(axisPlot&round==FALSE){
	axis(side=1, at=seq(from=xmin, to=xmax, by=xintbig), labels=seq(from=xmin, to=xmax, by=xintbig), las=1, tcl=-0.25)
	axis(side=2, at=seq(from=ymin, to=ymax, by=yintbig), labels=seq(from=ymin, to=ymax, by=yintbig), las=1, tcl=-0.25)	
} else if(axisPlot&round){
	axis(side=1, at=round(seq(from=xmin, to=xmax, by=xintbig), digit=digit[1]), labels=round(seq(from=xmin, to=xmax, by=xintbig), digit=digit[1]), las=1, tcl=-0.25)
	axis(side=2, at=round(seq(from=ymin, to=ymax, by=yintbig), digit=digit[2]), labels=round(seq(from=ymin, to=ymax, by=yintbig), digit=digit[2]), las=1, tcl=-0.25)	
}
}

#--


##---------------
# addCorner
##---------------

#To add segment corners to the map

#This functions add a background grid

addCorner <- function(
xmin, #min of abscisse
xmax, #max of abscisse
ymin, #min or ordinate
ymax, #max or ordinate
lengthSegment=1/10,#proportion of grid that should be covered by the segment
predefinedLength=NA,#will erase proportion, if not square it's to be used
col="black",#colour of the corners
lty=1,#lty of the corners
lwd=1#lwd of the corner
){

if(is.na(predefinedLength)){
#upper left
segments(
x0=xmin,
x1=xmin + (xmax-xmin)*lengthSegment,
y0=ymax ,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmin,
x1=xmin,
y0=ymax - (ymax-ymin)*lengthSegment,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

#upperright
segments(
x0=xmax,
x1=xmax - (xmax-xmin)*lengthSegment,
y0=ymax ,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmax,
x1=xmax,
y0=ymax - (ymax-ymin)*lengthSegment,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

#lowerleft
segments(
x0=xmin,
x1=xmin + (xmax-xmin)*lengthSegment,
y0=ymin,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmin,
x1=xmin,
y0=ymin + (ymax-ymin)*lengthSegment,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)

#lowerright
segments(
x0=xmax,
x1=xmax - (xmax-xmin)*lengthSegment,
y0=ymin,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmax,
x1=xmax,
y0=ymin + (ymax-ymin)*lengthSegment,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)
}
else{
#upper left
segments(
x0=xmin,
x1=xmin + predefinedLength,
y0=ymax ,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmin,
x1=xmin,
y0=ymax - predefinedLength,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

#upperright
segments(
x0=xmax,
x1=xmax - predefinedLength,
y0=ymax ,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmax,
x1=xmax,
y0=ymax - predefinedLength,
y1=ymax,
col=col,
lty=lty,
lwd=lwd)

#lowerleft
segments(
x0=xmin,
x1=xmin + predefinedLength,
y0=ymin,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmin,
x1=xmin,
y0=ymin + predefinedLength,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)

#lowerright
segments(
x0=xmax,
x1=xmax - predefinedLength,
y0=ymin,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)

segments(
x0=xmin,
x1=xmin,
y0=ymin + predefinedLength,
y1=ymin,
col=col,
lty=lty,
lwd=lwd)
}
}
 
##---------------
# Create the sd bars
##--------------
errorBars <- function(location, meanPt, barValue, refUnit, minValue, maxValue, upperBarValue=NA, lowerBarValue=NA, col="black", lty=1, horiz=FALSE, symmetrical=TRUE){
#Location indicates the loc on the x or y axis (if vert or horiz plot)
#bar value is a unique value.
#If not symmetrical, barValue should be set to a random value (won't be used) and you should use upperBarValue and lowerBarValue instead.
#ref unit should indicates the unit to have 10% to draw the little border lines of the bars.
#min and max value should indicate whether there are borders to the variable range (or plot) to adjust the bars
	if(symmetrical){
		if(meanPt - barValue < minValue){
			minPt=minValue
		} else{
			minPt=meanPt - barValue
		}
		if(meanPt + barValue > maxValue){
			maxPt=maxValue
		} else{
			maxPt=meanPt + barValue
		}
		if(horiz){
		  segments(y0=location, y1=location, x0=minPt, x1=maxPt, col=col, lty=lty)
		  segments(y0=location - 0.1*refUnit, y1=location + 0.1*refUnit, x0=maxPt, x1=maxPt, col=col)
		  segments(y0=location - 0.1*refUnit, y1=location + 0.1*refUnit, x0=minPt, x1=minPt, col=col)
		} else{
		  segments(x0=location, x1=location, y0=minPt, y1=maxPt, col=col, lty=lty)
		  segments(x0=location - 0.1*refUnit, x1=location + 0.1*refUnit, y0=maxPt, y1=maxPt, col=col)
		  segments(x0=location - 0.1*refUnit, x1=location + 0.1*refUnit, y0=minPt, y1=minPt, col=col)
		}
	} else{
		if(length(lowerBarValue < minValue)==0){
			minPt=lowerBarValue 
		} else{
			lowerBarValue[lowerBarValue < minValue] <- minValue
			minPt=lowerBarValue
		}
		if(length(upperBarValue > maxValue)==0){
			maxPt=upperBarValue
		} else{
			upperBarValue[upperBarValue > maxValue]=maxValue
			maxPt=upperBarValue
		}
		if(horiz){
		  segments(y0=location, y1=location, x0=minPt, x1=maxPt, col=col, lty=lty)
		  segments(y0=location - 0.1*refUnit, y1=location + 0.1*refUnit, x0=maxPt, x1=maxPt, col=col)
		  segments(y0=location - 0.1*refUnit, y1=location + 0.1*refUnit, x0=minPt, x1=minPt, col=col)
		} else{
		  segments(x0=location, x1=location, y0=minPt, y1=maxPt, col=col, lty=lty)
		  segments(x0=location - 0.1*refUnit, x1=location + 0.1*refUnit, y0=maxPt, y1=maxPt, col=col)
		  segments(x0=location - 0.1*refUnit, x1=location + 0.1*refUnit, y0=minPt, y1=minPt, col=col)
		}
	}
}

#--	
 
 
  
##---------------
# Labelling a panel in a regular manner based on par definition
##---------------

#modified from: https://seananderson.ca/2013/10/21/panel-letters/
  
#' @param xfrac The fraction over from the left side.
#' @param yfrac The fraction down from the top.
#' @param label The text to label with.
#' @param pos Position to pass to text()
#' @param ... Anything extra to pass to text(), e.g. cex, col.

#circle: plot circle or not
#radiuscircle, size of the radius, as in draw.circle, in x units
#circle.border/bg, colour of border/bg of circle if plotted
#font.col #font colour
#font.size #font size
addLabel <- function(xfrac, yfrac, label, circle = FALSE, radiuscircle=NA, circle.bg=NA, circle.border="black", font.col="black", font.size=1, ...) {
  if(circle){
  	  u <- par("usr")
	  x <- u[1] + xfrac * (u[2] - u[1])
	  y <- u[4] - yfrac * (u[4] - u[3])
	  
	  library(plotrix)
	  draw.circle(x=x, y=y, radius=radiuscircle, border=circle.border, col=circle.bg)
	  text(x, y, label, col=font.col, cex=font.size, ...)

  }else{
	  u <- par("usr")
	  x <- u[1] + xfrac * (u[2] - u[1])
	  y <- u[4] - yfrac * (u[4] - u[3])
	  text(x, y, label, col=font.col, cex=font.size, adj=c(0.5, 0.5),...)
  }
}
  
  
#--	
    
  
  
##---------------
# Create filled pattern for polygons/Necessitates the splancs package
##---------------

# Please note:
# the tmap package (tm_fill, option pattern) and cartography package (patternLayer function) for more advanced patterns (hexagonal or wave for instance !!)


patternLayer <- function(x, y, pch=19, lwd=1, cex=1, col="black", xspacing, yspacing){

	#Initial polygon coordinates
	polygon_sp <- cbind(as.numeric(as.character(x)),as.numeric(as.character(y)))
	#Extract border points
	minX <- min(polygon_sp[,1])
	maxX <- max(polygon_sp[,1])
	minY <- min(polygon_sp[,2])
	maxY <- max(polygon_sp[,2])
	#Create a grid of the points
	library(splancs)
	pointsToPlot_init <- gridpts(polygon_sp,xs=xspacing, ys=yspacing)
	pointsToPlot_spshifted <- pointsToPlot_init
	pointsToPlot_spshifted[,1] <- pointsToPlot_spshifted[,1] + xspacing/2
	pointsToPlot_spshifted[,2] <- pointsToPlot_spshifted[,2] + yspacing/2 
	#Remove points out of polygon
	library(sp)
	isInPolygon <- sp::point.in.polygon(pointsToPlot_spshifted[,1], pointsToPlot_spshifted[,2] , polygon_sp[,1], polygon_sp[,2], mode.checked=FALSE)
	pointsToPlot_spshifted <- pointsToPlot_spshifted[isInPolygon==1,]
	pointsToPlot <- rbind(pointsToPlot_init, pointsToPlot_spshifted)
	
	#Plot
	points(x=pointsToPlot[,1], y=pointsToPlot[,2], pch=pch, cex=cex, lwd=lwd, col=col)
}
			
#--			

##---------------
# Create a spatial dataframe of points for a filled pattern (usable with ggplot with geom_points)
##---------------

#argProjection should meet the requirement for the CRS() function to work properly
patternLayerCoordinates <- function(x, y, xspacing, yspacing, argProjection){
	#Initial polygon coordinates
	polygon_sp <- cbind(as.numeric(as.character(x)),as.numeric(as.character(y)))
	#Extract border points
	minX <- min(polygon_sp[,1])
	maxX <- max(polygon_sp[,1])
	minY <- min(polygon_sp[,2])
	maxY <- max(polygon_sp[,2])
	#Create a grid of the points
	library(splancs)
	pointsToPlot_init <- gridpts(polygon_sp,xs=xspacing, ys=yspacing)
	pointsToPlot_spshifted <- pointsToPlot_init
	pointsToPlot_spshifted[,1] <- pointsToPlot_spshifted[,1] + xspacing/2
	pointsToPlot_spshifted[,2] <- pointsToPlot_spshifted[,2] + yspacing/2 
	#Remove points out of polygon
	library(sp)
	isInPolygon <- sp::point.in.polygon(pointsToPlot_spshifted[,1], pointsToPlot_spshifted[,2] , polygon_sp[,1], polygon_sp[,2], mode.checked=FALSE)
	pointsToPlot_spshifted <- 	pointsToPlot_spshifted[isInPolygon==1,]
	pointsToPlot <- rbind(pointsToPlot_init, pointsToPlot_spshifted)
	
	#Get the spatialPoints dataframe
	pointsToPlot_sp<- SpatialPointsDataFrame(pointsToPlot, data=as.data.frame(seq(from=1, to=nrow(pointsToPlot), by=1)))
	proj4string(pointsToPlot_sp) <-  CRS(argProjection)
	return(pointsToPlot_sp)
}
			
#--		

##################
## R markdown related
##################

##---------------
# To count and display word and character number
##---------------

#from: https://stackoverflow.com/questions/46317934/exclude-sections-from-word-count-in-r-markdown

#The file is the RMd file. All parts between
#<!---TC:ignore--->
#
#will not be counted. References are not counted either.
#Note that you rmd file need to be up to date. Therefore you have to knit twice.

RmdWords <- function(file) {

	library(stringr)
	library(tidyverse)

  # Creates a string of text
  file_string <- file %>%
    readLines() %>%
    paste0(collapse = " ") %>%
    # Remove YAML header
    str_replace_all("^<--- .*?--- ", "") %>%    
    str_replace_all("^--- .*?--- ", "") %>%
    # Remove code
    str_replace_all("```.*?```", "") %>%
    str_replace_all("`.*?`", "") %>%
    # Remove LaTeX
    str_replace_all("[^\\\\]\\$\\$.*?[^\\\\]\\$\\$", "") %>%
    str_replace_all("[^\\\\]\\$.*?[^\\\\]\\$", "") %>%
    # Deletes text between tags
    str_replace_all("TC:ignore.*?TC:endignore", "") %>%
    str_replace_all("[[:punct:]]", " ") %>%
    str_replace_all("  ", "") %>%
    str_replace_all("<", "") %>%
    str_replace_all(">", "")

  # Save several different results
  word_count <- str_count(file_string, "\\S+")
  char_count <- str_replace_all(string = file_string, " ", "") %>% str_count()

   return(list(num_words = word_count, num_char = char_count, word_list = file_string))
}

#--	

##---------------
# Merging the .bib file with articles with a newly created .bib with used packages
##---------------

#To cite R package directly, from: https://stackoverflow.com/questions/60026015/use-citation-in-r-markdown-to-automatically-generate-a-bibliography-of-r-packa

#bibliographyArticle: complete (absolute) path towards the .bib with all cited articles
#bibliographyArticle: complete (absolute) path towards the .bib with all cited articles + added R packages (i.e. OUTPUT)
#then, the different packages to use, separated by commas: e.g. readr, lubridate, Rcpp and so on. You can then cite them within the Rmarkdown file using @packageName. 
#But be careful, in doing so, some packages might be referenced with multiple citation (i.e. no "Rcpp" found because it creates @Rcpp1, @Rcpp2 etc...


citeR <- function(bibliographyArticle, bibliographyOutput, ...)
{
  packages <- unlist(lapply(as.list(match.call()), deparse))[c(-1, -2, -3)]
  #Rbibs <- ""
  library(bibtex)
  write.bib(packages, file = "packages.bib", append = FALSE, verbose = TRUE)
  
  big_bib <- c(readLines(bibliographyArticle), "\n", readLines("packages.bib"))
  writeLines(big_bib, bibliographyOutput)
 
}