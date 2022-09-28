rm(list=ls())

library(ggplot2)
library(tidyr)
library(faux)

sampleSize=10000
# b=0.001

diagnostics.plot.dharma <-function(mod.res, col=grey(level=0.25, alpha=0.5), breaks.histo=20, quantreg=TRUE){
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
  plotResiduals(simulationOutput, quantreg=quantreg)
  #Old way without dharma, from Roger Mundry
  # qqnorm(residuals(mod.res), main="", pch=19)
  # qqline(residuals(mod.res))
  # mtext(text="qq-plot of residuals", side=3, line=0)
  # plot(fitted(mod.res), residuals(mod.res), pch=19, col=col)
  # abline(h=0, lty=2)
  # mtext(text="residuals against fitted values", side=3, line=0)
  par(old.par)
}



# The function to test allometric relationships

testWeightBM <- function(sampleSize, #Number of species
                         b,
                         c,
                         d
                         #noiseBMtoBS, #Variance in a Gaussian distribution added to brainsize
                         #noiseZBS, #Variance in a Gaussian distribution added to brainsize
                         #noiseZrBS #Variance in a Gaussian distribution added to brainsize
){
  
  library(dplyr)
  library(lme4)
  
  dataToTest <- data.frame(
    id = 1:sampleSize#brainSize = rnorm(sampleSize, 10, 1)
  )

  dataToTest <- dataToTest %>%
    group_by(id) %>%
    mutate(
      
      #  simulate BS, BM, trait linked to BS
 
      dat = rnorm_multi(n = 1, 
                        mu = c(10, 100, 100, 100),
                        sd = c(1, 5, 5, 5),
                        r = c(b, max(c(1-b,0)), 0, 0, 0, 0), 
                        varnames = c("brainSize", "bodyMass", "traitLinkedToBrainSize", "traitNotLinkedToBrainSize"),
                        empirical = FALSE),
      brainSize=dat$brainSize,
      bodyMass=dat$bodyMass,
      traitLinkedToBrainSize=dat$traitLinkedToBrainSize,
      traitNotLinkedToBrainSize=dat$traitNotLinkedToBrainSize,
      
      relativeBrainSize = brainSize/bodyMass,
      

        
    )
  
  dataToTest <- dataToTest %>%
    filter(brainSize > 0 & bodyMass > 0 & relativeBrainSize < 1)

  # sample from conditional probability:
  # https://online.stat.psu.edu/stat414/lesson/21/21.1
  
  # relativeBrainSize: mean of 0.1, sd of sd_relativeBrainSize
  # traitLinkedToRelativeBrainSize: mean of 1, sd of 0.1
  # covariance: 1-b
  
  mean_sd_relativeBrainSize = mean(dataToTest$relativeBrainSize)
  sd_relativeBrainSize=sd(dataToTest$relativeBrainSize)
  
  dataToTest$traitLinkedToRelativeBrainSize = rnorm(nrow(dataToTest), 1 + (1-b)*0.1/sd_relativeBrainSize*(dataToTest$relativeBrainSize-mean_sd_relativeBrainSize), 0.1*sqrt(1-(1-b)^2) )
   
  dataToTest$traitLinkedToRelativeBrainSize100 = 100*dataToTest$traitLinkedToRelativeBrainSize

  
  # hist(dataToTest$relativeBrainSize)
  # mean(dataToTest$relativeBrainSize)
  # sd(dataToTest$relativeBrainSize)
  
  
  # plot(dataToTest$traitLinkedToRelativeBrainSize, dataToTest$relativeBrainSize) #, asp = 1)
  # abline(a = 0, b = 0.1, col = "red")
  # summary(lm(dataToTest$relativeBrainSize~dataToTest$traitLinkedToRelativeBrainSize))


  # plot(dataToTest$bodyMass, dataToTest$brainSize)
  
  # 
  # plot(dataToTest$brainSize, dataToTest$bodyMass)
  # abline(a = 0, b = 10, col = "red")
  # summary(lm(dataToTest$relativeBrainSize~dataToTest$bodyMass))
  # summary(lm(dataToTest$brainSize~dataToTest$bodyMass))
  # 
  # summary(lm(dataToTest$traitLinkedToRelativeBrainSize100~dataToTest$bodyMass))
  # 
  # summary(lm(dataToTest$brainSize~dataToTest$bodyMass+dataToTest$traitLinkedToRelativeBrainSize100))

  # plot(dataToTest$traitLinkedToBrainSize, dataToTest$bodyMass)
  # abline(a = 0, b = 10, col = "red")
  # 
  # plot(dataToTest$traitNotLinkedToBrainSize, dataToTest$brainSize)
  # abline(a = 0, b = 10, col = "red")
  
  # plot(dataToTest$traitLinkedToBrainSize, dataToTest$bodyMass)
  # 
  # hist(dataToTest$traitLinkedToBrainSize)
  # mean(dataToTest$traitLinkedToBrainSize)
  # var(dataToTest$traitLinkedToBrainSize)
  # 
  # hist(dataToTest$bodyMass)
  # mean(dataToTest$bodyMass)
  # var(dataToTest$bodyMass)
  # 
  # hist(dataToTest$traitLinkedToRelativeBrainSize)
  # mean(dataToTest$traitLinkedToRelativeBrainSize)
  # var(dataToTest$traitLinkedToRelativeBrainSize)
  # 
  # hist(dataToTest$relativeBrainSize)
  
  #Model 1:
  
  modelBMcovariatea <- lm(brainSize ~ traitNotLinkedToBrainSize + traitLinkedToBrainSize + bodyMass, data = dataToTest)
  modelBMcovariateb <- lm(brainSize ~ traitNotLinkedToBrainSize + traitLinkedToRelativeBrainSize100 + bodyMass, data = dataToTest)
  
  #Model 2:
  
  modelRelativeBSa <- lm(relativeBrainSize ~ traitNotLinkedToBrainSize + traitLinkedToBrainSize, data = dataToTest)
  modelRelativeBSb <- NA #lm(brainSize ~ (traitNotLinkedToBrainSize + traitLinkedToBrainSize)*bodyMass, data = dataToTest)
  modelRelativeBSc <- lm(relativeBrainSize ~ traitNotLinkedToBrainSize + traitLinkedToRelativeBrainSize , data = dataToTest)
  modelRelativeBSd <- NA #lm(brainSize ~ (traitNotLinkedToBrainSize + traitLinkedToRelativeBrainSize)*bodyMass, data = dataToTest)
  
  return(
    list(
      dataToTest,
      modelBMcovariatea,
      modelBMcovariateb,
      modelRelativeBSa,
      modelRelativeBSb,
      modelRelativeBSc,
      modelRelativeBSd
    )
  )
}

toTest <- crossing(seq(0.1, 0.9, by = 0.1), seq(0.1, 0.9, by = 0.1)) %>% as.data.frame()
toTest <- toTest*10

# new
toTest <- cbind(seq(0.05, 0.95, by = 0.05),seq(0.05, 0.95, by = 0.05))

#~~~~~~~~~~~~~~
# BS as an output, BM as a covariate
#~~~~~~~~~~~~~~

library(tidyr)
library(car)

# Predictor linked to BS

resultsBS <- apply(toTest, 1, function(a){
  output <- testWeightBM(10000,a[1],a[2],5)
  return(c(summary(output[[2]])$coefficients[2:4, 1], max(vif(output[[2]]))))
}
)

resultsBS <- t(resultsBS)#do.call("rbind", resultsBS) %>% as.data.frame()
resultsBS <- as.data.frame(resultsBS)
colnames(resultsBS)[length(colnames(resultsBS))] <- "VIF"
resultsBS$ratioNoise <- toTest[,1]#/toTest[,2] # New

resultsBS_BSoutputBSPred <- resultsBS %>%
  pivot_longer(!c(ratioNoise, VIF), names_to = "trait", values_to = "effectSize") %>%
  mutate(
    what = "BSBMz"
  )

# Predictor linked to BS/BM

resultsBS <- apply(toTest, 1, function(a){
  output <- testWeightBM(10000,a[1],5,a[2])
  return(c(summary(output[[3]])$coefficients[2:4, 1], max(vif(output[[3]]))))
}
)

resultsBS <- t(resultsBS)#do.call("rbind", resultsBS) %>% as.data.frame()
resultsBS <- as.data.frame(resultsBS)
colnames(resultsBS)[length(colnames(resultsBS))] <- "VIF"
resultsBS$ratioNoise <- toTest[,1]#/toTest[,2] # New

resultsBS_BSoutputratioPred <- resultsBS %>%
  pivot_longer(!c(ratioNoise, VIF), names_to = "trait", values_to = "effectSize") %>%
  mutate(
    what = "BSBMzr"
  )

#~~~~~~~~~~~~~~
# BS/BM as an output
#~~~~~~~~~~~~~~

# Predictor linked to BS

# resultsBS <- apply(toTest, 1, function(a){
#   output <- testWeightBM(1000,a[1],a[2],5)
#   return(c(summary(output[[4]])$coefficients[2:3, 1], max(vif(output[[4]]))))
# }
# )
# 
# resultsBS <- t(resultsBS)#do.call("rbind", resultsBS) %>% as.data.frame()
# resultsBS <- as.data.frame(resultsBS)
# colnames(resultsBS)[length(colnames(resultsBS))] <- "VIF"
# resultsBS$ratioNoise <- toTest[,1]/toTest[,2]
# 
# resultsBS_ratioOutputBSpred <- resultsBS %>%
#   pivot_longer(!c(ratioNoise, VIF), names_to = "trait", values_to = "effectSize") %>%
#   mutate(
#     what = "BSrz"
#   )

# Predictor linked to BS/BM

resultsBS <- apply(toTest, 1, function(a){
  output <- testWeightBM(10000,a[1],5,a[2])
  return(c(summary(output[[6]])$coefficients[2:3, 1], max(vif(output[[6]]))))
}
)

resultsBS <- t(resultsBS)#do.call("rbind", resultsBS) %>% as.data.frame()
resultsBS <- as.data.frame(resultsBS)
colnames(resultsBS)[length(colnames(resultsBS))] <- "VIF"
resultsBS$ratioNoise <- toTest[,1]#/toTest[,2] # New

resultsBS_ratioOutputratioPred <- resultsBS %>%
  pivot_longer(!c(ratioNoise, VIF), names_to = "trait", values_to = "effectSize") %>%
  mutate(
    what = "BSrzr"
  )

# Plot results

completeResults <- rbind(
  resultsBS_BSoutputBSPred,
  resultsBS_BSoutputratioPred,
  #resultsBS_ratioOutputBSpred,
  resultsBS_ratioOutputratioPred
)

library(ggplot2)


completeResults$trait[which(completeResults$trait=="traitLinkedToRelativeBrainSize100")] <- "traitLinkedToRelativeBrainSize"

plot <- ggplot(completeResults %>%
                 filter(what != "BSrz" &
                          ratioNoise > 0.05 &
                          ratioNoise < 10) %>%
                 mutate(what = dplyr::recode(what, "BSBMz" = "BS ~ null + Z + BM", "BSBMzr" = "BS ~ null + Zr * 100+ BM", "BSrzr" = "BS/BM ~ null + Zr")),
               aes(x = ratioNoise, y = effectSize, colour = trait)) +
  #geom_line() +
  geom_point(shape = 21, bg = "white") +
  scale_colour_discrete(name = "Predictor", labels = c('Body mass (BM)', 'Variable related to BS', 'Variable related to BS/BM','Null variable' )) +
  facet_wrap(~ what) +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 15),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 10),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12,
                                  face = "bold"),
        panel.grid.minor=element_line(colour="grey93"),
        panel.grid.major=element_line(colour="grey93"),
        strip.text = element_text(face = "bold", size = rel(1.25)),
        strip.background = element_rect(fill = "white", colour = "black", size = 0)) +
  ylab("Linear regression estimate") + #scale_x_log10()+
  xlab("Covariance") #+ ylim(c(-0.1,1))
pdf("/Users/bperez/Nextcloud/Recherche/These/ENS/Papier/Papier_primates/PCI_revisions/relecture_3/simulations_covariances.pdf", width=10.5, height=6)
plot
dev.off()

