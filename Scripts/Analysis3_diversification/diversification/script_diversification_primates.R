
##################  Step 1: Prepare primate tree ###################################################

# from Dos Reis 2018
# set the MRCA at 74 Myr 

library(ape)

setwd("/Users/bperez/ownCloud/Recherche/These/ENS/Autres/Benjamin/diversification/")

tree <- read.tree("host_tree_primate_complete.tre")

tree$edge.length <- tree$edge.length/max(node.depth.edgelength(tree))*74

plot(tree)
tree <- ladderize(tree)

tree$edge.length
tree$tip.label

is.binary(tree)
is.rooted(tree)
is.ultrametric(tree)

write.tree(tree, "tree_primate_complete.tre")


# run using Julia on the cluster


##################  Step 2: Results ClaDS ###################################################

rm(list=ls())

library(RPANDA)
library(coda)

plot_ClaDS_phylo <- function (phylo, rates, rates2 = NULL, same.scale = T, main = NULL, lwd = 2, log = T, show.tip.label = F, ...) {
  Colors = colorRampPalette(c("steelblue2", "paleturquoise3", 
                              "palegreen2", "yellow2", "salmon1", "darkorange", "red", 
                              "red4"))(100)
  if (is.null(rates2)) {
    if (log) 
      rates = log(rates)
    if (isTRUE(all.equal(rep(as.numeric(rates[1]), length(rates)), 
                         as.numeric(rates)))) {
      col = rep(1, length(rates))
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
           main = main, edge.width = lwd, ...)
      if (log) {
        fields::image.plot(z = c(exp(rates[1]), 2 * exp(rates[1])), 
                           col = Colors, horizontal = T, legend.only = T)
      }
      else {
        fields::image.plot(z = c(rates[1], 2 * rates[1]), col = Colors, 
                           horizontal = T, legend.only = T)
      }
    }
    else {
      col = round((rates - min(rates))/diff(range(rates)) * 
                    99) + 1
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
           main = main, edge.width = lwd, ...)
      if (log) {
        min = min(rates)
        max = max(rates)
        m10 = floor(min/log(10))
        M10 = ceiling(max/log(10))
        if ((M10 - m10) < 4) {
          ticks = c(1, 2, 5)
        }
        else {
          ticks = 1
        }
        ticks = as.vector(sapply(m10:M10, function(k) {
          return(ticks * 10^k)
        }))
        lt = length(ticks[ticks > exp(min) & ticks < 
                            exp(max)])
        if (lt < 4) {
          ##### CHANGE THIS (03/09/2019)
          ticks = c(exp(min), max(0, -1 * m10 + 
                                    (lt < 2)), ticks, c(exp(max), max(0, 
                                                                      -1 * M10 + 1 + (lt < 2))))
        }
        fields::image.plot(z = c(min, max), col = Colors, horizontal = T, 
                           legend.only = T, axis.args = list(at = log(ticks), 
                                                             labels = ticks))
      }
      else {
        fields::image.plot(z = as.matrix(rates), col = Colors, 
                           horizontal = T, legend.only = T)
      }
    }
  }
  else {
    if (log) {
      rates = log(rates)
      rates2 = log(rates2)
    }
    if (same.scale) {
      min = min(min(rates), min(rates2))
      max = max(max(rates), max(rates2))
      par(mfrow = c(1, 2))
      col = round(((rates - min)/(max - min)) * 99) + 1
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
           edge.width = lwd, ...)
      col = round(((rates2 - min)/(max - min)) * 99) + 
        1
      plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
           edge.width = lwd, ...)
      par(mfrow = c(1, 1))
      if (log) {
        m10 = floor(min/log(10))
        M10 = ceiling(max/log(10))
        if ((M10 - m10) < 4) {
          ticks = c(1, 2, 5)
        }
        else {
          ticks = 1
        }
        ticks = as.vector(sapply(m10:M10, function(k) {
          return(ticks * 10^k)
        }))
        lt = length(ticks[ticks > exp(min) & ticks < 
                            exp(max)])
        if (lt < 4) {
          ticks = c(round(exp(min), max(0, -1 * m10 + 
                                          (lt < 2))), ticks, round(exp(max), max(0, 
                                                                                 -1 * M10 + 1 + (lt < 2))))
        }
        fields::image.plot(z = c(min, max), col = Colors, horizontal = T, 
                           legend.only = T, axis.args = list(at = log(ticks), 
                                                             labels = ticks))
      }
      else {
        fields::image.plot(z = c(min, max), col = Colors, horizontal = T, 
                           legend.only = T)
      }
    }
    else {
      par(mfrow = c(1, 2))
      if (isTRUE(all.equal(rep(rates[1], length(rates)), 
                           rates))) {
        col = rep(1, length(rates))
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
             edge.width = lwd, ...)
        if (log) {
          fields::image.plot(z = c(exp(rates[1]), 2 * exp(rates[1])), 
                             col = Colors, horizontal = T, legend.only = T)
        }
        else {
          fields::image.plot(z = c(rates[1], 2 * rates[1]), col = Colors, 
                             horizontal = T, legend.only = T)
        }
      }
      else {
        col = round(((rates - min(rates))/(max(rates) - 
                                             min(rates))) * 99) + 1
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
             edge.width = lwd, ...)
        if (log) {
          min = min(rates)
          max = max(rates)
          m10 = floor(min/log(10))
          M10 = ceiling(max/log(10))
          if ((M10 - m10) < 4) {
            ticks = c(1, 2, 5)
          }
          else {
            ticks = 1
          }
          ticks = as.vector(sapply(m10:M10, function(k) {
            return(ticks * 10^k)
          }))
          lt = length(ticks[ticks > exp(min) & ticks < 
                              exp(max)])
          if (lt < 4) {
            ticks = c(round(exp(min), max(0, -1 * m10 + 
                                            (lt < 2))), ticks, round(exp(max), max(0, 
                                                                                   -1 * M10 + 1 + (lt < 2))))
          }
          fields::image.plot(z = c(min, max), col = Colors, horizontal = T, 
                             legend.only = T, axis.args = list(at = log(ticks), 
                                                               labels = ticks))
        }
        else {
          fields::image.plot(z = as.matrix(rates), col = Colors, 
                             horizontal = T, legend.only = T)
        }
      }
      if (isTRUE(all.equal(rep(rates2[1], length(rates2)), 
                           rates2))) {
        col = rep(1, length(rates2))
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
             edge.width = lwd, ...)
        if (log) {
          fields::image.plot(z = c(exp(rates2[1]), 2 * exp(rates2[1])), 
                             col = Colors, horizontal = T, legend.only = T)
        }
        else {
          fields::image.plot(z = c(rates2[1], 2 * rates2[1]), 
                             col = Colors, horizontal = T, legend.only = T)
        }
      }
      else {
        col = round(((rates2 - min(rates2))/(max(rates2) - 
                                               min(rates2))) * 99) + 1
        plot(phylo, edge.color = Colors[col], show.tip.label = show.tip.label, 
             edge.width = lwd, ...)
        if (log) {
          min = min(rates2)
          max = max(rates2)
          m10 = floor(min/log(10))
          M10 = ceiling(max/log(10))
          if ((M10 - m10) < 4) {
            ticks = c(1, 2, 5)
          }
          else {
            ticks = 1
          }
          ticks = as.vector(sapply(m10:M10, function(k) {
            return(ticks * 10^k)
          }))
          lt = length(ticks[ticks > exp(min) & ticks < 
                              exp(max)])
          if (lt < 4) {
            ticks = c(round(exp(min), max(0, -1 * m10 + 
                                            (lt < 2))), ticks, round(exp(max), max(0, 
                                                                                   -1 * M10 + 1 + (lt < 2))))
          }
          fields::image.plot(z = c(min, max), col = Colors, horizontal = T, 
                             legend.only = T, axis.args = list(at = log(ticks), 
                                                               labels = ticks))
        }
        else {
          fields::image.plot(z = as.matrix(rates2), col = Colors, 
                             horizontal = T, legend.only = T)
        }
      }
    }
    par(mfrow = c(1, 1))
  }
}


setwd("/users/biodiv/bperez/data/others/Benji/diversification/")

f="60"

for (f in c( "95", "90", "80", "70", "60", "67")){
  
  print(f)
  name=paste0("ClaDS2_tree_primate_complete_f",f)
  
  load(paste0("ClaDS2_tree_primate_complete_f",f,".Rdata"))
  


  pdf(paste0("MAPS_rates_tree_",name,"_ClaDS2.pdf"), width=8,height=10)
  plot_ClaDS_phylo(tree, MAPS[5:npar], log = T, lwd = 3)
  title(main = paste0("sigma = ", round(MAPS[1],2), " ; alpha = ", round(MAPS[2],2), " ; epsilon = ", round(MAPS[3],2), " ; lambda_0 = ", round(MAPS[4],2)))
  dev.off()
  
  
  speciation_rates <- data.frame(cbind(tree$tip.label,(MAPS[(1:Ntip(tree)) + nrow(tree$edge) + 4])))
  colnames(speciation_rates) <- c("Species","Speciation_rate")
  speciation_rates$Speciation_rate <- as.numeric(as.character(speciation_rates$Speciation_rate))
  
  write.table(speciation_rates, paste0("MAPS_speciation_rates_tips_",name,".csv"), sep=";",row.names=F,quote=F)
  
  # density
  MAPS <- c()
  for (i in 1:length(mr[[1]][[1]])){
    D=density(sapply(1:length(mr), function(k) sapply(floor(length(mr[[1]])*0.3):length(mr[[1]]), function(j) mr[[k]][[j]][i])))
    MAPS <- c(MAPS, D[[1]][which.max(D[[2]])])
  }
  plot(MAPS)
  table_MAPS_rates <- cbind(seq(-74,0,length.out = length(mr[[1]][[1]])), MAPS)
  colnames(table_MAPS_rates) <- c("time","speciation_rates")
  write.table(table_MAPS_rates,paste0("MAPS_speciation_rates_trought_time_",name,".csv"),sep=";",quote=F,row.names = F)
  
  # mean speciation rates
  mean_rates <- vector("numeric",length = length(mr[[1]][[1]]))
  count <- 0
  for (i in c(1:3)){
    for (j in (length(mr[[i]])*0.3):length(mr[[i]])){ # 30% burn-in
      count <- count +1
      mean_rates <- mean_rates + mr[[i]][[j]]
    }
  }
  mean_rates <- mean_rates/count
  print(plot(seq(-505,0,length.out = length(mr[[1]][[1]])), mean_rates))
  
  table_mean_rates <- cbind(seq(-74,0,length.out = length(mr[[1]][[1]])),mean_rates)
  colnames(table_mean_rates) <- c("time","speciation_rates")
  write.table(table_mean_rates,paste0("mean_speciation_rates_trought_time_",name,".csv"),sep=";",quote=F,row.names = F)
  
  # index <- 1:12
  # 
  # for (rep in index){
  #   load(paste0("ClaDS2_tree_GMYC_units_LSU_replicate_",rep,"_f",f,".Rdata"))
  #   
  #   # density
  #   MAPS <- c()
  #   for (i in 1:length(mr[[1]][[1]])){
  #     D=density(sapply(1:length(mr), function(k) sapply(floor(length(mr[[1]])*0.3):length(mr[[1]]), function(j) mr[[k]][[j]][i])))
  #     MAPS <- c(MAPS, D[[1]][which.max(D[[2]])])
  #   }
  #   plot(MAPS)
  #   table_MAPS_rates <- cbind(seq(-505,0,length.out = length(mr[[1]][[1]])), MAPS)
  #   colnames(table_MAPS_rates) <- c("time","speciation_rates")
  #   write.table(table_MAPS_rates,paste0("MAPS_speciation_rates_trought_time_ClaDS2_tree_GMYC_units_LSU_replicate_",rep,"_f",f,".csv"),sep=";",quote=F,row.names = F)
  #   
  # }
  
  
  # Plot 
  
  pdf(paste0("MAPS_rates_throught_time_", name,".pdf"), width=5, height = 5)
  
  #max_y=0.019
  #max_y=0.05
  
  table_mean_rates <- read.table(paste0("MAPS_speciation_rates_trought_time_",name,".csv"), sep=";",header=T)
  
  print(plot(table_mean_rates$time,table_mean_rates$speciation_rates, type='l', lwd=4, pch=19, col="#212f3c", xlab="Time (in Myr)", ylab="MAPS speciation rates (/Myr)"))
  #ylim=c(0.007,max_y), axes = F))
  
  # for (i in index){
  #   table_mean_rates_rep <- read.table(paste0("MAPS_speciation_rates_trought_time_ClaDS2_tree_GMYC_units_LSU_replicate_",i,"_f",f,".csv"), sep=";",header=T)
  #   par(new=T)
  #   print(max(table_mean_rates_rep$speciation_rates))
  #   print(plot(table_mean_rates_rep$time,table_mean_rates_rep$speciation_rates, type='l', lwd=2.5, pch=19, col="#e67e22", xlab=" ", ylab=" ", ylim=c(0.007,max_y), axes = F))
  # }
  # 
  # par(new=T)
  # print(plot(table_mean_rates$time,table_mean_rates$speciation_rates, type='l', lwd=4, pch=19, col="#212f3c", xlab=" ", ylab=" ", ylim=c(0.007,max_y)))
  # 
  dev.off()
  
}

# scp bperez@jord.biologie.ens.fr:/users/biodiv/bperez/data/others/Benji/diversification/MAPS*  /Users/bperez/ownCloud/Recherche/These/ENS/Autres/Benjamin/diversification/





##################  Step 3: Look for correlations  ###################################################

rm(list=ls())


# setwd("/users/biodiv/bperez/data/others/Benji/Evolutionary_history/")
# load("geography_traits_biogeobears_2.RData")

setwd("/Users/bperez/ownCloud/Recherche/These/ENS/Autres/Benjamin/diversification/")


tree <- read.tree("tree_primate_complete.tre")

summaryData <- read.table("../results/Evolutionary_history/Dataplot.txt", header=TRUE, sep="\t")

nrow(summaryData) # 428 species
Ntip(tree) # 367 species 
table(summaryData$SpeciesForPhylogeny %in% tree$tip.label)

summaryData$SpeciesForPhylogeny[!summaryData$SpeciesForPhylogeny %in% tree$tip.label]
tree$tip.label[!tree$tip.label %in% summaryData$SpeciesForPhylogeny]

# replace the species name in summaryData by the long name 
for (i in 1:nrow(summaryData)){
  if (!summaryData$SpeciesForPhylogeny[i] %in% tree$tip.label){
    if (length(grep(summaryData$SpeciesForPhylogeny[i], tree$tip.label))>0){
      summaryData$SpeciesForPhylogeny[i] <- tree$tip.label[grep(summaryData$SpeciesForPhylogeny[i], tree$tip.label)]
    }
  }
}

table(summaryData$SpeciesForPhylogeny %in% tree$tip.label)

summaryData$SpeciesForPhylogeny[summaryData$SpeciesForPhylogeny=="Pan_troglodytes_schweinfurthii"] <- "Pan_troglodytes"
summaryData$SpeciesForPhylogeny[summaryData$SpeciesForPhylogeny=="Aotus_azarai"] <- "Aotus_azarae"
summaryData$SpeciesForPhylogeny[summaryData$SpeciesForPhylogeny=="Eulemur_fulvus_albifrons"] <- "Eulemur_fulvus"
summaryData$SpeciesForPhylogeny[summaryData$SpeciesForPhylogeny=="Eulemur_macaco_flavifrons"] <- "Eulemur_macaco"
summaryData$SpeciesForPhylogeny[summaryData$SpeciesForPhylogeny=="Gorilla_gorilla_gorilla"] <- "Gorilla_gorilla"
summaryData$SpeciesForPhylogeny[summaryData$SpeciesForPhylogeny=="Varecia_variegata_variegata"] <- "Varecia_variegata"

sort(summaryData$SpeciesForPhylogeny[!summaryData$SpeciesForPhylogeny %in% tree$tip.label])
sort(tree$tip.label[!tree$tip.label %in% summaryData$SpeciesForPhylogeny])
table(summaryData$SpeciesForPhylogeny %in% tree$tip.label)


#Create variable of rinterest
summaryData$Brain.log <- log(summaryData$Brain)
summaryData$ratioBrain <- summaryData$Brain*1.036*(10**-3)/summaryData$Bodymass #Following decasien for multiplication by 1.036
summaryData$EQ <- summaryData$Brain*1.036*(10**-3)/(0.085*summaryData$Bodymass**0.775) #Following decasien, according to #Jerison, H. J. Evolution of the Brain and Intelligence (Academic, 1973).
summaryData$ratioNeocortex <- summaryData$Neocortex/ summaryData$Brain
summaryData$ratioHippocampus <- summaryData$Hippocampus/ summaryData$Brain
summaryData$ratioCerebellum <- summaryData$Cerebellum/ summaryData$Brain
summaryData$ratioStriatum <- summaryData$Striatum/ summaryData$Brain
summaryData$ratioMOB <- summaryData$MOB/ summaryData$Brain
summaryData$ratioMOB.log <- log(summaryData$ratioMOB)


f=67

for (f in c(60, 67, 70, 80, 90, 95)){
  
  print(f)
  
  summaryData_div <- summaryData[which(summaryData$SpeciesForPhylogeny %in% tree$tip.label),]
  
  table_MAPS_rates <- read.table(paste0("MAPS_speciation_rates_tips_ClaDS2_tree_primate_complete_f",f,".csv"), sep=";",header=T)
  
  summaryData_div$speciation_rate <- NA
  for (i in 1:nrow(summaryData_div)){
    summaryData_div$speciation_rate[i] <- table_MAPS_rates$Speciation_rate[which(table_MAPS_rates$Species==summaryData_div$SpeciesForPhylogeny[i])]
  }
  
  hist(summaryData_div$speciation_rate)
  
  pdf(paste0("Correlation_speciation_rates_traits_f",f,".pdf"), width=4, height=3.5)
  plot(summaryData_div$Brain.log, summaryData_div$speciation_rate, xlab="Brain.log", ylab="Speciation rate (/Myr)", col="#e59866", pch=19)
  plot(summaryData_div$EQ, summaryData_div$speciation_rate, xlab="EQ", ylab="Speciation rate (/Myr)", col="#e59866", pch=19)
  plot(summaryData_div$ratioNeocortex, summaryData_div$speciation_rate, xlab="ratioNeocortex", ylab="Speciation rate (/Myr)", col="#e59866", pch=19)
  plot(summaryData_div$ratioHippocampus, summaryData_div$speciation_rate, xlab="ratioHippocampus", ylab="Speciation rate (/Myr)", col="#e59866", pch=19)
  plot(summaryData_div$ratioCerebellum, summaryData_div$speciation_rate, xlab="ratioCerebellum", ylab="Speciation rate (/Myr)", col="#e59866", pch=19)
  plot(summaryData_div$ratioStriatum, summaryData_div$speciation_rate, xlab="ratioStriatum", ylab="Speciation rate (/Myr)", col="#e59866", pch=19)
  plot(summaryData_div$ratioMOB.log, summaryData_div$speciation_rate, xlab="ratioMOB.log", ylab="Speciation rate (/Myr)", col="#e59866", pch=19)
  dev.off()
  
  
  # PGLS
  
  # Brain
  summaryData_omit<- summaryData_div[-which(is.na(summaryData_div$Brain.log)),]
  rownames(summaryData_omit) <- summaryData_omit$SpeciesForPhylogeny
  BM_tree <- corBrownian(1, drop.tip(tree,tip=tree$tip.label[!tree$tip.label %in% summaryData_omit$SpeciesForPhylogeny]), form = ~SpeciesForPhylogeny)
  bm.gls<-nlme::gls(speciation_rate ~ Brain.log, correlation=BM_tree, data=summaryData_omit)
  print(summary(bm.gls))
  
  # EQ
  summaryData_omit<- summaryData_div[-which(is.na(summaryData_div$EQ)),]
  rownames(summaryData_omit) <- summaryData_omit$SpeciesForPhylogeny
  BM_tree <- corBrownian(1, drop.tip(tree,tip=tree$tip.label[!tree$tip.label %in% summaryData_omit$SpeciesForPhylogeny]), form = ~SpeciesForPhylogeny)
  bm.gls<-nlme::gls(speciation_rate ~ EQ, correlation=BM_tree, data=summaryData_omit)
  print(summary(bm.gls))
  
  # ratioNeocortex
  summaryData_omit<- summaryData_div[-which(is.na(summaryData_div$ratioNeocortex)),]
  rownames(summaryData_omit) <- summaryData_omit$SpeciesForPhylogeny
  BM_tree <- corBrownian(1, drop.tip(tree,tip=tree$tip.label[!tree$tip.label %in% summaryData_omit$SpeciesForPhylogeny]), form = ~SpeciesForPhylogeny)
  bm.gls<-nlme::gls(speciation_rate ~ ratioNeocortex, correlation=BM_tree, data=summaryData_omit)
  print(summary(bm.gls))
  
  # ratioHippocampus
  summaryData_omit<- summaryData_div[-which(is.na(summaryData_div$ratioHippocampus)),]
  rownames(summaryData_omit) <- summaryData_omit$SpeciesForPhylogeny
  BM_tree <- corBrownian(1, drop.tip(tree,tip=tree$tip.label[!tree$tip.label %in% summaryData_omit$SpeciesForPhylogeny]), form = ~SpeciesForPhylogeny)
  bm.gls<-nlme::gls(speciation_rate ~ ratioHippocampus, correlation=BM_tree, data=summaryData_omit)
  print(summary(bm.gls))
  
  # ratioCerebellum
  summaryData_omit<- summaryData_div[-which(is.na(summaryData_div$ratioCerebellum)),]
  rownames(summaryData_omit) <- summaryData_omit$SpeciesForPhylogeny
  BM_tree <- corBrownian(1, drop.tip(tree,tip=tree$tip.label[!tree$tip.label %in% summaryData_omit$SpeciesForPhylogeny]), form = ~SpeciesForPhylogeny)
  bm.gls<-nlme::gls(speciation_rate ~ ratioCerebellum, correlation=BM_tree, data=summaryData_omit)
  print(summary(bm.gls))
  
  # ratioStriatum
  summaryData_omit<- summaryData_div[-which(is.na(summaryData_div$ratioStriatum)),]
  rownames(summaryData_omit) <- summaryData_omit$SpeciesForPhylogeny
  BM_tree <- corBrownian(1, drop.tip(tree,tip=tree$tip.label[!tree$tip.label %in% summaryData_omit$SpeciesForPhylogeny]), form = ~SpeciesForPhylogeny)
  bm.gls<-nlme::gls(speciation_rate ~ ratioStriatum, correlation=BM_tree, data=summaryData_omit)
  print(summary(bm.gls)) 
  
  # ratioMOB.log
  summaryData_omit<- summaryData_div[-which(is.na(summaryData_div$ratioMOB.log)),]
  rownames(summaryData_omit) <- summaryData_omit$SpeciesForPhylogeny
  BM_tree <- corBrownian(1, drop.tip(tree,tip=tree$tip.label[!tree$tip.label %in% summaryData_omit$SpeciesForPhylogeny]), form = ~SpeciesForPhylogeny)
  bm.gls<-nlme::gls(speciation_rate ~ ratioMOB.log, correlation=BM_tree, data=summaryData_omit)
  print(summary(bm.gls))

}





