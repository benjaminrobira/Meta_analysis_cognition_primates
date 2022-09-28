# What does this script?
# This script plots the phylogenetic history of diet as reconstructed based on our analysis.

#reminder of vector used:
frugivoryThresholdVector <- seq(from=20, to=40, by=20)
folivoryThresholdVector <- seq(from=40, to=60, by=20)
geographicThresholdVector <- c(10,30)/100
randomSampling=10
numberSimulations=10
numberTrees=1


# Loading data  --------
# Selecting the case where frugivory is the most stringent and
a = 2
b = 1
c = 1
d = 1
setwd("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/")
load(file=paste("Processed_data/Simmap/Output_simmap_transition", a, "_", b, "_", c, "_", d, ".Rdata", sep=""))


#Parameter for temporal legend
alpha = 0.2
show.lines = TRUE
cex = 0.5

colors <- setNames(
  c(rgb(255, 242, 127, 255, maxColorValue = 255), 
                     rgb(255, 230, 25, 255, maxColorValue = 255), rgb(253, 
                                                                      154, 82, 255, maxColorValue = 255), rgb(127, 
                                                                                                              198, 78, 255, maxColorValue = 255), rgb(52, 178, 
                                                                                                                                                      201, 255, maxColorValue = 255), rgb(129, 43, 
                                                                                                                                                                                          146, 255, maxColorValue = 255), rgb(240, 64, 
                                                                                                                                                                                                                              40, 255, maxColorValue = 255), rgb(103, 165, 
                                                                                                                                                                                                                                                                 153, 255, maxColorValue = 255), rgb(203, 140, 
                                                                                                                                                                                                                                                                                                     55, 255, maxColorValue = 255), rgb(179, 225, 
                                                                                                                                                                                                                                                                                                                                        182, 255, maxColorValue = 255), rgb(0, 146, 112, 
                                                                                                                                                                                                                                                                                                                                                                            255, maxColorValue = 255), rgb(127, 160, 86, 
                                                                                                                                                                                                                                                                                                                                                                                                           255, maxColorValue = 255), rgb(247, 67, 112, 
                                                                                                                                                                                                                                                                                                                                                                                                                                          255, maxColorValue = 255)), c("Quaternary", "Neogene", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "Paleogene", "Cretaceous", "Jurassic", "Triassic", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "Permian", "Carboniferous", "Devonian", "Silurian", 
                                                                                                                                                                                                                                                                                                                                                                                                                                                                        "Ordovician", "Cambrian", "Precambrian"))
leg <- rbind(c(2.588, 0), c(23.03, 2.588), c(66, 23.03), 
             c(145, 66), c(201.3, 145), c(252.17, 201.3), c(298.9, 
                                                            252.17), c(358.9, 298.9), c(419.2, 358.9), c(443.8, 
                                                                                                         419.2), c(485.4, 443.8), c(541, 485.4), c(4600, 
                                                                                                                                                   541))
rownames(leg) <- c("Quaternary", "Neogene", "Paleogene", 
                   "Cretaceous", "Jurassic", "Triassic", "Permian", 
                   "Carboniferous", "Devonian", "Silurian", "Ordovician", 
                   "Cambrian", "Precambrian")

plot(summarySimmapToPlot,
     colors = setNames(c("darkgreen","salmon"),
                       c("Leaf", "Fruit")),
     asp = 1, 
     cex = 0.3,
     fsize = 0.4,
     offset = 0.5,
     ylim = c(5,153),
     hold = FALSE)
obj <- get("last_plot.phylo", envir = .PlotPhyloEnv)
max(obj$xx)
min(obj$xx)
t.max <- max(obj$xx)

ii <- which(leg[, 2] <= t.max)
leg <- leg[ii, ]
leg[max(ii), 1] <- t.max
colors2 <- sapply(colors, make.transparent, alpha = alpha)
ylabelmin = 0
y <- c(rep(0, 2), rep(par()$usr[4], 2))
ylabel <- par()$usr[4]

# Computing the plot

##Main plot
library(phytools)
summarySimmapToPlot <- describe.simmap(simmapdiet1,plot=FALSE)
par(ask=F)

pdf("Plots/plotReconstructionDiet.pdf", width = 7, height = 14)

plot(summarySimmapToPlot,
     colors = setNames(c("darkgreen","salmon"),
                       c("Leaf", "Fruit")),
     asp = 1, 
     cex = 0.3,
     fsize = 0.4,
     offset = 0.5,
     ylim = c(5,153),
     hold = FALSE)

for(i in 1:nrow(leg)){
  polygon(abs(c(leg[i, 1:2], leg[i, 2:1]) - max(leg[,1])), y, col = colors2[rownames(leg)[i]], 
          border = NA)
  # text(x = abs(leg[i, 1] - max(leg[,1])), y = ylabel, labels = rownames(leg)[i], 
  #   #srt = 90,
  #   adj = c(0, 0.5),
  #   cex = cex)
}
legend(
  x = -1,#-0.001*par()$usr[1],
  y = 0.1*par()$usr[4],#par()$usr[4]/2,
  fill = colors[1:nrow(leg)],
  legend = names(colors[1:nrow(leg)]),
  bty = "n",
  cex = 0.8,
  adj = 0
)

add.simmap.legend(colors=setNames(c("darkgreen","salmon"),
                                  c("Folivore", "Frugivore")), x=0.75,#-0.1*par()$usr[1],
                  y=0.035*par()$usr[4],prompt=FALSE, shape = "circle", adj = 0, cex = 0.7, font = 1)

library(ape)
add.scale.bar(0, 0.01*par()$usr[4], cex = 0.8)
text(x = 14.5, y = 0.01*par()$usr[4], labels = "Ma", cex = 0.8)
dev.off()

##Then working with inkscape to:
# - add legend on phylogeny
# - add info on fossils

# Figure caption:
# 
# (https://www.cambridge.org/core/journals/paleobiology/article/abs/reconstructing-dietary-ecology-of-extinct-strepsirrhines-primates-mammalia-with-new-approaches-for-characterizing-and-analyzing-tooth-shape/563B3D7FBB4EEA9930C354FE5C073D74).
# 
# ##Adding chronological timeline
# 
# https://journals.plos.org/plosone/article?id=10.1371%2Fjournal.pone.0049521&fbclid=IwAR3Jv0pRcA8ghwwWWU5huZQcU_bEz0EnkzMRZ9QMSCcy48NsaytCMyiSQ-A
# https://books.google.it/books?hl=en&lr=&id=_-MOaksvLQ8C&oi=fnd&pg=PA133&ots=ZXp6hsMMYF&sig=W2osCli8HGqNrNScfqhq9WajARo&redir_esc=y#v=onepage&q&f=false
# https://www.sciencedirect.com/science/article/pii/S004724840900133X?fbclid=IwAR2iIpx6VkU-_6DCVIsAaBsVbBpvcOC1BvYc1gBGzQD7gA6SQX20AhlX-uo
# https://www.sciencedirect.com/science/article/pii/S0047248417303822?fbclid=IwAR16UVcUDJaNqdG-nGsLjpN2IICDzkBLiQjG1q5fckp6Mxk92Kn-vIvsvDs
# https://onlinelibrary.wiley.com/doi/full/10.1002/ajpa.21638?fbclid=IwAR3tlpYZcpefTT67J7ZlcD0IqqRxA1DzS9i9VUl2IeA5dvaQRZXYZwcPB5U
