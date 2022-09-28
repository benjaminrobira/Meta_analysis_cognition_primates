#Load environments
load("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/REnvironments/Data_spatial_primate.RData")
load("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/REnvironments/geography_traits_biogeobears.RData")

load("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/REnvironments/PGLSdiversification_withautocorr.RData")
load("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/REnvironments/PGLSdirectionSelection.RData")
load("C:/Users/robira/Documents/PhD/Meta_analysis/Meta_analysis_cognition_primates/REnvironments/PGLSdiversificationAndSympatry.RData")

head(summaryDataForPlot)


summaryDataWeighting <- data.frame(
  species = summaryData$SpeciesForPhylogeny,
  Neocortex_BS = summaryData$Neocortex/summaryData$Brain,
  Cerebellum_BS = summaryData$Cerebellum/summaryData$Brain,
  Hippocampus_BS = summaryData$Hippocampus/summaryData$Brain,
  Striatum_BS = summaryData$Striatum/summaryData$Brain,
  MOB_BS = summaryData$MOB/summaryData$Brain,
  Neocortex_BM = summaryData$Neocortex/summaryData$Bodymass,
  Cerebellum_BM = summaryData$Cerebellum/summaryData$Bodymass,
  Hippocampus_BM = summaryData$Hippocampus/summaryData$Bodymass,
  Striatum_BM = summaryData$Striatum/summaryData$Bodymass,
  MOB_BM = summaryData$MOB/summaryData$Bodymass
)
# colnamesDf <- colnames(summaryDataWeighting)
# summaryDataWeighting  <- t(apply(summaryDataWeighting , 1, as.numeric)) %>%  as.data.frame()
# colnames(summaryDataWeighting) <- colnamesDf

library(tidyr)
library(dplyr)

summaryDataWeighting <- summaryDataWeighting %>% 
  pivot_longer(!species, names_to = "area", values_to = "ratio") %>% 
  separate(
    col = area, 
    into = c("area", "weight"),
    sep = "_"
  ) %>% 
  # filter(!is.na(ratio)) %>% 
  pivot_wider(names_from = weight, values_from = ratio)



library(ape)
library(phylolm)
library(phytools)

resultsLinearRegressioNPerArea <- vector(mode = "list", length = length(unique(summaryDataWeighting$area)))
resultsRdcLinearRegressionPerArea <- data.frame(
  x = rep(NA, times = length(resultsLinearRegressioNPerArea)),
  y = rep(NA, times = length(resultsLinearRegressioNPerArea)),
  area = unique(summaryDataWeighting$area),
  intercept = rep(NA, times = length(resultsLinearRegressioNPerArea)),
  est = rep(NA, times = length(resultsLinearRegressioNPerArea)),
  r = rep(NA, times = length(resultsLinearRegressioNPerArea))
)

for(i in 1:length(resultsLinearRegressioNPerArea)){
  
  data = summaryDataWeighting[summaryDataWeighting$area == unique(summaryDataWeighting$area)[i],]
  rownames(data) <- data$species
  model <- phylolm(formula = BS ~ BM, 
                     data=data, phy = drop.tip(phylo,
                                               phylo$tip.label[
                                                 which(phylo$tip.label
                                                       %nin% data$species)]),
                     model = "lambda", measurement_error = FALSE, boot = repetitionBootstrap)
  resultsLinearRegressioNPerArea[[i]] <- summary(model)
  resultsRdcLinearRegressionPerArea$intercept[i] <- summary(model)$coefficients[1,1]
  resultsRdcLinearRegressionPerArea$est[i] <- summary(model)$coefficients[2,1]
  resultsRdcLinearRegressionPerArea$r[i] <- summary(model)$r.squared
}


library(ggplot2)

plotCorrelatioNBSandBM <- ggplot(summaryDataWeighting, aes(x = BM, y = BS, group = area)) +
  facet_grid(~ area, scales = "free_x") +
  #Not phylogenetic regressions
  # geom_smooth(method = "lm", colour = "black", fill = "gray") +
  # stat_cor(aes(label = ..rr.label..), color = "black", geom = "label", label.y.npc="top", label.x.npc = "left") +
  # Phylogenetic results
  geom_label(data = resultsRdcLinearRegressionPerArea, aes(x = 0, y = 0.85, group = area, label = paste("R² =", round(r, digit = 2))), hjust = 0) +
  geom_abline(data = resultsRdcLinearRegressionPerArea,aes(intercept = intercept, slope = est)) +
  geom_point() +
  theme_bw() +
  theme(axis.title = element_text(face = "bold", size = 16),
        axis.text.x = element_text(size = 12),
        axis.text.y = element_text(size = 13),
        legend.text=element_text(size=10),
        legend.title=element_text(size=12, face = "bold"),
        panel.grid.minor=element_line(colour="grey93"),
        panel.grid.major=element_line(colour="grey93"),
        strip.background = element_blank(),
        strip.text.x = element_text(face = "bold", size = 16)) +
  xlab("Body mass") +
  ylab("Body size")
plotCorrelatioNBSandBM

saveRDS(plotCorrelatioNBSandBM, "REnvironments/plotCorrelationBSandBM.rds")
