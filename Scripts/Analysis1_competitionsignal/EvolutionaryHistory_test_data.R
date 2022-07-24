rm(list=ls())

table <- read.table("/Users/bperez/Documents/GitHub/Meta_analysis_cognition_primates/Processed_data/Dataplot.txt", header=TRUE)

sort(colnames(table))

for (i in 1:nrow(table)){
  
  for (area in c("Body_mass_g", "Brain_volume_mm3", "Cerebellum_volume_mm3", "Hippocampus_volume_mm3",
                 "Neocortex_volume_mm3", "Striatum_volume_mm3"))
  data <- unlist(table[i,grep(area,colnames(table))])
  data <- data[!is.na(data)]
  if (length(data)!=length(unique(data))){
    print(i)
    print(area)
  }
}




