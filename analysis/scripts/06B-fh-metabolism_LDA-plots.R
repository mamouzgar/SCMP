

rm(list = ls())
library(tidyverse)
library(ggpubr)
library(viridis)
library(gridExtra)
setwd("/Users/mamouzgar/phd-projects/")




files.list <- list.files("SCMP/data/summary-tables/", full.names= TRUE) %>% .[grepl("continuous-axes-predictions",.)]
# df <- data.table::fread("SCMP/data/summary-tables/CD8naivefh-metabolism-continuous-axes-predictions.csv") %>%
#   mutate(day= factor(day, levels = c(0,1,2,3,4,5)))
 

for(filepath in files.list) { 
  print(filepath)
  df <- data.table::fread(filepath) %>%
    mutate(day= factor(day, levels = c(0,1,2,3,4,5)))
  
  head(df)
  
  ## plots
  ## full dataset
  df.plot <- filter(df, LDA.analysis == "full-balanced-dataset")
  output.name <- paste0(unique(df.plot$cell.type.subtype), "_", unique(df.plot$subject), "_", unique(df.plot$LDA.analysis))
  
  plot.full <- ggscatter(df.plot, 
                         x = "ld1", y= "ld2", color = "day",
                         title = output.name,
                         subtitle = paste(min(table(df.plot$day)), "cells per timepoint"),
                         legend = "right",
                         size = 0.05) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    scale_color_viridis(discrete =TRUE, option = "viridis")
  
  plot.facet <- ggscatter(df.plot, 
                          x = "ld1", y= "ld2", color = "day",
                          facet.by = "day",
                          nrow = 1, 
                          legend = "none",
                          size = 0.2) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    scale_color_viridis(discrete =TRUE, option = "viridis")
  
  
  pdf(paste0("SCMP/analysis/plots/fh-metabolism/LDA-timepoint-analysis/", output.name, ".pdf"),
      height = 7,
      width = 6)
  grid.arrange(plot.full,plot.facet, nrow = 2)
  dev.off()
  
  
  
  ## predicting dataset
  df.plot <- filter(df, LDA.analysis == "predicting-timepoint-balanced-subset")
  output.name <- paste0(unique(df.plot$cell.type.subtype), "_", unique(df.plot$subject), "_", unique(df.plot$LDA.analysis))
  
  plot.full <- ggscatter(df.plot, 
                         x = "ld1", y= "ld2", color = "day",
                         title =output.name,
                         legend = "right",
                         # subtitle = paste(min(table(df.plot$day)), "cells per timepoint"),
                         subtitle = "predicting-timepoint: train on day 0,2,4 + test on day 1,3,5",
                         size = 0.05) +
    guides(color = guide_legend(override.aes = list(size=5))) +
    scale_color_viridis(discrete =TRUE, option = "viridis")
  
  plot.facet <- ggscatter(df.plot, 
                          x = "ld1", y= "ld2", color = "day",
                          facet.by = "day",
                          nrow = 1, 
                          legend = "none",
                          size = 0.2) +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank()) +
    scale_color_viridis(discrete =TRUE, option = "viridis")
  
  
  pdf(paste0("SCMP/analysis/plots/fh-metabolism/LDA-timepoint-analysis/", output.name, ".pdf"),
      height = 7,
      width = 6)
  grid.arrange(plot.full,plot.facet, nrow = 2)
  dev.off()
}

