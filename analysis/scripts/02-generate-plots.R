#' this script generates plots from the 01-generate-features.R script output
#' 

library(ggpubr)
library(RColorBrewer)
set.seed(0)
setwd("/Users/mamouzgar/phd-projects")

###############
## variables ##
###############
result.files <- list.files("SCMP/data/analysis-ready", full.names = TRUE)

pca.files <- result.files[grepl("pca.csv",result.files)]
tsne.files <- result.files[grepl("tsne.csv",result.files)]
umap.files <- result.files[grepl("umap.csv",result.files)]
phate.files <- result.files[grepl("phate.csv",result.files)]

plot.output.path <- "SCMP/analysis/plots/"

#########
## pca ##
#########
pca.plots <- lapply(pca.files, function(filepath) { 
  df <- data.table::fread(filepath)
  df <- df %>%
    mutate(gate = factor(gate, levels = c("monocyte","lymphocyte","erythroid","blast", "neutrophil")))
  
  sp <-ggscatter(df, x= "PC1", y = "PC2", color = "gate", size = 1, legend = "none", main = "pca", subtitle = basename(filepath)) +
    scale_color_brewer(palette = "Set1")
  
  # sp + stat_density_2d(aes(fill = ..level..)) +
  #   gradient_fill("Set1")
  
  pdf(paste0(plot.output.path, basename(filepath), ".pdf"), height = 3, width =4)
  print(sp)
  dev.off()
  return(sp)
})


##########
## tsne ##
##########
tsne.plots <- lapply(tsne.files, function(filepath) { 
  df <- data.table::fread(filepath)
  df <- df %>%
    mutate(gate = factor(gate, levels = c("monocyte","lymphocyte","erythroid","blast", "neutrophil")))
  
  sp <-ggscatter(df, x= "X1", y = "X2", color = "gate", size = 1, legend = "none", main = "tsne", subtitle = basename(filepath)) +
    scale_color_brewer(palette = "Set1")
  # sp + stat_density_2d(aes(fill = ..level..)) +
  #   gradient_fill("Set1")
  
  pdf(paste0(plot.output.path, basename(filepath), ".pdf"), height = 3, width =3)
  print(sp)
  dev.off()
  return(sp)
})


##########
## umap ##
##########
umap.plots <- lapply(umap.files, function(filepath) { 
  df <- data.table::fread(filepath)
  df <- df %>%
    mutate(gate = factor(gate, levels = c("monocyte","lymphocyte","erythroid","blast", "neutrophil")))
  
  sp <-ggscatter(df, x= "X1", y = "X2", color = "gate", size = 1, legend = "none", main = "umap", subtitle = basename(filepath)) +
    scale_color_brewer(palette = "Set1")
  # sp + stat_density_2d(aes(fill = ..level..)) +
  #   gradient_fill("Set1")
  
  pdf(paste0(plot.output.path, basename(filepath), ".pdf"), height = 3, width =3)
  print(sp)
  dev.off()
  return(sp)
})



##########
## umap ##
##########
phate.plots <- lapply(phate.files, function(filepath) { 
  df <- data.table::fread(filepath)
  df <- df %>%
    mutate(gate = factor(gate, levels = c("monocyte","lymphocyte","erythroid","blast", "neutrophil")))
  
  sp <-ggscatter(df, x= "PHATE1", y = "PHATE2", color = "gate", size = 1, legend = "none", main = "phate", subtitle = basename(filepath)) +
    scale_color_brewer(palette = "Set1")
  # sp + stat_density_2d(aes(fill = ..level..)) +
  #   gradient_fill("Set1")
  
  pdf(paste0(plot.output.path, basename(filepath), ".pdf"), height = 3, width =3)
  print(sp)
  dev.off()
  return(sp)
})


##############
## LDA-plot ##
##############
dat <- data.table::fread("SCMP/data/analysis-ready/processed_LDaxes.csv")
dat <- dat %>%
  mutate(gate = factor(gate, levels = c("monocyte","lymphocyte","erythroid","blast", "neutrophil")))
dat[ , "cell.id"] <- 1:nrow(dat)

analyzed.cells <- data.table::fread("SCMP/data/analysis-ready/2020-09-18_analyzed-cells_50000-cells.csv")
df.ld <- dat[ dat$cell.id %in% analyzed.cells$x, ]

sp <-ggscatter(df.ld, x= "ld1", y = "ld2", color = "gate", size = 1, legend = "right", main = "LDA", subtitle = "subsetted-cells") +
  scale_color_brewer(palette = "Set1")
pdf(paste0(plot.output.path, "LD-results_subsetted-cells.pdf"), height = 3, width =3)
print(sp)
dev.off()

legend.plot <- get_legend(sp)
# Convert to a ggplot and print
legend.plot <- as_ggplot(legend.plot)

pdf(paste0(plot.output.path, "legend.pdf"), height = 1.5, width =1.5)
print(legend.plot)
dev.off()

###############
## aggregate ##
###############
library(gridExtra)
# library(grid)

do.call("grid.arrange", c(pca.plots, tsne.plots, umap.plots, phate.plots, ncol = 3, nrow = 4),  legend.plot)

# do.call("grid.arrange", c(pca.plots, tsne.plots, umap.plots, ncol = 3, nrow = 2))








