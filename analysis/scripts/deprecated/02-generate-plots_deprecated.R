#' this script generates plots from the 01-generate-features.R script output
#' 

library(ggpubr)
library(RColorBrewer)
set.seed(0)
setwd("/Users/mamouzgar/phd-projects")


###############
## functions ##
###############

sp_plot <-function(input.data, cluster.method = NULL, markers = NULL) { 
  ggscatter(input.data, x= "axis1", y = "axis2", 
            color = "gate", size = 1, legend = "none", main = cluster.method, subtitle = markers ) +
    scale_color_brewer(palette = "Set1")
} 

###############
## variables ##
###############
result.files <- list.files("SCMP/data/analysis-ready", full.names = TRUE)
result.files <- result.files[grepl("pca.csv|tsne.csv|umap.csv|phate.csv|LDaxes.csv",result.files)]
# pca.files <- result.files[grepl("pca.csv",result.files)]
# tsne.files <- result.files[grepl("tsne.csv",result.files)]
# umap.files <- result.files[grepl("umap.csv",result.files)]
# phate.files <- result.files[grepl("phate.csv",result.files)]

plot.output.path <- "SCMP/analysis/plots/"

##########
## MAIN ##
##########
plots.list <- lapply(result.files, function(filepath) { 
  cat(filepath)
  df <- data.table::fread(filepath) %>%
    mutate(gate = factor(gate, levels = c("monocyte","lymphocyte","erythroid","blast", "neutrophil")))
  
  if (grepl("LDaxes.csv", filepath)) { 
    df[ , "cell.id"] <- 1:nrow(dat)
    analyzed.cells <- data.table::fread("SCMP/data/analysis-ready/2020-09-18_analyzed-cells_50000-cells.csv")
    df <- df[ df$cell.id %in% analyzed.cells$x, c("ld1","ld2","gate") ]
  }
  
  df <- df[, c(1,2, ncol(df))]
  colnames(df) <- c("axis1","axis2","gate")
  
  sp <- sp_plot(df, cluster.method = gsub(".*_|.csv", "", filepath), markers = gsub(".*_cells_","",filepath) )
  
  pdf(paste0(plot.output.path, basename(filepath), ".pdf"), height = 3, width =4)
  print(sp)
  dev.off()
  return(sp)
})

###############
## aggregate ##
###############
library(gridExtra)
# library(grid)

do.call("grid.arrange", c(plots.list, ncol = 3, nrow = 4))

# do.call("grid.arrange", c(pca.plots, tsne.plots, umap.plots, ncol = 3, nrow = 2))



##############################################################################################################################
##################
## facet-method ##
##################

result.files <- result.files[!grepl("LDaxes.csv", result.files)]
data.compiled <- lapply(result.files, function(filepath) { 
  print(filepath)
  df <- data.table::fread(filepath) %>%
    mutate(gate = factor(gate, levels = c("monocyte","lymphocyte","erythroid","blast", "neutrophil")))
  
  if (grepl("LDaxes.csv", filepath)) { 
    df[ , "cell.id"] <- 1:nrow(dat)
    analyzed.cells <- data.table::fread("SCMP/data/analysis-ready/2020-09-18_analyzed-cells_50000-cells.csv")
    df <- df[ df$cell.id %in% analyzed.cells$x, c("ld1","ld2","gate") ]
  }
  
  df <- df[, c(1,2, ncol(df))]
  colnames(df) <- c("axis1","axis2","gate")
  df[ , "method"] <- gsub(".*_|.csv", "", filepath)
  df[ , "marker"] <- gsub("_.*","",gsub(".*-cells_","",filepath))
  return(df)
})

data.compiled <- do.call("bind_rows", data.compiled)

pdf(paste0(plot.output.path, "plotting-figures.pdf"), height = 16, width = 16)
ggscatter(data.compiled, x= "axis1", y = "axis2", 
          color = "gate", size = 1, legend = "top",
          facet.by = c("method","marker")) +
  scale_color_brewer(palette = "Set1")
dev.off()











##

#' this script generates plots from the 01-generate-features.R script output
#' 

library(ggpubr)
library(RColorBrewer)
set.seed(0)
setwd("/Users/mamouzgar/phd-projects")


###############
## functions ##
###############
sp_plot <- function(input.data, cluster.method = NULL, markers = NULL) { 
  ggscatter(input.data, x= "axis1", y = "axis2", 
            color = "gate", size = 1, legend = "none", main = cluster.method, subtitle = markers ) +
    scale_color_brewer(palette = "Set1")
} 

###############
## variables ##
###############
result.files <- list.files("SCMP/data/analysis-ready", full.names = TRUE)
result.files <- result.files[grepl("pca.csv|tsne.csv|umap.csv|phate.csv|LDaxes.csv",result.files)]
result.files <- result.files[grepl("5000-",result.files)]
# pca.files <- result.files[grepl("pca.csv",result.files)]
# tsne.files <- result.files[grepl("tsne.csv",result.files)]
# umap.files <- result.files[grepl("umap.csv",result.files)]
# phate.files <- result.files[grepl("phate.csv",result.files)]

plot.output.path <- "SCMP/analysis/plots/"

##########
## MAIN ##
##########
plots.list <- lapply(result.files, function(filepath) { 
  print(filepath)
  df <- data.table::fread(filepath) %>%
    mutate(gate = factor(gate, levels = c("monocyte","lymphocyte","erythroid","blast", "neutrophil")))
  
  # if (grepl("LDaxes.csv", filepath)) { 
  #   df[ , "cell.id"] <- 1:nrow(dat)
  #   analyzed.cells <- data.table::fread("SCMP/data/analysis-ready/2020-09-18_analyzed-cells_50000-cells.csv")
  #   df <- df[ df$cell.id %in% analyzed.cells$x, c("ld1","ld2","gate") ]
  # }
  
  df <- df[, c(1,2, ncol(df))]
  colnames(df) <- c("axis1","axis2","gate")
  
  sp <- sp_plot(df, cluster.method = gsub(".*_|.csv", "", filepath), markers = gsub(".*_cells_","",filepath) )
  
  pdf(paste0(plot.output.path, basename(filepath), ".pdf"), height = 3, width =4)
  print(sp)
  dev.off()
  return(sp)
})

###############
## aggregate ##
###############
library(gridExtra)
# library(grid)

do.call("grid.arrange", c(plots.list, ncol = 4, nrow = 1))

# do.call("grid.arrange", c(pca.plots, tsne.plots, umap.plots, ncol = 3, nrow = 2))



##############################################################################################################################
##################
## facet-method ##
##################

result.files <- result.files[!grepl("LDaxes.csv", result.files)]
data.compiled <- lapply(result.files, function(filepath) { 
  print(filepath)
  df <- data.table::fread(filepath) %>%
    mutate(gate = factor(gate, levels = c("monocyte","lymphocyte","erythroid","blast", "neutrophil")))
  
  if (grepl("LDaxes.csv", filepath)) { 
    df[ , "cell.id"] <- 1:nrow(dat)
    analyzed.cells <- data.table::fread("SCMP/data/analysis-ready/2020-09-18_analyzed-cells_50000-cells.csv")
    df <- df[ df$cell.id %in% analyzed.cells$x, c("ld1","ld2","gate") ]
  }
  
  df <- df[, c(1,2, ncol(df))]
  colnames(df) <- c("axis1","axis2","gate")
  df[ , "method"] <- gsub(".*_|.csv", "", filepath)
  df[ , "marker"] <- gsub("_.*","",gsub(".*-cells_","",filepath))
  return(df)
})

data.compiled <- do.call("bind_rows", data.compiled)

pdf(paste0(plot.output.path, "plotting-figures.pdf"), height = 8, width = 8)
ggscatter(data.compiled, x= "axis1", y = "axis2", 
          color = "gate", size = 1, legend = "top",
          scales = "free",
          facet.by = c("method")) +
  scale_color_brewer(palette = "Set1") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
dev.off()









#' this script generates plots from the 01-generate-features.R script output
#' 

library(ggpubr)
library(RColorBrewer)
set.seed(0)
setwd("/Users/mamouzgar/phd-projects")


###############
## functions ##

sp_plot <-function(input.data, cluster.method = NULL) { 
  ggscatter(input.data, x= input.data[,1], y = input.data[,2], 
            color = "gate", size = 1, legend = "none", main = cluster.method) +
  scale_color_brewer(palette = "Set1")
} 

###############
## variables ##
###############
result.files <- list.files("SCMP/data/analysis-ready", full.names = TRUE)
result.files <- result.files[grepl("pca.csv|tsne.csv|umap.csv|phate.csv|LDaxes.csv",result.files)]
# pca.files <- result.files[grepl("pca.csv",result.files)]
# tsne.files <- result.files[grepl("tsne.csv",result.files)]
# umap.files <- result.files[grepl("umap.csv",result.files)]
# phate.files <- result.files[grepl("phate.csv",result.files)]

plot.output.path <- "SCMP/analysis/plots/"

##########
## MAIN ##
##########
plots.list <- lapply(result.files, function(filepath) { 
  df <- data.table::fread(filepath) %>%
    mutate(gate = factor(gate, levels = c("monocyte","lymphocyte","erythroid","blast", "neutrophil"))) %>%
    
  if (grepl("LDaxes.csv", filepath)) { 
  } else if (grepl("tsne.csv")) { 
      
    
    }
  
  })

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

do.call("grid.arrange", c(pca.plots, tsne.plots, umap.plots, phate.plots, legend.plot, ncol = 3, nrow = 4))

# do.call("grid.arrange", c(pca.plots, tsne.plots, umap.plots, ncol = 3, nrow = 2))








