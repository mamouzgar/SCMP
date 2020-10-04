#' this script generates plots from the 01-generate-features.R script output
#' 

library(ggpubr)
library(RColorBrewer)
set.seed(0)
setwd("/Users/mamouzgar/phd-projects")


###############
## functions ##
###############
# sp_plot <- function(input.data, cluster.method = NULL, markers = NULL) { 
  # ggscatter(data.subset, x= "axis1", y = "axis2",
  #           color = "gate", size = 1 ) +
  # scale_color_brewer(palette = "Set1")
# } 

density_2d_plot <- function(input.data, cluster.method = NULL, cell_counts = NULL) { 
  d2d.plot <- ggplot(input.data, aes(x=axis1,y = axis2)) +
    geom_density_2d(aes(color = gate)) + 
    scale_color_brewer(palette = "Set1") +
    theme_pubr() + 
    theme(legend.position = "none") +
    ggtitle(paste(cluster.method, cell_counts,sep = "-")) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),

          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  return(d2d.plot)
  }

sp_plot  <- function(input.data, cluster.method = NULL, cell_counts = NULL) { 
  
  sp.plot <- ggscatter(input.data, x= "axis1", y = "axis2", 
                          color = "gate", size = 0.1, legend = "top",
                          title = paste(cluster.method, cell_counts,sep = "-"),
                          shape = 1) +
    scale_color_brewer(palette = "Set1") +
    theme_pubr() +
    theme(legend.position = "none") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),

          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  return(sp.plot)
  
}

###############
## variables ##
###############
result.files <- list.files("SCMP/data/analysis-ready", full.names = TRUE)
result.files <- result.files[grepl("pca.csv|tsne.csv|umap.csv|phate.csv",result.files)]
# result.files <- result.files[!grepl("176664",result.files)]
result.files 
result.files <- data.frame(result.files, stringsAsFactors = FALSE) %>%
  dplyr::mutate(cell.counts = gsub(".*processed_|-cells.*", "", result.files),
                cell.counts = as.numeric(cell.counts),
                method = gsub(".*_|.csv", "", result.files),
                method = factor(method, levels = c("pca","tsne","umap","phate"))) %>%
  arrange(cell.counts, method)
# pca.files <- result.files[grepl("pca.csv",result.files)]
# tsne.files <- result.files[grepl("tsne.csv",result.files)]
# umap.files <- result.files[grepl("umap.csv",result.files)]
# phate.files <- result.files[grepl("phate.csv",result.files)]

plot.output.path <- "SCMP/analysis/plots/"

##########
## MAIN ##
##########
plots.list <- lapply(result.files$result.files, function(filepath) { 
  print(filepath)
  df <- data.table::fread(filepath) %>%
    mutate(gate = factor(gate, levels = c("monocyte","lymphocyte","erythroid","blast", "neutrophil")))
  
  df <- df[, c(1,2, ncol(df))]
  colnames(df) <- c("axis1","axis2","gate")
  
  if (grepl("phate", filepath) | any(table(df$gate) <= 1)) { 
    sp_plot(df, cluster.method = gsub(".*_|.csv", "", filepath), cell_counts = gsub(".*processed_|-cells.*", "", filepath))
    
  } else { 
    density_2d_plot(df, cluster.method = gsub(".*_|.csv", "", filepath), cell_counts = gsub(".*processed_|-cells.*", "", filepath))
    }
  })

###############
## aggregate ##
###############
library(gridExtra)
# library(grid)
pdf(paste0(plot.output.path, "plotting-figures.pdf"), height = 18, width = 18)
do.call("grid.arrange", c(plots.list, ncol = 4))
dev.off()
  
  
# pdf(paste0(plot.output.path, "legend-labels.pdf"), height = 1.5, width = 1.5)
# plot(get_legend(plots.list[[8]] + theme(legend.position="right")))
# sp_plot(dat.clean, cluster.method = "LDA", cell_counts = 176664) + theme(legend.position="right"))
# dev.off()



#########
## LDA ##
#########

cat("reading global variables")
dat <- data.table::fread("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/processed_LDaxes.csv", sep = ",", stringsAsFactors = FALSE)

dat.clean <- dat %>%
  dplyr::filter(gate != "ungated") %>%
  dplyr::select( axis1=ld1, axis2=ld2, gate) ## remove irrelevant markers - confirm with David

density_2d_plot(dat.clean, cluster.method = "LDA", cell_counts = 176664)


pdf(paste0(plot.output.path, "LDA-plot_all-cells.pdf"), height = 2, width = 3)
density_2d_plot(dat.clean, cluster.method = "LDA", cell_counts = 176664)+theme(legend.position="right")
dev.off()






