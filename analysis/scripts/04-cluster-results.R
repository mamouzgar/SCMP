


library(tidyverse)
library(ggpubr)
setwd("/Users/mamouzgar/phd-projects")


plot.output.path <- "SCMP/analysis/plots/k-medoids-clustering-analysis/"

k.medoid.results <- data.table::fread("SCMP/data/summary-tables/k-medoid-accuracy-results.csv")




ggbarplot(filter(k.medoid.results, balance == "unbalanced", accuracy != 0 ) , 
          x = "labels", y = "accuracy", facet.by = c("cell.count","method")) + 
  scale_y_continuous( limits=c(0, 1)) +
  coord_flip()


pdf(paste0(plot.output.path, "line-plot_unbalanced-data_accuracy-vs-cell.count.pdf"), height = 8, width = 8)
ggline(filter(k.medoid.results, balance == "unbalanced", accuracy != 0 ),
       x = "cell.count", y = "accuracy", group= "lables", color = "labels", facet.by = c("method"),
       title = "unbalanced-data (distribution preserved)",
       subtitle = "k-medoids euclidean distance clustering accuracy")
dev.off()


pdf(paste0(plot.output.path, "line-plot_balanced-data_accuracy-vs-cell.count.pdf"), height = 8, width = 8)
ggline(filter(k.medoid.results, balance == "balanced", accuracy != 0 ),
       x = "cell.count", y = "accuracy", group= "lables", color = "labels", facet.by = c("method"),
       title = "balanced-data (no resampling-blast cells are limiting)",
       subtitle = "k-medoids euclidean distance clustering accuracy")
dev.off()



kmed.input.files <- list.files("SCMP/data/analysis-ready/unsupervised-analysis-ready", full.names = TRUE) %>%
  .[grepl("100-|200-|300-|400-|500-|548-|1000-|5000-|10000-",.)]

result.files <- data.frame(kmed.input.files, stringsAsFactors = FALSE) %>%
  dplyr::mutate(cell.counts = gsub(".*processed_|-cells.*", "", kmed.input.files),
                cell.counts = as.numeric(cell.counts),
                method = gsub(".*_|.csv", "", kmed.input.files),
                method = factor(method, levels = c("pca","tsne","umap","phate", "lda")),
                balanced = gsub(".*scatterbodies_|_.*","",kmed.input.files)) %>%
  arrange(cell.counts, method, balanced)


plots.list <- lapply(result.files$kmed.input.files, function(filepath) {
  
  k.medoid.data <- data.table::fread(filepath) 
  k.medoid.data.melt <- k.medoid.data %>%
    gather(key = "method", value="label", -axis1,-axis2,-cell.id,-k.med.cluster,-k.med.id)
  ggscatter(k.medoid.data.melt, x= "axis1",y="axis2", color = "label", facet.by = "method", size = 0.1, 
            subtitle = basename(filepath)) 
})
names(plots.list) <- result.files$kmed.input.files
library(gridExtra)
# library(grid)
pdf(paste0(plot.output.path, "scatterplots_balanced-data_visuals.pdf"), height = 28, width = 30)
do.call("grid.arrange", c(plots.list[grepl("_balanced", names(plots.list))], ncol = 5))
dev.off()

pdf(paste0(plot.output.path, "scatterplots_unbalanced-data_visuals.pdf"), height = 20, width = 32)
do.call("grid.arrange", c(plots.list[grepl("_unbalanced", names(plots.list))], ncol = 5))
dev.off()




