
rm(list = ls())
library(cluster) 
require(tidyverse)
library(viridis)
setwd("/Users/mamouzgar/phd-projects")



###################
## scatterbodies ##
###################
scatterbodies.files <- list.files("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready", full.names = TRUE) %>%
  .[grepl(".csv", .)] %>%
  .[!grepl("processed_LDaxes.csv",.)] %>%
  .[!grepl("176664",.)] %>%
  .[!grepl("50000",.)] %>%
  .[!grepl("_1e",.)]

silh.result<- list()
silh.results <- lapply(scatterbodies.files, function(filepath) {
  print(filepath)
  filename <- basename(filepath)


  df <- data.table::fread(filepath)
  if (grepl("lda.csv", filepath)) {
    df <- df %>% rename(X1=ld1,X2=ld2, labels = gate)
  } else if ((grepl("pca.csv", filepath))) {
    df <- df %>% rename(X1=PC1,X2=PC2)
  } else if ((grepl("phate.csv", filepath))) {
    df <- df %>% rename(X1=PHATE1,X2=PHATE2)
  } else if ((grepl("tsne.csv", filepath))) {
    ## OK
  } else if ((grepl("umap.csv", filepath))) {
    ## OK
  }

  class.labels <- data.frame(labels = df$labels) %>%
    dplyr::mutate(labels = factor(labels , levels = c("lymphocyte","neutrophil","monocyte","erythroid","blast")),
                  clusters = as.numeric(labels))

  dist.matrix <- dplyr::select(df, X1,X2) %>%
    dist(., method = "euclidean")%>%
    as.matrix(.)
  silhoutte.output <- silhouette(x = class.labels$clusters, dmatrix = dist.matrix, do.clus.stat = TRUE, do.n.k=TRUE )


  silh.result.all<- df %>% dplyr::select(X1, X2, cell.id, labels ) %>%
    bind_cols(., data.frame(cluster = silhoutte.output[,1], neighbor = silhoutte.output[,2], sil_width = silhoutte.output[,3], filename = filename ) )

  sum.sil <- summary(silhoutte.output)
  silh.result.sum  <- data.frame(unique(class.labels), silh.width = sum.sil$clus.avg.widths, filename=filename)


  silh.result[["all"]] <- silh.result.all
  silh.result[["summary"]] <- silh.result.sum
  return(silh.result)

})

silh.results <- do.call("bind_rows",silh.results)
silh.results <- silh.results %>% mutate(cell.count = gsub(".*processed_|-cells.*","", filename),
                                          method = gsub(".*_|.csv","",filename),
                                          data_balance = gsub(".*scatterbodies_|_.*","", filename))

silh.results.all <- silh.results %>% dplyr::select(X1, X2, cell.id, labels, cluster, neighbor, sil_width, filename, cell.count, method, data_balance) %>% na.omit() %>%
  dplyr::mutate(cell.count = as.numeric(cell.count))
write.table(silh.results.all, "SCMP/data/summary-tables/scatterbodies-silhouette-metric-results-individual-cells.csv",sep=",",row.names = FALSE, col.names = TRUE)

silh.results.summary <- silh.results %>% dplyr::select(labels, clusters, silh.width, filename, cell.count, method, data_balance) %>% na.omit()
write.table(silh.results.summary, "SCMP/data/summary-tables/scatterbodies-silhouette-metric-results-summary.csv",sep=",",row.names = FALSE, col.names = TRUE)


silh.results.all <- data.table::fread("SCMP/data/summary-tables/scatterbodies-silhouette-metric-results-individual-cells.csv")




## balanced data
pdf("/Users/mamouzgar/phd-projects/SCMP/analysis/plots/silhouette-scores/scatterbodies-silhouette-results-balanced.pdf", width= 5,height = 5)
ggpubr::ggviolin(filter(silh.results.all, data_balance == "balanced"),
                  x = "cell.count", y = "sil_width",
                 subtitle = "scatterbodies: balanced",
                 add = "mean",
                 xlab = "cell count for each class",
                 add.params = list(size = 0.2),
                 # color = "labels",
                 facet.by = "method")+
  ggpubr::rotate_x_text(90)
dev.off()

pdf("/Users/mamouzgar/phd-projects/SCMP/analysis/plots/silhouette-scores/scatterbodies-silhouette-results-balanced-celltypes.pdf", width= 6,height = 6)
ggpubr::ggline(filter(silh.results.all, data_balance == "balanced"),
               x = "cell.count", y = "sil_width",
               legend = "right",
               subtitle = "scatterbodies: balanced: mean_sd",
               xlab = "cell count for each class",
               add = "mean_sd",
               color = "labels",
               facet.by = "method")+
  ggpubr::rotate_x_text(90)
dev.off()

## unbalanced data
pdf("/Users/mamouzgar/phd-projects/SCMP/analysis/plots/silhouette-scores/scatterbodies-silhouette-results-unbalanced.pdf", width= 5,height = 5)
ggpubr::ggviolin(filter(silh.results.all, data_balance == "unbalanced"),
                 x = "cell.count", y = "sil_width",
                 xlab = "total cell count",
                 subtitle = "scatterbodies: unbalanced",
                 add = "mean",
                 add.params = list(size = 0.2),
                 # color = "labels",
                 facet.by = "method") +
  ggpubr::rotate_x_text(90)
dev.off()

pdf("/Users/mamouzgar/phd-projects/SCMP/analysis/plots/silhouette-scores/scatterbodies-silhouette-results-unbalanced-celltypes.pdf", width= 6,height = 6)
ggpubr::ggline(filter(silh.results.all, data_balance == "unbalanced"),
               x = "cell.count", y = "sil_width",
               subtitle = "scatterbodies: unbalanced: mean+sd",
               xlab = "total cell count",
               add = "mean_sd",
               legend = "right",
               color = "labels",
               facet.by = "method") +
  ggpubr::rotate_x_text(90)
dev.off()



## silhouette width plots
pdf("/Users/mamouzgar/phd-projects/SCMP/analysis/plots/silhouette-scores/scatterbodies-silhouette-results-balanced-scatter-plots-548cells.pdf", width= 12,height = 3)
ggpubr::ggscatter(filter(silh.results.all, cell.count == 548),
                  x = "X1", y = "X2", color = "labels", 
                  subtitle = "scatterbodies: balanced 548 cells each",
                  facet.by = "method", scales = "free", 
                  size = 0.1,
                  nrow = 1,
                  legend = "right"
) +
  guides(color = guide_legend(override.aes = list(size=5))) +
  scale_color_viridis(discrete= TRUE, option = "viridis") 
dev.off()

pdf("/Users/mamouzgar/phd-projects/SCMP/analysis/plots/silhouette-scores/scatterbodies-silhouette-results-balanced-scatter-silhouette.values-548cells.pdf", width= 12,height = 3)
ggpubr::ggscatter(filter(silh.results.all, cell.count == 548),
                  x = "X1", y = "X2",
                  subtitle = "scatterbodies: balanced 548 cells each",
                  # add = "mean",
                  # xlab = "cell count for each class",
                  # add.params = list(size = 0.2),
                  color = "sil_width",
                  scales = "free",
                  size = 0.1,
                  nrow = 1,
                  facet.by = "method",
                  legend = "right"
) +
  scale_color_viridis(option = "cividis")
dev.off()

pdf("/Users/mamouzgar/phd-projects/SCMP/analysis/plots/silhouette-scores/scatterbodies-silhouette-results-unbalanced-scatter-plots-20000cells.pdf", width= 12,height = 3)
ggpubr::ggscatter(filter(silh.results.all, cell.count == "20000"),
                  x = "X1", y = "X2", color = "labels", 
                  subtitle = "scatterbodies: unbalanced 20000 cells total",
                  facet.by = "method", scales = "free", 
                  size = 0.1,
                  nrow = 1,
                  legend = "right"
) +
  guides(color = guide_legend(override.aes = list(size=5))) +
  scale_color_viridis(discrete= TRUE, option = "viridis") 
dev.off()

pdf("/Users/mamouzgar/phd-projects/SCMP/analysis/plots/silhouette-scores/scatterbodies-silhouette-results-unbalanced-scatter-silhouette.values-20000cells.pdf.pdf", width= 12,height = 3)
ggpubr::ggscatter(filter(silh.results.all, cell.count == "20000"),
                  x = "X1", y = "X2",
                  subtitle = "scatterbodies: unbalanced 20000 cells total",
                  # add = "mean",
                  # xlab = "cell count for each class",
                  # add.params = list(size = 0.2),
                  color = "sil_width",
                  scales = "free",
                  size = 0.1,
                  nrow = 1,
                  facet.by = "method",
                  legend = "right"
) +
  scale_color_viridis(option = "cividis")
dev.off()


###################
## fh-metabolism ##
###################

metabolism.files <- list.files("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/fh-metabolism", full.names = TRUE)  %>%
  .[grepl(".csv", .)] %>%
  .[grepl("c393_CD8naive", .)] %>%
  .[grepl("1000-", .)] 
  
  
silh.result<- list()
silh.results <- lapply(metabolism.files, function(filepath) { 
  print(filepath)
  filename <- basename(filepath)
  
  
  df <- data.table::fread(filepath)
  if (grepl("lda.csv", filepath)) { 
    df <- df %>% rename(X1=ld1,X2=ld2)
  } else if ((grepl("pca.csv", filepath))) {
    df <- df %>% rename(X1=PC1,X2=PC2)
  } else if ((grepl("phate.csv", filepath))) {
    df <- df %>% rename(X1=PHATE1,X2=PHATE2)
  } else if ((grepl("tsne.csv", filepath))) {
    ## OK
  } else if ((grepl("umap.csv", filepath))) {
    ## OK
  }
  
  class.labels <- data.frame(labels = df$labels) %>%
    dplyr::mutate(labels = factor(labels , levels = c( "CD8_naive_0", "CD8_naive_1", "CD8_naive_2", "CD8_naive_3", "CD8_naive_4","CD8_naive_5")),
                  clusters = as.numeric(labels))
  
  dist.matrix <- dplyr::select(df, X1,X2) %>%
    dist(., method = "euclidean")%>%
    as.matrix(.)
  silhoutte.output <- cluster::silhouette(x = class.labels$clusters, dmatrix = dist.matrix, do.clus.stat = TRUE, do.n.k=TRUE )
  
  
  silh.result.all<- df %>% dplyr::select(X1, X2, cell.id, labels ) %>%
    bind_cols(., data.frame(cluster = silhoutte.output[,1], neighbor = silhoutte.output[,2], sil_width = silhoutte.output[,3], filename = filename ) )
  
  sum.sil <- summary(silhoutte.output)
  silh.result.sum  <- data.frame(unique(class.labels), silh.width = sum.sil$clus.avg.widths, filename=filename)
  
  
  silh.result[["all"]] <- silh.result.all
  silh.result[["summary"]] <- silh.result.sum
  return(silh.result)
  
})

silh.results <- do.call("bind_rows",silh.results)
silh.results <- silh.results %>% mutate(cell.count = gsub(".*processed_|-cells.*","", filename),
                                        method = gsub(".*_|.csv","",filename),
                                        data_balance = gsub(".*metabolism_|_.*","", filename),
                                        neighbor = factor(neighbor, levels = c(1,2,3,4,5,6)))

silh.results.all <- silh.results %>% dplyr::select(X1, X2, cell.id, labels, cluster, neighbor, sil_width, filename, cell.count, method, data_balance) %>% na.omit() %>%
  dplyr::mutate(cell.count = as.numeric(cell.count))
write.table(silh.results.all, "SCMP/data/summary-tables/fh-metabolism-CD8naive-silhouette-metric-results-individual-cells.csv",sep=",",row.names = FALSE, col.names = TRUE)

silh.results.summary <- silh.results %>% dplyr::select(labels, clusters, silh.width, filename, cell.count, method, data_balance) %>% na.omit()
write.table(silh.results.summary, "SCMP/data/summary-tables/fh-metabolism-CD8naive-silhouette-metric-results-summary.csv",sep=",",row.names = FALSE, col.names = TRUE)


silh.results.all <- data.table::fread("SCMP/data/summary-tables/fh-metabolism-CD8naive-silhouette-metric-results-individual-cells.csv")
## balanced data
pdf("/Users/mamouzgar/phd-projects/SCMP/analysis/plots/silhouette-scores/fh-metabolism-silhouette-results-balanced-violin.pdf", width= 5,height = 5)
ggpubr::ggviolin(filter(silh.results.all, data_balance == "balanced"),
                 x = "labels", y = "sil_width",
                 subtitle = "fh-metabolism: balanced",
                 add = "mean",
                 # xlab = "cell count for each class",
                 add.params = list(size = 0.2),
                 color = "labels",
                 facet.by = "method"
                 )+
  ggpubr::rotate_x_text(90) +
  scale_color_viridis(discrete= TRUE, option = "viridis") 
dev.off()


pdf("/Users/mamouzgar/phd-projects/SCMP/analysis/plots/silhouette-scores/fh-metabolism-silhouette-results-balanced-scatter-plots.pdf", width= 12,height = 3)
ggpubr::ggscatter(filter(silh.results.all, data_balance == "balanced"),
                  x = "X1", y = "X2", color = "labels", 
                  subtitle = "CD8 Naive: 1000 cells per label",
                  facet.by = "method", scales = "free", 
                  size = 0.1,
                  nrow = 1,
                  legend = "right"
                  ) +
  guides(color = guide_legend(override.aes = list(size=5))) +
  scale_color_viridis(discrete= TRUE, option = "viridis") 
dev.off()

pdf("/Users/mamouzgar/phd-projects/SCMP/analysis/plots/silhouette-scores/fh-metabolism-silhouette-results-balanced-scatter-silhouette.values.pdf", width= 12,height = 3)
ggpubr::ggscatter(filter(silh.results.all, data_balance == "balanced"),
                 x = "X1", y = "X2",
                 subtitle = "CD8 Naive: 1000 cells per label",
                 # add = "mean",
                 # xlab = "cell count for each class",
                 # add.params = list(size = 0.2),
                 color = "sil_width",
                 scales = "free",
                 size = 0.1,
                 nrow = 1,
                 facet.by = "method",
                 legend = "right"
) +
  scale_color_viridis(option = "cividis")
dev.off()



pdf("/Users/mamouzgar/phd-projects/SCMP/analysis/plots/silhouette-scores/fh-metabolism-silhouette-results-balanced-scatter-plots-neighbor.pdf", width= 12,height = 3)
ggpubr::ggscatter(filter(silh.results.all, data_balance == "balanced"),
                  x = "X1", y = "X2", color = "neighbor", 
                  subtitle = "CD8 Naive: 1000 cells per label",
                  facet.by = "method", scales = "free", 
                  size = 0.1,
                  nrow = 1,
                  legend = "right",
                  # font.legend = list(size = 20)
                  ) +
  guides(color = guide_legend(override.aes = list(size=5))) +
  scale_color_viridis(discrete= TRUE, option = "viridis") 
dev.off()



