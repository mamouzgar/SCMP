

library(tidyverse)
library(caret)

library(e1071)
library(LiblineaR)

setwd("/Users/mamouzgar/phd-projects")



## pulls relevant columns from processed_... table for k-means clustering
## z-score scaling of x and y axis to ensure fair euclidean comparison
extract_plot_columns <- function(filepath) {
  print(filepath)
  processed.dat <- data.table::fread(filepath)
  
  if (grepl("lda", filepath)) { 
    processed.dat <- dplyr::select(processed.dat, axis1=ld1, axis2=ld2, labels= gate, cell.id)
    data.output <- scale(processed.dat[ , c("axis1","axis2")]) %>%
      data.frame() %>%
      dplyr::mutate(labels = processed.dat$labels,
                    cell.id = processed.dat$cell.id)
  } else if (grepl("pca", filepath)) {
    processed.dat <- dplyr::select(processed.dat, axis1=PC1, axis2=PC2, labels, cell.id)
    data.output <- scale(processed.dat[ , c("axis1","axis2")]) %>%
      data.frame() %>%
      dplyr::mutate(labels = processed.dat$labels,
                    cell.id = processed.dat$cell.id)
  } else if (grepl("phate", filepath)) {
    processed.dat <- dplyr::select(processed.dat, axis1=PHATE1, axis2=PHATE2, labels, cell.id)
    data.output <- scale(processed.dat[ , c("axis1","axis2")]) %>%
      data.frame() %>%
      dplyr::mutate(labels = processed.dat$labels,
                    cell.id = processed.dat$cell.id)
  } else if (grepl("tsne", filepath)) {
    processed.dat <- dplyr::select(processed.dat, axis1=X1, axis2=X2, labels, cell.id)
    data.output <- scale(processed.dat[ , c("axis1","axis2")]) %>%
      data.frame() %>%
      dplyr::mutate(labels = processed.dat$labels,
                    cell.id = processed.dat$cell.id)
  } else if (grepl("umap", filepath)) {
    processed.dat <- dplyr::select(processed.dat, axis1=X1, axis2=X2, labels, cell.id)
    data.output <- scale(processed.dat[ , c("axis1","axis2")]) %>%
      data.frame() %>%
      dplyr::mutate(labels = processed.dat$labels,
                    cell.id = processed.dat$cell.id)
  }
  write.table(data.output, paste0("SCMP/data/analysis-ready/unsupervised-analysis-ready/", basename(filepath)),sep=",",row.names=FALSE,col.names=TRUE)
}




##
## DEPRECATED: DOES NOT ALLOW MORE THAN 65536 OBSERVATIONS 
## run k-mediods clustering
## input data, k : # of clusters
# kmedoids_function <- function(input.data, k = 5, metric = "euclidean", stand=FALSE) { 
#   kmedoids.model <- cluster::pam(input.data[ , c("axis1","axis2")], k, metric = metric, stand )
#   kmedoids.data <- data.frame(cbind(input.data, k.med.cluster = kmedoids.model$clustering), k.med.id = 1:length(kmedoids.model$clustering))
#   kmedoids.results <- data.frame(cbind(k.med.id=kmedoids.model$id.med)) %>%
#     left_join(.,kmedoids.data, by = "k.med.id") %>%
#     dplyr::select(k.med.cluster, k.med.results = labels) %>%
#     left_join(kmedoids.data,., by = "k.med.cluster")
#   
#   return(kmedoids.results)
# }


## run k-mediods clustering
## input data, k : # of clusters
## dist.method : calcualte distance matrix using XXX method (default euclidean)
kmedoids_function <- function(input.data, k = 5, dist.method = "euclidean", k.method = "PAM") { 
  dist.object <- dist(input.data[ , c("axis1","axis2")], method = dist.method)
  kmedoids.model <-WeightedCluster::wcKMedoids(as.matrix(dist.object), k=k, method = k.method)
  kmedoids.data <- data.frame(cbind(input.data, k.med.cluster = kmedoids.model$clustering), k.med.id = 1:length(kmedoids.model$clustering))
  kmedoids.results <- data.frame(cbind(k.med.id=as.numeric(rownames(kmedoids.model$ASW)))) %>%
    left_join(.,kmedoids.data, by = "k.med.id") %>%
    dplyr::select(k.med.cluster, k.med.results = labels) %>%
    left_join(kmedoids.data,., by = "k.med.cluster")
  
  return(kmedoids.results)
}

accuracy_summary <- function(input.data) { 
  
  tab <- with(input.data, table(labels = labels, CORRECT = labels == k.med.results))
  
  accuracy.table <- data.frame(tab) %>%
    tidyr::spread(key = "CORRECT", value = "Freq") %>%
    dplyr::rename(positive=`TRUE`, negative=`FALSE`) %>%
    dplyr::group_by(labels) %>%
    dplyr::mutate(accuracy = positive/(positive+negative)) %>%
    dplyr::ungroup() 
  
  summary.accuracy <- with(accuracy.table, 
                           data.frame(labels = "all.cells", 
                                      negative = sum(negative),
                                      positive=sum(positive)))
  summary.accuracy$accuracy = with(summary.accuracy, positive/(positive+negative))
  
  accuracy.table <- rbind(accuracy.table, summary.accuracy)
  
  return(accuracy.table)
  }



# ggpubr::ggscatter(kmdoids.data %>% group_by(labels) %>% gather(key = "method", value = "group", -cell.id,-axis1,-axis2), x = "axis1", y ="axis2", facet.by = "method",color = "group")


########################
## prepare data model ##
########################
data.sets <- list.files("SCMP/data/analysis-ready", full.names = TRUE) %>% .[grepl("scatterbodies", .)]
lapply(data.sets, function(filepath) { 
  extract_plot_columns(filepath)})

kmed.input.files <- list.files("SCMP/data/analysis-ready/unsupervised-analysis-ready", full.names = TRUE) %>%
  .[grepl("100-|200-|300-|400-|500-|548-|1000-|5000-|10000-",.)]

# kmed.input.files <- list.files("SCMP/data/analysis-ready/unsupervised-analysis-ready", full.names = TRUE) %>%
#   .[grepl("548-",.)]
accuracy.outputs <- lapply(kmed.input.files, function(filepath) {
  print(filepath)
  input.data <- data.table::fread(filepath)
  kmedoids_output<-kmedoids_function(input.data, k = 5, dist.method = "euclidean")
  
  write.table(kmedoids_output, filepath, sep = ",",row.names=FALSE, col.names = TRUE )
  final_results <- accuracy_summary(kmedoids_output)
  final_results$file <- basename(filepath)
  final_results$method = gsub(".*_|.csv","",final_results$file)
  final_results$cell.count = gsub("processed_|-cells.*","", final_results$file)
  final_results$balance = gsub(".*scatterbodies_|_.*.csv", "",final_results$file)
  
  return(final_results)
})

df.accuracy.output <- do.call("bind_rows",accuracy.outputs)
write.table(df.accuracy.output, "SCMP/data/summary-tables/k-medoid-accuracy-results.csv",sep=",",row.names=FALSE, col.names = TRUE)










