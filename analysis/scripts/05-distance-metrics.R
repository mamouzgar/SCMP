
library(gridExtra)
library(ggpubr)
setwd("/Users/mamouzgar/phd-projects")

#### function to calculate the distance from centroid for files outputted from the function "extract_plot_columns" in 03-cluster-analysis script.
## calculates the INTRACLUSTER distance from its centroid
calculate_dist_from_centroid <- function(filepath, writeTable = FALSE) { 
  data.input <- data.table::fread(filepath)
  centroids <- data.input %>% group_by(labels) %>% summarize(x.centroid = mean(axis1), y.centroid = mean(axis2)) %>% ## find centroids
    assign("object.list", value = . , envir = .GlobalEnv)  %>% ## save minimum subsetted table 
    right_join(., data.input, by = "labels") %>% rowwise() %>% mutate(dist.from.centroid = sqrt( (x.centroid-axis1)^2 + (y.centroid-axis2)^2 ),
                                                                      abs.dist.from.centroid = abs(dist.from.centroid)) %>%
    mutate(filepath = basename(filepath))
  
  if (writeTable == TRUE) {
    write.table(centroids, paste0("SCMP/data/summary-tables/dist-from-centroid/dist-from-centroid_", basename(filepath)), sep = ",", row.names = FALSE, col.names = TRUE)
  }
  
  return(centroids)
}


#### function to summarize and aggregate data
summarize_dist_from_centroid <- function(centroids) { 
  
  centroids.summary <- centroids %>% group_by(labels) %>% summarize(mean.dist.from.centroid = mean(abs.dist.from.centroid),
                                                                    median.dist.from.centroid = median(abs.dist.from.centroid),
                                                                    sd.dist.from.centroid = sd(abs.dist.from.centroid),
                                                                    coef.of.variation.from.centroid = sd.dist.from.centroid/mean.dist.from.centroid)  %>%
    mutate(filepath = basename(filepath)) 
  return(centroids.summary)
}



#### function to calculate mahalanobis distance of a set of points with "label A" to a centroid of another set of points "label B"
## tells you how many standard deviations P is (a point) is from the mean of D ( a distribution)
## in other words, it can be a metric for the global relationship between a centroid of "label A" to the distribution of "label B"
## this is often used for outlier detection, and is better than just calculating euclidean distance to centroid because it takes into account the distribution/covariance fo the data
## but the larger the mahalanobis distance for a centroid to a set of points, the more separated it is from that set of points as well
calculate_mahalanobis_dist <- function(filepath) { 
  data.input <- data.table::fread(filepath)
  centroids <- data.input %>% group_by(labels) %>% summarize(x.centroid = mean(axis1), y.centroid = mean(axis2))
  label.combinations <- data.frame(tidyr::crossing(points = centroids$labels, centroid = centroids$labels), stringsAsFactors = FALSE) ## X1 is the point, CENTROID, and X2 are the cells

  mahalanobis_outputs <- apply(label.combinations , MARGIN = 1, function(label.combination) { 
    class.labels = c(label.combination[[1]],label.combination[[2]]) ## 1 =  class label for POINTS, 2 = class label for CENTROID
    
    data.output <- data.input %>% dplyr::filter(., labels %in% class.labels[1]) %>%
      mutate(mahalanobis_dist  = mahalanobis(x = dplyr::select(., axis1, axis2),
                                             center = c(centroids[centroids$labels == class.labels[2], ]$x.centroid, centroids[centroids$labels == class.labels[2], ]$y.centroid),
                                             cov = cov(dplyr::select(., axis1, axis2)) ),
             ## mahalanobis function outputs SQUARED MAHALANOBIS DISTANCE,
             ## take the sqrt of mahalanovis output
             mahalanobis_dist = sqrt(mahalanobis_dist),
             centroid.identity = class.labels[2]) %>%
      spread(key = "centroid.identity", value = "mahalanobis_dist") 
    
    return(data.output)
  })
  
  mahalanobis_summary <- plyr::join_all(mahalanobis_outputs, by = "cell.id", type = "full")
  
  
  ## filename is a global variable
  write.table(mahalanobis_summary, paste0("SCMP/data/summary-tables/mahalanobis-distance/mahalanobis-dist_", basename(filepath)), sep = ",", row.names = FALSE, col.names = TRUE)
}


##########
## MAIN ##
##########

analysis.ready.tables <- list.files("SCMP/data/analysis-ready/unsupervised-analysis-ready/", full.names = TRUE)

### get table of all cells for distances from centroid
dist_from_centroid <- lapply(analysis.ready.tables, function(filepath) { 
  filepath <<- filepath
  print(filepath)
  summary.table <- calculate_dist_from_centroid(filepath, writeTable = TRUE) 
  return(summary.table)
}) 
dist_from_centroid <- do.call("bind_rows", dist_from_centroid)

dist_from_centroid <- dist_from_centroid %>%
  mutate(cell.count = gsub("processed_|-cells.*","", filepath),
         method = gsub(".*_|.csv","",filepath),
         data_balance = gsub(".*scatterbodies_|_.*","", filepath)) %>%
  arrange(data_balance, cell.count, method,labels)
write.table(dist_from_centroid, "SCMP/data/summary-tables/dist-from-centroid_full-analysis.csv",sep=",",row.names = F,col.names = T)





#### calculates distances for every cell from centroid  and summarizes
summary_dist_from_centroid <- lapply(analysis.ready.tables, function(filepath) { 
  filepath <<- filepath
  print(filepath)
  summary.table <- summarize_dist_from_centroid(calculate_dist_from_centroid(filepath) )
  return(summary.table)
}) 
summary_dist_from_centroid <- do.call("bind_rows", summary_dist_from_centroid)

summary_dist_from_centroid <- summary_dist_from_centroid %>%
  mutate(cell.count = gsub("processed_|-cells.*","", filepath),
         method = gsub(".*_|.csv","",filepath),
         data_balance = gsub(".*scatterbodies_|_.*","", filepath)) %>%
  arrange(data_balance, cell.count, method,labels)
  

write.table(summary_dist_from_centroid, "SCMP/data/summary-tables/dist-from-centroid_summary-analysis.csv",sep=",",row.names = F,col.names = T)



## calculates mahalanobis distances for every cluster algorithm 
lapply(analysis.ready.tables, function(filepath) { 
  print(filepath)
  calculate_mahalanobis_dist(filepath)  
}) 



##############
## analysis ##
##############

## summary centroid data 
centroid.dists <- data.table::fread("SCMP/data/summary-tables/dist-from-centroid_summary-analysis.csv")

line_plot <- function(data , x , y, group, color, method, data.filter) { 
  
  p<- ggline(filter(centroid.dists, data_balance == data.filter), 
             x=x, y=y, group=group, color=color,
             facet.by=method, nrow=1, title = data.filter) +
    rotate_x_text(90)

  # p <- ggplot(filter(centroid.dists, data_balance == data.filter),
  #             aes(x=x, y=y, group=group, color=color)) + 
  #   geom_line() +
  #   geom_point() +
  #   # geom_errorbar(aes(ymin=median.dist.from.centroid-sd.dist.from.centroid, ymax=median.dist.from.centroid+sd.dist.from.centroid), width=.2,
  #   #               position=position_dodge(0.05)) +
  #   facet_wrap(~method , nrow = 1) +
  #   theme_pubr()
  return(p)
}
plot.list <- list()
plot.list <- setNames(vector("list", 6), paste0("plot", c(1:6)))

plot.list[["plot1"]] <- line_plot(data=centroid.dists,x="cell.count", y="mean.dist.from.centroid", group = "labels", color ="labels", method = "method", data.filter="balanced")
plot.list[["plot2"]] <- line_plot(data=centroid.dists,x="cell.count", y="median.dist.from.centroid", group = "labels", color ="labels", method = "method", data.filter="balanced")
plot.list[["plot3"]] <- line_plot(data=centroid.dists,x="cell.count", y="coef.of.variation.from.centroid", group = "labels", color ="labels", method = "method", data.filter="balanced")

plot.list[["plot4"]] <- line_plot(data=centroid.dists,x="cell.count", y="mean.dist.from.centroid", group = "labels", color ="labels", method = "method", data.filter="unbalanced")
plot.list[["plot5"]] <- line_plot(data=centroid.dists,x="cell.count", y="median.dist.from.centroid", group = "labels", color ="labels", method = "method", data.filter="unbalanced")
plot.list[["plot6"]] <- line_plot(data=centroid.dists,x="cell.count", y="coef.of.variation.from.centroid", group = "labels", color ="labels", method = "method", data.filter="unbalanced")

pdf("SCMP/analysis/plots/distance-metrics/centroid-distance-summaries.pdf", width = 8, height =18)
do.call("grid.arrange", c(plot.list, nrow = 6))
dev.off()
# p<- ggplot(filter(centroid.dists, data_balance == "balanced"), aes(x=cell.count, y=median.dist.from.centroid, group=labels, color=labels)) + 
#   geom_line() +
#   geom_point() +
#   # geom_errorbar(aes(ymin=median.dist.from.centroid-sd.dist.from.centroid, ymax=median.dist.from.centroid+sd.dist.from.centroid), width=.2,
#   #               position=position_dodge(0.05)) +
#   facet_wrap(~method , nrow = 1, scales = "free_x") +
#   theme_pubr()
# print(p)




## all aggregated centroid data (very large data)
centroid.dists <- data.table::fread("SCMP/data/summary-tables/dist-from-centroid_full-analysis.csv")
ggboxplot(filter(centroid.dists, data_balance == "balanced"), 
          x = "cell.count", y = "abs.dist.from.centroid", 
          color = "labels",
          yscale = "log2",scales = "free_x",
          facet.by = c("method", "labels"), nrow = 2  )

ggboxplot(filter(centroid.dists, data_balance == "unbalanced"), 
         x = "cell.count", y = "abs.dist.from.centroid", 
         add = "median",
         color = "labels",
         # yscale = "log2",
         scales = "free_x",
         facet.by = c("method", "labels"), nrow = 2  )

p<- ggplot(filter(centroid.dists, data_balance == "balanced"), aes(x=cell.count, y=median.dist.from.centroid, group=labels, color=labels)) + 
  geom_line() +
  geom_point() +
  # geom_errorbar(aes(ymin=median.dist.from.centroid-sd.dist.from.centroid, ymax=median.dist.from.centroid+sd.dist.from.centroid), width=.2,
  #               position=position_dodge(0.05)) +
  facet_wrap(~method , nrow = 1, scales = "free_x") +
  theme_pubr()
print(p)












