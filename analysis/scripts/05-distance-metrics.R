


## function to calculate the distance from centroid for files outputted from the function "extract_plot_columns" in 03-cluster-analysis script.
## calculates the INTRACLUSTER distance from its centroid
calculate_dist_from_centroid <- function(filepath) { 
  data.input <- data.table::fread(filepath)
  centroids <- data.input %>% group_by(labels) %>% summarize(x.centroid = mean(axis1), y.centroid = mean(axis2)) %>% ## find centroids
    assign("object.list", value = . , envir = .GlobalEnv)  %>% ## save minimum subsetted table 
    right_join(., data.input, by = "labels") %>% rowwise() %>% mutate(dist.from.centroid = sqrt( (x.centroid-axis1)^2 + (y.centroid-axis2)^2 ),
                                                                      abs.dist.from.centroid = abs(dist.from.centroid))
  return(centroids)
}


## function to calculate mahalanobis distance of a set of points with "label A" to a centroid of another set of points "label B"
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




analysis.ready.tables <- list.files("SCMP/data/analysis-ready/unsupervised-analysis-ready/", full.names = TRUE)
lapply(analysis.ready.tables, function(filepath) { 
  print(filepath)
  calculate_mahalanobis_dist(filepath)  
  }) 
