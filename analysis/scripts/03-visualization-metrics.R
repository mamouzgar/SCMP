
#' this script tests different metrics of visualization quality of clusters
#' required package: clv
#' 

#############
## methods ##
#############
#' a) Davies–Bouldin index : a general metric for determining the best separated clusters. Assesses intra-cluster and inter-cluster distances. 
#' b) Dunn index : aims to identify dense and well-separated clusters. Ratio between minimal inter-cluster distance ton maximal intra-cluster distance. 
#' https://www.rdocumentation.org/packages/clv/versions/0.3-2.1
#' c) Compare how well k-means of 4 clusters performs relative to the visual at hand. Assess performance based on: 
#' Purity : a measure of the extent to which clusters contain a single class (doesn't work well with imbalanced data)
#' Rand index : a benchmark classification to see how well the k-means clustering matches the true clusters.
#' rand-index = (TP + TN) / (TP+FP+FN+TN)
#' 

library(clv)
setwd("/Users/mamouzgar/phd-projects")

###########
final_tables_filepath <- "SCMP/data/analysis-ready/final-methods-tables" ## output directory 
final_tables <- list.files(final_tables_filepath, full.names = TRUE)
  

data_subset <- data.table::fread(final_tables[5])
head(data_subset)
data_metric <- data_subset %>%
  select(PC1,PC2)

km.res <- kmeans(data_metric, 4, iter.max = 10, nstart = 4)
km.res <- bind_cols(data_subset, cluster = km.res$cluster) %>%
  select(PC1, PC2, gate, cluster) %>%
  gather(key = "method", value = "value", -PC1,-PC2)



ggplot(km.res, aes(x=PC1,y = PC2)) +
  geom_density_2d(aes(color = value)) + 
  scale_color_brewer(palette = "Set1") +
  facet_grid(~method) +
  theme_pubr() + 
  # theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
        
        


data_subset <- data.table::fread("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/processed_LDaxes.csv")


data_metric <- data_subset %>%
  filter(gate != "ungated") %>%
  select(ld1,ld2,gate)

km.res <- kmeans(select(data_metric, -gate), 5, iter.max = 10000, nstart = 5)
comb.res <- bind_cols(data_metric, cluster = km.res$cluster) %>%
  select(axis1=ld1, axis2=ld2, gate, cluster) %>%
  gather(key = "method", value = "value", -axis1,-axis2)



ggplot(comb.res, aes(x=axis1,y = axis2)) +
  # geom_density_2d(aes(color = value)) +
  geom_point(aes(color = value)) +
  # scale_color_brewer(palette = "Set1") +
  facet_grid(~method) +
  theme_pubr() + 
  # theme(legend.position = "none") +
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())

data_metric$gate
km.res$cluster
(sum(pred==true, na.rm=T) + sum(is.na(pred) & is.na(true))) / length(pred)



