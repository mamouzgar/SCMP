


library(KernSmooth)
library(raster)
library(tidyverse)
setwd("/Users/mamouzgar/phd-projects")

###############
## FUNCTIONS ##
###############
##' creates a 100 x 100 pixel grid
create_pixel_grid <- function(xbreaks =100, ybreaks = 100) {
  xbreaks <-xbreaks
  ybreaks <-ybreaks
  pixel.grid = expand.grid(1:xbreaks,1:ybreaks) %>% data.frame(.) %>% 
    rename(x=Var1, y = Var2) %>%
    mutate(pixel = paste(x, y,sep= "."))
  return(pixel.grid)
  }

##' generates a pixel density map for each class label
generate_density_map <- function(data, pixel.grid = pixel.grid, xbreaks = 100, ybreaks = 100) { 

  data_pixel_counts <- lapply(unique(data$labels), function(class.label) { 
    print(class.label)
    data.label <-  data %>% 
      filter(labels == class.label)
    xbin <- cut(data.label$x, xbreaks, include.lowest = TRUE)
    ybin <- cut(data.label$y, ybreaks, include.lowest = TRUE) 
    
    data.label <- data.label %>%
      ungroup() %>%
      mutate(xout = as.numeric(xbin), 
             yout = as.numeric(ybin),
             pixel = paste(xout,yout,sep=".")) %>%
      group_by(labels, pixel) %>%
      summarize(count = n())  %>%
      right_join(.,pixel.grid,by="pixel") %>%
      ungroup() %>%
      mutate(count.0 = ifelse(is.na(count), 0, count),
             percent.0 = count.0 / sum(count.0),
             labels = class.label)
    
  })
  data_pixel_counts <- data_pixel_counts %>% bind_rows(.)
  return(data_pixel_counts)
}

##' calculates the pixel density overlap metric for a class label relative to all other class labels
calculate_pixel_density_metric <- function(data = density_metric_output) { 
  print("calculate_pixel_density_metric")
  data_pixel_intensity <- lapply(unique(data$labels), function(class.label) { 
    print(class.label)
    pixel.label.density <- data %>% 
      ungroup() %>%
      mutate(binary.labels = ifelse(labels == class.label, "class.of.interest", "other")) %>%
      dplyr::select(pixel, x, y, count.0, percent.0, binary.labels) %>%
      group_by(pixel, x, y, binary.labels) %>%
      summarize(count = sum(count.0),
                percent = sum(percent.0)) %>%
      gather(key = "approach", value = "quantity", -pixel,-x,-y,-binary.labels) %>%
      spread(key = "binary.labels", value = "quantity") %>%
      ungroup() %>%
      mutate(density.metric =   (class.of.interest / (other+class.of.interest)),
             # density.metric = ifelse()
             labels = class.label)  %>%
      dplyr::select(pixel,x,y, approach, density.metric, labels)
    return(pixel.label.density)
  })
  print("binding rows for: calculate_pixel_density_metric")
  
  pixel.label.density <- data_pixel_intensity %>% bind_rows(.) %>% na.omit() %>%
    ungroup() %>%
    mutate(density.metric = case_when(is.infinite(density.metric) ~ 1,
                                      TRUE ~ density.metric))
  return(pixel.label.density)
}

aggregate_pixel_density_metric <- function(data = pixel.label.density) { 
  print("aggregating output from calculate_pixel_density_metric")
  
  pixel.label.density <- data %>%
    group_by(approach, labels) %>% 
    summarize(density.summary = sum(density.metric)) %>%
    mutate(density.summary.normalized = density.summary / 10000) ## normalize to the total # of pixels in the pixel grid
  return(pixel.label.density)
}


################################
## calculate density metrics ###
################################
scatterbodies_analysis.ready.tables <- list.files("SCMP/data/analysis-ready/unsupervised-analysis-ready/", full.names = TRUE)
metabolism_analysis.ready.tables <- list.files("SCMP/data/analysis-ready/fh-metabolism", full.names = TRUE) %>% .[grepl("CD8",.)]

analysis.ready.tables <- c(scatterbodies_analysis.ready.tables,metabolism_analysis.ready.tables)


coordinate.density <- list()
pixel.label.density_sc <- list()
pixel.density.metric.output <- lapply(scatterbodies_analysis.ready.tables, function(filepath){
  print(filepath)
  file.name<-basename(filepath)
  
  df <- data.table::fread(filepath) 
  if (!any(colnames(df) %in% c("axis1","PHATE1","PC1"))) { ## skips if the data is missing the appropriate columns
    print("skipped file")
    return() 
  }
  print("calculating density")
  df.density <- df[ , c(1, 2)] 
  df.density$labels <- df$labels
  colnames(df.density) <- c("x","y","labels")
  
  pixel.grid <- create_pixel_grid()
  density_metric_output <- generate_density_map(data = df.density, pixel.grid = pixel.grid)
  coordinate.density <<- dplyr::bind_rows(coordinate.density, density_metric_output %>% mutate(filename = file.name))
  density_metric_output <- calculate_pixel_density_metric(data=density_metric_output)
  pixel.label.density_sc <<-  dplyr::bind_rows(pixel.label.density_sc, density_metric_output %>% mutate(filename = file.name))
  density_metric_output <- aggregate_pixel_density_metric(data=density_metric_output)
  density_metric_output$filename = file.name
  return(density_metric_output)
})
  
  
pixel.counts.metric.bind <- coordinate.density %>%
  dplyr::mutate(cell.count = gsub(".*processed_|-cells_.*","",filename),
                cell.count = as.numeric(cell.count),
                method = gsub(".*_|.csv","", filename),
                distribution = case_when(grepl("_balanced",filename)~"balanced",
                                         grepl("_unbalanced",filename)~"unbalanced"),
                dataset = case_when(grepl("scatterbodies",filename)~"scatterbodies",
                                    grepl("metabolism",filename)~"metabolism")) 
write.table(pixel.counts.metric.bind, "/Users/mamouzgar/phd-projects/SCMP/data/summary-tables/density-metric/scatterbodies-x-y-coordinates_density-metric.csv", row.names = FALSE, sep=",",col.names = TRUE)

sc_pixel.density.metric.bind <- pixel.label.density_sc %>%
  dplyr::mutate(cell.count = gsub(".*processed_|-cells_.*","",filename),
                cell.count = as.numeric(cell.count),
                method = gsub(".*_|.csv","", filename),
                distribution = case_when(grepl("_balanced",filename)~"balanced",
                                         grepl("_unbalanced",filename)~"unbalanced"),
                dataset = case_when(grepl("scatterbodies",filename)~"scatterbodies",
                                    grepl("metabolism",filename)~"metabolism")) 
write.table(sc_pixel.density.metric.bind, "/Users/mamouzgar/phd-projects/SCMP/data/summary-tables/density-metric/scatterbodies-x-y-coordinates_sc_density-metric.csv", row.names = FALSE, sep=",",col.names = TRUE)

pixel.density.metric.bind<- pixel.density.metric.output %>% bind_rows() %>%
  dplyr::mutate(cell.count = gsub("processed_|-cells_.*","",filename),
                cell.count = as.numeric(cell.count),
                method = gsub(".*_|.csv","", filename),
                distribution = case_when(grepl("_balanced",filename) ~"balanced",
                                         grepl("_unbalanced",filename) ~"unbalanced"),
                dataset = case_when(grepl("scatterbodies",filename) ~"scatterbodies",
                                    grepl("metabolism",filename) ~"metabolism"))
write.table(pixel.density.metric.bind, "/Users/mamouzgar/phd-projects/SCMP/data/summary-tables/density-metric/scatterbodies_density-metric.csv",sep=",", row.names = FALSE, col.names = TRUE)


coordinate.density <- list()
pixel.density.metric.output <- lapply(metabolism_analysis.ready.tables, function(filepath){
  print(filepath)
  file.name<-basename(filepath)
  
  df <- data.table::fread(filepath) 
  # if (!any(colnames(df) %in% c("axis1","PHATE1","PC1", "ld1"))) { ## skips if the data is missing the appropriate columns
  #   print("skipped file")
  #   return() 
  # } 
  
  if (any(colnames(df) %in% c("ld1"))) { 
    df <- df %>% dplyr::select(ld1, ld2, labels)
    } 
  print("calculating density")
  df.density <- df[ , c(1, 2)] 
  df.density$labels <- df$labels
  colnames(df.density) <- c("x","y","labels")
  
  pixel.grid <- create_pixel_grid()
  density_metric_output <- generate_density_map(data = df.density, pixel.grid = pixel.grid)
  coordinate.density[[file.name]] <<- density_metric_output %>% mutate(filename = file.name)
  
  density_metric_output <- calculate_pixel_density_metric(density_metric_output)
  density_metric_output$filename = file.name
  return(density_metric_output)
})

pixel.counts.metric.bind <- coordinate.density %>% bind_rows()  %>%
  dplyr::mutate(cell.count = gsub(".*processed_|-cells_.*","",filename),
                cell.count = as.numeric(cell.count),
                method = gsub(".*_|.csv","", filename),
                distribution = case_when(grepl("_balanced",filename)~"balanced",
                                         grepl("_unbalanced",filename)~"unbalanced"),
                dataset = case_when(grepl("scatterbodies",filename)~"scatterbodies",
                                    grepl("metabolism",filename)~"metabolism")) 
write.table(pixel.counts.metric.bind, "/Users/mamouzgar/phd-projects/SCMP/data/summary-tables/density-metric/fh-metabolism-x-y-coordinates_density-metric.csv",sep=",", row.names = FALSE, col.names = TRUE)


pixel.density.metric.bind <- pixel.density.metric.output %>% bind_rows() %>%
  dplyr::mutate(cell.count = gsub(".*processed_|-cells_.*","",filename),
                cell.count = as.numeric(cell.count),
                method = gsub(".*_|.csv","", filename),
                distribution = case_when(grepl("_balanced",filename)~"balanced",
                                         grepl("_unbalanced",filename)~"unbalanced"),
                dataset = case_when(grepl("scatterbodies",filename)~"scatterbodies",
                                    grepl("metabolism",filename)~"metabolism")) 

write.table(pixel.density.metric.bind, "/Users/mamouzgar/phd-projects/SCMP/data/summary-tables/density-metric/fh-metabolism_density-metric.csv",sep=",", row.names = FALSE,col.names = TRUE)

library(patchwork)

###########
## plots ##
###########

density.files <- list.files("SCMP/data/summary-tables/density-metric/", full.names = TRUE)

## SCATTERBODIES

scatter <- data.table::fread(scatterbodies_analysis.ready.tables[21])
scatter.density.sc <- data.table::fread("SCMP/data/summary-tables/density-metric/scatterbodies-x-y-coordinates_sc_density-metric.csv")

scatter.density.sc2 <- scatter.density.sc %>%
  mutate(pixel = as.character(pixel))
  # right_join(.,pixel.grid, by = "pixel") %>%
  # mutate(pixel = as.numeric(pixel),)
df <- data.table::fread( "SCMP/data/summary-tables/density-metric/scatterbodies-x-y-coordinates_density-metric.csv")  
df.plot <- df %>%
  dplyr::mutate(cell.count = gsub(".*processed_|-cells_.*","",filename),
                percent.0 = na_if(percent.0, 0),
                cell.count = as.numeric(cell.count),
                method = gsub(".*_|.csv","", filename),
                distribution = case_when(grepl("_balanced",filename)~"balanced",
                                         grepl("_unbalanced",filename)~"unbalanced"),
                dataset = case_when(grepl("scatterbodies",filename)~"scatterbodies",
                                    grepl("metabolism",filename)~"metabolism")) 

ggscatter(scatter,x = "axis1",y="axis2", color = "labels", nrow=1) 
t <- ggscatter(scatter,x = "axis1",y="axis2", color = "labels", nrow=1) 
t1 <- ggscatter(scatter,x = "axis1",y="axis2", color = "labels", facet.by ="labels", nrow=1) 
t2 <- ggplot(filter(df.plot, filename == "processed_176664-cells_LDA-scatterbodies_unbalanced_minimum_lda.csv"), aes(x, y, fill= count)) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis", na.value = "white") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()
t3 <- ggplot(filter(scatter.density.sc2, filename ==  "processed_176664-cells_LDA-scatterbodies_unbalanced_minimum_lda.csv"), aes(x, y, fill= density.metric )) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis", na.value = "white") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()
pdf("SCMP/data/summary-tables/density-metric/visualization-LDA-example.pdf", width = 14,height= 8)
t + (t1/t2/t3)
dev.off()



scatter <- data.table::fread(scatterbodies_analysis.ready.tables[22])
t <- ggscatter(scatter,x = "axis1",y="axis2", color = "labels", nrow=1) 
t1 <- ggscatter(scatter,x = "axis1",y="axis2", color = "labels", facet.by ="labels", nrow=1) 
t2 <- ggplot(filter(df.plot, filename == "processed_176664-cells_LDA-scatterbodies_unbalanced_minimum_pca.csv"), aes(x, y, fill= count)) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis", na.value = "white") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()
t3 <- ggplot(filter(scatter.density.sc2, filename ==  "processed_176664-cells_LDA-scatterbodies_unbalanced_minimum_pca.csv"), aes(x, y, fill= density.metric )) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis", na.value = "white") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()
pdf("SCMP/data/summary-tables/density-metric/visualization-pca-example.pdf", width = 14,height= 8)
t + (t1/t2/t3)
dev.off()



scatter <- data.table::fread(scatterbodies_analysis.ready.tables[23])
t <- ggscatter(scatter,x = "axis1",y="axis2", color = "labels", nrow=1) 
t1 <- ggscatter(scatter,x = "axis1",y="axis2", color = "labels", facet.by ="labels", nrow=1) 
t2 <- ggplot(filter(df.plot, filename == "processed_176664-cells_LDA-scatterbodies_unbalanced_minimum_phate.csv"), aes(x, y, fill= count)) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis", na.value = "white") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()
t3 <- ggplot(filter(scatter.density.sc2, filename ==  "processed_176664-cells_LDA-scatterbodies_unbalanced_minimum_phate.csv"), aes(x, y, fill= density.metric )) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis", na.value = "white") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()
pdf("SCMP/data/summary-tables/density-metric/visualization-phate-example.pdf", width = 14,height= 8)
t + (t1/t2/t3)
dev.off()



scatter <- data.table::fread(scatterbodies_analysis.ready.tables[24])
t <- ggscatter(scatter,x = "axis1",y="axis2", color = "labels", nrow=1) 
t1 <- ggscatter(scatter,x = "axis1",y="axis2", color = "labels", facet.by ="labels", nrow=1) 
t2 <- ggplot(filter(df.plot, filename == "processed_176664-cells_LDA-scatterbodies_unbalanced_minimum_tsne.csv"), aes(x, y, fill= count)) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis", na.value = "white") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()
t3 <- ggplot(filter(scatter.density.sc2, filename ==  "processed_176664-cells_LDA-scatterbodies_unbalanced_minimum_tsne.csv"), aes(x, y, fill= density.metric )) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis", na.value = "white") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()
pdf("SCMP/data/summary-tables/density-metric/visualization-tsne-example.pdf", width = 14,height= 8)
t + (t1/t2/t3)
dev.off()

scatter <- data.table::fread(scatterbodies_analysis.ready.tables[25])
t <- ggscatter(scatter,x = "axis1",y="axis2", color = "labels", nrow=1) 
t1 <- ggscatter(scatter,x = "axis1",y="axis2", color = "labels", facet.by ="labels", nrow=1) 
t2 <- ggplot(filter(df.plot, filename == "processed_176664-cells_LDA-scatterbodies_unbalanced_minimum_umap.csv"), aes(x, y, fill= count)) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis", na.value = "white") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()
t3 <- ggplot(filter(scatter.density.sc2, filename ==  "processed_176664-cells_LDA-scatterbodies_unbalanced_minimum_umap.csv"), aes(x, y, fill= density.metric )) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis", na.value = "white") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()
pdf("SCMP/data/summary-tables/density-metric/visualization-umap-example.pdf", width = 14,height= 8)
t + (t1/t2/t3)
dev.off()


pixel.density.metric.bind <- data.table::fread( "SCMP/data/summary-tables/density-metric/scatterbodies_density-metric.csv")  

ggline(filter(pixel.density.metric.bind, approach == "count", distribution =="unbalanced"), x = "cell.count", y = "density.summary", 
       color = "labels") +
  facet_grid(method ~ distribution, scales = "free_x") 

pdf("SCMP/data/summary-tables/density-metric/density-metric-summary-example.pdf", width = 5,height= 5)
ggline(filter(pixel.density.metric.bind, approach == "count", distribution =="unbalanced"), x = "cell.count", y = "density.summary.normalized", 
       color = "labels", ylab = "overlap") +
  rotate_x_text(90)+
  # theme_pubclean()+
  facet_grid(method ~ distribution, scales = "free_x") 
dev.off()

ggline(filter(pixel.density.metric.bind, approach == "percent", distribution =="unbalanced"), x = "cell.count", y = "density.summary", 
       color = "labels") +
  facet_grid(method ~ distribution, scales = "free_x") 
ggline(filter(pixel.density.metric.bind, approach == "percent", distribution =="unbalanced"), x = "cell.count", y = "density.summary.normalized", 
       color = "labels") +
  facet_grid(method ~ distribution, scales = "free_x") 

scatter.phate <- data.table::fread(scatterbodies_analysis.ready.tables[23])
t <- ggscatter(scatter.phate,x = "axis1",y="axis2", color = "labels", nrow=1) 
t1 <- ggscatter(scatter.phate,x = "axis1",y="axis2", color = "labels", facet.by ="labels", nrow=1) 
t2 <- ggplot(filter(df.plot, filename == "processed_176664-cells_LDA-scatterbodies_unbalanced_minimum_phate.csv"), aes(x, y, fill= count)) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis", na.value = "white") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()
t3 <- ggplot(filter(df.plot, filename ==  "processed_176664-cells_LDA-scatterbodies_unbalanced_minimum_phate.csv"), aes(x, y, fill= percent.0)) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis", na.value = "white") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()


t + (t1/t2/t3)




## METABOLISM
scatter <- data.table::fread(metabolism_analysis.ready.tables[6])
df <- data.table::fread( "SCMP/data/summary-tables/density-metric/fh-metabolism-x-y-coordinates_density-metric.csv")  
df.plot <- df %>%
  dplyr::mutate(cell.count = gsub(".*processed_|-cells_.*","",filename),
                percent.0 = na_if(percent.0, 0),
                cell.count = as.numeric(cell.count),
                method = gsub(".*_|.csv","", filename),
                distribution = case_when(grepl("_balanced",filename)~"balanced",
                                         grepl("_unbalanced",filename)~"unbalanced"),
                dataset = case_when(grepl("scatterbodies",filename)~"scatterbodies",
                                    grepl("metabolism",filename)~"metabolism")) 

ggscatter(scatter,x = "ld1",y="ld2", color = "labels", nrow=1) 
t1 <- ggscatter(scatter,x = "ld1",y="ld2", color = "labels", facet.by ="labels", nrow=1) 
t2 <- ggplot(filter(df.plot, filename == "c393_CD8naive_processed_10000-cells_metabolism_balanced_c393_lda.csv"), aes(x, y, fill= count)) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis", na.value = "white") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()
t3 <- ggplot(filter(df.plot, filename == "c393_CD8naive_processed_10000-cells_metabolism_balanced_c393_lda.csv"), aes(x, y, fill= percent.0)) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis", na.value = "white") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()
t1/t2/t3


ggline(filter(pixel.density.metric.bind, approach == "count"), x = "cell.count", y = "density.summary", 
       color = "labels") +
  facet_grid(method ~ distribution, scales = "free_x") 
ggline(filter(pixel.density.metric.bind, approach == "count"), x = "cell.count", y = "density.summary.normalized", 
       color = "labels") +
  facet_grid(method ~ distribution, scales = "free_x") 

ggline(filter(pixel.density.metric.bind, approach == "percent"), x = "cell.count", y = "density.summary", 
       color = "labels") +
  facet_grid(method ~ distribution, scales = "free_x") 
ggline(filter(pixel.density.metric.bind, approach == "percent"), x = "cell.count", y = "density.summary.normalized", 
       color = "labels") +
  facet_grid(method ~ distribution, scales = "free_x") 


ggline(filter(pixel.density.metric.bind, approach == "percent"), x = "cell.count", y = "density.summary.normalized", 
       color = "method") +
  facet_grid(~ distribution, scales = "free_x") 
t1 / t2 


ggplot(density_metric_output, aes(x, y, fill= count.0)) +
  geom_tile() +
  viridis::scale_fill_viridis(discrete= FALSE, option = "cividis") +
  facet_wrap(~labels, nrow = 1) +
  ggpubr::theme_pubr()
