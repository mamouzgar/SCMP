
library(tidyverse)
library(factoextra)
library(Rtsne)
library(umap)

reticulate::use_python("/Users/mamouzgar/opt/anaconda3/bin/python")
library(phateR)
# library(diffusionMap)
# library(caret)

setwd("/Users/mamouzgar/phd-projects")
source("SCMP/analysis/scripts/00-feature-gen-functions.R") 


cat("reading global variables")
dat <- data.table::fread("SCMP/data/intermediate-data/scATAC/Satpathy_Nature/Satpathy_scATAC_TME_row-cell_col-peak_df_top15-percent.csv", sep = ",", stringsAsFactors = FALSE)
dat$V1 <- NULL
output_filepath <- "SCMP/data/analysis-ready/scATAC-Satpathy-NatBiotech/" ## output directory , include / at the end

colnames(dat)[3:ncol(dat)] <- paste0("peak",colnames(dat)[3:ncol(dat)])
features <- colnames(dat[3:ncol(dat)])
features_summary <- "scATAC_top15percent-var-peaks"
label.levels <- c("CD4.naive", "CD4.Th1","CD4.Th17","CD4.Tfh", "CD4.Treg", "CD8.naive", "CD8.effector","CD8.memory","CD8.ex" )

cell.labels <- data.frame(sample.cell.id = dat$cell.id, cell.id = 1:nrow(dat))
dat$cell.id <- cell.labels$cell.id
# table(dat %>% group_by(labels) %>% dplyr:::sample_n(size = 5) %>% .$labels)
##########
## MAIN ##
##########

balanced_data <- "balanced" ## options are: "balanced" or "unbalanced" or "unbalanced_minimum" you can specify the minimum count
subset_numbers <- c(nrow(dat) ) ## for unbalanced or unbalanced_minimum data
subset_number <- subset_numbers
dat_preproc <<-  dat %>%
  dplyr::select(labels, cell.id, all_of(features)) %>%
  dplyr::filter(labels != "ungated")

# run_algorithms(data.input = data.input)

generate_subset <- function(data, subset_number, label.levels = NULL, balanced_data = "unbalanced", unbalanced.min.count = 50, features, features_summary, output_filepath = output_filepath) { 
  
    dat.clean <- data
    dat.clean$labels <- factor(dat.clean$labels, levels = c(label.levels)) ## genereates factor variable for balanced subsetting
    output_filename <<- paste0(output_filepath, "processed_", subset_number, "-cells_", features_summary,"_", balanced_data, "_")
    data.input <- dat.clean %>%
      ungroup() %>%
      dplyr::select(labels, cell.id, all_of(features))
  return(data.input)
}

data.input <- generate_subset(data = dat_preproc, label.levels = label.levels, balanced_data = balanced_data, subset_number = subset_number, features = features, features_summary = features_summary, output_filepath=output_filepath)
subset_rows <<- caret::createDataPartition(dat.clean[ , "labels"], p = subset_fraction, list =FALSE)

caret::createDataPartition()
subset.cells <- data.input  %>% group_by(labels) %>% sample_n(341) %>% ungroup()
pca.model <- prcomp(dplyr::select(subset.cells, -labels, -cell.id), center = TRUE, scale = TRUE)

pca.data <- data.frame(pca.model$x) %>%
  dplyr::mutate(labels = subset.cells$labels,
                cell.id = subset.cells$cell.id)
output_filename <<- paste0(output_filepath, "processed_", subset_number, "-cells_", features_summary,"_", balanced_data, "_")
algorithm_filename <<-  paste0(output_filename, "pca.csv")
write.table(pca.data, algorithm_filename, sep = ",", row.names = FALSE, col.names = TRUE)
saveRDS(pca.model, paste0(output_filepath, "RDS-files/", "processed_", subset_number, "-cells_", features_summary, "_", balanced_data, "_", "pca.RDS"))
ggpubr::ggscatter(pca.data, x = "PC2", y ="PC2", color = "labels")

umap.config.custom <- umap.defaults
umap.config.custom$a <- 0.3
umap.config.custom$b <- 0.82
umap.model <- umap(dplyr::select(pca.data, -labels, -cell.id)[ , 1:100], config = umap.config.custom)
umap.data <- data.frame(umap.model$layout) %>%
  dplyr::mutate(labels = subset.cells$labels,
                cell.id = subset.cells$cell.id)

algorithm_filename <<-  paste0(output_filename, "umap.csv")
write.table(umap.data, algorithm_filename, sep = ",", row.names = FALSE, col.names = TRUE)
saveRDS(umap.model, paste0(output_filepath, "RDS-files/", "processed_", subset_number, "-cells_", features_summary, "_", balanced_data, "_", "pca.RDS"))
ggpubr::ggscatter(umap.data, x = "X1", y ="X2", color = "labels")


phate.model <- phate(dplyr::select(pca.data, -labels, -cell.id), gamma= 1 )
phate.data <- data.frame(phate.model$embedding) %>%
  dplyr::mutate(labels = subset.cells$labels,
                cell.id = subset.cells$cell.id)
algorithm_filename <<- paste0(output_filename, "phate.csv")
write.table(phate.data, algorithm_filename, sep = ",", row.names = FALSE, col.names = TRUE)
saveRDS(phate.model, paste0(output_filepath, "RDS-files/", "processed_", subset_number, "-cells_", features_summary, "_", balanced_data, "_", "phate.RDS"))
ggpubr::ggscatter(phate.data, x = "PHATE1", y ="PHATE2", color = "labels")


lda.model <- MASS::lda(x = dplyr::select(pca.data, -labels, -cell.id)[ , 1:1000],grouping =pca.data$labels  )
lda.data <- data.frame(lda.model$means %*% lda.model$scaling[ , 1:2])
lda.data$labels = rownames(lda.data)

co = lda.model$scaling[, 1:2, drop=F]
dt <- as.data.table(dplyr::select(pca.data, -labels, -cell.id)[ , 1:1000])
x <- as.matrix(dplyr::select(pca.data, -labels, -cell.id)[ , 1:1000])
for (i in 1:ncol(co)) {
  dt[, eval(paste0("ld", i)):=x %*% co[, i]]
}
dt$labels = pca.data$labels
ggpubr::ggscatter(dt, x = "ld1", y ="ld2", color = "labels", ellipse = TRUE)

d2d.plot <- ggplot(dt, aes(x=ld1,y = ld2)) +
  geom_density_2d(aes(color = labels)) + 
  scale_color_brewer(palette = "Set1") +
  ggpubr::theme_pubr() + 
  # theme(legend.position = "none") +
  # ggtitle(paste(cluster.method, cell_counts,sep = "-")) + 
  theme(axis.title.x=element_blank(),
        axis.text.x=element_blank(),
        axis.ticks.x=element_blank(),
        
        axis.title.y=element_blank(),
        axis.text.y=element_blank(),
        axis.ticks.y=element_blank())
d2d.plot
# balanced_data <- "balanced" ## options are: "balanced" or "unbalanced"
# subset_numbers <- c(100, 200, 300, 400, 548) ## for balanced data, limiter is blast cells

use_previously_subsetted_data <- FALSE ## TRUE : use subsetted data (a training set) previously produced, FALSE : generates new subset of data
data.files_subsetted <- list.files(output_filepath, full.names = TRUE) 

lapply(subset_numbers, function(subset_number) { 
  print(subset_number)
  subset_number <<- subset_number ## assign to global environment to pass to other functions easily
  
  if (use_previously_subsetted_data == TRUE ) { 
    data.subset.filepath <- data.files_subsetted[grepl(paste0(subset_number,"-"), data.files_subsetted)]
    data.subset.reference <- data.table::fread(data.subset.filepath, sep = ",", stringsAsFactors = FALSE)
    
    dat_preproc <<- dat %>%
      dplyr::filter(cell.id %in% data.subset.reference[[1]]) %>%
      dplyr::select("labels", "cell.id", all_of(features)) %>%
      dplyr::filter(labels != "ungated")
  } else if (use_previously_subsetted_data == FALSE) { 
    dat_preproc <<-  dat %>%
      dplyr::select("labels", "cell.id", all_of(features)) %>%
      dplyr::filter(labels != "ungated")
  }
  
  data.input <- generate_subset(data = dat_preproc, label.levels = label.levels, balanced_data = balanced_data, subset_number = subset_number, features = features, features_summary = features_summary, output_filepath=output_filepath)
  
  run_algorithms(data.input = data.input)
})



## get list of sampled cells used in the analysis
samples <- list.files("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready", full.names = TRUE) %>%
  .[grepl("pca", .)]
lapply(samples, function(filepath){ 
  
  # filename<<-gsub("*balanced.*", "balanced.csv", basename(filepath))
  filename<<-gsub("_pca","", basename(filepath))
  output_cellspath <- paste0("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/sampled-cells/", filename)
  write.table(data.table::fread(filepath) %>%
                .[,"cell.id"],output_cellspath, sep=",", row.names = FALSE, col.names = TRUE)
}) 

