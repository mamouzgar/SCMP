
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

output_filepath <- "SCMP/data/analysis-ready/scATAC-Satpathy-NatBiotech/" ## output directory , include / at the end
features <- colnames(dat[3:ncol(dat)])
features_summary <- "scATAC_top15percent-var-peaks"
label.levels <- c("CD4.naive", "CD4.Th1","CD4.Th17","CD4.Tfh", "CD4.Treg", "CD8.naive", "CD8.effector","CD8.memory","CD8.ex" )

cell.labels <- data.frame(sample.cell.id = dat$cell.id, cell.id = 1:nrow(dat))
dat$cell.id <- cell.labels$cell.id
# table(dat %>% group_by(labels) %>% dplyr:::sample_n(size = 5) %>% .$labels)
##########
## MAIN ##
##########

balanced_data <- "unbalanced" ## options are: "balanced" or "unbalanced" or "unbalanced_minimum" you can specify the minimum count
subset_numbers <- c(nrow(dat) ) ## for unbalanced or unbalanced_minimum data
subset_number <- subset_numbers
dat_preproc <<-  dat %>%
  dplyr::select("labels", "cell.id", all_of(features)) %>%
  dplyr::filter(labels != "ungated")
data.input <- generate_subset(data = dat_preproc, label.levels = label.levels, balanced_data = balanced_data, subset_number = subset_number, features = features, features_summary = features_summary, output_filepath=output_filepath)
# run_algorithms(data.input = data.input)

pca_output <- pca_function(data.input)
generate_final_datatables(pca_output)

tsne_output <- tsne_function(data.input, pca.prior = TRUE)
generate_final_datatables(tsne_output)

phate_output <- uwot_umap_function(data.input, gamma.value = 0)
generate_final_datatables(phate_output)

umap_output <- umap_function(data.input)
generate_final_datatables(umap_output)

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

