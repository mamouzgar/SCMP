#' Master script that Runs various dimensionality reduction or clustering methods on the inputted data and outputs the data.
#' Functions are edited in: 00-feature-gen-functions.R
#' Input variables are edited in: 01-generate-features.R 
#' 01-generate-features.R requires the following inputs:
#' 
#' Requires cell-types to be specified by a column named "gate"
#' Requires a filepath to save outputs
#' Requires a # to subset the data on
#' Requires you to provide a name for the features included in this analysis
rm(list = ls())
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
# set.seed(0)


##########
##########
##########
##########

#############################################################
## LDA-scatterbody channels - different clustering methods ##
#############################################################
## analysis-specific variables ##
#############################################################

cat("reading global variables")
dat <- data.table::fread("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/processed_LDaxes.csv", sep = ",", stringsAsFactors = FALSE)
dat[ , "cell.id"] <- 1:nrow(dat)

## rename label of interest
dat$labels <- dat$gate 
dat$gate <- NULL

output_filepath <- "/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/" ## output directory 
features <- c("WGA_106", "beta_actin", "HP1b", "rRNA", "lamin_A_C", "lamin_B", "lysozyme","VAMP_7", "lactoferrin", "MPO", "serpin_B1", "CD45")
features_summary <- "LDA-scatterbodies"
label.levels <- c("lymphocyte", "neutrophil", "monocyte", "erythroid", "blast" )

# table(dat %>% group_by(labels) %>% dplyr:::sample_n(size = 5) %>% .$labels)
##########
## MAIN ##
##########

balanced_data <- "unbalanced_minimum" ## options are: "balanced" or "unbalanced" or "unbalanced_minimum" you can specify the minimum count
subset_numbers <- c(1000, 5000, 10000, 15000, 20000,50000, 100000, 176664 ) ## for unbalanced or unbalanced_minimum data

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
  
  filename<-gsub("*balanced.*", "balanced.csv", basename(filepath))
  output_cellspath <- paste0("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/sampled-cells/", filename)
  write.table(data.table::fread(filepath) %>%
                .[,"cell.id"],output_cellspath, sep=",", row.names = FALSE, col.names = TRUE)
  }) 















