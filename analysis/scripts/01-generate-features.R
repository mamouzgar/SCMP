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
set.seed(0)


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

output_filepath <- "/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/" ## output directory 

# subset_number <- 500 ## # of cells (rows) to include in analysis
features <- c("WGA_106", "beta_actin", "HP1b", "rRNA", "lamin_A_C", "lamin_B",
              "lysozyme","VAMP_7", "lactoferrin", "MPO", "serpin_B1", "CD45")
features_summary <- "LDA-scatterbodies"

##########
## MAIN ##
##########

## functions
subset_numbers <- c(1000, 5000, 10000, 50000, 100000, 176664 )
use_previously_subsetted_data <- FALSE ## TRUE : use subsetted data (a training set) previously produced, FALSE : generates new subset of data
data.files_subsetted <- list.files(output_filepath, full.names = TRUE) %>%
  .[grepl("2020-09-21", . )]

lapply(subset_numbers, function(subset_number) { 
  subset_number <<- subset_number ## assign to global environment to pass to other functions easily
  
  if (use_previously_subsetted_data == TRUE ) { 
    
    data.subset.filepath <- data.files_subsetted[grepl(subset_number, data.files_subsetted)]
    data.subset.reference <- data.table::fread(data.subset.filepath, sep = ",", stringsAsFactors = FALSE)
    
    data.input <- dat %>%
      dplyr::filter(cell.id %in% data.subset.reference[[1]]) %>%
      column_to_rownames("cell.id") %>%
      dplyr::select("gate", all_of(features))
    
  } else { 
    data.input <- generate_subset(dat = dat, subset_number = subset_number, features = features, features_summary = features_summary, output_filepath=output_filepath)
    }
  
  
  # run_algorithms(data.input = data.input)
  })











