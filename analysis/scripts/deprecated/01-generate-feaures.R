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
library(diffusionMap)
set.seed(0)

setwd("/Users/mamouzgar/phd-projects")
source("SCMP/analysis/scripts/00-feature-gen-functions.R") 


###############################################################
## (1) LDA-matched channels - the key channels of interest
###############################################################

######################
## GLOBAL VARIABLES ##
######################
cat("reading global variables")

dat <- data.table::fread("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/processed_LDaxes.csv", sep = ",", stringsAsFactors = FALSE)
dat[ , "cell.id"] <- 1:nrow(dat)

output_filepath <- "/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/" ## output directory 
subset_number <- 50000 ## # of cells (rows) to include in analysis

## relevant features in analysis 
# CD.Ig.markers <- c("IgL", "IgK", "CD235ab", "CD71", "CD61", "CD3", "CD8", "CD2", "CD5", "CD4", "CD7", "CD11c", "CD23", "CD123", "CD56", "CD45", "CD10", "CD13", "CD117", "CD34", "CD20", "CD19", "CD22", "CD79a", "CD15", "CD33", "CD14", "CD64", "CD16" , "CD38" )
CD.Ig.markers <- colnames(dat)[grepl("CD|Ig|HLA", colnames(dat))]
scatterbodies <- c("MPO","WGA_102", "WGA_104", "WGA_105", "WGA_106", "WGA_108", "WGA_110", "barium", "lactoferrin", "rRNA", "HP1b", "lamin_B", "VAMP_7", "lysozyme", "serpin_B1", "lamin_A_C", "beta_actin" )
all.markers <- c(CD.Ig.markers, scatterbodies)
# LDA-analysis matched channels
channels <- c("WGA_106", "beta_actin", "HP1b", "rRNA", "lamin_A_C", "lamin_B",
              "lysozyme","VAMP_7", "lactoferrin", "MPO", "serpin_B1", "CD45")

## generate NEW subset
dat.clean <- dat %>%
  dplyr::filter(gate != "ungated") %>%
  dplyr::sample_n(subset_number)  %>%
  tibble::column_to_rownames("cell.id") %>%
  dplyr::select( -Time, -Event_length, -Bead_1, -DNA_1, -DNA_2, -Viability, -ld1,-ld2,-beadDist,-file) ## remove irrelevant markers - confirm with David

## write cell.ids used in this analysis
write.table(rownames(dat.clean), paste0(output_filepath, Sys.Date(), "_analyzed-cells_", subset_number, "-cells.csv"), sep=",",row.names = FALSE, col.names = TRUE)

## use previous subset
analyzed.cells <- data.table::fread("SCMP/data/analysis-ready/2020-09-18_analyzed-cells_50000-cells.csv")
dat.clean <- dat %>%
  dplyr::filter(gate != "ungated") %>%
  dplyr::filter(cell.id %in% analyzed.cells$x) %>%
  tibble::column_to_rownames("cell.id") %>%
  dplyr::select( -Time, -Event_length, -Bead_1, -DNA_1, -DNA_2, -Viability, -ld1,-ld2,-beadDist,-file) ## remove irrelevant markers - confirm with David

cat("beginning analysis")

#################################
## (1) all markers
#################################
## analysis-specific variables ##
#################################
## data characteristics and output
features_included <- "all-markers"
relevant_features <- all.markers
output_filename <- paste0(output_filepath, "processed_", subset_number, "-cells_", features_included, "_")

data.input <- dat.clean %>%
  dplyr::select("gate", all_of(relevant_features))
##########
## MAIN ##
##########
pca_function(data.input)
tsne_function(data.input)
umap_function(data.input)
phate_function(data.input, gamma.value = 0) ## reticulate::use_python("/Users/mamouzgar/opt/anaconda3/bin/python") ## solution if phate does not work
# dmap_function(data.input) ## not ready




#################################
## (2) CD (and Ig) markers
#################################
## analysis-specific variables ##
#################################
## data characteristics and output
features_included <- "CD-markers"
relevant_features <- CD.Ig.markers
output_filename <- paste0(output_filepath, "processed_", subset_number, "-cells_", features_included, "_")

data.input <- dat.clean %>%
  dplyr::select("gate", all_of(relevant_features))
##########
## MAIN ##
##########
pca_function(data.input)
tsne_function(data.input)
umap_function(data.input)
phate_function(data.input, gamma.value = 0) ## reticulate::use_python("/Users/mamouzgar/opt/anaconda3/bin/python") ## solution if phate does not work
# dmap_function(data.input) ## not ready


#################################
## (3) scatterbodies
#################################
## analysis-specific variables ##
#################################
## data characteristics and output
features_included <- "scatterbodies"
relevant_features <- scatterbodies
output_filename <- paste0(output_filepath, "processed_", subset_number, "-cells_", features_included, "_")

data.input <- dat.clean %>%
  dplyr::select("gate", all_of(relevant_features))
##########
## MAIN ##
##########
pca_function(data.input)
tsne_function(data.input)
umap_function(data.input)
phate_function(data.input, gamma.value = 0) ## reticulate::use_python("/Users/mamouzgar/opt/anaconda3/bin/python") ## solution if phate does not work
# dmap_function(data.input) ## not ready


#################################
## (4) LDA-matched channels - the key channels of interest
#################################
## analysis-specific variables ##
#################################
## data characteristics and output
features_included <- "LDA-matched-channels"
relevant_features <- channels
output_filename <- paste0(output_filepath, "processed_", subset_number, "-cells_", features_included, "_")

data.input <- dat.clean %>%
  dplyr::select("gate", all_of(relevant_features))
##########
## MAIN ##
##########
pca_function(data.input)
tsne_function(data.input)
umap_function(data.input)
phate_function(data.input, gamma.value = 0) ## reticulate::use_python("/Users/mamouzgar/opt/anaconda3/bin/python") ## solution if phate does not work
# dmap_function(data.input) ## not ready


