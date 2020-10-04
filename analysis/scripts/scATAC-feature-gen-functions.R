#' Author: Meelad Amouzgar
#' Goal:
#' Runs various dimensionality reduction or clustering methods on the inputted data and outputs the data.
#' This script does not need to be edited unless changes to the fundamental methods want to be made.
#' 
#' Input variables are edited in: 01-generate-features.R and this script requires the following inputs:
#' Requires cell-types to be specified by a column named "labels"
#' Requires a filepath to save outputs
#' Requires a # to subset the data on
# library(tidyverse)
# library(factoextra)
# library(Rtsne)
# library(umap)
# library(phateR)
# library(diffusionMap)
# set.seed(0)
##################################
##################################
##################################
## FEATURE-GENERATING FUNCTIONS ##
##################################
##################################
##################################

#########################################
## CLUSTERING OR DIM-REDUCTION METHODS ## 
#########################################

## PCA
## data is centered and scaled
pca_function <- function(input.data) { 
  pca.model <- prcomp(dplyr::select(input.data, -labels, -cell.id), center = TRUE, scale = TRUE)
  pca.data <- data.frame(pca.model$x) %>%
    dplyr::mutate(labels = input.data$labels,
                  cell.id = input.data$cell.id)
  
  algorithm_filename <<-  paste0(output_filename, "pca.csv")
  write.table(pca.data, algorithm_filename, sep = ",", row.names = FALSE, col.names = TRUE)
  saveRDS(pca.model, paste0(output_filepath, "RDS-files/", "processed_", subset_number, "-cells_", features_summary, "_", balanced_data, "_", "pca.RDS"))
  
  print("pca is complete")
  return(pca.data)
}


## tsne (with PCA)
## note: does NOT use Rtsne's internal PCA by default, generates 2-dimensions. 
## Without PCA: Scales and normalized the data 
## With PCA: scales and centers the data prior to PCA, uses to
## perplexity is defaulted to 25 if arguement is not provided
## potential future added functionality: perplexity is defaulted to 5% of data size if arguement is not provided
## other features set to default 
tsne_function <- function(input.data, pca.prior = FALSE, perplexity.value = NULL, initial_dims = 50)  { 
  if(is.null(perplexity.value)) { 
    # perplexity.value  <- nrow(input.data)*0.05
    perplexity.value  <- 25
  }
  
  if (pca.prior == FALSE) { 
    tsne.model <- Rtsne(dplyr::select(input.data, -labels, -cell.id), dims=2, 
                        perplexity= perplexity.value , 
                        pca = FALSE, pca_center = FALSE, pca_scale = TRUE, normalize = TRUE, 
                        theta=0.5,  max_iter = 1000, check_duplicates = FALSE)
  } else if (pca.prior == TRUE) { 
    tsne.model <- Rtsne(dplyr::select(input.data, -labels, -cell.id), dims=2, 
                        initial_dims = initial_dims,
                        perplexity= perplexity.value , 
                        pca = TRUE, pca_center = TRUE, pca_scale = TRUE, normalize = TRUE, 
                        theta=0.5,  max_iter = 1000, check_duplicates = FALSE)
  }
  
  tsne.model <- Rtsne(dplyr::select(input.data, -labels), dims=2, 
                      perplexity= perplexity.value , 
                      pca = TRUE, pca_center = TRUE, pca_scale = TRUE, normalize = TRUE, 
                      theta=0.5,  max_iter = 1000, check_duplicates = FALSE)
  tsne.data <- data.frame(tsne.model$Y) %>%
    dplyr::mutate(labels = input.data$labels,
                  cell.id = input.data$cell.id)
  
  algorithm_filename <<- paste0(output_filename, "tsne.csv")
  write.table(tsne.data, algorithm_filename, sep = ",", row.names = FALSE, col.names = TRUE)
  saveRDS(tsne.model, paste0(output_filepath, "RDS-files/", "processed_", subset_number, "-cells_", features_summary, "_", balanced_data, "_", "tsne.RDS"))
  
  print("tsne is complete")
  return(tsne.data)
}


## UMAP
umap_function <- function(input.data, custom.config = NULL )  {
  
  ## umap parameter tuning
  if (is.null(custom.config)) { 
    custom.config <- umap.defaults
    custom.config$a <- 0.3
    custom.config$b <- 0.82
  }
  
  umap.model <- umap(dplyr::select(input.data, -labels, -cell.id), config = custom.config)
  umap.data <- data.frame(umap.model$layout) %>%
    dplyr::mutate(labels = input.data$labels,
                  cell.id = input.data$cell.id)
  
  algorithm_filename <<- paste0(output_filename, "umap.csv")
  write.table(umap.data, algorithm_filename,sep = ",", row.names = FALSE, col.names = TRUE)
  saveRDS(umap.model, paste0(output_filepath, "RDS-files/", "processed_", subset_number, "-cells_", features_summary, "_", balanced_data, "_", "umap.RDS"))
  
  print("umap is complete")
  return(umap.data)
}



## PHATE
## knn is selected using k=sqrt(n) unless otherwise stated
## if a decay.value is supplied, it uses alpha for estimating graph connectivity instead of knn
## gamma is an informational distance metric is defaulted to 1. Change to 0 (or tune to another #) if points are concentrated in one spot
phate_function <-  function(input.data, knn.value = NULL, decay.value = NULL, gamma.value = NULL)  { 
  if(is.null(knn.value)) { 
    knn.value  <- round(sqrt(nrow(input.data)))
  }
  
  if(is.null(decay.value)) { 
    decay.value  <- NULL
  } else { 
    knn.value <- NULL
  }
  
  if(is.null(gamma.value)) { 
    gamma.value  <- 1
  } 
  
  
  phate.model <- phate(dplyr::select(input.data, -labels, -cell.id), knn = knn.value, decay = decay.value, )
  phate.data <- data.frame(phate.model$embedding) %>%
    dplyr::mutate(labels = input.data$labels,
                  cell.id = input.data$cell.id)
  
  algorithm_filename <<- paste0(output_filename, "phate.csv")
  write.table(phate.data, algorithm_filename, sep = ",", row.names = FALSE, col.names = TRUE)
  saveRDS(phate.model, paste0(output_filepath, "RDS-files/", "processed_", subset_number, "-cells_", features_summary, "_", balanced_data, "_", "phate.RDS"))
  
  print("phate is complete")
  return(phate.data)
}

## INCOMPLETE
## diffusion map 
# Uses the pair-wise distance matrix for a data set to compute the diffusion map coefficients.
# Computes the Markov transition probability matrix, and its eigenvalues and left & right eigenvectors.
# dmap_function <- function(input.data)  {
# 
#   dists.data <- dist(dplyr::select(input.data, -labels)) ## compute pairwise distance matrix
# 
#   dmap.model <- diffuse(dists.data)
#   dmap.data <- data.frame(dmap.model$layout) %>%
#     dplyr::mutate(labels = input.data[ , "labels"])
# 
#   write.table(dmap.data, paste0(output_filename, "umap.csv"),sep = ",", row.names = FALSE, col.names = TRUE)
# saveRDS(dmap.model, paste0(output_filename, "dmap.RDS"))
# 
#   return(dmap.data)
# }
#####################################################################################################################



##########################################################
## DATA PRE-PROCESSING FUNCTIONS FOR CLUSTERING METHODS ##
##########################################################



## generates the subset of data used for the clustering
## dat = dataset with the classification feature of interest column named as "label" 
## subset_number : # of cells to subset
## can do balanced or unbalanced data...BUT if unbalanced, it won't preserve all labels so you can end up with no observations sampled for labels with few total cells 
## features : the channels  that are included in the run
## features_summary : a summary of the features (eg, "scatterbodies" or "CD-markers" )
generate_subset <- function(data, subset_number, label.levels = NULL, balanced_data = "unbalanced", unbalanced.min.count = 50, features, features_summary, output_filepath = output_filepath) { 
  
  ## generates NEW subset of data
  
  if (balanced_data == "unbalanced") { 
    print("subsetting unbalanced data")
    dat.clean <- data
    dat.clean$labels <- factor(dat.clean$labels, levels = c(label.levels)) ## genereates factor variable for balanced subsetting
    subset_fraction <- round(subset_number/nrow(dat.clean), digits = 30)
    subset_rows <<- caret::createDataPartition(dat.clean[ , "labels"], p = subset_fraction, list =FALSE)
    
    dat.clean <- dat.clean[ dat.clean$cell.id %in% subset_rows, ] 
    # write.table(rownames(dat.clean), paste0(output_filepath, Sys.Date(), "_analyzed-cells_", subset_number, "-cells.csv"), sep=",",row.names = FALSE, col.names = TRUE)
    
    ## data characteristics and output
    output_filename <<- paste0(output_filepath, "processed_", subset_number, "-cells_", features_summary,"_", balanced_data, "_")
    data.input <- dat.clean %>%
      dplyr::select("labels", cell.id, all_of(features))
    
  } else if (balanced_data == "balanced") {
    print("subsetting balanced data")
    dat.clean <- data
    dat.clean$labels <- factor(dat.clean$labels, levels = c(label.levels)) ## genereates factor variable for balanced subsetting
    dat.clean <- dat.clean %>% dplyr::group_by(labels) %>% dplyr::sample_n(subset_number)
    subset_rows <<- dat.clean$cell.id
    # write.table(rownames(dat.clean), paste0(output_filepath, Sys.Date(), "_analyzed-cells_", subset_number, "-cells.csv"), sep=",",row.names = FALSE, col.names = TRUE)
    
    ## data characteristics and output
    output_filename <<- paste0(output_filepath, "processed_", subset_number, "-cells_", features_summary,"_", balanced_data, "_")
    
    data.input <- dat.clean %>%
      ungroup() %>%
      dplyr::select("labels", cell.id, all_of(features))
    
  } else if (balanced_data == "unbalanced_minimum") { ## subsets samples while maintaining distribution but extracts a MINIMUM # of cells for each cell first
    print("subsetting unbalanced data with a minimum # of cells")
    cat(unbalanced.min.count)
    
    
    try(if(subset_number < (unbalanced.min.count*length(unique(data$label)) )) stop("total requested subset of observations is smaller than the minimum # of observaations required for each class label"))
    
    remaining.subset_number = subset_number - (unbalanced.min.count*length(unique(data$label)) )
    
    dat.clean <- data
    dat.clean$labels <- factor(dat.clean$labels, levels = c(label.levels)) ## genereates factor variable for balanced subsetting
    
    ## get minimum cells first then anti_join with original data
    dat.clean <- dat.clean %>% dplyr::group_by(labels) %>% dplyr::sample_n(unbalanced.min.count) %>% ## select minimum cells first
      assign("df.min.subset",.,envir = .GlobalEnv)  %>% ## save minimum subsetted table 
      anti_join(dat.clean, ., by = "cell.id" ) 
    
    ## extract remaining cells  from new dataframe not containing the minimum cells
    subset_fraction <- round(remaining.subset_number/nrow(dat.clean), digits = 30)
    subset_rows <<- caret::createDataPartition(dat.clean[ , "labels"], p = subset_fraction, list =FALSE)
    dat.clean <- dat.clean[ dat.clean$cell.id %in% subset_rows, ] %>%
      bind_rows(., df.min.subset)
    # write.table(rownames(dat.clean), paste0(output_filepath, Sys.Date(), "_analyzed-cells_", subset_number, "-cells.csv"), sep=",",row.names = FALSE, col.names = TRUE)
    ## data characteristics and output
    output_filename <<- paste0(output_filepath, "processed_", subset_number, "-cells_", features_summary,"_", balanced_data, "_")
    data.input <- dat.clean %>%
      dplyr::select("labels", cell.id ,all_of(features))
    
  }
  
  return(data.input)
}

## runs various clustering functions
## creates final tables from clustering function outputs
generate_final_datatables <- function(algorithm_input) {
  
  # algorithm_input$cell.id <- as.vector(subset_rows) ## subset_rows is assigned as a gobal variable from the generate_subset function
  algorithm_input <- algorithm_input %>%
    left_join(., dplyr::select(dat_preproc, -labels), by = "cell.id") ## from main script function
  
  write.table(algorithm_input, algorithm_filename, sep=",",row.names = FALSE, col.names = TRUE) ## algorithm_filename is assigned as a global varialbe in each algorithm method
  
}


## runs various clustering functions and calls generate_final_datatables to generate final summary table
run_algorithms <- function(data.input) { 
  
  pca_output <- pca_function(data.input)
  generate_final_datatables(pca_output)
  
  tsne_output <- tsne_function(data.input, pca.prior = TRUE)
  generate_final_datatables(tsne_output)
  
  phate_output <- phate_function(data.input, gamma.value = 0)
  generate_final_datatables(phate_output)
  
  umap_output <- umap_function(data.input)
  generate_final_datatables(umap_output)
  
}






# generate_final_datatables <- function(data.files_subsetted, channels) {
# 
#   lapply(data.files_subsetted, function(cells_filepath) {
#     algorithm_outputs <- list.files(dirname(cells_filepath), full.names = TRUE) %>%
#       .[grepl(gsub(".*_|cells.*","", cells_filepath), .)] %>%
#       .[grepl("processed", .) ] %>%
#       .[grepl(".csv", .) ]
# 
#     # print(algorithm_outputs)
#     cells_subset <- data.table::fread(cells_filepath) %>%
#       rename(cell.id = x)
#     print(cells_filepath)
#     print(dim(cells_subset))
# 
#     lapply(algorithm_outputs, function(algorithm_file) {
# 
#       print(algorithm_file)
#       processed_data <- data.table::fread(algorithm_file) %>%
#       bind_cols(., cells_subset) %>%
#         left_join(., dat.channels, by = "cell.id")
# 
#       write.table(processed_data, paste0("SCMP/data/analysis-ready/final-methods-tables/", basename(algorithm_file)), sep=",",row.names = FALSE, col.names = TRUE)
#     })
#   })
# 
#     processed_data <- data.table::fread(algorithm_file) %>%
#       bind_cols(., )
# 
# 
#   subsetted_cells <- data.table::fread(cells_filepath)
# 
#   algorithm_outputs <- list.files(basename(cells_filepath), full.names = TRUE) %>%
#     .[grepl(gsub(".*_|-cells.*","",.))]
# 
#   lapply(algorithm_outputs, function(algorithm_file) {
#     processed_data <- data.table::fread(algorithm_file) %>%
#       bind_cols(., )
# 
# 
#   })
# }











