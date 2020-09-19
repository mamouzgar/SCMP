#' Author: Meelad Amouzgar
#' Goal:
#' Runs various dimensionality reduction or clustering methods on the inputted data and outputs the data.
#' This script does not need to be edited unless changes to the fundamental methods want to be made.
#' 
#' Input variables are edited in: 01-generate-features.R and this script requires the following inputs:
#' Requires cell-types to be specified by a column named "gate"
#' Requires a filepath to save outputs
#' Requires a # to subset the data on
# library(tidyverse)
# library(factoextra)
# library(Rtsne)
# library(umap)
# library(phateR)
# library(diffusionMap)
set.seed(0)

##################################
## FEATURE-GENERATING FUNCTIONS ##
##################################
## PCA
## data is centered and scaled
pca_function <- function(input.data) { 
  pca.model <- prcomp(dplyr::select(input.data, -gate), center = TRUE, scale = TRUE)
  pca.data <- data.frame(pca.model$x) %>%
    dplyr::mutate(gate = input.data[ , "gate"])
  
  write.table(pca.data, paste0(output_filename, "pca.csv"),sep = ",", row.names = FALSE, col.names = TRUE)
  saveRDS(pca.model, paste0(output_filename, "pca.RDS"))
  
  return(pca.data)
}


## tsne (with PCA)
## note: uses Rtsne's internal PCA, generates 2-dimensions
## perplexity is defaulted to 25 if arguement is not provided
## potential future added functionality: perplexity is defaulted to 5% of data size if arguement is not provided
## other features set to default 
tsne_function <- function(input.data, perplexity.value = NULL)  { 
  if(is.null(perplexity.value)) { 
    # perplexity.value  <- nrow(input.data)*0.05
    perplexity.value  <- 25
  }
  
  tsne.model <- Rtsne(dplyr::select(input.data, -gate), dims=2, 
                     perplexity= perplexity.value , 
                     pca = TRUE, pca_center = TRUE, pca_scale = TRUE, normalize = TRUE, 
                     theta=0.5,  max_iter = 1000, check_duplicates = FALSE)
  tsne.data <- data.frame(tsne.model$Y) %>%
    dplyr::mutate(gate = input.data[ , "gate"])
  
  write.table(tsne.data, paste0(output_filename, "tsne.csv"),sep = ",", row.names = FALSE, col.names = TRUE)
  saveRDS(tsne.model, paste0(output_filename, "tsne.RDS"))
  
  return(tsne.data)
}


## UMAP
umap_function <- function(input.data)  {
  umap.model <- umap(dplyr::select(input.data, -gate))
  umap.data <- data.frame(umap.model$layout) %>%
    dplyr::mutate(gate = input.data[ , "gate"])
  
  write.table(umap.data, paste0(output_filename, "umap.csv"),sep = ",", row.names = FALSE, col.names = TRUE)
  saveRDS(umap.model, paste0(output_filename, "umap.RDS"))
  
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
  
  
  phate.model <- phate(dplyr::select(input.data, -gate), knn = knn.value, decay = decay.value, )
  phate.data <- data.frame(phate.model$embedding) %>%
    dplyr::mutate(gate = input.data[ , "gate"])
  
  write.table(phate.data, paste0(output_filename, "phate.csv"),sep = ",", row.names = FALSE, col.names = TRUE)
  saveRDS(phate.model, paste0(output_filename, "phate.RDS"))
  
  return(phate.data)
}




## INCOMPLETE
## diffusion map 
# Uses the pair-wise distance matrix for a data set to compute the diffusion map coefficients.
# Computes the Markov transition probability matrix, and its eigenvalues and left & right eigenvectors.
# dmap_function <- function(input.data)  {
#   
#   dists.data <- dist(dplyr::select(input.data, -gate),) ## compute pairwise distance matrix
#   
#   dmap.model <- diffuse(dists.data)
#   dmap.data <- data.frame(dmap.model$layout) %>%
#     dplyr::mutate(gate = input.data[ , "gate"])
#   
#   write.table(dmap.data, paste0(output_filename, "umap.csv"),sep = ",", row.names = FALSE, col.names = TRUE)
# saveRDS(dmap.model, paste0(output_filename, "dmap.RDS"))

#   return(dmap.data)
# }


## diffusion map














