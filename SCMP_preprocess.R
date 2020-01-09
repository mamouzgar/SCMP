####################################################################################################################
#
# Script: SCMP_preprocess.R
# Project: Multiplexed Single Cell Morphometry for Hematopathology Diagnostics
# Author: David Glass
# Date: 1-8-20
#
# Purpose: Preprocess fcs files from Figure 4 so they may be used as an example to generate MM axes in SCMP_LDA.R
#
# FCS files:
#   http://flowrepository.org/id/FR-FCM-Z2DX
#
# Pseudocode:
#   Read in all fcs files
#   Combine data into a single data.table with metadata
#   Asinh transform and scale
#   Write csv
#
# Instructions:
# Install all required packages (see LIBRARIES and/or Versions below)
# Download the dataset linked above
# Put the dataset into a unique directory
# In the USER INPUTS section below, assign path variable to the path of the fcs directory
# In the MAIN section at the bottom of the script, if you do not wish to generate a csv of processed data,
#   set csv=FALSE in the processData function (a data.table will be returned)
# By default, the csv is written to the fcs file path, you can change that by assigning mp in the processData function
#   to a different directory
# Run script
# The csv file or data.table can be used in SCMP_axes.R to recreated the MM axes
#
### NOTE: Processing often takes a few minutes to run
#
# Versions:
# R 3.5.1
# RStudio 1.1.463
# flowCore_1.48.1
# dplyr_0.8.3
# data.table_1.12.8
#
#######################################################################################################################



##### USER INPUTS #####

### Path to folder containing fcs files
path <- "~/example/"



###### LIBRARIES ######

require(flowCore)
require(dplyr)
require(data.table)



##### METADATA #####

files <- paste0(c("all-cells-9878",
              "all-cells-10874",
              "AML-4Cn5",
              "AML-APL-1An7", 
              "AMML-1An6",
              "B-ALL-4Cr9",
              "MDS-EB2-1Ar1",
              "MDS-RS-1Ar0",
              "MPAL-3Cn3",
              "PCM-2Ar",
              "SM-CMML-1Ar0",
              "T-ALL-3An1",
              "TCL-1An7",
              "training_MM_all-cells-10874_B_cells",
              "training_MM_all-cells-10874_blasts",
              "training_MM_all-cells-10874_erythroids",
              "training_MM_all-cells-10874_hematogones",
              "training_MM_all-cells-10874_monocytes",
              "training_MM_all-cells-10874_neutrophils",
              "training_MM_all-cells-10874_T_and_NK_cells"), ".fcs")
gates <- c(rep("ungated", 13), c("lymphocyte", "blast", "erythroid", "lymphocyte", "monocyte", "neutrophil", "lymphocyte"))
meta.dat <- data.table(file=files, gate=gates)
factors <- c("file", "gate")
rm(files, gates)



##### FUNCTIONS #####

readFile <- function(p=path) {
  # takes in a path with fcs files and returns a list of data.tables of the expression matrix
  # Inputs:
  #   p - character vector with directory storing fcs files
  # Outputs:
  #   frames - a list of data.tables
  files <- list.files(path=p, pattern=".fcs", full.names=FALSE, recursive=FALSE)
  frames <- setNames(vector("list", length(files)), files)
  for (i in 1:length(files)) {
    fcs <- read.FCS(paste0(path, files[i]), transformtablation = FALSE, emptyValue = FALSE)
    frames[[i]] <- data.table(exprs(fcs))
  }
  return(frames)
}


combineFiles <- function(frames, fa=factors, meta=meta.dat) {
  # combines a list of data.tables into a single data.table with factor columns added
  # Inputs:
  #   frames - a list of data.tables
  #   fa - vector of factor column names
  #   meta - data.table with metadata for the dataset being used
  # Outputs:
  #   dt - a single data.table with all data for that dataset
  all.cols <- unique(unlist(lapply(frames, colnames))) %>% c(fa)
  dt <- data.table(matrix(nrow=0, ncol=length(all.cols))) %>% setnames(all.cols)
  for (f in meta$file) {
    dt <- frames[[f]] %>%
      cbind(meta[file==f]) %>%
      rbind(dt, fill=T)
  }
  return(dt[, setdiff(all.cols, "Scatter"), with=F])
}


asinTransform <- function(dt, exceptions=c(factors, "ld1", "ld2")) {
  # asinh transforms data
  # Inputs:
  #   dt - data.table
  #   exceptions - character vector of channels not to transform
  # Outputs:
  #   dat - data.table
  to.transform <- colnames(dt)[!colnames(dt) %in% exceptions]
  dt[, (to.transform) := asinh(dt[, to.transform, with=F]/5)]
  return(dt)
}


scaleData <- function(dt, exceptions=c(factors, "ld1", "ld2"), ref="all-cells-9878.fcs") {
  # scales data to 99.5th percentile
  # Inputs:
  #   dt - data.table
  #   exceptions - character vector of channels not to scale
  #   ref - character vector with patient to use as reference for scaling
  # Outputs:
  #   dt - data.table
  findRefs <- function(x) quantile(x, probs=c(0.995), na.rm=T)
  channels <- colnames(dt)[!colnames(dt) %in% exceptions]
  mins <- apply(dt[, (channels), with=F], 2, min)
  refs <- apply(dt[file==ref, (channels), with=F], 2, findRefs) - mins
  for (channel in channels) {
    dt[, eval(channel) := (dt[, eval(channel), with=F] - mins[channel]) / refs[channel]]
  }
  return(dt)
}


processData <- function(csv, mp=paste0(path, "processed_mm.csv")) {
  # Processes data for a dataset
  # Inputs:
  #   csv - logical, true if csv should be written
  #   mp - character vector of location to write csv file
  # Outputs:
  #   dt - data.table
  dt <- readFile() %>%
    combineFiles() %>%
    asinTransform() %>%
    scaleData()
  if (csv) fwrite(dt, mp)
  return(dt)
}



##### MAIN #####

dat <- processData(csv=TRUE)