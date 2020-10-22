

rm(list = ls())
require(MASS)
require(dplyr)
require(data.table)
setwd("/Users/mamouzgar/phd-projects/")

# source("SCMP/analysis/scripts/00B_SCMP_axes.R") 

##### USER INPUTS #####

# path to your processed csv file
# path <- "/Users/mamouzgar/phd-projects/SCMP/FlowRepository_FR-FCM-Z2DX_files/processed_mm.csv"

# list of channels you want explored for dimensionality reduction
# channels <- c("WGA_106", "beta_actin", "HP1b", "rRNA", "lamin_A_C", "lamin_B",
#               "lysozyme","VAMP_7", "lactoferrin", "MPO", "serpin_B1", "CD45")


##### FUNCTIONS ##### 

hybridSubsetSelection <- function(x, y, two.d=TRUE) {
  # performs hybrid stepwise subset selection on LDA reduced dimensions
  # Inputs:
  #   x - data.table to evaluate, rows are cells, columns are columns to evaluate
  #   y - vector of observation classes
  #   two.d - logical if true creates two new axes, if false only one
  # Outputs:
  #   matrix of coefficients where rows are markers and columns are axes
  
  ### global data structures ###
  keep <- NULL
  channels <- colnames(x)
  n.channels <- length(channels)
  current.score <- 0
  continue <- TRUE
  results <- setNames(data.frame(matrix(nrow=1, ncol=n.channels)), channels)
  results[1,] <- as.list(rep(F, n.channels))
  subtract.log <- data.table(results[0,]) # record of keep values inputted into subtractOne
  results$score <- 0
  results <- data.table(results)
  
  
  ### functions ###
  addOne <- function() {
    # Evaluates the addition of each channel not in keep to keep. Adds best and updates current.score
    temp.results <- results[0,]
    for (channel in channels[!channels %in% keep]) {
      temp.keep <- c(keep, channel)
      temp.score <- getScore(temp.keep)
      temp.results <- rbind(temp.results, as.list(channels %in% temp.keep) %>% append(temp.score))
    }
    current.score <<- max(temp.results$score)
    new.keep <- temp.results[score==current.score, channels, with=F]
    if (nrow(new.keep) > 1) new.keep <- new.keep[sample(.N,1)]
    keep <<- channels[as.logical(new.keep)]
    results <<- unique(rbind(results, temp.results))
  }
  
  subtractOne <- function() {
    # Evaluates the subtraction of each channel from keep. Removes worst if it improves score and updates current.score
    # If a better subset is found, it calls itself.
    # If this keep has been evaluted before, exits
    subtract.log <<- rbind(subtract.log, as.list(channels %in% keep))
    if (anyDuplicated(subtract.log) > 0) {
      subtract.log <<- unique(subtract.log)
      return()
    }
    temp.results <- results[0,]
    for (channel in keep) {
      temp.keep <- keep[!keep %in% channel]
      temp.score <- getScore(temp.keep)
      temp.results <- rbind(temp.results, as.list(channels %in% temp.keep) %>% append(temp.score))
    }
    if (max(temp.results$score) > current.score) {
      current.score <<- max(temp.results$score)
      new.keep <- temp.results[score==current.score, channels, with=F]
      if (nrow(new.keep) > 1) new.keep <- new.keep[sample(.N,1)]
      keep <<- channels[as.logical(new.keep)]
      results <<- unique(rbind(results, temp.results))
      subtractOne()
    } else results <<- unique(rbind(results, temp.results))
  }
  
  getScore <- function(cols) {
    # performs LDA using columns provided and returns lowest euclidean distance between pop means
    lda.out <- lda(y~., data=x[, cols, with=F])
    if (two.d) return(min(dist(lda.out$means %*% lda.out$scaling[,1:2])))
    return(min(dist(lda.out$means %*% lda.out$scaling[,1])))
  }
  
  initializeKeep <- function() {
    # chooses the best scoring pair of markers to initialize keep
    temp.results <- results[0,]
    for (channel.1 in channels) {
      for (channel.2 in channels) {
        if (channel.1==channel.2) next
        temp.keep <- c(channel.1, channel.2)
        temp.score <- getScore(temp.keep)
        temp.results <- rbind(temp.results, as.list(channels %in% temp.keep) %>% append(temp.score))
      }
    }
    current.score <<- max(temp.results$score)
    new.keep <- temp.results[score==current.score, channels, with=F]
    if (nrow(new.keep) > 1) new.keep <- new.keep[sample(.N,1)]
    keep <<- channels[as.logical(new.keep)]
    results <<- unique(rbind(results, temp.results))
  }
  
  getElbow <- function(res) {
    # takes results and returns the elbow point
    res[, no.markers:=apply(res[, channels, with=F], 1, sum)]
    res.lite <- res[, max(score), by=no.markers][!1]
    slope <- (res.lite$V1[nrow(res.lite)] - res.lite$V1[1]) / (res.lite$no.markers[nrow(res.lite)] - res.lite$no.markers[1])
    intercept <- res.lite$V1[1] - slope * res.lite$no.markers[1]
    perp.slope <- -1 / slope
    perp.int <- res.lite$V1 - (perp.slope * res.lite$no.markers)
    xcross <- (intercept - perp.int) / (perp.slope - slope)
    ycross <- slope * xcross + intercept
    dists <- sqrt((res.lite$no.markers - xcross)^2 + (res.lite$V1 - ycross)^2) %>% round(1)
    elbowi <- max(which(dists==max(dists))) # if dists are tie, take the largest number of channels
    return(elbowi+1)
  }
  
  ### main ###
  initializeKeep()
  while(continue) {
    print(paste("Number of markers:", length(keep)))
    addOne()
    print(paste("Number of markers:", length(keep)))
    if (length(keep) > 3) subtractOne()
    if (length(keep)==length(channels)) continue <- FALSE
  }
  elbow <- getElbow(results)
  markers <- results[no.markers==elbow] %>%
    .[score==max(.[, score]), channels, with=F] %>%
    unlist() %>%
    names(.)[.]
  lda.out <- lda(y~., data=x[, markers, with=F])
  assign("lda.model.output", lda.out, envir = .GlobalEnv)
  if(two.d) return(lda.out$scaling[, 1:2])
  return(lda.out$scaling[, 1, drop=F])
}


makeAxes <- function(dt=dat, co=coefficients, axis.name="ld") {
  # makes new axes based on coefficients
  # Inputs:
  #   dt - data.table of data
  #   co - matrix of coefficients
  #   axis.name - character vector of new axis name (e.g. "ld" results in "ld1" and "ld2")
  # Outputs:
  #   dt - data.table
  x <- as.matrix(dt[, rownames(co), with=F])
  for (i in 1:ncol(co)) {
    dt[, eval(paste0(axis.name, i)):=x %*% co[, i]]
  }
  return(dt)
}




# subsetted.files <- list.files("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/fh-metabolism/sampled-cells/", full.names = TRUE) %>%
#   .[grepl("3ea6", .) ]
# # subsetted.files <- list.files("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/fh-metabolism/sampled-cells/", full.names = TRUE)
# # methodready.df = data.table::fread(subsetted.files[1])
# output_path <- "/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/fh-metabolism/"
# filepath=subsetted.files[3]



metadata <- data.table::fread( "fh-metabolism/data/analysis-ready/metadata-invitro-samples.csv", stringsAsFactors = FALSE) %>%
  mutate(sample.id = paste(cell.type,cell.subtype,day , sep  = "_"),
         cell.type.subtype = paste(cell.type, cell.subtype, sep = "_"),
         day = factor(day)) 

## METHOD 
##' 1) apply all of metabolism dataset to all timepoints, obtain the optimal features, plot
##' SKIP - does not make sense for LDA utility and it's extra work. not a battle to choose 2) extract those features,  apply LDA to the 0 and 3rd time points, and predict on the remaining timepoints, plot
##' SKIP - this produces 1 linear discriminant (since there are only 2 classes) see how it performs
##' 3) Apply LDA to timepoints 0, 2, and 4 (this generates 2 linear discriminants) and see how it performs on timepoints 1,3,and 5
##' 4) Apply with cross-validation to the dataset (using the same extracted features as above) and estimate performance. (maybe, this is overkill)
# apply(subsetted.files, function(filepath) {

  print(filepath)
  # lda.ready_all <-  data.table::fread(filepath)
  lda.ready_all <-  data.table::fread("fh-metabolism/data/analysis-ready/invitro-tcells-3ea6.csv")
  
  comparison.lists<- list(c("naive", "CD4"),
                          c("naive", "CD8"),
                          c("memory","CD4"),
                          c("memory","CD8"))

  ## start loop
  lapply(comparison.lists, function(cell.types.list){
    
    
    dat <- lda.ready_all %>%
      dplyr::filter(grepl(cell.types.list[1], gate.source) & grepl(cell.types.list[2], gate.source))  %>%
      mutate(cell.id = 1:nrow(.),
             labels = gate.source) %>%
      dplyr::select(-Time,-Event_length, -DNA,-dead,-beadDist,-bc_separation_dist,-mahalanobis_dist,-gate.source,
                    -BC102, -BC104, -BC105, -BC106, -BC108, -BC110, -barium, -IdU, H3_p)
    
    minimum.balanced.count <- min(table(dat$labels))
    
    dat <- dat %>%
      group_by(labels) %>%
      # sample_n(minimum.balanced.count) %>%
      sample_n(minimum.balanced.count) %>%
      as.data.table(.)
    
    rownames(dat) <- dat$cell.id
    dat$cell.id <- NULL
    # dat$labels <- NULL
    # dat$BC110 <- NULL
    # dat$CD4 <-NULL
    # dat$CD8 <- NULL
    # dat <- dat[, -c(2,3)]
    channels <- colnames(dat)[!colnames(dat) %in% c("labels")]
    
    
    #################################
    ## APPLY LDA TO ENTIRE DATASET ##
    #################################
    
    # training data with appropriate channels
    train.x <- dat[ , channels, with=F]
    # training class labels - must be 3+ unique classes
    train.y <- dat[, labels]
    
    coefficients <- hybridSubsetSelection(x=train.x, y=train.y, two.d=TRUE)
    dat.output <- makeAxes(dt=dat, co=coefficients, axis.name="ld")
    
    full.dataset <- dat.output %>% mutate(LDA.analysis = "full-balanced-dataset", days= "all", cell.id = rownames(dat))
    final.analysis.features <- rownames(coefficients)
    
    #######################################
    ## TRAIN LDA TO TIMEPOINTS 0,2,and 4 ##
    ## TEST LDA ON TIMEPOINTS  1,3,and 5 ##
    #######################################
    
    
    
    # training data with appropriate channels
    train.x <- dat[grepl("day0|day2|day4", dat$labels), final.analysis.features, with=F]
    # training class labels - must be 3+ unique classes
    train.y <- dat[grepl("day0|day2|day4", dat$labels), labels]
    
    test.x <- dat[grepl("day1|day3|day5", dat$labels), final.analysis.features, with=F]
    test.y <- dat[grepl("day1|day3|day5", dat$labels), labels]
    # test.x$labels = test.y
    
    
    ## run training
    lda.model <- MASS::lda(x = train.x, grouping =train.y  )
    lda.data <- data.frame(lda.model$means %*% lda.model$scaling[ , 1:2])
    lda.data$labels = rownames(lda.data)
    
    co = lda.model$scaling[, 1:2, drop=F]
    dt <- data.table::as.data.table(train.x)
    x <- as.matrix(train.x)
    for (i in 1:ncol(co)) {
      dt[, eval(paste0("ld", i)):=x %*% co[, i]]
    }
    dt$labels = train.y
    dt$cell.id <- rownames(dat)[grepl("day0|day2|day4", dat$labels)]
    
    
    test.predict <- predict(lda.model, test.x)
    predict.axes <- data.frame(test.predict$x) %>%
      bind_cols(.,test.x, data.frame(labels = test.y)) %>%
      mutate(days = "days-1-3-5-test",
             cell.id = rownames(dat)[grepl("day1|day3|day5", dat$labels)]) %>%
      rename(ld1=LD1, ld2 = LD2) %>% 
      bind_rows(.,dt) %>%
      mutate(days = case_when(is.na(days) ~ "days-0-2-4-train",
                              TRUE ~ days),
             LDA.analysis = "predicting-timepoint-balanced-subset") 
    
    full.dataset <- bind_rows(full.dataset, predict.axes) %>%
      left_join(.,metadata, by = c("labels"= "gate.source"))
    
    
    
    write.table(full.dataset, paste0("SCMP/data/summary-tables/",cell.types.list[2],cell.types.list[1], "fh-metabolism-continuous-axes-predictions.csv",sep=",",row.names = FALSE,col.names = TRUE ))
    
  })
  
  # 
  # coefficients <- hybridSubsetSelection(x=train.x, y=train.y, two.d=tryNew())
  # dat.output <- makeAxes()
  # 
  # test.predict <- predict(lda.model.output, test.x)
  # predict.axes <- data.frame(test.predict$x) %>%
  #   bind_cols(.,test.x) %>%
  #   rename(ld1=LD1, ld2 = LD2) %>%
  #   mutate(data.partition = "test") %>%
  #   bind_rows(., train.x)
  # 
  # 
  # # training data with appropriate channels
  # train.x <- dat[ , channels, with=F]
  # # training class labels - must be 3+ unique classes
  # train.y <- dat[ , labels]
  # 
  # 
  # # training data with appropriate channels
  # train.x <- dat[grepl("_0|_3", dat$labels), channels, with=F]
  # # training class labels - must be 3+ unique classes
  # train.y <- dat[grepl("_0|_3", dat$labels), labels]
  # 
  # test.x <- dat[!grepl("_0|_3", dat$labels), channels, with=F]
  # test.y <- dat[!grepl("_0|_3", dat$labels), labels]
  # test.x$labels = test.y
  # 
  # coefficients <- hybridSubsetSelection(x=train.x, y=train.y, two.d=FALSE)
  # dat.output <- makeAxes()
  # 
  # test.predict <- predict(lda.model.output, test.x)
  # predict.axes <- data.frame(test.predict$x) %>%
  #   bind_cols(.,test.x) %>%
  #   rename(ld1=LD1) %>%
  #   mutate(data.partition = "test") %>%
  #   bind_rows(., train.x)
  # 
  # 
  # print(filepath)
  # lda.ready_all <-  data.table::fread(filepath)
  # 
  # lda.ready <- lda.ready_all %>% dplyr::filter(grepl("_0|_2|_5", labels), grepl("naive", labels), grepl("CD8", labels))
  # # rownames(lda.ready) <- lda.ready$labels
  # metadata <- lda.ready %>% dplyr::select(cell.id,labels)
  # lda.ready <-lda.ready %>% dplyr::select(-cell.id,-labels)
  # 
  # lda.model <- MASS::lda(x = lda.ready,grouping =metadata$labels  )
  # lda.data <- data.frame(lda.model$means %*% lda.model$scaling[ , 1:2])
  # # lda.data <- data.frame(lda.model$means %*% lda.model$scaling[ , 1])
  # lda.data$labels = rownames(lda.data)
  # 
  # co = lda.model$scaling[, 1:2, drop=F]
  # dt <- data.table::as.data.table(lda.ready)
  # x <- as.matrix(lda.ready)
  # for (i in 1:ncol(co)) {
  #   dt[, eval(paste0("ld", i)):=x %*% co[, i]]
  # }
  # dt$labels = metadata$labels
  # 
  # write.table(dt, paste0(output_path, gsub(".csv", "_lda.csv", basename(filepath) )), sep = ",",row.names = FALSE, col.names = TRUE)
  # 
# } )
# 
# 
# ggpubr::ggscatter(total.output, x = "ld1", y= "ld2", color = "labels", subtitle = "trained on timepoints 0,1,and 5")
# 
# ggpubr::ggscatter(dat.output, x = "ld1", y= "ld2", color = "labels", subtitle = "trained on timepoints 0,1,and 5")
# 
# ggpubr::ggscatter(dat.output, x = "ld1", y= "ld2", color = "labels", subtitle = "trained on timepoints 0,1,and 5")




