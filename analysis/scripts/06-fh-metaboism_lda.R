

rm(list = ls())
require(MASS)
require(dplyr)
require(data.table)






subsetted.files <- list.files("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/fh-metabolism/sampled-cells/", full.names = TRUE) %>%
  .[grepl("c393", .) ]
subsetted.files <- list.files("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/fh-metabolism/sampled-cells/", full.names = TRUE)
# methodready.df = data.table::fread(subsetted.files[1])
output_path <- "/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/fh-metabolism/"

filepath=subsetted.files[1]

## apply LDA only to the first and 5th timepoints, impute the remaining 
apply(subsetted.files, function(filepath) {
  
  
  print(filepath)
  lda.ready_all <-  data.table::fread(filepath)
  
  dat <- lda.ready_all[
    # grepl("_0|_1|_5", lda.ready_all$labels) & grepl("naive", lda.ready_all$labels) &grepl("CD8", lda.ready_all$labels)]
     grepl("naive", lda.ready_all$labels) &grepl("CD8", lda.ready_all$labels)]

  rownames(dat) <- dat$cell.id
  dat$cell.id <- NULL
  # dat$labels <- NULL
  dat$BC110 <- NULL
  # dat$CD4 <-NULL
  # dat$CD8 <- NULL
  # dat <- dat[, -c(2,3)]
  channels <- colnames(dat)[!colnames(dat) %in% c("labels")]
  
  # training data with appropriate channels
  train.x <- dat[labels!="ungated", channels, with=F]
  # training class labels - must be 3+ unique classes
  train.y <- dat[labels!="ungated", labels]
  
  
  coefficients <- hybridSubsetSelection(x=train.x, y=train.y, two.d=TRUE)
  dat.output <- makeAxes()
  
  
  
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
  
  # write.table(dt, paste0(output_path, gsub(".csv", "_lda.csv", basename(filepath) )), sep = ",",row.names = FALSE, col.names = TRUE)
  
} )


lda.model.output

test <- lda.ready_all %>% dplyr::filter(!grepl("_0|_1|_5", labels), grepl("naive", labels), grepl("CD8", labels)) 
test.input <- lda.ready_all[!grepl("_0|_1|_5", lda.ready_all$labels) & grepl("naive", lda.ready_all$labels) &grepl("CD8", lda.ready_all$labels)]
test.input <- lda.ready_all[grepl("naive", lda.ready_all$labels) &grepl("CD8", lda.ready_all$labels)]
# test.input <- test %>% dplyr::select(-cell.id,-labels)
cl <- test.input$labels
test.predict <- predict(lda.model.output, test.input)
predict.axes <- data.frame(test.predict$x) %>%
  bind_cols(.,test) %>%
  rename(ld1=LD1,ld2=LD2)

total.output <- bind_rows(predict.axes, dat.output)

ggpubr::ggscatter(total.output, x = "ld1", y= "ld2", color = "labels", subtitle = "trained on timepoints 0,1,and 5")

ggpubr::ggscatter(dat.output, x = "ld1", y= "ld2", color = "labels", subtitle = "trained on timepoints 0,1,and 5")

ggpubr::ggscatter(dat.output, x = "ld1", y= "ld2", color = "labels", subtitle = "trained on timepoints 0,1,and 5")




