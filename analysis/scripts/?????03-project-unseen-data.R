



## test script to project unseen data onto the original embedding

## pca
# load pca RDS
pca.model <- readRDS("SCMP/data/analysis-ready/processed_10000-cells_LDA-scatterbodies_pca.RDS")

# load train data
train.set <- dat <- data.table::fread("SCMP/data/analysis-ready/processed_10000-cells_LDA-scatterbodies_pca.csv", sep = ",", stringsAsFactors = FALSE)
# trained cells
train.cells <- data.table::fread("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/2020-09-19_analyzed-cells_10000-cells.csv", sep = ",", stringsAsFactors = FALSE)

## generate new, unseen data (test set)
cat("reading global variables")
dat <- data.table::fread("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/processed_LDaxes.csv", sep = ",", stringsAsFactors = FALSE)
dat[ , "cell.id"] <- 1:nrow(dat)

data.test <- dat %>%
  dplyr::filter(!cell.id %in% train.cells$x)

output_filepath <- "/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/test-" ## output directory 
channels <- c("WGA_106", "beta_actin", "HP1b", "rRNA", "lamin_A_C", "lamin_B",
              "lysozyme","VAMP_7", "lactoferrin", "MPO", "serpin_B1", "CD45")
features_summary <- "LDA-scatterbodies"

data.test.input <- generate_subset(dat = data.test, subset_number = 10000, features = channels, features_summary = features_summary, output_filepath=output_filepath) 
gates <- data.test.input$gate
data.test.input$gate <- NULL

## project new data onto PCA-space
# predict(pca.model, data.test.input)
test.output <- scale(data.test.input, pca.model$center, pca.model$scale) %*% pca.model$rotation %>%
  data.frame(.) %>%
  dplyr::mutate(gate = gates)

## plot 
train.set$set <- "train"
test.output$set <- "test"
train.test.merge <- bind_rows(train.set, test.output)

ggscatter(train.test.merge, x = "PC1", y = "PC2", color = "gate", shape = "set", facet.by = "set")

TAB <- table(pred, winevalid$class) # table of preditions vs. original classifications
TAB
# pred         barolo grignolino barbera
#   barolo         29          1       0
#   grignolino      1         30       0
#   barbera         0          1      27

sum(diag(TAB)) / sum(TAB) # overall accuracy
# [1] 0.9662921

