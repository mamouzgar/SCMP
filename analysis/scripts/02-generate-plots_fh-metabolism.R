#' this script generates plots from the 01-generate-features.R script output
#' 
rm(list = ls())
library(ggpubr)
library(RColorBrewer)
set.seed(0)
setwd("/Users/mamouzgar/phd-projects")


###############
## functions ##
###############
# sp_plot <- function(input.data, cluster.method = NULL, markers = NULL) { 
  # ggscatter(data.subset, x= "axis1", y = "axis2",
  #           color = "gate", size = 1 ) +
  # scale_color_brewer(palette = "Set1")
# } 

density_2d_plot <- function(input.data, cluster.method = NULL, cell_counts = NULL) { 
  d2d.plot <- ggplot(input.data, aes(x=axis1,y = axis2)) +
    geom_density_2d(aes(color = gate)) + 
    scale_color_brewer(palette = "Set1") +
    theme_pubr() + 
    theme(legend.position = "right") +
    ggtitle(paste(cluster.method, cell_counts,sep = "-")) + 
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),

          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  return(d2d.plot)
  }

sp_plot  <- function(input.data, cluster.method = NULL, cell_counts = NULL) { 
  
  sp.plot <- ggscatter(input.data, x= "axis1", y = "axis2", 
                          color = "gate", size = 0.1, legend = "top",
                          title = paste(cluster.method, cell_counts,sep = "-"),
                          shape = 1) +
    scale_color_brewer(palette = "Set1") +
    theme_pubr() +
    theme(legend.position = "none") +
    theme(axis.title.x=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),

          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank())
  
  return(sp.plot)
  
}

###############
## variables ##
###############
result.files <- list.files("SCMP/data/analysis-ready", full.names = TRUE)
result.files <- result.files[grepl("pca.csv|tsne.csv|umap.csv|phate.csv|lda.csv",result.files)]
# result.files <- result.files[!grepl("176664",result.files)]
result.files 
result.files <- data.frame(result.files, stringsAsFactors = FALSE) %>%
  dplyr::mutate(cell.counts = gsub(".*processed_|-cells.*", "", result.files),
                cell.counts = as.numeric(cell.counts),
                method = gsub(".*_|.csv", "", result.files),
                method = factor(method, levels = c("pca","tsne","umap","phate", "lda"))) %>%
  arrange(cell.counts, method)
# pca.files <- result.files[grepl("pca.csv",result.files)]
# tsne.files <- result.files[grepl("tsne.csv",result.files)]
# umap.files <- result.files[grepl("umap.csv",result.files)]
# phate.files <- result.files[grepl("phate.csv",result.files)]

plot.output.path <- "SCMP/analysis/plots/"

##########
## MAIN ##
##########
plots.list <- lapply(result.files$result.files, function(filepath) { 
  print(filepath)
  
  if ( grepl("lda", filepath)) { 
    df <- data.table::fread(filepath) %>%
      mutate(gate = factor(gate, levels = c("monocyte","lymphocyte","erythroid","blast", "neutrophil")))
  } else { 
    df <- data.table::fread(filepath) %>%
      mutate(gate = factor(labels, levels = c("monocyte","lymphocyte","erythroid","blast", "neutrophil")))
  }
  
  if ( grepl("lda", filepath)) { 
    df <- df[, c("ld1", "ld2", "gate")]
  } else { 
    df <- df[, c(1,2, ncol(df))]
  }
  
  colnames(df) <- c("axis1","axis2","gate")
  
  if (grepl("phate", filepath) | any(table(df$gate) <= 1)) { 
    sp_plot(df, cluster.method = gsub(".*_|.csv", "", filepath), cell_counts = gsub(".*processed_|-cells.*", "", filepath))
    
  } else { 
    density_2d_plot(df, cluster.method = gsub(".*_|.csv", "", filepath), cell_counts = gsub(".*processed_|-cells.*", "", filepath))
    }
  })

###############
## aggregate ##
###############
library(gridExtra)
# library(grid)
pdf(paste0(plot.output.path, "all-data-plotting-figures.pdf"), height = 25, width = 22)
do.call("grid.arrange", c(plots.list, ncol = 5))
dev.off()
  
result.files <- list.files("SCMP/data/analysis-ready", full.names = TRUE)
result.files <- result.files[grepl("pca.csv|tsne.csv|umap.csv|phate.csv|lda.csv",result.files)]
# result.files <- result.files[!grepl("176664",result.files)]
result.files 
result.files <- data.frame(result.files, stringsAsFactors = FALSE) %>%
  dplyr::mutate(cell.counts = gsub(".*processed_|-cells.*", "", result.files),
                cell.counts = as.numeric(cell.counts),
                method = gsub(".*_|.csv", "", result.files),
                method = factor(method, levels = c("pca","tsne","umap","phate", "lda"))) %>%
  arrange(cell.counts, method)
plot.output.path <- "SCMP/analysis/plots/"




###########################################
## generate features for metabolism data ##
###########################################

cat("reading global variables")
metadata <- data.table::fread( "fh-metabolism/data/analysis-ready/metadata-invitro-samples.csv", stringsAsFactors = FALSE) %>%
  mutate(sample.id = paste(cell.type,cell.subtype,day , sep  = "_"))
dat <- data.table::fread("fh-metabolism/data/analysis-ready/all-invitro-samples-Tcells-LDA-ready.csv", sep = ",", stringsAsFactors = FALSE)
dat

head(dat)
all_na <- function(x) any(!is.na(x))
dat.clean <- dat %>%
  filter(grepl("3ea6", filename)) %>%
  select_if(all_na) %>%
  filter(H3_p < 0.1) %>%
  dplyr::select(-mahalanobis_dist, -Time, -Event_length, -H3_p, -IdU,-barium, -dead, -BC102,-DNA,-BC104,-BC105,-BC106,-BC108, -beadDist,-bc_separation_dist) %>%
  left_join(., dplyr::select(metadata,gate.source, sample.id), by= c("filename" = "gate.source")) %>%
  dplyr::select(-filename) %>%
  mutate(cell.id = paste0("cell", 1:nrow(.))) %>%
  ungroup()
dat.clean$labels = dat.clean$sample.id
dat.clean$sample.id <- NULL

dat <- dat.clean
output_filepath <- "/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/fh-metabolism/" ## output directory 
features <- colnames(dat.clean)[-c(46,47)]

features_summary <- "metabolism"
label.levels <- c("CD4_naive_0", "CD4_naive_1", "CD4_naive_2", "CD4_naive_3", "CD4_naive_4","CD4_naive_5", 
                  "CD4_memory_0", "CD4_memory_1", "CD4_memory_2", "CD4_memory_3", "CD4_memory_4","CD4_memory_5",
                  "CD8_naive_0", "CD8_naive_1", "CD8_naive_2", "CD8_naive_3", "CD8_naive_4","CD8_naive_5", 
                  "CD8_memory_0", "CD8_memory_1", "CD8_memory_2", "CD8_memory_3", "CD8_memory_4","CD8_memory_5")
balanced_data <- "balanced" ## options are: "balanced" or "unbalanced" or "unbalanced_minimum" you can specify the minimum count
subset_numbers <- c(1000, 4000,  8336 ) ## for unbalanced or unbalanced_minimum data


# balanced_data <- "balanced" ## options are: "balanced" or "unbalanced"
# subset_numbers <- c(100, 200, 300, 400, 548) ## for balanced data, limiter is blast cells

use_previously_subsetted_data <- FALSE ## TRUE : use subsetted data (a training set) previously produced, FALSE : generates new subset of data
# data.files_subsetted <- list.files(output_filepath, full.names = TRUE) 

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



density_2d_plot(dat.clean, cluster.method = "LDA", cell_counts = 176664)


pdf(paste0(plot.output.path, "LDA-plot_all-cells.pdf"), height = 2, width = 3)
density_2d_plot(dat.clean, cluster.method = "LDA", cell_counts = 176664)+theme(legend.position="right")
dev.off()



#########
## LDA ##
#########
## metabolism data


## extract cell.id for LDA 
dat <- data.table::fread("fh-metabolism/data/analysis-ready/all-invitro-samples-Tcells-LDA-ready.csv", sep = ",", stringsAsFactors = FALSE)
dat
head(dat)
all_na <- function(x) any(!is.na(x))
dat.clean <- dat %>%
  filter(grepl("3ea6", filename)) %>%
  select_if(all_na) %>%
  filter(H3_p < 0.1) %>%
  dplyr::select(-mahalanobis_dist, -Time, -Event_length, -H3_p, -IdU,-barium, -dead, -BC102,-DNA,-BC104,-BC105,-BC106,-BC108, -beadDist,-bc_separation_dist) %>%
  left_join(., dplyr::select(metadata,gate.source, sample.id), by= c("filename" = "gate.source")) %>%
  dplyr::select(-filename) %>%
  mutate(cell.id = paste0("cell", 1:nrow(.))) %>%
  ungroup()
dat.clean$labels = dat.clean$sample.id
dat.clean$sample.id <- NULL

## get sampled-cells for lda
output_path <- "/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/fh-metabolism/sampled-cells/"
subsetted.files <- list.files("SCMP/data/analysis-ready/fh-metabolism/", full.names = TRUE) %>%
  .[grepl("umap",.)]

lapply(subsetted.files, function(filepath) {
  print(filepath)
  methodready.df = data.table::fread(filepath)

  lda.ready.df <- dat.clean %>%
    filter(cell.id %in% methodready.df$cell.id)
  
  file.name <- gsub("balanced_.*","balanced",basename(filepath))
  write.table(lda.ready.df, paste0(output_path, file.name, "_3ea6.csv"  ),sep=",",row.names = FALSE, col.names = TRUE)
  
  
})

subsetted.files <- list.files("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/fh-metabolism/sampled-cells/", full.names = TRUE)
# methodready.df = data.table::fread(subsetted.files[1])
output_path <- "/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/fh-metabolism/"

lapply(subsetted.files, function(filepath) {
  print(filepath)
  lda.ready <-  data.table::fread(filepath)
  # rownames(lda.ready) <- lda.ready$labels
  metadata <- lda.ready %>% dplyr::select(cell.id,labels)
  lda.ready <-lda.ready %>% dplyr::select(-cell.id,-labels)
  
  lda.model <- MASS::lda(x = lda.ready,grouping =metadata$labels  )
  lda.data <- data.frame(lda.model$means %*% lda.model$scaling[ , 1:2])
  lda.data$labels = rownames(lda.data)
  
  co = lda.model$scaling[, 1:2, drop=F]
  dt <- data.table::as.data.table(lda.ready)
  x <- as.matrix(lda.ready)
  for (i in 1:ncol(co)) {
    dt[, eval(paste0("ld", i)):=x %*% co[, i]]
  }
  dt$labels = metadata$labels
  
  write.table(dt, paste0(output_path, gsub(".csv", "_lda.csv", basename(filepath) )), sep = ",",row.names = FALSE, col.names = TRUE)
  
  } )
  




##########
## PLOT ##
##########
subsetted.files <- list.files("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/fh-metabolism/", full.names = TRUE) %>%
  .[grepl("csv",.)]
metadata <- data.table::fread( "fh-metabolism/data/analysis-ready/metadata-invitro-samples.csv", stringsAsFactors = FALSE) %>%
  mutate(sample.id = paste(cell.type,cell.subtype,day , sep  = "_"),
         cell.type.subtype = paste(cell.type, cell.subtype, sep = "_"),
         day = factor(day)) 
## lda only
lapply(subsetted.files, function(file.path) { 
  print(file.path)
  lda.ready <-  data.table::fread(file.path) %>%
    left_join(., metadata, by = c("labels"= "sample.id")) %>%
    filter(day != 4)
  
  filename = gsub("processed_|-cells_metabolism|.csv","",basename(file.path))
  outputpath = paste0("SCMP/analysis/plots/fh-metabolism/T-cells_invitro", filename, ".pdf")
  
  
  if (grepl("lda", filename)) { 
    lda.plot <- ggpubr::ggscatter(lda.ready, x = "ld1",y="ld2", color = "day",facet.by = "cell.type.subtype", size = 0.1,
                                  subtitle = basename(outputpath)) +
      viridis::scale_color_viridis(discrete=TRUE) +
      guides(color = guide_legend(override.aes = list(size=5))) 
    pdf(outputpath, width = 6, height = 6)
    print(lda.plot)
    dev.off()
    
  lda.plot <- ggpubr::ggscatter(lda.ready, x = "ld1",y="ld2", color = "day", size = 0.1,
                                subtitle = basename(outputpath)) +
    viridis::scale_color_viridis(discrete=TRUE) +
    guides(color = guide_legend(override.aes = list(size=5))) 
  pdf(paste0(outputpath, "no-facet.pdf"), width = 6, height = 6)
   print(lda.plot)
  dev.off()
  
  } else if (grepl("pca", filename)) { 
    lda.plot <- ggpubr::ggscatter(lda.ready, x = "PC1",y="PC2", color = "day",facet.by = "cell.type.subtype", size = 0.1,
                                  subtitle = basename(outputpath)) +
      viridis::scale_color_viridis(discrete=TRUE) +
      guides(color = guide_legend(override.aes = list(size=5))) 
    pdf(outputpath, width = 6, height = 6)
    print(lda.plot)
    dev.off()
    
    lda.plot <- ggpubr::ggscatter(lda.ready, x = "PC1",y="PC2", color = "day", size = 0.1,
                                  subtitle = basename(outputpath)) +
      viridis::scale_color_viridis(discrete=TRUE) +
      guides(color = guide_legend(override.aes = list(size=5))) 
    pdf(paste0(outputpath, "no-facet.pdf"), width = 6, height = 6)
    print(lda.plot)
    dev.off()
  } else if (grepl("umap", filename)) { 
    lda.plot <- ggpubr::ggscatter(lda.ready, x = "X1",y="X2", color = "day",facet.by = "cell.type.subtype", size = 0.1,
                                  subtitle = basename(outputpath)) +
      viridis::scale_color_viridis(discrete=TRUE) +
      guides(color = guide_legend(override.aes = list(size=5))) 
    pdf(outputpath, width = 6, height = 6)
    print(lda.plot)
    dev.off()
    
    lda.plot <- ggpubr::ggscatter(lda.ready, x = "X1",y="X2", color = "day", size = 0.1,
                                  subtitle = basename(outputpath)) +
      viridis::scale_color_viridis(discrete=TRUE) +
      guides(color = guide_legend(override.aes = list(size=5))) 
    
    pdf(paste0(outputpath, "no-facet.pdf"), width = 6, height = 6)
    print(lda.plot)
    dev.off()
  } else if (grepl("phate", filename)) { 
    lda.plot <- ggpubr::ggscatter(lda.ready, x = "PHATE1",y="PHATE2", color = "day",facet.by = "cell.type.subtype", size = 0.1,
                                  subtitle = basename(outputpath)) +
      viridis::scale_color_viridis(discrete=TRUE) +
      guides(color = guide_legend(override.aes = list(size=5))) 
    pdf(outputpath, width = 6, height = 6)
    print(lda.plot)
    dev.off()
    
    lda.plot <- ggpubr::ggscatter(lda.ready, x = "PHATE1",y="PHATE2", color = "day", size = 0.1,
                                  subtitle = basename(outputpath)) +
      viridis::scale_color_viridis(discrete=TRUE) +
      guides(color = guide_legend(override.aes = list(size=5))) 
    pdf(paste0(outputpath, "no-facet.pdf"), width = 6, height = 6)
    print(lda.plot)
    dev.off()
  }
  
  })

  

# subsetted.files <- list.files("/Users/mamouzgar/phd-projects/SCMP/data/analysis-ready/fh-metabolism/", full.names = TRUE) %>%
#   .[grepl("csv",.)] %>%
#   .[!grepl("lda",.)] 
# metadata <- data.table::fread( "fh-metabolism/data/analysis-ready/metadata-invitro-samples.csv", stringsAsFactors = FALSE) %>%
#   mutate(sample.id = paste(cell.type,cell.subtype,day , sep  = "_"),
#          cell.type.subtype = paste(cell.type, cell.subtype, sep = "_"),
#          day = factor(day)) 
# 
# ## NOT lda 
# lapply(subsetted.files, function(file.path) { 
#   lda.ready <-  data.table::fread(file.path) %>%
#     left_join(., metadata, by = c("labels"= "sample.id")) %>%
#     filter(day != 4)
#   
#   filename = gsub("processed_|-cells_metabolism|.csv","",basename(file.path))
#   outputpath = paste0("SCMP/analysis/plots/fh-metabolism/T-cells_invitro", filename, ".pdf")
#   lda.plot <- ggpubr::ggscatter(lda.ready, x = "ld1",y="ld2", color = "day",facet.by = "cell.type.subtype", size = 0.1,
#                                 subtitle = basename(outputpath)) +
#     viridis::scale_color_viridis(discrete=TRUE) +
#     guides(color = guide_legend(override.aes = list(size=5))) 
#   
#   pdf(outputpath, width = 6, height = 6)
#   print(lda.plot)
#   dev.off()
#   
# })


