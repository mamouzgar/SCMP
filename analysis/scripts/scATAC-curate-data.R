###########
## scATC ##
###########

scdata <- readRDS("/Users/mamouzgar/data/pd1-epigenetics/data/Satpathy2019NatBiotech/Satpathy_TME_scATAC_MeeladPeaks.rds")


se <- SummarizedExperiment(t(assay(scdata, "counts")))
# Too large, let's just use the most variable peaks

colVars_sparse <- function(sum.exp) { 
  
  lapply(1:nrow(sum.exp), function(row.num) { 
    cat(row.num)
    var(sum.exp[ row.num, ])
    
  })
  
}


variance.features <- colVars_sparse(assay(scdata))


## compute variance for each peak
colVars_spm <- function( spm ) {
  stopifnot( methods::is( spm, "dgCMatrix" ) )
  ans <- sapply( base::seq.int(spm@Dim[2]), function(j) {
    print(j)
    if( spm@p[j+1] == spm@p[j] ) { return(0) } # all entries are 0: var is 0
    mean <- base::sum( spm@x[ (spm@p[j]+1):spm@p[j+1] ] ) / spm@Dim[1]
    sum( ( spm@x[ (spm@p[j]+1):spm@p[j+1] ] - mean )^2 ) +
      mean^2 * ( spm@Dim[1] - ( spm@p[j+1] - spm@p[j] ) ) } ) / ( spm@Dim[1] - 1 )
  names(ans) <- spm@Dimnames[[2]]
  ans
}

variance.features <- colVars_spm(assay(se))
sum.cells <- Matrix::rowSums(assay(se))

# top 25% most variable peaks
keep.peaks <- variance.features >= quantile(variance.features, c(0.85))[1]
keep.peaks.index <- grep("TRUE", keep.peaks)
se.filtered <- se[ , keep.peaks]


se.matrix <- as.matrix(assay(se.filtered))
colnames(se.matrix) <- keep.peaks.index
write.table(se.matrix, "SCMP/data/intermediate-data/scATAC/Satpathy_Nature/Satpathy_scATAC_TME_row-cell_col-peak_matrix_top15-percent.csv",sep=",",row.names = TRUE, col.names = TRUE)



sc_metadata <- utils::read.csv("/Users/mamouzgar/phd-projects/SCMP/data/intermediate-data/scATAC/Satpathy_Nature/scATAC_clusters_labeled.txt",sep ="\t",stringsAsFactors = FALSE)
sc_metadata <- sc_metadata %>%
  dplyr::select(Group_Barcode, cell.type)

se.df <- as.data.frame(se.matrix)
se.df$cell.id <-rownames(se.df)
se.df <- se.df %>%
  left_join(.,sc_metadata, by  = c("cell.id" = "Group_Barcode")) %>%
  dplyr::rename(labels=cell.type) %>%
  dplyr::select(cell.id, labels, everything())


write.table(se.df, "SCMP/data/intermediate-data/scATAC/Satpathy_Nature/Satpathy_scATAC_TME_row-cell_col-peak_df_top15-percent.csv",sep=",",row.names = TRUE, col.names = TRUE)




