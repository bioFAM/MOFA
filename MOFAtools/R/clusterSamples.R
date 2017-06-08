#' @title Cluster samples in a MOFA model based on latent factors
#' @name clusterSamples
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factors factors to use for the clustering
#' @details asd
#' @return A hclust object containing the clustering of the samples
#' @reference fill this
#' @export
#' 

clusterSamples <- function(object, factors="all", anno_df=NULL, main=NULL){
  
  # Sanity checks
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  Z <- model@Expectations$Z$E
  N <- model@Dimensions[["N"]]
  
  if(is.null(main)) main <- "Clustering based on latent factors"
  
  # Define factors
  if (factors=="all") { 
    factors <- factorNames(model) 
    if(is.null(factors)) factors <- 1:ncol(Z)
  } else {
    stopifnot(all(factors %in% factorNames(model)))  
  }
  
  hc.out <- hclust(dist(Z[, factors]))

  if(!is.null(anno_df)) {
  if(!is.data.frame(anno_df)) stop("anno_df should be a dataframe containing covaraites on samples")
  if(nrow(anno_df) != N) stop("The number of rows in anno_df has to match the number of samples in the model")
  if(!setequal(rownames(df_anno), rownames(Z))) stop("df_anno needs to have rownames matching the samples names in the MOFA object")
  } 
  
  pheatmap::pheatmap(Z[,factors, drop=F], annotation_row = anno_df, 
                       cluster_cols = length(factors)>1, show_rownames = T,
                     main = main)
  
  return(hc.out)

}
