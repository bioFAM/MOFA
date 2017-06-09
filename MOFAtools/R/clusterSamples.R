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
  
  Z <- object@Expectations$Z$E
  N <- object@Dimensions[["N"]]
  
  if(is.null(main)) main <- "Clustering based on latent factors"
  
  # Define factors
  if (factors=="all") { 
    factors <- factorNames(object) 
    if(is.null(factors)) factors <- 1:ncol(Z)
  } else {
    stopifnot(all(factors %in% factorNames(object)))  
  }
  
  hc.out <- hclust(dist(Z[, factors]))

  if(!is.null(anno_df)) {
  if(!is.data.frame(anno_df)) stop("anno_df should be a dataframe containing covaraites on samples")
  if(nrow(anno_df) != N) stop("The number of rows in anno_df has to match the number of samples in the model")
  if(!setequal(rownames(anno_df), rownames(Z))) stop("anno_df needs to have rownames matching the samples names in the MOFA object")
  } 
  
  pheatmap::pheatmap(t(Z[,factors, drop=F]), annotation_col = anno_df, 
                       cluster_rows = length(factors)>1, show_colnames = T,
                     main = main)
  
  return(hc.out)

}
