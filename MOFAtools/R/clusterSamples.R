
##########################################################
## Functions to cluster sampels based on latent factors ##
##########################################################


#' @title Hierarchical clustering on samples based on latent factors
#' @name clusterSamples
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factors factors to use for the clustering
#' @param anno_df annotation data frame that will be passed to pheatmap
#' @param main title of the plot
#' @details If only one factor (...). Plots heatmap
#' @return A hclust object containing the clustering of the samples
#' @import pheatmap
#' @export
#' 

clusterSamples <- function(object, factors="all", anno_df=NULL, main=NULL, ...){
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # Collect relevant data
  Z <- object@Expectations$Z$E
  N <- object@Dimensions[["N"]]
  
  
  # Define factors
  if (paste(factors, collapse="")=="all") { 
    factors <- factorNames(object) 
    if(is.null(factors)) factors <- 1:ncol(Z)
  } else {
    stopifnot(all(factors %in% factorNames(object)))  
  }
  
  # Perform hierarchical clustering
  hc.out <- hclust(dist(Z[, factors]))

  # Add annotation data frame
  if (!is.null(anno_df)) {
    if (!is.data.frame(anno_df)) stop("anno_df should be a dataframe containing information on samples")
    if (nrow(anno_df) != N) stop("The number of rows in anno_df has to match the number of samples in the model")
    if (!setequal(rownames(anno_df), rownames(Z))) stop("anno_df needs to have rownames matching the samples names in the MOFA object")
  } 
  
  # Plot heatmap
  if(is.null(main)) main <- "Clustering based on latent factors"
  pheatmap::pheatmap(t(Z[,factors, drop=F]), annotation_col = anno_df, 
                     cluster_rows = length(factors)>1, show_colnames = T,
                     main = main, ...)
  
  return(hc.out)

}
