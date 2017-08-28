
##########################################################
## Functions to cluster samples based on latent factors ##
##########################################################

#' @title Hierarchical clustering on samples based on latent factors
#' @name clusterSamples
#' @description latent factors are continuous in nature but they can be used to predict clusters of samples, similar to what the iCluster model does (Shen, 2009). \cr
#' The clustering can be performed in a single factor, which is equivalent to setting a manual threshold or using multiple factors. \cr
#' The result of the clustering is displayed using a heatmap of the latent factors.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factors vector with the factors indices (numeric) or factor names (character) to use for clustering. Default is "all"
#' @param ... further arguments to be passed to pheatmap
#' @return Plots a heatmap and returns a \link{hclust} object containing the clusters
#' @importFrom pheatmap pheatmap
#' @export
#' 

clusterSamples <- function(object, factors = "all", ...) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # Define factors
  if (paste(factors, collapse="")=="all") { 
    factors <- factorNames(object) 
  } else {
    stopifnot(all(factors %in% factorNames(object)))  
  }
  
  # Collect relevant data
  Z <- getFactors(object, factors=factors)
  N <- getDimensions(object)[["N"]]
  
  # Perform hierarchical clustering
  hc.out <- hclust(dist(Z))

  # Plot heatmap
  pheatmap::pheatmap(t(Z), cluster_rows = length(factors)>1, show_colnames = T, ...)
  
  return(hc.out)

}
