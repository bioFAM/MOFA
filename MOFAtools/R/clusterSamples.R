
##########################################################
## Functions to cluster samples based on latent factors ##
##########################################################

#' @title K-means clustering on samples based on latent factors
#' @name clusterSamples
#' @description latent factors are continuous in nature but they can be used to predict clusters of samples, similar to what the iCluster model does (Shen, 2009). \cr
#' The clustering can be performed in a single factor, which is equivalent to setting a manual threshold, or using multiple factors. \cr
#' The result of the clustering is displayed using a heatmap of the latent factors.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param k number of clusters
#' @param factors vector with the factors indices (numeric) or factor names (character) to use for clustering. Default is "all"
#' @return Plots a heatmap and returns a \link{hclust} object containing the clusters
#' @importFrom pheatmap pheatmap
#' @export
#' 

clusterSamples <- function(object, k, factors = "all") {
  
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
  
  # For now remove sample with missing values on factors
  # TO-DO incorporate a clustering functions able to cope with missing values
  Z <- Z[apply(Z,1, function(x) all(!is.na(x))),]

  # Perform k-means clustering
  kmeans.out <- kmeans(Z, k)

  return(kmeans.out)  

}
