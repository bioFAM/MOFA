
##########################################################
## Functions to cluster samples based on latent factors ##
##########################################################

#' @title clusterSamples: K-means clustering on samples based on latent factors
#' @name clusterSamples
#' @description MOFA factors are continuous in nature but they can be used to predict discrete clusters of samples, 
#' similar to the iCluster model (Shen, 2009). \cr
#' The clustering can be performed in a single factor, which is equivalent to setting a manual threshold; 
#' or using multiple factors, where multiple sources of variation are aggregated. \cr
#' Importantly, this type of clustering is not weighted and does not take into account the different importance of the latent factors. 
#' @param object a trained \code{\link{MOFAmodel}} object.
#' @param k number of clusters
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. 
#' Default is 'all'
#' @param ... extra arguments  passed to \code{\link{kmeans}}
#' @details In some cases, samples can have missing values in the latent factor space. 
#' This happens if a factors is active only in views where these samples have no data.
#' In such a case, these samples are currently ignored in the clustering procedure and NAs are returned.
#' @return output from \code{\link{kmeans}} function
#' @export
#' @examples
#' # Example on the CLL data
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' MOFA_CLL <- loadModel(filepath)
#' # cluster samples based into 3 groups based on all factors
#' clusterSamples(MOFA_CLL, k=3, factors="all")
#' # cluster samples based into 2 groups based on factor 1
#' clusters <- clusterSamples(MOFA_CLL, k=2, factors=1)
#' # cluster can be visualized for example on the factors values:
#' plotFactorBeeswarm(MOFA_CLL, factor=1, color_by=clusters)
#'
#' # Example on the scMT data
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAtools")
#' MOFA_scMT <- loadModel(filepath)
#' # cluster samples based into 2 groups based on all factor 1 and 2
#' clusters <- clusterSamples(MOFA_CLL, k=2, factors=1:2)
#' # cluster can be visualized for example on the factors values:
#' plotFactorScatter(MOFA_CLL, factors=1:2, color_by=clusters)


clusterSamples <- function(object, k, factors = "all", ...) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  

  # Define factors
  if (paste0(factors,collapse="") == "all") { factors <- factorNames(object) } 
    else if(is.numeric(factors)) {
      if (object@ModelOptions$learnIntercept == T) factors <- factorNames(object)[factors+1]
      else factors <- factorNames(object)[factors]
    }
      else{ stopifnot(all(factors %in% factorNames(object))) }
  
  # Collect relevant data
  Z <- getFactors(object, factors=factors, include_intercept=F)
  N <- getDimensions(object)[["N"]]
  
  # For now remove sample with missing values on factors
  # (TO-DO) incorporate a clustering function that is able to cope with missing values
  haveAllZ <- apply(Z,1, function(x) all(!is.na(x)))
  if(!all(haveAllZ)) warning(paste("Removing", sum(!haveAllZ), "samples with missing values on at least one factor"))
  Z_sub <- Z[haveAllZ,]

  # Perform k-means clustering
  kmeans.out <- kmeans(Z_sub, centers=k,  ...)
  clusters <- rep(NA, length(sampleNames(object)))
  names(clusters) <- sampleNames(object)
  clusters[haveAllZ] <- kmeans.out$cluster
  
  return(clusters)  

}
