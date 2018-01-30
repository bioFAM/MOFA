
##########################################################
## Functions to cluster samples based on latent factors ##
##########################################################

#' @title K-means clustering on samples based on latent factors
#' @name clusterSamples
#' @description latent factors are continuous in nature but they can be used to predict clusters of samples, similar to what the iCluster model does (Shen, 2009). \cr
#' The clustering can be performed in a single factor, which is equivalent to setting a manual threshold; or using multiple factors, where you combine multiple sources of variation. \cr
#' @param object a \code{\link{MOFAmodel}} object.
#' @param k number of clusters
#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. Default is 'all'\#' @param factors character vector with the factor name(s), or numeric vector with the index of the factor(s) to use. Default is 'all'
#' @param ... other arguments that can be passed to the kmeans function
#' @details In some cases, due to model technicalities, samples can have missing values in the latent factor space. In such a case, these samples are currently removed.
#' @return output from \code{\link{kmeans}} function
#' @export
#' 

clusterSamples <- function(object, k, factors = "all", ...) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  

  # Define factors
  if (paste0(factors,collapse="") == "all") { factors <- factorNames(object) } 
    else if(is.numeric(factors)) {
      if (object@ModelOpts$learnIntercept == T) factors <- factorNames(object)[factors+1]
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
  Z <- Z[haveAllZ,]

  # Perform k-means clustering
  kmeans.out <- kmeans(Z, centers=k,  ...)

  return(kmeans.out)  

}
