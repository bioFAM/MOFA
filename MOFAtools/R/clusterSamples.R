
##########################################################
## Functions to cluster samples based on latent factors ##
##########################################################


#' @title Hierarchical clustering on samples based on latent factors
#' @name clusterSamples
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factors factors to use for the clustering (default is "all")
#' @param anno_df annotation data frame that will be passed to pheatmap
#' @param main title of the plot
#' @details if only one factor (...). Plots heatmap
#' @return Plots a heatmap and returns a hclust object containing the clusters
#' @import pheatmap
#' @export
#' 

clusterSamples <- function(object, factors="all", main=NULL, ...){
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # Collect relevant data
  Z <- object@Expectations$Z$E
  N <- object@Dimensions[["N"]]
  
  # Define factors to use
  if (paste(factors, collapse="")=="all") { 
    factors <- factorNames(object) 
    if(is.null(factors)) factors <- 1:ncol(Z)
  } else {
    stopifnot(all(factors %in% factorNames(object)))  
  }
  
  # Perform hierarchical clustering
  hc.out <- hclust(dist(Z[, factors]))

  # Plot heatmap
  if(is.null(main)) main <- "Clustering based on latent factors"
  pheatmap::pheatmap(t(Z[,factors, drop=F]), 
                     cluster_rows = length(factors)>1, show_colnames = T,
                     main = main, ...)
  
  return(hc.out)

}
