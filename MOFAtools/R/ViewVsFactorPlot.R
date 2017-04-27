
#' @title View versus Factor plot
#' @name ViewFactorPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @details asd
#' @return fill this
#' @reference fill this
#' @import pheatmap
#' @export
ViewFactorPlot <- function(object, views="all", factors="all") {
  
  # Define views
  if (paste0(views,sep="",collapse="") =="all") { 
    views <- viewNames(object) 
  } else {
    stopifnot(all(views %in% viewNames(object)))  
  }
  M <- length(views)
  
  # Define factors
  if (paste0(factors,sep="",collapse="") == "all") { 
    factors <- factorNames(object) 
  } else {
    stopifnot(all(factors %in% factorNames(object)))  
  }
  K <- length(factors)
  
  # Sanity checks
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  # Calculate proportion of residual variation explained by each factor in each view
  fvar_mk <- CalculateVariance_Views(object, views, factors)
    
  # Generate plot
  color <- colorRampPalette(c("grey100", "grey0"))(100)
  pheatmap::pheatmap(fvar_mk, color=color)
}