
#' @title View versus Factor plot
#' @name ViewFactorPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @details asd
#' @return fill this
#' @reference fill this
#' @import pheatmap
#' @export
ViewFactorPlot <- function(model, ...) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  # Calculate proportion of residual variation explained by each factor in each view
  prvar_mk <- CalculateVariance_Views(model)
  
  # Generate plot
  color <- colorRampPalette(c("grey100", "grey0"))(100)
  pheatmap::pheatmap(prvar_mk, color=color)
}