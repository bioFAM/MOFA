
#' @title View versus Factor plot
#' @name FactorsCorPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param method a character string indicating which correlation coefficient is to be computed: pearson (default), kendall, or spearman.
#' @param ... arguments passed to \code{corrplot}
#' @details asd
#' @return fill this
#' @reference fill this
#' @import RColorBrewer, pheatmap
#' @export
ViewFactorPlot <- function(model, ...) {
  # Calculate proportion of residual variation explained by each factor in each view
  prvar_mk <- CalculateProportionResidualVariance(model, plot=F)
  
  # Generate plot
  color <- RColorBrewer::colorRampPalette(c("grey100", "grey0"))(100)
  p <- pheatmap::pheatmap(prvar_mk, color=color, ...)
  
  return(p)
}