
library(corrplot)
library(dplyr)
library(tidyr)

#' @title Visualize correlation matrix between the latent variables
#' @name FactorsCorPlot
#' @description: fill this ...
#' @param object a \code{\link{MOFAmodel}} object.
#' @param method a character string indicating which correlation coefficient is to be computed: pearson (default), kendall, or spearman.
#' @details asd
#' @return fill this...
#' @export
#' @reference fill this...
#' @examples fill this...

FactorsCorPlot <- function(object, method="pearson", ...) {
  Z <- model@Expectations$Z$E
  h <- cor(x=Z, y=Z, method=method)
  p <- corrplot(h, ...)
  return(p)
}
