#' 
#' library(corrplot)
#' library(dplyr)
#' library(tidyr)
#' 
#' #' @title FactorsCorPlot: plot a correlation matrix between latent variables
#' #' @name FactorsCorPlot
#' #' @rdname FactorsCorPlot
#' #' @description: fill this ...
#' #' @param object a \code{\link{GFATrainedModel}} object. 
#' #' @param cormethod a character string indicating which correlation coefficient is to be computed.  One of pearson (default), kendall, or spearman.
#' #' @param diag Logical, whether display the correlation coefficients on the principal diagonal.
#' #' @param display_type Character, "full" (default), "upper" or "lower", display full matrix, lower triangular or upper triangular matrix.
#' #' @param visualisation_method Character, "full" (default), "upper" or "lower", display full matrix, lower triangular or upper triangular matrix.
#' #' @details asd
#' #' @return a \code{\link{GFATrainedModel}} object.
#' #' @export
#' #' @import rhdf5
#' 
#' FactorsCorPlot <- function(object, cormethod="pearson", diag=T, display_type="upper", visualisation_method="circle",
#'                      tl.cex=1.8, cl.cex=0.9, tl.col="black") {
#'   # tl.cex: size of text label for variable names
#'   # cl.cex: size of color label
#'   # tl.col: the colour of the text label
#'   Z <- model@Expectations$Z$E
#'   h <- cor(x=Z, y=Z, method=cormethod)
#'   corrplot(h, method=visualisation_method, order="original", diag=diag, type=display_type,
#'            tl.col=tl.col, tl.cex=tl.cex, cl.cex=cl.cex)
#' }
#' 
#' # corrPlot(gfa, cormethod="pearson", diag=T, display_type="upper", visualisation_method="circle",
#'          # title="", tl.cex=1.8, cl.cex=0.9, tl.col="black")
#' 
#' # Faced scatterplot between latent variables
#' FactorsFacedScatterPlot <- function(object, z_order=NULL, title="") {
#'   
#'   # Convert latent variable matrix into dataframe
#'   if (is.null(z_order))
#'     z_order <- 1:ncol(object@Expectations$Z)
#'   Z <- Z[,z_order]; colnames(Z) <- z_order
#'   Z_long <- as.data.frame(t(Z)) %>% tbl_df %>% mutate(Z=factor(1:n())) %>%
#'     gather(sample, value, -Z)
#'   
#'   # Concate PCs
#'   joined <- Z_long %>% inner_join(Z_long, by='sample') %>%
#'     dplyr::rename(Z1=Z.x, Z2=Z.y, value1=value.x, value2=value.y)
#'   
#'   # Plot
#'   p <- ggplot(joined, aes(x=value1, y=value2)) +
#'     ggtitle(title) +
#'     stat_smooth(method=lm, color='black') +
#'     geom_point(aes(color=sample), size=0.5, alpha=0.5) +
#'     xlab('') + ylab('') +
#'     facet_grid(Z1~Z2) +
#'     guides(color=F)
#'   return(p)
#' }
#' 
#' # FactorsFacedScatterPlot(gfa, z_order=c(1,2,3), title="")
