
#' @title Visualize correlation matrix between the latent variables
#' @name FactorsCorPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param method a character string indicating which correlation coefficient is to be computed: pearson (default), kendall, or spearman.
#' @param ... arguments passed to \code{corrplot}
#' @details asd
#' @return fill this
#' @reference fill this
#' @import corrplot
#' @export
FactorsCorPlot <- function(object, method="pearson", ...) {
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  Z <- model@Expectations$Z$E
  h <- cor(x=Z, y=Z, method=method)
  p <- corrplot::corrplot(h, ...)
  return(p)
}


#' @title Faced scatterplot between latent variables
#' @name FactorsScatterPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param z_order ...
#' @param title fill this
#' @details asd
#' @return fill this
#' @import ggplot2 dplyr
#' @export
FactorsScatterPlot <- function(object, z_order=NULL, title="") {
  
  if (is.null(z_order))
    z_order <- 1:object@Dimensions[["K"]]
  
  # Convert latent variable matrix into dataframe
  Z <- Z[,z_order]; colnames(Z) <- z_order
  Z_long <- as.data.frame(t(Z)) %>% dplyr::tbl_df %>% dplyr::mutate(Z=factor(1:n())) %>%
    dplyr::gather(sample,value,-Z)
  
  # Concate PCs
  joined <- Z_long %>% inner_join(Z_long, by='sample') %>%
    dplyr::rename(Z1=Z.x, Z2=Z.y, value1=value.x, value2=value.y)
  
  # Plot
  p <- ggplot2::ggplot(joined, aes(x=value1, y=value2)) +
    ggtitle(title) +
    stat_smooth(method=lm, color='black') +
    geom_point(aes(color=sample), size=0.5, alpha=0.5) +
    xlab('') + ylab('') +
    facet_grid(Z1~Z2) +
    guides(color=F)
  
  return(p)
}
