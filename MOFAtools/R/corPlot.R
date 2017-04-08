
#' @title Visualize correlation matrix between the features
#' @name FeaturesCorPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @details asd
#' @return fill this
#' @reference fill this
#' @import pheatmap
#' @export
FeaturesCorPlot <- function(object, view, method="pearson", regress_factors=NULL, top=500, ...) {
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  # Select 'top' most variable features
  data <- model@TrainData[[view]]
  top_features <- names(tail(sort(apply(data,2,var)),n=top))
  data <- data[,top_features]
  
  # Regress out latent variables
  if (!is.null(regress_factors)) {
    if (regress_factors=="all")
      regress_factors <- 1:model@Dimensions$K
    SW <- getExpectations(model,"SW","E")[[view]]; rownames(SW) <- colnames(model@TrainData[[view]])
    Z <- getExpectations(model,"Z","E")
    data <- data - t(SW[top_features,] %*% t(Z[,regress_factors]))
  }
  
  # Compute correlation
  r <- cor(data, method=method)
  
  # Draw heatmap
  p <- pheatmap::pheatmap(r, cluster_cols=F, cluster_rows=F, show_colnames=F, show_rownames=F, ...)
  
  return(p)
}


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

#' @title Correlate latent factors to principal components on single views
#' @name CorrplotLFvsPC
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views fill
#' @param noPCs fill
#' @param method fill, can be svd, ppca, bpca
#' @details fill
#' @return Correlation matrix of latent factors versus principal components
#' @reference fill this
#' @import corrplot pcaMethods
#' @export

CorrplotLFvsPC<-function(model, views="all", noPCs=5, method="svd"){
  
  if (views=="all") {
    views <- viewNames(model)  
  } else {
    stopifnot(all(views%in%viewNames(model)))
  }
  
  # Collect expectations
  Z <- getExpectations(model,"Z","E")
  
  # Perform PCAs
  listPCs <- lapply(views, function(m) {
    pc.out <- pcaMethods::pca(model@TrainData[[m]], method=method, center=TRUE, nPcs=noPCs)
    tmp <- pc.out@scores
    colnames(tmp) <- paste(m, colnames(tmp), sep="_")
    tmp
  })
  
  # Calculate correlation matrix between latent factors and PCs
  matPCs <- do.call(cbind,listPCs)
  corrmatrix <- cor(matPCs,Z)
  
  # Plot correlation matrix
  corrplot::corrplot(corrmatrix, order="original", title=viewname4PC,mar = c(1, 1, 3, 1))
  
  return(corrmatrix)
}
