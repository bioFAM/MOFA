
##########################################################################################
## Functions to analyse the correlation between latent factors and principal components ##
##########################################################################################

#' @title Visualize correlation matrix between the features
#' @name FeaturesCorPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @details asd
#' @return fill this
#' @references fill this
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @export
FeaturesCorPlot <- function(object, view, method="pearson", regress_factors=NULL, top=500, ...) {
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  # Select 'top' most variable features
  data <- object@TrainData[[view]]
  top_features <- names(tail(sort(apply(data,2,var)),n=top))
  data <- data[,top_features]
  
  # Regress out latent variables
  if (!is.null(regress_factors)) {
    if (regress_factors=="all") { regress_factors <- 1:object@Dimensions$K }
    SW <- getExpectations(object,"SW","E")[[view]]#; rownames(SW) <- colnames(object@TrainData[[view]])
    Z <- getExpectations(object,"Z","E")
    data <- data - t(SW[top_features,] %*% t(Z[,regress_factors]))
  }
  
  # Compute correlation
  r <- cor(data, method=method)
  
  # Draw heatmap
  breaksList = seq(-1,1, by=0.01)
  pheatmap(r, 
    color=colorRampPalette(rev(brewer.pal(n=10, name="RdYlBu")))(length(breaksList)),
    breaks=breaksList,
    cluster_cols=F, cluster_rows=F, 
    show_colnames=F, show_rownames=F, ...
  )
}


#' @title Visualize correlation matrix between the latent variables
#' @name FactorsCorPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param method a character string indicating which correlation coefficient is to be computed: pearson (default), kendall, or spearman.
#' @param ... arguments passed to \code{corrplot}
#' @details asd
#' @return fill this
#' @references fill this
#' @import corrplot
#' @export
FactorsCorPlot <- function(object, method="pearson", ...) {
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  Z <- getExpectations(object,"Z","E")
  if(object@ModelOpts$learnMean==T) Z <- Z[,-1]
  h <- cor(x=Z, y=Z, method=method)
  p <- corrplot::corrplot(h, tl.col="black", ...)
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
#' @import ggplot2 dplyr tidyr
#' @export
FactorsScatterPlot <- function(object, z_order=NULL, title="") {
  # THIS HAS TO BE FINISHED, WE SHOULDNT USE PIPES OR DPLYR 
  if (is.null(z_order))
    z_order <- 1:object@Dimensions[["K"]]
  
  # Convert latent variable matrix into dataframe
  Z <- Z[,z_order]; colnames(Z) <- z_order
  Z_long <- as.data.frame(t(Z)) %>% dplyr::tbl_df %>% dplyr::mutate(Z=factor(1:n())) %>%
    tidyr::gather(sample,value,-Z)
  
  # Concate PCs
  joined <- Z_long %>% inner_join(Z_long, by='sample') %>%
    dplyr::rename(Z1=Z.x, Z2=Z.y, value1=value.x, value2=value.y)
  
  # Plot
  p <- ggplot(joined, aes(x=value1, y=value2)) +
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
#' @references fill this
#' @import corrplot pcaMethods
#' @importFrom stats cor
#' @importFrom plyr rbind.fill.matrix
#' @export

CorrplotLFvsPC <- function(model, views="all", noPCs=5, method="svd") {

  if (class(model) != "MOFAmodel") stop("'model' has to be an instance of MOFAmodel")
  if (views[1]=="all") { views <- viewNames(model) } else { stopifnot(all(views%in%viewNames(model))) }
  
  # Collect expectations
  # Z <- getExpectations(model,"Z","E")
  Z <- getFactors(model)
  if (model@ModelOpts$learnMean) Z <- Z[,-1]
  
  
  # Perform PCAs
  listPCs <- lapply(views, function(m) {
    data <- getTrainData(model,m)
    # Replace Nan by NA
    data[is.nan(data)] <- NA
    # Remove samples with missing views
    data <- data[,apply(data,2, function(x) mean(is.na(x))) < 1]
    # Perform PCA
    pc.out <- pcaMethods::pca(t(data), method=method, center=TRUE, scale="none", nPcs=noPCs)
    # Extract principal components
    tmp <- t(pc.out@scores)
    rownames(tmp) <- paste(m, rownames(tmp), sep="_")
    tmp
  })
  
  # Calculate correlation matrix between latent factors and PCs
  # matPCs <- do.call(cbind,listPCs)
  matPCs <- t(plyr::rbind.fill.matrix(listPCs))
  colnames(matPCs) <- unlist(lapply(listPCs,rownames))
  corrmatrix <- cor(matPCs,Z, use="complete.obs")
  
  # Plot correlation matrix
  # corrplot::corrplot(t(corrmatrix), order="original", title="", tl.col="black", mar=c(1,1,3,1))
  corrplot::corrplot(t(abs(corrmatrix)), order="original", title="", tl.col="black")
  
  return(corrmatrix)
}
