
#################################################################
## Functions to analyse the correlation between latent factors ##
#################################################################

#' @title Visualize correlation matrix between the features
#' @name FeaturesCorPlot
#' @description the function plots the correlation between all features in a given view. 
#' This is useful to see if the factors learnt by the model capture the correlations between features, as the residuals should be uncorrelated.
#' To check this, run the function twice, first without regressing out any factor (regress_factors=NULL) and then regress all factors (regress_factors="all").
#' The first plot should be enriched by correlations whereas the second should be fairly uncorrelated. 
#' If not, it suggests that there is more signal to be learnt and you should increase the number of factors.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param view view name
#' @param features character vector with the feature names or numeric vector with the feature indices (default is "all").
#' @param method a character string indicating which correlation coefficient is to be computed: pearson (default), kendall, or spearman.
#' @param regress_factors character or numeric vector with the factors to regress (default is NULL)
#' @param ... arguments to be passed to \link{pheatmap} function
#' @return correlation matrix between features
#' @importFrom pheatmap pheatmap
#' @importFrom RColorBrewer brewer.pal
#' @importFrom grDevices colorRampPalette
#' @export
FeaturesCorPlot <- function(object, view, features = "all", method = "pearson", regress_factors = NULL, ...) {
  
  # Sanity check
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # Fetch data
  data <- object@TrainData[[view]]
  
  # Feature selection
  if (paste0(features,collapse="") == "all") { 
    features <- featureNames(object)[[view]]
  } else {
    stopifnot(all(features %in% featureNames(model)[[view]]))  
  }
  # top_features <- names(tail(sort(apply(data,2,var)),n=top))
  data <- data[,features]
  
  # Regress out latent variables
  if (!is.null(regress_factors)) {
    if (regress_factors=="all") { regress_factors <- 1:getDimensions(object)[["K"]] }
    SW <- getWeights(views=view)
    # SW <- getExpectations(object,"SW","E")[[view]]
    # Z <- getExpectations(object,"Z","E")
    Z <- getFactors(object)
    data <- data - t(SW[features,] %*% t(Z[,regress_factors]))
  }
  
  # Compute correlation
  r <- cor(data, method=method)
  
  # Draw heatmap
  breaksList = seq(-1,1, by=0.01)
  pheatmap::pheatmap(r, 
    color=colorRampPalette(rev(brewer.pal(n=10, name="RdYlBu")))(length(breaksList)),
    breaks=breaksList,
    cluster_cols=F, cluster_rows=F, 
    show_colnames=F, show_rownames=F, ...
  )
  
  return(r)
}


#' @title Visualize correlation matrix between the latent variables
#' @name FactorsCorPlot
#' @description plot the correlation matrix between the latent variables. The model encourages factors to be uncorrelated, but it is nos hard constraint such as in PCA.
#' Ideally all learnt factors should be uncorrelated to make interpretation easier, but correlations can happen, particularly with large number factors or with non-linear relationship.
#' If you have too many correlated factors try to train the model again reducing the number of factors.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param method a character string indicating which correlation coefficient is to be computed: pearson (default), kendall, or spearman.
#' @param ... arguments passed to \code{corrplot}
#' @return symmetric matrix with the correlation coefficient between factors
#' @import corrplot
#' @export
FactorsCorPlot <- function(object, method = "pearson", ...) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # Fetch factors
  Z <- getFactors(object)
  
  # Remove intercept
  if(object@ModelOpts$learnIntercept==TRUE) Z <- Z[,-1]
  
  # Compute and plot correlation
  r <- cor(x=Z, y=Z, method=method, use = "complete.obs")
  p <- corrplot::corrplot(r, tl.col="black", ...)
  
  return(r)
}


#' #' @title Faced scatterplot between latent variables
#' #' @name FactorsScatterPlot
#' #' @description plot a faced scatterplot with combinations of latent variables
#' #' @param object a \code{\link{MOFAmodel}} object.
#' #' @param factors factors to plot (default="all")
#' #' @details asd
#' #' @return fill this
#' #' @import ggplot2 dplyr tidyr
#' #' @export
#' FactorsScatterPlot <- function(object, factors="all") {
#'   # THIS HAS TO BE FINISHED, WE SHOULDNT USE PIPES OR DPLYR
#'   if (is.null(factors))
#'     factors <- 1:object@Dimensions[["K"]]
#' 
#'   # Convert latent variable matrix into dataframe
#'   Z <- Z[,factors]; colnames(Z) <- factors
#'   Z_long <- as.data.frame(t(Z)) %>% dplyr::tbl_df %>% dplyr::mutate(Z=factor(1:n())) %>%
#'     tidyr::gather(sample,value,-Z)
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
#' 
#'   return(p)
#' }

#' @title Correlation plot of latent factors to principal components of single views
#' @name CorrplotLFvsPC
#' @description This function is used to identify the relationship between latent factors and principal components.
#' Usually the first latent factors also correspond to the main principal components,
#' We only recommend to run this function with complete data or with small number of missing values (with method="nipals"). 
#' PCA is very inaccurate with large number of missing values.
#' This function uses the package \link{pcaMethods} to calculate the PCA solution
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views view names (default is "all")
#' @param noPCs Number of principal components to calculate
#' @param method method to compute PCA, we recommend svd (traditional PCA no missing values) or nipals (small number of missing values)
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
  if (model@ModelOpts$learnIntercept==TRUE) { Z <- Z[,-1] }
  
  
  # Perform PCAs
  listPCs <- lapply(views, function(m) {
    data <- getTrainData(model,m)
    # Replace NaN by NA
    data[is.nan(data)] <- NA
    # Remove missing samples
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
