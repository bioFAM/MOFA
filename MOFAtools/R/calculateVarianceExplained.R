
#' @title Calculate variance explained by the model
#' @name calculateVarianceExplained
#' @description Method to calculate variance explained by the MOFA model for each view and latent factor. \cr
#' As a measure of variance explained for gaussian data we adopt the coefficient of determination (R2). \cr
#' For non-gaussian views the calculations are based on the normally-distributed pseudo-data 
#' (for more information on the non-gaussian model see Supplementary Methods of the MOFA paper or Seeger & Bouchard, 2012).
#' @param object a trained \code{\link{MOFAmodel}} object.
#' @param views character vector with the view names, or numeric vector with view indexes. Default is 'all'
#' @param factors character vector with the factor names, or numeric vector with the factor indexes. Default is 'all'
#' @param include_intercept include the intercept factor for calculation of variance explained (only used when an intercept was learned)
#' @details This function takes a trained MOFA model as input and calculates for each view the coefficient of determination (R2),
#' i.e. the proportion of variance in the data explained by the MOFA factor(s) (both jointly and for each individual factor). 
#' In case of non-Gaussian data the variance explained on the Gaussian pseudo-data is calculated. 
#' @return a list with matrices with the amount of variation explained per factor and view and the total variance explained per view.
#' @export
#' @examples
#' # Using an existing trained model on the CLL data
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' MOFA_CLL <- loadModel(filepath)
#' plotVarianceExplained(MOFA_CLL)
#'
#' # Using an existing trained model on the scMT data
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAtools")
#' MOFA_scMT <- loadModel(filepath)
#' plotVarianceExplained(MOFA_scMT)

calculateVarianceExplained <- function(object, views = "all", factors = "all", include_intercept = TRUE) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # check whether the intercept was learned
  if(!object@ModelOptions$learnIntercept==TRUE & include_intercept==TRUE) {
    include_intercept <- FALSE
    # warning("No intercept was learned in MOFA.\n Intercept is not included in the model prediction.")
  }
  
  # Define views
  if (paste0(views,sep="",collapse="") =="all") { 
    views <- viewNames(object) 
  } else {
    stopifnot(all(views %in% viewNames(object)))  
  }
  M <- length(views)
  
  # Define factors
  if (paste0(factors,collapse="") == "all") { 
    factors <- factorNames(object) 
  } else if (is.numeric(factors)) {
    if (include_intercept == T) {
      factors <- factorNames(object)[factors+1] 
    } else {
      factors <- factorNames(object)[factors]
    }
  } else { 
    stopifnot(all(factors %in% factorNames(object))) 
  }
  factors <- factors[factors!="intercept"]
  K <- length(factors)
  
  # Collect relevant expectations
  W <- getWeights(object,views,factors)
  Z <- getFactors(object, factors)
  Y <- getExpectations(object,"Y") # for non-Gaussian likelihoods the pseudodata is considered
  
  # replace masked values on Z by 0 (so that they do not contribute to predictions)
  Z[is.na(Z)] <- 0 
  
  # If an intercept is included, regress out the intercept from the data
  if (include_intercept) {
    intercept <- getWeights(object,views,"intercept")
    Y <- lapply(views, function(m) sweep(Y[[m]],2,intercept[[m]],"-"))
    names(Y) <- views
  }
  
  # Calculate coefficient of determination per view
  tmp <- sapply(views, function(m) sum(scale(Y[[m]],center=T, scale=F)**2, na.rm=T))
  fvar_m <- sapply(views, function(m) 1 - sum((Y[[m]]-tcrossprod(Z,W[[m]]))**2, na.rm=T) / tmp[m])
  names(fvar_m) <- views
  
  # Calculate coefficient of determination per factor and view
  fvar_mk <- matrix(99, ncol=length(views), nrow=length(factors))
  colnames(fvar_mk) <- views
  rownames(fvar_mk) <- factors
  for (m in views) {
    for (k in factors) {
      fvar_mk[k,m] <- 1 - sum( (Y[[m]]-tcrossprod(Z[,k],W[[m]][,k]))**2, na.rm=T) / tmp[m]
    }
  }
  
  # Replace negative values by zero
  fvar_m[fvar_m<0] <- 0
  fvar_mk[fvar_mk<0] <- 0
  
  # Store results
  R2_list <- list(R2Total = fvar_m, R2PerFactor = fvar_mk)
  
  return(R2_list)
  
}

#' @title Plot variance explained by the model
#' @name plotVarianceExplained
#' @description Method to plot variance explained (R-squared) by the MOFA model for each view and latent factor. \cr
#' As a measure of variance explained for gaussian data we adopt the coefficient of determination (R2). \cr
#' For details on the computation see the help of the \code{\link{calculateVarianceExplained}} function
#' @param object a \code{\link{MOFAmodel}} object.
#' @param cluster logical indicating whether to do hierarchical clustering on the plot
#' @param ... extra arguments to be passed to \code{\link{calculateVarianceExplained}}
#' @return ggplot object
#' @import pheatmap ggplot2 reshape2
#' @importFrom cowplot plot_grid
#' @export
#' @examples
#' # Using an existing trained model on the CLL data
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' MOFA_CLL <- loadModel(filepath)
#' plotVarianceExplained(MOFA_CLL)
#'
#' # Using an existing trained model on the scMT data
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAtools")
#' MOFA_scMT <- loadModel(filepath)
#' plotVarianceExplained(MOFA_scMT)

plotVarianceExplained <- function(object, cluster = TRUE, ...) {
  
  # Calculate Variance Explained
  R2_list <- calculateVarianceExplained(object, ...)
  fvar_m <- R2_list$R2Total
  fvar_mk <- R2_list$R2PerFactor
  
  ## Plot variance explained by factor ##
  
  # convert matrix to data frame for ggplot2  
  fvar_mk_df <- reshape2::melt(fvar_mk, varnames=c("factor","view"))
  fvar_mk_df$factor <- factor(fvar_mk_df$factor)
  
  # If multiple views, sort factors according to hierarchical clustering
  if (cluster & ncol(fvar_mk)>1) {
    hc <- hclust(dist(t(fvar_mk)))
    fvar_mk_df$view <- factor(fvar_mk_df$view, levels = colnames(fvar_mk)[hc$order])
  }
  
  # Grid plot with the variance explained per factor and view
  hm <- ggplot(fvar_mk_df, aes(view,factor)) + 
    geom_tile(aes(fill=value), color="black") +
    guides(fill=guide_colorbar("R2")) +
    scale_fill_gradientn(colors=c("gray97","darkblue"), guide="colorbar") +
    ylab("Latent factor") +
    theme(
      # plot.margin = margin(5,5,5,5),
      plot.title = element_text(size=17, hjust=0.5),
      axis.title.x = element_blank(),
      axis.text.x = element_text(size=11, angle=60, hjust=1, vjust=1, color="black"),
      axis.text.y = element_text(size=12, color="black"),
      axis.title.y = element_text(size=15),
      axis.line = element_blank(),
      axis.ticks =  element_blank(),
      panel.background = element_blank()
    )
  hm <- hm + ggtitle("Variance explained per factor")  + 
    guides(fill=guide_colorbar("R2"))
  
  ## Plot variance explained per view ##
  
  # Create data.frame for ggplot
  fvar_m_df <- data.frame(view=factor(names(fvar_m), levels=names(fvar_m)), R2=fvar_m)
  
  # If multiple views, sort factors according to hierarchical clustering
  if (cluster==TRUE & ncol(fvar_mk)>1) {
    fvar_m_df$view <- factor(fvar_m_df$view, levels = colnames(fvar_mk)[hc$order])
  }
  
  # Barplot with variance explained per view
  bplt <- ggplot(fvar_m_df, aes(x=view, y=R2)) + 
    ggtitle("Total variance explained per view") +
    geom_bar(stat="identity", fill="deepskyblue4", width=0.9) +
    xlab("") + ylab("R2") +
    scale_y_continuous(expand=c(0.01,0.01)) +
    theme(
      plot.margin = unit(c(1,2.4,0,0), "cm"),
      panel.background = element_blank(),
      plot.title = element_text(size=17, hjust=0.5),
      axis.ticks.x = element_blank(),
      # axis.text.x = element_text(angle = 90, hjust = 1, vjust=0.5)
      axis.text.x = element_blank(),
      axis.text.y = element_text(size=12, color="black"),
      axis.title.y = element_text(size=13, color="black"),
      axis.line = element_line(size=rel(1.0), color="black")
    )
  
  # Join the two plots
  p <- plot_grid(bplt, hm, align="v", nrow=2, rel_heights=c(1/3,2/3), axis="l")
  
  return(p)
}


