
#####################################
## Functions to visualise the data ##
#####################################

#' @title showDataHeatmap: heatmap of the observed data of the most relevant features for a given factor
#' @name showDataHeatmap
#' @description Function to plot a heatmap of the original data using the most relevant features of a given latent variables
#' @param object a fitted MOFA model
#' @param view name of view
#' @param factor factor
#' @param nfeatures total number of features, sorted by weight (default=???)
#' @param ... further arguments that can be passed to pheatmap
#' @details fill this
#' @return fill this
#' @import pheatmap
#' @export
showDataHeatmap <- function(object, view, factor, nfeatures=50, manual_features=NULL, include_weights=F, transpose=F, ...) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(view %in% viewNames(object))
  stopifnot(factor %in% factorNames(object)) 
  
  # Collect relevant expectations
  W <- getExpectations(object,"SW","E")[[view]][,factor]
  Z <- getFactors(object)[,factor]
  
  # Define features
  features <- c()
  if (nfeatures>0) {
    tmp <- names(tail(sort(abs(W)), n=nfeatures))
    stopifnot(all(tmp %in% featureNames(object)[[view]]))
    features <- c(features,tmp)
  }
  if (length(manual_features)>0) {
    stopifnot(all(manual_features %in% featureNames(object)[[view]]))
    features <- c(features,manual_features)
  }
  
  # Get train data
  data <- getTrainData(object, view)
  
  # Ignore samples with full missing views
  data <- data[,apply(data, 2, function(x) !all(is.na(x)))]
  
  # Sort samples according to latent factors
  order_samples <- names(sort(Z, decreasing=T))
  order_samples <- order_samples[order_samples %in% colnames(data)]
  data <- data[features,order_samples]
  
  if (transpose==T) {
    data <- t(data)
  }
  
  # Plot heatmap
  # if(is.null(main)) main <- paste(view, "observations for the top weighted features of factor", factor)
  if (include_weights) { 
    anno <- data.frame(row.names=names(W[features]), weight=W[features]) 
    if (transpose) {
      pheatmap::pheatmap(t(data), annotation_col=anno, ...)
    } else {
      pheatmap::pheatmap(t(data), annotation_row=anno, ...)
    }
  } else {
    pheatmap::pheatmap(t(data), ...)
  }
  
}



#' @title showDataScatter: scatterplot of the observed data for the most relevant features for a given factor
#' @name showDataScatter
#' @description Function to plot a scatterplot of the observed values for the most relevant features for a given factor
#' @param object a fitted MOFA model
#' @param view name of view
#' @param factor factor
#' @param nfeatures total number of features (default=???)
#' @param color_by vector
#' @param shape_by vector
#' @details fill this
#' @return ggplot2 model
#' @import ggplot2
#' @import dplyr
#' @export

showDataScatter <- function(object, view, factor, nfeatures=50, color_by=NULL, shape_by=NULL,
                            xlabel="", ylabel="") {
  
  # Sanity checcks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(length(factor)==1)
  stopifnot(length(view)==1)
  if (!factor %in% factorNames(object)) stop(sprintf("The factor %s is not present in the object",factor))
  if (!view %in% viewNames(object)) stop(sprintf("The view %s is not present in the object",view))
  
  # Collect data
  N <- object@Dimensions[["N"]]
  z <- getExpectations(object, "Z", "E")[,factor]
  w <- getExpectations(object, "SW", "E")[[view]][,factor]
  Y <- object@TrainData[[view]]
  
  # Select top features (by absolute value)
  features <- names(tail(sort(abs(w)), n=nfeatures))
  w <- w[features]
  Y <- Y[features,]
  
  
  # Set color
  if (!is.null(color_by)) {
    if (length(color_by) == 1 & is.character(color_by)) { # It is the name of a covariate 
      color_by <- as.factor(getCovariates(object, color_by))
    } else if (length(color_by) > 1) { # It is a vector of length N
      stopifnot(length(color_by) == N)
    } else {
      stop("'color_by' was specified but it was not recognised, please read the documentation")
    }
    colorLegend <- T
  } else {
    color_by <- rep(TRUE,N)
    colorLegend <- F
  }
  
  # Set shape
  if (!is.null(shape_by)) {
    if (length(shape_by) == 1 & is.character(shape_by)) { # It is the name of a covariate 
      shape_by <- as.factor(getCovariates(object, shape_by))
    } else if (length(shape_by) > 1) { # It is a vector of length N
      stopifnot(length(shape_by) == N)
      shape_by <- as.factor(shape_by)
    } else {
      stop("'shape_by' was specified but it was not recognised, please read the documentation")
    }
    shapeLegend <- T
  } else {
    shape_by <- rep(TRUE,N)
    shapeLegend <- F
  }
  
  
  # Create data frame 
  df1 <- data.frame(sample=names(z), x = z, shape_by = shape_by, color_by = color_by, stringsAsFactors=F)
  df2 <- getTrainData(object, views=view, features = list(features), as.data.frame=T)
  df <- left_join(df1,df2, by="sample")
  
  #remove values missing color or shape annotation
  # if(!showMissing) df <- df[!(is.nan(df$shape_by) & !(is.nan(df$color_by))]
  
  p <- ggplot(df, aes(x, value, color = color_by, shape = shape_by)) + 
    geom_point(color="black") + 
    ggtitle("") + xlab(xlabel) + ylab(ylabel) + 
    stat_smooth(method="lm", color="blue", alpha=0.5) +
    facet_wrap(~feature, scales="free_y") +
    scale_shape_manual(values=c(19,1,2:18)[1:length(unique(shape_by))]) +
    theme(plot.margin = margin(20, 20, 10, 10), 
          axis.text = element_text(size = rel(1), color = "black"), 
          axis.title = element_text(size = 16), 
          axis.title.y = element_text(size = rel(1.1), margin = margin(0, 15, 0, 0)), 
          axis.title.x = element_text(size = rel(1.1), margin = margin(15, 0, 0, 0)), 
          axis.line = element_line(color = "black", size = 0.5), 
          axis.ticks = element_line(color = "black", size = 0.5),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.key = element_rect(fill = "white")
          # legend.text = element_text(size = titlesize),
          # legend.title = element_text(size =titlesize)
          )
  if (colorLegend) { p <- p + labs(color = name_color) } else { p <- p + guides(color = FALSE) }
  if (shapeLegend) { p <- p + labs(shape = name_shape) }  else { p <- p + guides(shape = FALSE) }
  
  return(p)
}


