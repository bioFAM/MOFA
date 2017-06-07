
#####################################
## Functions to visualise the data ##
#####################################

#' @title showDataHeatmap: heatmap of the observed values of the most relevant features for a given factor
#' @name showDataHeatmap
#' @description Function to plot a heatmap of the original data using the most relevant features of a given latent variables
#' @param model a fitted MOFA model
#' @param view name of view
#' @param factor factor
#' @param nfeatures total number of features, sorted by weight (default=???)
#' @param ... further arguments that can be passed to pheatmap
#' @details fill this
#' @return pheatmap model
#' @import pheatmap
#' @export

showDataHeatmap <- function(model, view, factor, nfeatures=50, main=NULL, ...) {
  
  # Sanity checks
  if (class(model) != "MOFAmodel")
    stop("'model' has to be an instance of MOFAmodel")
  stopifnot(view %in% viewNames(model))
  stopifnot(factor %in% factorNames(model)) 
  
  # Collect relevant expectations
  W <- getExpectations(model,"SW","E")[[view]][,factor]
  
  # Define features
  features <- names(tail(sort(abs(W)), n=nfeatures))
  stopifnot(all(features %in% featureNames(model)[[view]]))
  
  # Ignore samples with full missing views
  tmp <- model@TrainData[[view]][,apply(model@TrainData[[view]],2, function(x) !all(is.na(x)))]
  
  # Plot heatmap
  if(is.null(main)) main <- paste(view, "values on top weighted feautres for factor", factor)
  pheatmap::pheatmap(t(tmp[features,]), fontsize = 8, show_rownames = F, main=main,...)
}



#' @title showDataScatter: scatterplot of the observed values for the most relevant features for a given factor
#' @name showDataScatter
#' @description Function to plot a scatterplot of the observed values for the most relevant features for a given factor
#' @param model a fitted MOFA model
#' @param view name of view
#' @param factor factor
#' @param nfeatures total number of features (default=???)
#' @param color_by vector
#' @param shape_by vector
#' @param ... further arguments that can be passed to ...?
#' @details fill this
#' @return ggplot2 model
#' @import ggplot2
#' @export

showDataScatter <- function(model, view, factor, nfeatures=50, colour_by=NULL, shape_by=NULL, ...) {
  
  # Sanity checcks
  if (class(model) != "MOFAmodel") 
    stop("'model' has to be an instance of MOFAmodel")
  if (!factor %in% factorNames(model)) 
    stop(sprintf("The factor %s is not present in the model",factor))
  if (!view %in% viewNames(model)) 
    stop(sprintf("The view %s is not present in the model",view))
  
  # Collect data
  N <- model@Dimensions[["N"]]
  z <- getExpectations(model, "Z", "E")[,factor]
  w <- getExpectations(model, "SW", "E")[[view]][,factor]
  Y <- model@TrainData[[view]]
  
  # Selecct top features
  features <- names(tail(sort(abs(w)), n=nfeatures))
  w <- w[features]
  Y <- Y[,features]
  
  
  # collect color groups
  if (!is.null(colour_by)) {
    if (length(colour_by) != N) 
      stop("'colour_by' has to be a vector of length N")
    colorLegend<-T
  } else {
    colour_by <- rep(TRUE, N)
    colorLegend<-F
  }
  
  # collect shape groups
  if (!is.null(shape_by)) {
    if (length(shape_by) != N) 
      stop("'shape_by' has to be a vector of length N")
    shapeLegend<-T
  } else {
    shape_by <- rep(TRUE, N)
    shapeLegend<-F
  }
  
  # TO-DO: REARRANGE DATA AND PLOT
  
  # Create data frame 
  df = data.frame(
    x = z, y = Y, w = w,
    shape_by = shape_by, colour_by = colour_by
    )
  
  #remove values missing color or shape annotation
  # if(!showMissing) df <- df[!(is.nan(df$shape_by) & !(is.nan(df$colour_by))]
  
  p <- ggplot(df, aes(x, y, color = colour_by, shape = shape_by)) + 
    geom_point(size = dotsize) + 
    ggtitle(title) + 
    xlab(xlabel) + 
    ylab(ylabel) + 
    scale_y_continuous(limits = c(ylim_down, ylim_up)) + 
    scale_x_continuous(limits = c(xlim_down, xlim_up)) +
    theme(plot.margin = margin(40, 40, 20, 20), 
          axis.text = element_text(size = rel(1), color = "black"), 
          axis.title = element_text(size = titlesize), 
          axis.title.y = element_text(size = rel(1.1), margin = margin(0, 15, 0, 0)), 
          axis.title.x = element_text(size = rel(1.1), margin = margin(15, 0, 0, 0)), 
          axis.line = element_line(colour = "black", size = 0.5), 
          axis.ticks = element_line(colour = "black", size = 0.5),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.key = element_rect(fill = "white"),
          legend.text = element_text(size = titlesize),
          legend.title = element_text(size =titlesize))+
    guides(color=guide_legend(title=name_colour), shape=guide_legend(title=name_shape))+
    scale_shape_manual(values=c(19,1,2:18)[1:length(unique(shape_by))])+
    if(!colorLegend) p <- p + guides(color = FALSE)
  if(!shapeLegend) p <- p + guides(shape = FALSE)
  return(p)
}


