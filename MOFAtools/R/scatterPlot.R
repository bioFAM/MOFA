
###########################################################
## Functions to visualise scatterplots of latent factors ##
###########################################################


#' @title Visualize histogram of one latent variable
#' @name histPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param id name of latent variable to plot
#' @details asd
#' @return fill this
#' @references fill this
#' @import ggplot2
#' @export
histPlot <- function(object, factor, xlabel = NULL, fill=NULL, name_fill="", alpha=0.6, binwidth=NULL, showNA=F) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") 
    stop("'object' has to be an instance of MOFAmodel")
  if ((!is.null(fill_discrete)) & (!is.null(fill_continuous)))
    stop("do not specify both fill_discrete and fill_continuous")
  if(!factor %in% factorNames(object)) { stop("factor not recognised") }
  
  # Collect relevant data
  N <- object@Dimensions[["N"]]
  Z <- getExpectations(object, "Z", "E")
  factor <- as.character(factor)
  
  # Create data frame
  if (is.null(fill)) { 
    fill <- 1 
  } else {
    stopifnot(is.character(fill) | is.factor(fill))
    fill <- as.factor(fill)
  }
  
  df = data.frame(x=Z[,factor], fill=fill)
  
  if(!showNA) df <- df[!is.na(fill),]

  if (is.null(xlabel)) { xlabel <- paste("Latent factor", factor) }
  
  p <- ggplot(df, aes(x=x)) + 
    geom_histogram(aes(fill=fill), alpha=alpha, binwidth=binwidth, position="identity") + 
    xlab(xlabel) + 
    ylab("Count") + 
    scale_y_continuous(expand=c(0,0)) +
    guides(fill=guide_legend(title=name_fill)) +
    theme(plot.margin = margin(40,40,20,20), 
          axis.text = element_text(size=rel(1.3), color = "black"), 
          axis.title.y = element_text(size=rel(1.5), margin=margin(0,15,0,0)), 
          axis.title.x = element_text(size=rel(1.5), margin=margin(15,0,0,0)), 
          axis.line = element_line(colour="black", size=rel(1.0)),
          # axis.ticks = element_line(colour="black", size=0.5),
          panel.border = element_blank(), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(), 
          panel.background = element_blank(),
          legend.key = element_rect(fill = "white")
    )
  
  if (is.null(fill)) { p <- p + guides(fill = FALSE) }
  
  return(p)
}


#' @title Visualize beeswarm plot of one latent variable
#' @name beeswarmPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param id latent variable
#' @details asd
#' @return fill this
#' @references fill this
#' @import beeswarm
#' @export
beeswarmPlot <- function(object, factor, label = NULL, color = NULL, color_name="", legend=F) {
  # TO-DO: BETTER USE GGPLOT AND ALLOW FOR CONTINUOUS AS WELL AS DISCRETE COLORUING
  
  if (class(object) != "MOFAmodel") 
    stop("'object' has to be an instance of MOFAmodel")
  
  N <- object@Dimensions[["N"]]
  Z <- getExpectations(object, "Z", "E")
  
  if (is.null(label)) { label <- paste("Latent factor", factor) }
  if (is.null(color)) { 
    color <- rep("1",N) 
  } 
  beeswarm::beeswarm(Z[,factor], pwcol = color, pch = 16, ylab = label, xlab = "")
  
  if(!is.null(color) & legend==T) {
    legend("topright", legend = levels(col), title = color_name, pch = 16, col = 1:length(levels(color)))
  }
  
}

#' @title Visualize scatterplot of two latent variables
#' @name scatterPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factorx latent variable on the x axis
#' @param factory latent variable on the y axis
#' @details asd
#' @return fill this
#' @references fill this
#' @import ggplot2
#' @export
scatterPlot <- function (object, factors, title = "", titlesize = 16, xlabel = NULL, 
      ylabel = NULL, xlim_down = NA, xlim_up = NA, ylim_down = NA, ylim_up = NA, 
      dotsize = 2.5, colour_by = NULL, shape_by = NULL, name_colour="", name_shape="", showMissing = T) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(length(factors)==2)
  stopifnot(all(factors %in% factorNames(object)))
  
  # Collect relevant data  
  N <- object@Dimensions[["N"]]
  Z <- getExpectations(object, "Z", "E")
  factors <- as.character(factors)
  
  # Set color
  if (!is.null(colour_by)) {
    # It is the name of a covariate 
    if (length(colour_by) == 1 & is.character(colour_by)) {
      colour_by <- as.factor(getCovariates(object, colour_by))
    # It is a vector of length N
    } else if (length(colour_by) > 1) {
      stopifnot(length(colour_by) == N)
      colour_by <- as.factor(colour_by)
    } else {
      stop("'colour_by' was specified but it was not recognised, please read the documentation")
    }
    colorLegend <- T
  } else {
    colour_by <- rep(TRUE,N)
    colorLegend <- F
  }
  
  # Set shape
  if (!is.null(shape_by)) {
    # It is the name of a covariate 
    if (length(shape_by) == 1 & is.character(shape_by)) {
      shape_by <- as.factor(getCovariates(object, shape_by))
    # It is a vector of length N
    } else if (length(shape_by) > 1) {
      stopifnot(length(shape_by) == N)
      shape_by <- as.factor(shape_by)
    } else {
      stop("'shape_by' was specified but it was not recognised, please read the documentation")
    }
    colorLegend <- T
  } else {
    shape_by <- rep(TRUE,N)
    shapeLegend <- F
  }
  
  # Create data frame to plot
  df = data.frame(x = Z[, factors[1]], y = Z[, factors[2]], shape_by = shape_by, colour_by = colour_by)
  
  # remove values missing color or shape annotation
  if (!showMissing) df <- df[!((df$shape_by=="NaN") | (df$colour_by=="NaN"))]
  
  if (is.null(xlabel)) xlabel <- paste("Latent factor", factors[1])
  if (is.null(ylabel)) ylabel <- paste("Latent factor", factors[2])
                                
  p <- ggplot(df, aes(x, y, color = colour_by, shape = shape_by)) + 
      geom_point(size = dotsize) + 
      ggtitle(title) + xlab(xlabel) + ylab(ylabel) + 
      scale_y_continuous(limits = c(ylim_down, ylim_up)) + 
      scale_x_continuous(limits = c(xlim_down, xlim_up)) +
      scale_shape_manual(values=c(19,1,2:18)[1:length(unique(shape_by))]) +
      theme(plot.margin = margin(20, 20, 10, 10), 
            axis.text = element_text(size = rel(1), color = "black"), 
            axis.title = element_text(size = titlesize), 
            axis.title.y = element_text(size = rel(1.1), margin = margin(0, 10, 0, 0)), 
            axis.title.x = element_text(size = rel(1.1), margin = margin(10, 0, 0, 0)), 
            axis.line = element_line(colour = "black", size = 0.5), 
            axis.ticks = element_line(colour = "black", size = 0.5),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            legend.key = element_rect(fill = "white"),
            legend.text = element_text(size = titlesize),
            legend.title = element_text(size =titlesize))
  if (colorLegend) { p <- p + labs(color = name_colour) } else { p <- p + guides(color = FALSE) }
  if (shapeLegend) { p <- p + labs(shape = name_shape) }  else { p <- p + guides(shape = FALSE) }
  return(p)
}
  
  
#' @title Visualize scatterplot of all latent variables in a pair-wise grid
#' @name scatterPairs
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param colour_by vector of length N with discrete sample groups to colour by
#' @param colour_name name of the factor used to colour_by (only relevant if colour_by is not NULL)
#' @param showMissing boolean wether to include missing values in colour_by
#' @details asd
#' @return fill this
#' @references fill this
#' @import ggplot2 GGally
#' @export
#' 
scatterPairs <- function(object, factors = "all", showMissing=T, 
                         colour_by=NULL, colour_name=NULL, colourLegend=F, 
                         shape_by=NULL, shape_name=NULL, shapeLegend=F) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") 
    stop("'object' has to be an instance of MOFAmodel")

  # Collect relevant data
  N <- object@Dimensions[["N"]]
  Z <- getExpectations(object, "Z", "E")
  factors <- as.character(factors)
  
  # Get factors
  if (paste0(factors,collapse="") == "all") { 
    factors <- factorNames(object) 
    # if(is.null(factors)) factors <- 1:ncol(Z) # old object are not compatible with factro names
  } else {
    stopifnot(all(factors %in% factorNames(object)))  
  }
  Z <- Z[,factors]

  # Remove constant factors 
  tmp <- apply(Z,2,var)
  if (any(tmp==0)) {
    message(paste0("Removing constant factors: ", paste(which(tmp==0), collapse="")))
    Z <- Z[,!tmp==0]
    factors <- factors[!tmp==0]
  }
  
  # Set color
  if (!is.null(colour_by)) {
    # It is the name of a covariate 
    if (length(colour_by) == 1 & is.character(colour_by)) {
      colour_by <- as.factor(getCovariates(object, colour_by))
      # It is a vector of length N
    } else if (length(colour_by) > 1) {
      stopifnot(length(colour_by) == N)
      colour_by <- as.factor(colour_by)
    } else {
      stop("'colour_by' was specified but it was not recognised, please read the documentation")
    }
    colorLegend <- T
  } else {
    colour_by <- rep(TRUE,N)
    colorLegend <- F
  }
  
  # Set shape
  if (!is.null(shape_by)) {
    # It is the name of a covariate 
    if (length(shape_by) == 1 & is.character(shape_by)) {
      shape_by <- as.factor(getCovariates(object, shape_by))
      # It is a vector of length N
    } else if (length(shape_by) > 1) {
      stopifnot(length(shape_by) == N)
      shape_by <- as.factor(shape_by)
    } else {
      stop("'shape_by' was specified but it was not recognised, please read the documentation")
    }
    colorLegend <- T
  } else {
    shape_by <- rep(TRUE,N)
    shapeLegend <- F
  }

  # Remove missing values
  if(!showMissing) {
    Z <- Z[!is.na(colour_by),]
    colour_by <- colour_by[!is.na(colour_by)]
  }
  
  # Crete data.frame
  df <- as.data.frame(Z); colnames(df) <- paste0("LF_",colnames(df))
  df <- cbind(df, colour_by=colour_by, shape_by=shape_by)
  
  # Define title and legend of the plot
  main <- "" 
  if (!is.null(colour_name)) { colourLegend <- T }
  if (!is.null(shape_name)) { shapeLegend <- T }
  p <- ggplot(df, aes_string(x=colnames(df)[1], y=colnames(df)[2], color="colour_by", shape="shape_by")) + geom_point()
  if (colourLegend) { p <- p + labs(color=colour_name) } else { p <- p + guides(color = FALSE) }
  if (shapeLegend) { p <- p + labs(shape=shape_name) } else { p <- p + guides(shape = FALSE) }
  if (colourLegend | shapeLegend) { 
    p <- p +
      theme(
        legend.title=element_text(size=15, hjust=0.5, color="black"),
        legend.position = "right", 
        legend.direction = "vertical",
        legend.key = element_blank()
      )
    legend <- GGally::grab_legend(p)
  } else {
    legend <- NULL
  }
  
  # Generate plot
  GGally::ggpairs(df, columns = colnames(df[,!colnames(df) %in% c("colour_by","shape_by")]), 
                  lower=list(continuous="points"), diag=list(continuous='blankDiag'), upper=list(continuous='points'),
          mapping=aes(colour=colour_by, shape=shape_by), title=main, legend=legend) +
    theme_bw() +
    theme(plot.title = element_text(size = 16, hjust=0.5, color="black"), 
          axis.title = element_text(size = 10, color="black"), 
          axis.text = element_text(size = 9, color="black"),
          legend.position = "right", 
          legend.direction = "vertical"
          )
}
  