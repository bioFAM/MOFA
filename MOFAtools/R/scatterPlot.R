
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
          axis.line = element_line(color="black", size=rel(1.0)),
          # axis.ticks = element_line(color="black", size=0.5),
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
#' @export
beeswarmPlot <- function(object, factors, color_by = NULL, color_name="", colorLegend=T, showNA=F) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel")  stop("'object' has to be an instance of MOFAmodel")
  factors <- as.character(factors)
  stopifnot(all(factors %in% factorNames(object)))
  
  # Collect relevant data
  N <- object@Dimensions[["N"]]
  Z <- getExpectations(object, "Z", "E", as.data.frame=T)
  Z <- Z[Z$factor %in% factors,]
  
  # Z <- getExpectations(object, "Z", "E")
  
  # Set color
  if (!is.null(color_by)) {
    # It is the name of a covariate 
    if (length(color_by) == 1 & is.character(color_by)) {
      color_by <- getCovariates(object, color_by)
    # It is a vector of length N
    } else if (length(color_by) > 1) {
      stopifnot(length(color_by) == N)
    } else {
      stop("'color_by' was specified but it was not recognised, please read the documentation")
    }
  } else {
    color_by <- rep(TRUE,N)
    colorLegend <- F
  }
  
  # Remove samples with missing values
  if (!showNA) {
    Z <- Z[!is.nan(Z$value),]
  }
  
  
  # Generate plot
  p <- ggplot(Z, aes(factor, value)) + 
    ggbeeswarm::geom_quasirandom(aes(color=color_by)) +
    labs(x="Latent factor(s)", y="") +
    theme(
      # axis.text.x = element_blank(),
      axis.text.x = element_text(size = rel(1.5), color = "black", margin=margin(5,0,0,0)),
      axis.text.y = element_text(size = rel(1.5), color = "black"),
      axis.title.x = element_text(size = rel(1.4), color = "black"),
      axis.line = element_line(color = "black", size = 0.4),
      axis.ticks.length = unit(0.25,"cm"),
      axis.ticks = element_line(color = "black"),
      panel.border = element_blank(), 
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(), 
      panel.background = element_blank(),
      legend.title=element_text(size=15, hjust=0.5, color="black"),
      legend.position = "right", 
      legend.direction = "vertical",
      legend.key = element_blank()
      )
  
  # If color_by is numeric, define the default gradient
  if (is.numeric(color_by)) { p <- p + scale_color_gradientn(colors=terrain.colors(10)) }
  
  # Add legend
  if (colorLegend) { p <- p + labs(color=color_name) } else { p <- p + guides(color = FALSE) }
  
  return(p)
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
      dotsize = 2.5, color_by = NULL, shape_by = NULL, name_color="", name_shape="", showMissing = T) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(length(factors)==2)
  stopifnot(all(factors %in% factorNames(object)))
  
  # Collect relevant data  
  N <- object@Dimensions[["N"]]
  Z <- getExpectations(object, "Z", "E")
  factors <- as.character(factors)
  
  # Set color
  if (!is.null(color_by)) {
    # It is the name of a covariate 
    if (length(color_by) == 1 & is.character(color_by)) {
      color_by <- as.factor(getCovariates(object, color_by))
    # It is a vector of length N
    } else if (length(color_by) > 1) {
      stopifnot(length(color_by) == N)
      # color_by <- as.factor(color_by)
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
    shapeLegend <- T
  } else {
    shape_by <- rep(TRUE,N)
    shapeLegend <- F
  }
  
  # Create data frame to plot
  df = data.frame(x = Z[, factors[1]], y = Z[, factors[2]], shape_by = shape_by, color_by = color_by)
  
  # remove values missing color or shape annotation
  if (!showMissing) df <- df[!((df$color_by=="NaN") | (df$shape_by=="NaN")),]
  
  if (is.null(xlabel)) xlabel <- paste("Latent factor", factors[1])
  if (is.null(ylabel)) ylabel <- paste("Latent factor", factors[2])
                                
  p <- ggplot(df, aes(x, y, color = color_by, shape = shape_by)) + 
      geom_point(size = dotsize) + 
      ggtitle(title) + xlab(xlabel) + ylab(ylabel) + 
      scale_y_continuous(limits = c(ylim_down, ylim_up)) + 
      scale_x_continuous(limits = c(xlim_down, xlim_up)) +
      # scale_shape_manual(values=c(19,1,2:18)[1:length(unique(shape_by))]) +
      theme(plot.margin = margin(20, 20, 10, 10), 
            axis.text = element_text(size = rel(1), color = "black"), 
            axis.title = element_text(size = titlesize), 
            axis.title.y = element_text(size = rel(1.1), margin = margin(0, 10, 0, 0)), 
            axis.title.x = element_text(size = rel(1.1), margin = margin(10, 0, 0, 0)), 
            axis.line = element_line(color = "black", size = 0.5), 
            axis.ticks = element_line(color = "black", size = 0.5),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            legend.key = element_rect(fill = "white"),
            legend.text = element_text(size = titlesize),
            legend.title = element_text(size =titlesize)
            )
  if (colorLegend) { p <- p + labs(color = name_color) } else { p <- p + guides(color = FALSE) }
  if (shapeLegend) { p <- p + labs(shape = name_shape) }  else { p <- p + guides(shape = FALSE) }
  return(p)
}
  
  
#' @title Visualize scatterplot of all latent variables in a pair-wise grid
#' @name scatterPairs
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param color_by vector of length N with discrete sample groups to color by
#' @param color_name name of the factor used to color_by (only relevant if color_by is not NULL)
#' @param showMissing boolean wether to include missing values in color_by
#' @details asd
#' @return fill this
#' @references fill this
#' @import ggplot2 GGally
#' @export
#' 
scatterPairs <- function(object, factors = "all", showMissing=T, 
                         color_by=NULL, color_name=NULL, colorLegend=F, 
                         shape_by=NULL, shape_name=NULL, shapeLegend=F) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")

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
  tmp <- apply(Z,2,var,na.rm=T)
  if (any(tmp==0)) {
    message(paste0("Removing constant factors: ", paste(which(tmp==0), collapse="")))
    Z <- Z[,!tmp==0]
    factors <- factors[!tmp==0]
  }
  
  # Set color
  if (!is.null(color_by)) {
    # It is the name of a covariate 
    if (length(color_by) == 1 & is.character(color_by)) {
      color_by <- getCovariates(object, color_by)
      # It is a vector of length N
    } else if (length(color_by) > 1) {
      stopifnot(length(color_by) == N)
      # color_by <- as.factor(color_by)
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
    shapeLegend <- T
  } else {
    shape_by <- rep(TRUE,N)
    shapeLegend <- F
  }

  # Remove missing values
  if(!showMissing) {
    Z <- Z[!is.na(color_by),]
    color_by <- color_by[!is.na(color_by)]
    shape_by <- shape_by[!is.na(shape_by)]
  }
  
  # Crete data.frame
  df <- as.data.frame(Z); colnames(df) <- paste0("LF_",colnames(df))
  df <- cbind(df, color_by=color_by, shape_by=shape_by)
  
  # Define title and legend of the plot
  main <- "" 
  if (!is.null(color_name)) { colorLegend <- T }
  if (!is.null(shape_name)) { shapeLegend <- T }
  p <- ggplot(df, aes_string(x=colnames(df)[1], y=colnames(df)[2], color="color_by", shape="shape_by")) + geom_point()
  if (colorLegend) { p <- p + labs(color=color_name) } else { p <- p + guides(color = FALSE) }
  if (shapeLegend) { p <- p + labs(shape=shape_name) } else { p <- p + guides(shape = FALSE) }
  if (colorLegend | shapeLegend) { 
    p <- p +
      theme(
        legend.title=element_text(size=15, hjust=0.5, color="black"),
        legend.position = "right", 
        legend.direction = "vertical",
        legend.key = element_blank()
      )
    
    # If color_by is numeric, define the default gradient
    if (is.numeric(color_by)) { p <- p + scale_color_gradientn(colors=terrain.colors(10)) }
    
    # Extract the legend
    legend <- GGally::grab_legend(p)
  } else {
    legend <- NULL
  }
  
  # Generate plot
  p <- GGally::ggpairs(df, columns = colnames(df[,!colnames(df) %in% c("color_by","shape_by")]), 
                  lower=list(continuous="points"), diag=list(continuous='blankDiag'), upper=list(continuous='points'),
          mapping=aes(color=color_by, shape=shape_by), title=main, legend=legend) +
    theme_bw() +
    theme(plot.title = element_text(size = 16, hjust=0.5, color="black"), 
          axis.title = element_text(size = 10, color="black"), 
          axis.text = element_text(size = 9, color="black"),
          legend.position = "right", 
          legend.direction = "vertical"
          )
  
  # If color_by is numeric, define the default gradient
  if (is.numeric(color_by)) { 
    for(i in 1:p$nrow) {
      for(j in 1:p$ncol){
        p[i,j] <- p[i,j] + scale_color_gradientn(colors=terrain.colors(10)) 
      }
    }
  }
  
  return(p)
}
  