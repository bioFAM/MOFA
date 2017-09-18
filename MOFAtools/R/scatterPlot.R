
###########################################################
## Functions to visualise scatterplots of latent factors ##
###########################################################


#' @title Visualize histogram of one latent variable
#' @name histPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factor one factors to plot
#' @param groups specifies groups or values used to split histogram into groups and color. This can be either a character giving the name of a feature or covariate or a vector of same length as number of samples specifying a group or value for each sample.
#' @param name_groups name for groups (usually only used if groups is not a character itself)
#' @param alpha transparency parameter
#' @param binwidth binwidth for histogram
#' @param showMissing boolean, if false, removes sample for which shape_by or color_by is missing
#' @details asd
#' @return fill this
#' @references fill this
#' @import ggplot2
#' @export
histPlot <- function(object, factor, groups=NULL, name_groups="", alpha=0.6, binwidth=NULL, showMissing=FALSE) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") 
    stop("'object' has to be an instance of MOFAmodel")

  if(!factor %in% factorNames(object)) { stop("factor not recognised") }
  
  # Collect relevant data
  N <- object@Dimensions[["N"]]
  Z <- getFactors(object, factors = factor, as.data.frame = TRUE)
  
  # get groups
  groupLegend <- T
  if (!is.null(groups)) {
    # It is the name of a covariate or a feature in the TrainData
    if (length(groups) == 1 & is.character(groups)) {
      if(name_groups=="") name_groups <- groups
      TrainData <- getTrainData(object)
      featureNames <- lapply(TrainData(object), rownames)
      if(groups %in% Reduce(union,featureNames)) {
        viewidx <- which(sapply(featureNames, function(vnm) groups %in% vnm))
        groups <- TrainData[[viewidx]][groups,]
      } else if(class(object@InputData) == "MultiAssayExperiment"){
        groups <- getCovariates(object, groups)
      }
      else stop("'groups' was specified but it was not recognised, please read the documentation")
      # It is a vector of length N
    } else if (length(groups) > 1) {
      stopifnot(length(groups) == N)
    } else {
      stop("'groups' was specified but it was not recognised, please read the documentation")
    }
  } else {
    groups <- rep(TRUE,N)
    groupLegend <- F
  }
  names(groups) <- sampleNames(object)
  Z$groups <- groups[Z$sample]

  if(!showMissing) Z <- Z[!is.na(Z$groups),]
  Z$groups <- as.factor(Z$groups)
  
  xlabel <- paste("Latent factor", factor)
  
  p <- ggplot(Z, aes(x=value, group=groups)) + 
    geom_histogram(aes(fill=groups), alpha=alpha, binwidth=binwidth, position="identity") + 
    xlab(xlabel) + 
    ylab("Count") + 
    scale_y_continuous(expand=c(0,0)) +
    guides(fill=guide_legend(title=name_groups)) +
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
  
  if (!groupLegend) { p <- p + guides(fill = FALSE) }
  
  return(p)
}


#' @title Visualize beeswarm plot of one latent variable
#' @name beeswarmPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factors factors to plot
#' @param color_by specifies groups or values used to color points. This can be either a character giving the name of a feature or covariate or a vector of same length as number of samples specifying a group or value for each sample.
#' @param name_color name for color legend (usually only used if color_by is not a character itself)
#' @param showMissing boolean, if false, removes sample for which shape_by or color_by is missing
#' @details asd
#' @return ggplot object containing the beeswarm plots
#' @references fill this
#' @import ggplot2
#' @export
beeswarmPlot <- function(object, factors, color_by = NULL, name_color="", showMissing=TRUE) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel")  stop("'object' has to be an instance of MOFAmodel")

  # Collect relevant data
  N <- object@Dimensions[["N"]]
  Z <- getFactors(object, factors=factors, include_intercept=FALSE, as.data.frame=T)
  Z$factor <- as.factor(Z$factor)
  
  # Z <- getExpectations(object, "Z", "E")
  
 # Set color
  colorLegend <- T
  if (!is.null(color_by)) {
    # It is the name of a covariate or a feature in the TrainData
    if (length(color_by) == 1 & is.character(color_by)) {
      if(name_color=="") name_color <- color_by
      TrainData <- getTrainData(object)
      featureNames <- lapply(TrainData(object), rownames)
      if(color_by %in% Reduce(union,featureNames)) {
        viewidx <- which(sapply(featureNames, function(vnm) color_by %in% vnm))
        color_by <- TrainData[[viewidx]][color_by,]
      } else if(class(object@InputData) == "MultiAssayExperiment"){
        color_by <- getCovariates(object, color_by)
    }
    else stop("'color_by' was specified but it was not recognised, please read the documentation")
    # It is a vector of length N
    } else if (length(color_by) > 1) {
      stopifnot(length(color_by) == N)
      # color_by <- as.factor(color_by)
    } else {
      stop("'color_by' was specified but it was not recognised, please read the documentation")
    }
  } else {
    color_by <- rep(TRUE,N)
    colorLegend <- F
  }
  names(color_by) <- sampleNames(object)
  if(length(unique(color_by)) < 5) color_by <- as.factor(color_by)

  # Remove samples with missing values
  if (!showMissing) {
    Z <- Z[!is.nan(Z$value),]
  }
  Z$color_by <- color_by[Z$sample]
  
  # Generate plot
  p <- ggplot(Z, aes(x=factor, y=value)) + 
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
      ) + facet_wrap(~factor)
  
  # If color_by is numeric, define the default gradient
  if (is.numeric(color_by)) { p <- p + scale_color_gradientn(colors=terrain.colors(10)) }
  
  # Add legend
  if (colorLegend) { p <- p + labs(color=name_color) } else { p <- p + guides(color = FALSE) }
  
  return(p)
}

#' @title Visualize scatterplot of two latent variables
#' @name scatterPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factors vector of two factor to plot
#' @param color_by specifies groups or values used to color points. This can be either a character giving the name of a feature or covariate or a vector of same length as number of samples specifying a group or value for each sample.
#' @param shape_by specifies groups or values used for point shapes. This can be either a character giving the name of a feature or covariate or a vector of same length as number of samples specifying a group or value for each sample.
#' @param name_color name for color legend (usually only used if color_by is not a character itself)
#' @param name_shape name for shape legend (usually only used if shape_by is not a character itself)
#' @param showMissing boolean, if false, removes sample for which shape_by or color_by is missing
#' @details asd
#' @return ggplot object containing the scatterplot
#' @references fill this
#' @import ggplot2
#' @export
scatterPlot <- function (object, factors, color_by = NULL, shape_by = NULL, name_color="",
                         name_shape="", showMissing = TRUE) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(length(factors)==2)
  stopifnot(all(factors %in% factorNames(object)))
  
  # Collect relevant data  
  N <- object@Dimensions[["N"]]
  Z <- getExpectations(object, "Z", "E")
  factors <- as.character(factors)
  
  # Set color
  colorLegend <- T
  if (!is.null(color_by)) {
    # It is the name of a covariate or a feature in the TrainData
    if (length(color_by) == 1 & is.character(color_by)) {
      if(name_color=="") name_color <- color_by
      TrainData <- getTrainData(object)
      featureNames <- lapply(TrainData(object), rownames)
      if(color_by %in% Reduce(union,featureNames)) {
        viewidx <- which(sapply(featureNames, function(vnm) color_by %in% vnm))
        color_by <- TrainData[[viewidx]][color_by,]
      } else if(class(object@InputData) == "MultiAssayExperiment"){
        color_by <- getCovariates(object, color_by)
    }
    else stop("'color_by' was specified but it was not recognised, please read the documentation")
    # It is a vector of length N
    } else if (length(color_by) > 1) {
      stopifnot(length(color_by) == N)
      # color_by <- as.factor(color_by)
    } else {
      stop("'color_by' was specified but it was not recognised, please read the documentation")
    }
  } else {
    color_by <- rep(TRUE,N)
    colorLegend <- F
  }

  # Set shape
  shapeLegend <- T
  if (!is.null(shape_by)) {
    # It is the name of a covariate 
    if (length(shape_by) == 1 & is.character(shape_by)) {
      if(name_shape=="") name_shape <- shape_by
      TrainData <- getTrainData(object)
      featureNames <- lapply(TrainData(object), rownames)
      if(shape_by %in% Reduce(union,featureNames)) {
        viewidx <- which(sapply(featureNames, function(vnm) shape_by %in% vnm))
        shape_by <- TrainData[[viewidx]][shape_by,]
      } else if(class(object@InputData) == "MultiAssayExperiment"){
        shape_by <- getCovariates(object, shape_by)
    }
    else stop("'shape_by' was specified but it was not recognised, please read the documentation")
    # It is a vector of length N
    # It is a vector of length N
    } else if (length(shape_by) > 1) {
      stopifnot(length(shape_by) == N)
    } else {
      stop("'shape_by' was specified but it was not recognised, please read the documentation")
    }
  } else {
    shape_by <- rep(TRUE,N)
    shapeLegend <- F
  }
  
  # Create data frame to plot
  df = data.frame(x = Z[, factors[1]], y = Z[, factors[2]], shape_by = shape_by, color_by = color_by)
  
  # remove values missing color or shape annotation
  if (!showMissing) df <- df[!is.na(df$shape_by) & !is.na(df$color_by),]

   #turn into factors
   df$shape_by[is.na(df$shape_by)] <- "NA"
   df$shape_by <- as.factor(df$shape_by)
   if(length(unique(df$color_by)) < 5) df$color_by <- as.factor(df$color_by)
 
  
  xlabel <- paste("Latent factor", factors[1])
  ylabel <- paste("Latent factor", factors[2])
                                
  p <- ggplot(df, aes(x, y, color = color_by, shape = shape_by)) + 
      geom_point() + xlab(xlabel) + ylab(ylabel) +
      # scale_shape_manual(values=c(19,1,2:18)[1:length(unique(shape_by))]) +
      theme(plot.margin = margin(20, 20, 10, 10), 
            axis.text = element_text(size = rel(1), color = "black"), 
            axis.title = element_text(size = 16), 
            axis.title.y = element_text(size = rel(1.1), margin = margin(0, 10, 0, 0)), 
            axis.title.x = element_text(size = rel(1.1), margin = margin(10, 0, 0, 0)), 
            axis.line = element_line(color = "black", size = 0.5), 
            axis.ticks = element_line(color = "black", size = 0.5),
            panel.border = element_blank(), 
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(), 
            panel.background = element_blank(),
            legend.key = element_rect(fill = "white"),
            legend.text = element_text(size = 16),
            legend.title = element_text(size =16)
            )
  if (colorLegend) { p <- p + labs(color = name_color) } else { p <- p + guides(color = FALSE) }
  if (shapeLegend) { p <- p + labs(shape = name_shape) }  else { p <- p + guides(shape = FALSE) }
  return(p)
}
  
  
#' @title Visualize scatterplot of all latent variables in a pair-wise grid
#' @name scatterPairs
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factors vector of factors to plot or "all"
#' @param color_by specifies groups or values used to color points. This can be either a character giving the name of a feature or covariate or a vector of same length as number of samples specifying a group or value for each sample.
#' @param shape_by specifies groups or values used for point shapes. This can be either a character giving the name of a feature or covariate or a vector of same length as number of samples specifying a group or value for each sample.
#' @param name_color name for color legend (usually only used if color_by is not a character itself)
#' @param name_shape name for shape legend (usually only used if shape_by is not a character itself)
#' @param showMissing boolean, if false, removes sample for which shape_by or color_by is missing
#' @details asd
#' @return fill this
#' @references fill this
#' @import ggplot2 GGally
#' @export
#' 
scatterPairs <- function(object, factors = "all", showMissing=TRUE, 
                         color_by=NULL, name_color="",  
                         shape_by=NULL, name_shape="") {
  
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
  colorLegend <- T
  if (!is.null(color_by)) {
    # It is the name of a covariate or a feature in the TrainData
    if (length(color_by) == 1 & is.character(color_by)) {
      if(name_color=="") name_color <- color_by
      TrainData <- getTrainData(object)
      featureNames <- lapply(TrainData(object), rownames)
      if(color_by %in% Reduce(union,featureNames)) {
        viewidx <- which(sapply(featureNames, function(vnm) color_by %in% vnm))
        color_by <- TrainData[[viewidx]][color_by,]
      } else if(class(object@InputData) == "MultiAssayExperiment"){
        color_by <- getCovariates(object, color_by)
    }
    else stop("'color_by' was specified but it was not recognised, please read the documentation")
    # It is a vector of length N
    } else if (length(color_by) > 1) {
      stopifnot(length(color_by) == N)
      # color_by <- as.factor(color_by)
    } else {
      stop("'color_by' was specified but it was not recognised, please read the documentation")
    }
  } else {
    color_by <- rep(TRUE,N)
    colorLegend <- F
  }

  # Set shape
  shapeLegend <- T
  if (!is.null(shape_by)) {
    # It is the name of a covariate 
    if (length(shape_by) == 1 & is.character(shape_by)) {
      if(name_shape=="") name_shape <- shape_by
      TrainData <- getTrainData(object)
      featureNames <- lapply(TrainData(object), rownames)
      if(shape_by %in% Reduce(union,featureNames)) {
        viewidx <- which(sapply(featureNames, function(vnm) shape_by %in% vnm))
        shape_by <- TrainData[[viewidx]][shape_by,]
      } else if(class(object@InputData) == "MultiAssayExperiment"){
        shape_by <- getCovariates(object, shape_by)
    }
    else stop("'shape_by' was specified but it was not recognised, please read the documentation")
    # It is a vector of length N
    # It is a vector of length N
    } else if (length(shape_by) > 1) {
      stopifnot(length(shape_by) == N)
    } else {
      stop("'shape_by' was specified but it was not recognised, please read the documentation")
    }
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

    #turn into factors
   df$shape_by[is.na(df$shape_by)] <- "NA"
   df$shape_by <- as.factor(df$shape_by)
   if(length(unique(df$color_by)) < 5) df$color_by <- as.factor(df$color_by)
  
  
  # Define title and legend of the plot
  main <- "" 
  p <- ggplot(df, aes_string(x=colnames(df)[1], y=colnames(df)[2], color="color_by", shape="shape_by")) + geom_point()
  if (colorLegend | shapeLegend) { 
    p <- p +
      theme(
        legend.title=element_text(size=15, hjust=0.5, color="black"),
        legend.position = "right", 
        legend.direction = "vertical",
        legend.key = element_blank()
      )
    
    # If color_by is numeric, define the default gradient
    if (is.numeric(df$color_by)) { p <- p + scale_color_gradientn(colors=terrain.colors(10)) }
    
    if (colorLegend) { p <- p + labs(color = name_color) } else { p <- p + guides(color = FALSE) }
    if (shapeLegend) { p <- p + labs(shape = name_shape) }  else { p <- p + guides(shape = FALSE) }
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
  if (is.numeric(df$color_by)) { 
    for(i in 1:p$nrow) {
      for(j in 1:p$ncol){
        p[i,j] <- p[i,j] + scale_color_gradientn(colors=terrain.colors(10)) 
      }
    }
  }
  
  return(p)
}
  