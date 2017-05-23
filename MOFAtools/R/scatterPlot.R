
#' @title Visualize histogram of one latent variable
#' @name histPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param id latent variable
#' @details asd
#' @return fill this
#' @reference fill this
#' @import ggplot2
#' @export
histPlot <- function(object, id, xlabel = NULL, groups = NULL, alpha=1.0, binwidth=NULL) {
  
  if (class(object) != "MOFAmodel") 
    stop("'object' has to be an instance of MOFAmodel")
  
  N <- object@Dimensions[["N"]]
  Z <- getExpectations(object, "Z", "E")
  
  if (is.null(xlabel)) { xlabel <- paste("Latent factor", id) }
  if (is.null(groups)) { groups <- rep("1",N) }
  
  df = data.frame(x=Z[,id], group=groups)
  
  p <- ggplot(df, aes(x=x)) + 
    geom_histogram(aes(fill=group), alpha=alpha, binwidth=binwidth) + 
    xlab(xlabel) + 
    ylab("Count") + 
    scale_y_continuous(expand=c(0,0)) +
    # guides(fill=guide_legend(title=name_colour))
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
  
  if (length(unique(df$group))==1) { p <- p + guides(fill = FALSE) }
  
  return(p)
}

#' @title Visualize scatterplot of two latent variables
#' @name scatterPlot
#' @description fill this
#' @param object a \code{\link{MOFAmodel}} object.
#' @param idx latent variable on the x axis
#' @param idy latent variable on the y axis
#' @details asd
#' @return fill this
#' @reference fill this
#' @import ggplot2
#' @export
  scatterPlot<-function (object, idx, idy, title = "", titlesize = 16, xlabel = NULL, 
      ylabel = NULL, xlim_down = NA, xlim_up = NA, ylim_down = NA, 
      ylim_up = NA, dotsize = 2.5, colour_by = NULL, shape_by = NULL, name_colour="", name_shape="",
      showMissing = T) 
{
    if (class(object) != "MOFAmodel") 
        stop("'object' has to be an instance of MOFAmodel")
  
    N <- object@Dimensions[["N"]]
    Z <- getExpectations(object, "Z", "E")
    
    if (!is.null(colour_by)) {
        if (length(colour_by) != N) 
            stop("'colour_by' has to be a vector of length N")
        colorLegend<-T
    } else {
        colour_by <- rep(TRUE, N)
        colorLegend<-F
    }
    if (!is.null(shape_by)) {
        if (length(shape_by) != N) 
            stop("'shape_by' has to be a vector of length N")
        shapeLegend<-T
    } else {
        shape_by <- rep(TRUE, N)
        shapeLegend<-F
    }
    if (is.null(xlabel)) xlabel <- paste("Latent factor", idx)
    if (is.null(ylabel)) ylabel <- paste("Latent factor", idy)

    df = data.frame(x = Z[, idx], y = Z[, idy], shape_by = shape_by, colour_by = colour_by)
    
    #remove values missing color or shape annotation
    if(!showMissing) df <- filter(df, !(shape_by=="NaN") & !(colour_by=="NaN"))
                                  
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
  
  
  #' @title Visualize scatterplot of all latent variables in a pair-wise grid
  #' @name scatterPairs
  #' @description fill this
  #' @param object a \code{\link{MOFAmodel}} object.
  #' @param colour_by discrete vector of length N (samples in object) to colour samples by
  #' @param showMissing boolean wether to include NA values of coulour_by
  #' @details asd
  #' @return fill this
  #' @reference fill this
  #' @import ggplot2, GGally
  #' @export
  #' 
  
  scatterPairs <- function (object, factors = "all", colour_by=NULL, showMissing=T, colour_name=NULL){
    
    if (class(object) != "MOFAmodel") 
      stop("'object' has to be an instance of MOFAmodel")
    
    
    N <- object@Dimensions[["N"]]
    Z <- getExpectations(object, "Z", "E")
    
    if (paste0(factors,sep="",collapse="") == "all") { 
      factors <- factorNames(object) 
      #old object are not compatible with factro names
      if(is.null(factors)) factors <- 1:ncol(Z)
    } else {
      stopifnot(all(factors %in% factorNames(object)))  
    }
  
    if(object@ModelOpts$learnMean==T & any(apply(Z[,factors],2,var)==0)) warning("The constant intercept factor is included in the plot")
    
    if (!is.null(colour_by)) {
      if (length(colour_by) != N) 
        stop("'colour_by' has to be a vector of length N")
    } else {
      colour_by <- rep(TRUE, N)
    }
    
    Z2include <- Z[,factors]
    colnames(Z2include) <- paste("LF", colnames(Z2include), sep="_")
  
    if(!showMissing) {
      Z2include <- Z2include[!is.na(colour_by),]
      colour_by <- colour_by[!is.na(colour_by)]
    }
    
    df_pairs <- as.data.frame(Z2include)
    df_pairs<- cbind(df_pairs, colour_by=as.factor(colour_by))
    if(!is.null(colour_name)) {
      main <- paste("Scatterplots of latent factors coloured by", colour_name)
      p <- ggplot(df_pairs, aes(x= colnames(Z2include)[1], y=colnames(Z2include)[2], colour=colour_by)) + guides(color=guide_legend(title=colour_name))+geom_point()
      legend <- grab_legend(p)
    }
    else {
      main <- "Scatterplots of latent factors" 
      legend <- NULL
    }
    ggpairs(df_pairs, columns = colnames(Z2include), lower=list(continuous="points"), upper=list(continuous="density"), 
            mapping= aes(colour=colour_by), title = main, legend =legend)
}
  