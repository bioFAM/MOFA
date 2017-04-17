
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
#' 
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
    ylim_up = NA, dotsize = 2.5, colour_by = NULL, shape_by = NULL, name_colour="", name_shape="") 
{
    if (class(object) != "MOFAmodel") 
        stop("'object' has to be an instance of MOFAmodel")
  
    N <- object@Dimensions[["N"]]
    Z <- getExpectations(object, "Z", "E")
    
    if (!is.null(colour_by)) {
        if (length(colour_by) != N) 
            stop("'colour_by' has to be a vector of length N")
        colorLegend<-T
    }
    else {
        colour_by <- rep(TRUE, N)
        colorLegend<-F
    }
    if (!is.null(shape_by)) {
        if (length(shape_by) != N) 
            stop("'shape_by' has to be a vector of length N")
        shapeLegend<-T
    }
    else {
        shape_by <- rep(TRUE, N)
        shapeLegend<-F
    }
    if(is.null(xlabel)) xlabel <- paste("Latent factor", idx)
    if(is.null(ylabel)) ylabel <- paste("Latent factor", idy)

    df = data.frame(x = Z[, idx], y = Z[, idy], shape = shape_by, colour = colour_by)
    
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
                legend.title = element_text(size =titlesize)) + 
         guides(color=guide_legend(title=name_colour), shape=guide_legend(title=name_shape))
    if(!colorLegend) p <- p + guides(color = FALSE)
    if(!shapeLegend) p <- p + guides(shape = FALSE)
    return(p)
}