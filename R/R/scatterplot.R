# Scatterplot of the latent variables


scatterPlot <- function(model, idx, idy, title="", titlesize=16, xlabel="", ylabel="", 
                        xlim_down=NA, xlim_up=NA, ylim_down=NA, ylim_up=NA,
                        dotsize=2.5, colour_by=NULL, shape_by=NULL) {
  
  N <- nrow(model@Data[[1]])
  
  if ( !is.null(colour_by) ) {
    if (length(colour_by) != N)
      stop("'colour_by' has to be a vector of length N")
  } else {
    colour_by <- rep(TRUE,N)
  }
  
  if ( !is.null(shape_by) ) {
    if (length(shape_by) != N)
      stop("'shape_by' has to be a vector of length N")
  } else {
    shape_by <- rep(TRUE,N)
  }
  
  df = data.frame(x=model@Expectations$Z$E[,idx], y=model@Expectations$Z$E[,idy], shape=shape_by, colour=colour_by)
  p <- ggplot(df,aes(x, y, color=colour_by, shape=shape_by)) + 
    geom_point(size=dotsize) +
    ggtitle(title) + xlab(xlabel) + ylab(ylabel) +
    scale_y_continuous(limits=c(ylim_down,ylim_up)) +
    scale_x_continuous(limits=c(xlim_down,xlim_up)) +
    theme(
      plot.margin = margin(40,40,20,20),
      axis.text = element_text(size=rel(1.3), color='black'),
      axis.title = element_text(size=titlesize),
      axis.title.y = element_text(size=rel(1.1), margin=margin(0,15,0,0)),
      axis.title.x = element_text(size=rel(1.1), margin=margin(15,0,0,0)),
      axis.line = element_line(colour="black", size=0.5),
      axis.ticks = element_line(colour="black", size=0.5),
      legend.position='none',
      panel.border=element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  return(p)
}

# colour_by <- c(rep(T,N/2),rep(F,N/2))
# shape_by <- c(rep(T,N/2),rep(F,N/2))
# scatterPlot(gfa, 1, 2, title="", titlesize=16, xlabel="Latent variable 1", ylabel="Latent variable 2",
            # dotsize=2.5, colour_by=colour_by, shape_by=shape_by, xlim_down=NULL, xlim_up=NULL, ylim_down=NULL, ylim_up=NULL)
  