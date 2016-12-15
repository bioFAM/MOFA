
##############################
## Plot training statistics ##
##############################

# Line plot of the number of active factors
# Line plot of the lower bound
# Bar plot with the different components of the lower bound 

library(ggplot2)
library(tidyr)
library(scales)

plotFn <- function(data, title="", titlesize=16, xlabel="", ylabel="", linesize=1.2, xlim_down=NA, xlim_up=NA,
                   ylim_down=NA, ylim_up=NA, scientific=F) {
  
  colnames(data) <- c("time","value")
  
  # if (scientific==T) {
  #   ylabels <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) }
  # } else {
  #   ylabels <- waiver()
  # }
  
  # if (is.null(x_nticks)) { 
  #   xbreaks = waiver() 
  # } else {
  #   stopifnot(x_nticks <= length(data$time))
  #   xbreaks <- number_ticks(x_nticks)
  # }
  
  # if (is.null(y_nticks)) { 
  #   ybreaks = waiver() 
  # } else {
  #   ybreaks <- number_ticks(y_nticks)
  # }
  
  p <- ggplot(data, aes(x=time, y=value)) +
    geom_line(size=linesize) +
    ggtitle(title) + xlab(xlabel) + ylab(ylabel) +
    theme(
      plot.title = element_text(size=titlesize),
      plot.margin = margin(10,10,10,10),
      axis.title.x=element_text(colour="black",size=rel(2.2), margin=margin(20,0,3,0)),
      axis.title.y=element_text(colour="black",size=rel(2.2), margin=margin(0,20,0,3)),
      axis.text.x=element_text(colour="black",size=rel(2.0)),
      axis.text.y=element_text(colour="black",size=rel(2.0)),
      axis.ticks.x = element_line(colour="black"),
      axis.ticks.y = element_line(colour="black"),
      axis.line.x = element_line(color="black"),
      axis.line.y = element_line(color="black"),
      legend.position="none",
      # legend.title = element_blank(),
      # legend.direction="horizontal",
      # legend.text = element_text(size=rel(2.5)),
      # legend.key = element_rect(fill="white"),
      # legend.key.size = unit(1.2, "cm"),
      # legend.title=element_text(size=rel(1.5)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  return(p)
}


activeK_trainCurve <- function(model, title="", titlesize=16, xlabel="Time", ylabel="Number of active latent variables",
                               linesize=1.2, xlim_down=NA, xlim_up=NA, ylim_down=NA, ylim_up=NA) {


  df <- data.frame(time=1:length(model@TrainStats$activeK), activeK=model@TrainStats$activeK)
  p <- plotFn(data=df, title=title, titlesize=titlesize, xlabel=xlabel, ylabel=ylabel, linesize=linesize,
              xlim_down=xlim_down, xlim_up=xlim_up, ylim_down=ylim_down, ylim_up=ylim_up)
  p <- p +
    scale_x_continuous(limits=c(xlim_down,xlim_up)) +
    scale_y_continuous(limits=c(ylim_down,ylim_up))
  
  return(p)
}


elbo_trainCurve <- function(model, title="", titlesize=16, xlabel="Time", ylabel="Evidence lower bound (ELBO)", 
                            linesize=1.2, xlim_down=NA, xlim_up=NA, ylim_down=NA, ylim_up=NA, nticks=NULL) {
  
  df <- data.frame(time=1:length(model@TrainStats$elbo), elbo=model@TrainStats$elbo)

  scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) }
  
  number_ticks <- function(n) { function(limits) pretty(limits, n) }
  if (is.null(nticks)) {
    breaks = waiver()
  } else {
    breaks <- number_ticks(nticks)
  }

  p <- plotFn(data=df, title=title, titlesize=titlesize, xlabel=xlabel, ylabel=ylabel, linesize=linesize,
              xlim_down=xlim_down, xlim_up=xlim_up, ylim_down=ylim_down, ylim_up=ylim_up)
  p <- p +
      scale_x_continuous(limits=c(xlim_down,xlim_up)) +
      scale_y_continuous(limits=c(ylim_down,ylim_up), labels=scientific_10, breaks=breaks)
  
  return(p)
}
