
###########################################
## Functions to plot training statistics ##
###########################################

#' @rdname trainCurveFactors
#' @name trainCurveFactors
#' @title Training Curve for number of active factors
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @import ggplot2 scales
#' @export

trainCurveFactors <- function(object, title="", titlesize=16, xlabel="Iteration", ylabel=NULL,
                       linesize=1.2, xlim_down=NULL, xlim_up=NULL) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") { stop("'object' has to be an instance of MOFAmodel") }
  
  if (is.null(ylabel)) { ylabel <- "Number of active latent varaibles" }
  
  # Collect training statistics
  idx = seq(1+object@TrainOpts$startdrop,length(object@TrainStats$activeK),object@TrainOpts$freqdrop)
  stat = object@TrainStats$activeK[idx] 
  data <- data.frame(time=idx, value=stat)
  
  # Plot
  p <- ggplot2::ggplot(data, aes(x=time, y=value)) +
    geom_line(size=linesize) +
    ggtitle(title) + xlab(xlabel) + ylab(ylabel) +
    scale_x_continuous(limits=c(xlim_down,xlim_up)) +
    # scale_y_discrete(limits=c(min(data$value)-1, max(data$value)+1)) +
    scale_y_continuous(limits=c(min(data$value)-1, max(data$value)+1), breaks=scales::pretty_breaks()) +
    theme(
      plot.title = element_text(size=titlesize),
      plot.margin = margin(10,10,10,10),
      axis.title.x=element_text(colour="black",size=rel(1.75), margin=margin(20,0,3,0)),
      axis.title.y=element_text(colour="black",size=rel(1.75), margin=margin(0,20,0,3)),
      axis.text.x=element_text(colour="black",size=rel(1.5)),
      axis.text.y=element_text(colour="black",size=rel(1.5)),
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


#' @rdname trainCurveELBO
#' @name trainCurveELBO
#' @title Training curve for Evidence Lower Bound (ELBO)
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @param statistic fill this
#' @import ggplot2 scales
#' @export

trainCurveELBO <- function(object, log=F, title="", titlesize=16, xlabel="Iteration", ylabel=NULL,
                       linesize=1.2, xlim_down=NULL, xlim_up=NULL) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") { stop("'object' has to be an instance of MOFAmodel") }
  if (is.null(ylabel)) { ylabel <- "Evidence lower bound" }
  
  idx = seq(1,length(object@TrainStats$elbo),object@TrainOpts$elbofreq)
  stat = object@TrainStats$elbo[idx]
  if (log) { stat <- -log(-stat); ylabel <- paste(ylabel,"(log)") }
  data <- data.frame(time=idx, value=stat)
  
  if (!is.null(xlim_down)) { data <- data[data$time>=xlim_down,] }
  if (!is.null(xlim_up)) { data <- data[data$time<=xlim_up,] }
  
  p <- ggplot2::ggplot(data, aes(x=time, y=value)) +
    geom_line(size=linesize) +
    ggtitle(title) + xlab(xlabel) + ylab(ylabel) +
    scale_x_continuous(limits=c(xlim_down,xlim_up)) +
    # scale_y_continuous(limits=c(ylim_down,ylim_up), breaks=scales::pretty_breaks()) +
    # scale_y_continuous(breaks=scales::pretty_breaks()) +
    theme(
      plot.title = element_text(size=titlesize),
      plot.margin = margin(10,10,10,10),
      axis.title.x=element_text(colour="black",size=rel(1.75), margin=margin(20,0,3,0)),
      axis.title.y=element_text(colour="black",size=rel(1.75), margin=margin(0,20,0,3)),
      axis.text.x=element_text(colour="black",size=rel(1.5)),
      axis.text.y=element_text(colour="black",size=rel(1.5)),
      axis.ticks.x = element_line(colour="black"),
      axis.ticks.y = element_line(colour="black"),
      axis.line.x = element_line(color="black"),
      axis.line.y = element_line(color="black"),
      legend.position="none",
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  
  if (log==F) {
      scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x))) }
      p <- p + scale_y_continuous(labels=scientific_10)
  }
  
  return(p)
}

