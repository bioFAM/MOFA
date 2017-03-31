
##############################
## Plot training statistics ##
##############################

#' @rdname trainCurve
#' @name trainCurve
#' @title Training Statistic
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @param statistic fill this
#' @import ggplot2 scales
#' @export
trainCurve <- function(object, statistic, title="", titlesize=16, xlabel="Time", ylabel="",
                               linesize=1.2, xlim_down=NA, xlim_up=NA, ylim_down=NA, ylim_up=NA) {
  if (class(object) != "MOFAmodel") { stop("'object' has to be an instance of MOFAmodel") }
  if (!statistic %in% c("activeK","elbo")) { stop("'statistic' have to be either 'activeK' or 'elbo'") }
  if (statistic=="activeK" & ylabel=="") { ylabel <- "Number of active latent varaibles" }
  if (statistic=="elbo" & ylabel=="") { ylabel <- "Evidence lower bound (log)" }
  
  if (statistic=="activeK") { 
    idx = seq(1+object@TrainOpts$startdrop,length(object@TrainStats[[statistic]]),object@TrainOpts$freqdrop)
    stat = object@TrainStats$activeK[idx] 
  }
  if (statistic=="elbo") { 
    idx = seq(1+object@TrainOpts$startdrop,length(object@TrainStats[[statistic]]),object@TrainOpts$elbofreq)
    stat = -log(object@TrainStats$activeK[idx] )
  }  
  
  data <- data.frame(time=idx, value=stat)
  p <- ggplot2::ggplot(data, aes(x=time, y=value)) +
    geom_line(size=linesize) +
    ggtitle(title) + xlab(xlabel) + ylab(ylabel) +
    scale_x_continuous(limits=c(xlim_down,xlim_up)) +
    scale_y_continuous(limits=c(ylim_down,ylim_up)) +
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
  
  # if (statistic=="elbo") {
  #   scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x))) }
  #   p <- p + scale_y_continuous(limits=c(ylim_down,ylim_up),labels=scientific_10)
  # }

  return(p)
}

# library(scales)
# number_ticks <- function(n) { function(limits) pretty(limits, n) }
# if (is.null(nticks)) {
#   breaks = waiver()
# } else {
#   breaks <- number_ticks(nticks)
# }

# scale_y_continuous(limits=c(ylim_down,ylim_up), labels=scientific_10, breaks=breaks)
  
