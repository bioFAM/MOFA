
# Number of active factors or latent variables
# Lower bound
# Likelihood

trainCurve <- function(data, title="", titlesize=16, xlabel="", ylabel="", linesize=1.2, xlim_down=NA, xlim_up=NA,
                        ylim_down=NA, ylim_up=NA, group.colors=NA, nticks=5) {
  
  # number_ticks <- function(n) {function(limits) pretty(limits, n)}
  # scientific_10 <- function(x) { parse(text=gsub("e", " %*% 10^", scientific_format()(x))) }
  
  p <- ggplot(data) +
    geom_line(aes(x=iter, y=value, group=interaction(model,trial), color=model), size=linesize) +
    ggtitle(title) + xlab(xlabel) + ylab(ylabel) +
    scale_y_continuous(limits = c(ylim_down, ylim_up), labels = comma, breaks=number_ticks(nticks)) +
    # scale_y_continuous(limits = c(ylim_down, ylim_up), labels=scientific_10) +
    # scale_y_continuous(limits = c(ylim_down, ylim_up), breaks=c(5,6,7,8,9)) +
    scale_x_continuous(limits = c(xlim_down, xlim_up)) +
    scale_colour_manual(values=group.colors) +
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
      legend.position="top",#c(0.8, 0.8),
      legend.title = element_blank(),
      legend.direction="horizontal",
      legend.text = element_text(size=rel(2.5)),
      legend.key = element_rect(fill="white"),
      legend.key.size = unit(1.2, "cm"), 
      legend.title=element_text(size=rel(1.5)),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.background = element_blank()
    )
  print(p)
}