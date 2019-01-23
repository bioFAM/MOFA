
###########################################
## Functions to plot training statistics ##
###########################################

#' @rdname trainCurveFactors
#' @name trainCurveFactors
#' @title Training curve for the number of active factors
#' @description the MOFA model starts with an initial number of factors and inactive factors
#'  can be dropped during training if they explain small amounts of variation 
#'  (as defined in \code{\link{getDefaultModelOptions}}. 
#' This allows the model to automatically infer the dimensionality of the latent space.
#' The corresponding hyperparameters are defined in \code{\link{prepareMOFA}}. \cr
#' All training statistics, including the number of active factors, can be fetch from the
#'  TrainStats slot of \code{\link{MOFAmodel}} .
#' @param object a \code{\link{MOFAmodel}} object.
#' @return plot of number of active factors during training
#' @import ggplot2 scales
#' @export
#' @examples 
#' # Example on the CLL data
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' MOFA_CLL <- loadModel(filepath)
#' trainCurveFactors(MOFA_CLL)

#' # Example on the scMT data
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAtools")
#' MOFA_scMT <- loadModel(filepath)
#' trainCurveFactors(MOFA_scMT)

trainCurveFactors <- function(object) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") { stop("'object' has to be an instance of MOFAmodel") }
  
  # Collect training statistics
  if(is.null(object@TrainOptions$freqdrop)) {
      freqdrop <- 1
  } else {
      freqdrop <- object@TrainOptions$freqdrop
  }
  idx = seq(1, length(object@TrainStats$activeK), freqdrop)
  stat = object@TrainStats$activeK[idx] 
  data <- data.frame(time=idx, value=stat)
  
  # Plot
  p <- ggplot(data, aes_string(x="time", y="value")) +
    geom_line() +
    labs(title="", x="Iteration", y="Number of active latent varaibles")
    # scale_y_discrete(limits=c(min(data$value)-1, max(data$value)+1)) +
    # scale_y_continuous(limits=c(min(data$value)-1, max(data$value)+1), breaks=scales::pretty_breaks()) +
    theme(
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
    
  return(p)
}


#' @title Training curve for Evidence Lower Bound (ELBO)
#' @name trainCurveELBO
#' @rdname trainCurveELBO
#' @param object a \code{\link{MOFAmodel}} object.
#' @param logScale boolean indicating whether to apply log transform
#' @return plot of ELBO values during training
#' @description MOFA inference is done using the variational Bayes algorithm,
#'  which maximises a quantity called the Evidence Lower Bound (ELBO).
#' The ELBO is supposed to increase monotonically up to convergence,
#'  but it can decrease substantially when dropping inactive factors.
#' For more details read the supplementary methods
#' The frequency of ELBO computation as well as the convergence criteria are defined
#'  as hyperparameters in \code{\link{prepareMOFA}}. \cr
#' All Training statistics, including the ELBO,
#'  can be fetch from the TrainStats slot of \code{\link{MOFAmodel}} .
#' @import ggplot2 scales
#' @export
#' @examples
#' # Example on the CLL data
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' MOFA_CLL <- loadModel(filepath)
#' trainCurveELBO(MOFA_CLL)
#' trainCurveELBO(MOFA_CLL, logScale= TRUE)
#'
#' # Example on the scMT data
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAtools")
#' MOFA_scMT <- loadModel(filepath)
#' trainCurveELBO(MOFA_scMT)

trainCurveELBO <- function(object, logScale = FALSE) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") { 
    stop("'object' has to be an instance of MOFAmodel") 
    }
  if (object@Status != "trained") { 
    stop("MOFAmodel is untrained. Use runMOFA to train it.") 
  }
  # Fetch ELBO from TrainStats  
  if(is.null(object@TrainOptions$elbofreq)) {
    elbofreq <- 1
    } else {
      elbofreq <- object@TrainOptions$elbofreq
  }
  idx = seq(1,length(object@TrainStats$elbo), elbofreq)
  stat = object@TrainStats$elbo[idx]
  
  # Apply log transform
  if (logScale==TRUE) { stat <- -log(-stat) }
  
  data <- data.frame(time=idx, value=stat)
  p <- ggplot2::ggplot(data, aes_string(x="time", y="value")) +
    geom_line() +
    labs(title="", x="Iteration", y="ELBO")
    theme(
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
  
  if (logScale==FALSE) {
      scientific_10 <- function(x) {
        parse(text=gsub("e", " %*% 10^", scales::scientific_format()(x))) 
        }
      p <- p + scale_y_continuous(labels=scientific_10)
  }
  
  return(p)
}

