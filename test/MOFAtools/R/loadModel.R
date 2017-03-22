
#' @title LoadModel: loading a MOFA model
#' @name loadModel
#' @description: function to load a trained MOFA model from a hdf5 file. \cr
#' Models are saved in hdf5 format using the 'save' function from utils.py
#' @param file an hd5f object with ...
#' @details fill this...
#' @return a \code{\link{MOFAmodel}} object.
#' @export
#' @import rhdf5

loadModel <- function(file) {
  
  # Load expectations and parameters
  expectations <- h5read(file,"expectations")
  parameters <- h5read(file,"parameters")
  
  # Load training statistics
  training_stats <- h5read(file,"training_stats")
  colnames(training_stats$elbo_terms) <- attr(h5read(file,"training_stats/elbo_terms", read.attributes=T),"colnames")
  
  # Load training options
  # training_opts <- as.list(h5read(file,"training_opts", read.attributes=T))
  
  # Load training data
  data <- h5read(file,"data")
  featuredata <- h5read(file,"features")
  sampledata <- h5read(file,"samples")
  for (m in names(data)) {
    rownames(data[[m]]) <- sampledata
    colnames(data[[m]]) <- featuredata[[m]]
  }
  
  # Specify dimensions 
  M=length(data)
  N=nrow(data[[1]])
  D=sapply(data,ncol)
  K=tail(training_stats$activeK,n=1)
  dim=list("M"=M, "N"=N, "D"=D, "K"=K)
  
  # Create object
  mofa <- new("MOFAmodel", 
             TrainData=data, TrainStats=training_stats, 
             Expectations=expectations, Parameters=parameters,
             Dimensions=dim)
  
  return(mofa)
}

