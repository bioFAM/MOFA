
#' @title LoadModel: loading a trained GFA model
#' @name loadModel
#' @rdname loadModel
#' @description: function to load a trained GFA model from a hdf5 file. \cr
#' Models are saved in hdf5 format using the 'save' function from utils.py
#' @param file an hd5f object with ...
#' @details asd
#' @return a \code{\link{GFATrainedModel}} object.
#' @export
#' @import rhdf5
loadModel <- function(file) {
  
  # Load expectations and parameters
  expectations <- rhdf5::h5read(file,"expectations")
  parameters <- rhdf5::h5read(file,"parameters")
  
  # Load training statistics
  training_stats <- rhdf5::h5read(file,"training_stats")
  colnames(training_stats$elbo_terms) <- attr(rhdf5::h5read(file,"training_stats/elbo_terms", read.attributes=T),"colnames")
  
  # Load training options
  # training_opts <- as.list(h5read(file,"training_opts", read.attributes=T))
  
  # Load training data
  # training_data <- h5read(file,"training_data", read.attributes=T)
  # for (view in names(training_data)) {
  #   tmp <- h5read(file,str_c("training_data/",view), read.attributes=T)
  #   rownames(training_data[[view]]) <- attr(tmp,"rownames")
  #   colnames(training_data[[view]]) <- attr(tmp,"colnames")
  # }
  
  # Load training data
  data <- rhdf5::h5read(file,"data")
  featuredata <- rhdf5::h5read(file,"features")
  sampledata <- rhdf5::h5read(file,"samples")
  for (m in names(data)) {
    rownames(data[[m]]) <- sampledata
    colnames(data[[m]]) <- featuredata[[m]]
  }
  
  # Specify dimensions 
  M=length(data)
  N=nrow(data[[1]])
  D=sapply(data,ncol)
  K=tail(training_stats$activeK,n=1)
  dim=list("M"=M, "N"=N, "D"=D, K=K)
  
  # Create object
  gfa <- new("GFATrainedModel", 
             TrainData=data, TrainStats=training_stats, 
             Expectations=expectations, Parameters=parameters,
             Dimensions=dim)
  
  return(gfa)
}

