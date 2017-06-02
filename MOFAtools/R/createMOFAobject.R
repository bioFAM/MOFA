
#' @title Initialize a \code{\link{MOFAmodel}} object with multi-omic input data
#' @name createMOFAobject
#' @description Method to intilize the InputData slot of the MOFA object starting from a MutliAssayExperiment object
#' @param mae A MultiAssayExperiment containing the input data for MOFA
#' @param minViews Number of views a sample has to be observed to be included (default 1, i.e. including all samples)
#' @details fill this
#' @return an untrained MOFA model object
#' @export


createMOFAobject <- function(mae, minViews = 1) {
  
  if(class(mae) != "MultiAssayExperiment") stop ("mae has to be a MultiAssayExperiment object, for matrices use createMOFAobjectFromMatrixList")
  
  ##=== initialize MOFA object and set inputData====
  object <- new("MOFAmodel")
  object@Status <- "untrained"
  object@InputData <- mae 
  message(paste("Creating MOFA object from", length(mae), "experiments"))
  
  
  ##=== arrange to matrices in TrainData====
  # find samples which have been onserved in less than minViews
  if(!minViews %in% 1: length(mae)) 
    stop("'minViews' needs to be in 1 : number of experiments in 'mae'")
  sub <- combn(1:length(mae), minViews)
  samples2include <- apply(apply(sub, 2, function(s) complete.cases(mae[,,s])), 1, any)

  # drop samples which have been onserved in less than minViews
  mae_sub <- mae[,samples2include,]
  if(minViews>1) 
    message(paste("Removing ", sum(!samples2include), " sample not present in at least 'minViews'= ", minViews, ". Remaing samples for training: ", sum(samples2include), sep=""))
  
  # re-arrange data for training in MOFA to matrices, fill in NAs and store in TrainData slot
  data4train <- lapply(names(mae), function(exp) subset_augment(mae@ExperimentList[[exp]], sampleMap(mae_sub)$colname))
  object@TrainData <- data4train

  # set dimensionalities
  object@Dimensions[["M"]] <- length(object@TrainData)
  object@Dimensions[["N"]] <- ncol(object@TrainData[[1]])
  object@Dimensions[["D"]] <- sapply(object@TrainData,nrow)
  object@Dimensions[["K"]] <- NULL
  
  # set view names
  viewNames(object) <- names(mae)
  
  message("Returning the following MOFA object:")
  print(object)
  
  return(object)
}


createMOFAobjectFromList <- function(ExpList, pData = S4Vectors::DataFrame(), minViews=1){

  mae <- MultiAssayExperiment(experiments = ExpList, 
                             pData = pData)
  
  object <- createMOFAobject(mae, minViews=minViews)
  
  return(object)
}
