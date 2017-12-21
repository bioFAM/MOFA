
#' @title Prepare an untrained MOFA object for training
#' @name prepareMOFA
#' @description Function to prepare a MOFA model for training by defining data, model and training options.
#' @param object an untrained MOFA object
#' @param DirOptions list with I/O options, it must contain a 'dataDir' element where temporary text files will be stored and a 'outFile' where the final model will be stored as an hdf5 object.
#' @param DataOptions list of DataOptions (see getDefaultDataOpts for what options can be set here). If NULL, default data options are used.
#' @param ModelOptions list of ModelOptions (see getDefaultModelOpts for what options can be set here). If NULL, default model options are used.
#' @param TrainOptions list of TrainOptions (see getDefaultTrainOpts for what options can be set here). If NULL, default training options are used.
#' @return a untrained MOFA object with specified DataOpts, ModelOpts and TrainOpts 
#' @export

prepareMOFA <- function(object, DirOptions, DataOptions = NULL, ModelOptions = NULL, TrainOptions = NULL) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  # Create temporary folder to store data
  dir.create(DirOptions$dataDir, showWarnings = F)
  
  # Get data options
  message("Checking data options...")
  if (is.null(DataOptions)) {
    message("No data options specified, using default...")
    object@DataOpts <- getDefaultDataOpts()
  } else {
    # (To-do) Check that DataOpts is correct
    object@DataOpts <- DataOptions
  }
  
  # Get training options
  message("Checking training options...")
  if (is.null(TrainOptions)) {
    message("No training options specified, using default...")
    object@TrainOpts <- getDefaultTrainOpts()
  } else {
    # (To-do) Check that TrainOpts is correct
    # if(!class(TrainOptions) == "list" & !all(names(TrainOptions) == names(getDefaultTrainOpts()))) 
    #   stop("'TrainOptions' are misspecified, use the list format provided by getDefaultTrainOpts()")
    object@TrainOpts <- TrainOptions
  }
  
  # Get model options
  message("Checking model options...")
  if(is.null(ModelOptions)) {
    message("No model options specified, using default...")
    object@ModelOpts <- getDefaultModelOpts(object)
  } else {
    # (To-do) Check that ModelOpts is correct
    # if(!class(ModelOptions) == "list" & !all(names(ModelOptions) == names(getDefaultModelOpts(object, silent=T)))) 
      # stop("'TrainOpts' are misspecified, use the list format provided by getDefaultModelOpts()")
    object@ModelOpts <- ModelOptions
  }
  
  # Store views as matrices in .txt files
  message(sprintf("Storing input views in %s...", DirOptions$dataDir))
  for(view in viewNames(object)) {
    write.table(t(object@TrainData[[view]]), file=file.path(DirOptions$dataDir, paste0(view,".txt")),
                sep=object@DataOpts$delimiter, row.names=T, col.names=T, quote=F)
  }
  
  # Store covariates as a .txt file
  if (!is.null(ModelOptions$covariates)) {
    write.table(ModelOptions$covariates, file=file.path(DirOptions$dataDir, "covariates.txt"),
                sep=object@DataOpts$delimiter, row.names=F, col.names=F, quote=F)
  }
  
  # If output already exists, remove it
  if (file.exists(DirOptions$outFile)) {
    file.remove(DirOptions$outFile)
  }
  
  return(object)
}



#' @title Get default training options
#' @name getDefaultTrainOpts
#' @description Function to obtain default training options
#' @details The training options are the following: \cr
#' maxiter: Maximum number of iterations. \cr
#' tolerance: Convergence threshold based on the change in Evidence Lower Bound, we recommend this be around 0.01. \cr
#' DropFactorThreshold: Threshold on fraction of variance explained to drop a factor. That is, factors explaining less than 'DropFactorThreshold' fraction of variance (in all views) will be dropped. \cr
#' @return list with default training options
#' @export
getDefaultTrainOpts <- function() {
  TrainOpts <- list(
    maxiter = 10000,              # Maximum number of iterations
    tolerance = 0.01,            # Convergence threshold based on change in the evidence lower bound
    DropFactorThreshold = 0.03,   # Threshold on fraction of variance explained to drop a factor
    verbose = F
  )
  return(TrainOpts)
}


#' @title Get default data options
#' @name getDefaultDataOpts
#' @description Function to obtain default data options
#' @details The data options are the following: \cr
#' centerFeatures: boolean indicating whether to center the features to zero mean. This is not required as long as the option learnIntercept is set to TRUE in the model options. Default is FALSE.\cr
#' scaleViews: boolean indicating whether to scale the views to unit variance. This is optional and recommended, but not required. Default is FALSE. \cr
#' removeIncompleteSamples: boolean indicating whether to remove incomplete samples that are not profiled in all omics. Default is FALSE
#' @return list with default data options
#' @export
getDefaultDataOpts <- function() {
  DataOpts <- list(
    delimiter = "\t",   # Delimiter for the data
    centerFeatures = F,   # Center features to zero mean (does not apply to binary or count views)
    scaleViews = F,        # Scale views to unit variance (does not apply to binary or count views)
    removeIncompleteSamples = F # Remove incomplete samples that are not profiled in all omics?
  )
  return(DataOpts)
}

#' @title Get default model options
#' @name getDefaultModelOpts
#' @param object untrained MOFA object to get model options for
#' @description Function to obtain default model options
#' @details The model options are the following: \cr
#' likelihood: character vector with data likelihoods per view, 'gaussian' for (roughly) normally distributed data, 'bernoulli' for binary data and 'poisson' for count data
#' numFactors: initial number of factors, we recommend this to be large enough, larger than 10. \cr
#' learnIntercept: boolean indicating whether to learn the intercept (the means) per feature. This is always recommended, particularly if you have non-gaussian likelihoods. \cr
#' sparsity: boolean indicating whether to use sparsity. This is always recommended, as it will make the loadings more interpretable. \cr
#' covariates: (TO-DEFINE).
#' @return  list with default model options
#' @export
getDefaultModelOpts <- function(object) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  if (!.hasSlot(object,"Dimensions") | length(object@Dimensions) == 0) stop("Dimensions of object need to be defined before getting ModelOpts")
  if (!.hasSlot(object,"InputData")) stop("InputData slot needs to be specified before getting ModelOpts")
  if (!.hasSlot(object,"TrainData")) stop("TrainData slot needs to be specified before getting ModelOpts")
  
  # Guess likelihood type
  likelihood = rep("gaussian", object@Dimensions[["M"]]); names(likelihood) <- viewNames(object)
  # for (view in viewNames(object)) {
  #   data <- getTrainData(object, view)
  #   if (all(data %in% c(0,1,NA))) {
  #     likelihood[view] <- "bernoulli"
  #   } else if (all(data%%1==0)) {
  #     likelihood[view] <- "poisson"
  #   }
  # }
  
  # Define default model options
  ModelOptions <- list(
    likelihood = likelihood,    # (character vector) likelihood per view [gaussian/bernoulli/poisson]
    learnIntercept = TRUE,      # (bool) include a constant factor of 1s to learn the mean of features (intercept)? If not, you need to center the data
    numFactors = 25,            # (numeric) initial number of latent factors
    sparsity = T,               # use feature-wise sparsity?
    covariates = NULL           # no covariates by default
  )
  
  return(ModelOptions)
}
