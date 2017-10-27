
#' @title prepareMOFA: Prepare an untrained MOFA object for training
#' @name prepareMOFA
#' @description Function to set the training and model options, produces .txt files that are used for python as input and 
#' creates an .sh file for calling MOFA with the specified options from the command line. These files are all stored in the specified directory.
#' @param object an untrained MOFA object
#' @param dir directory to store .txt and .sh files in
#' @param DataOptions list of DataOptions (see getDefaultDataOpts for what options can be set here). If NULL, default options are used.
#' @param ModelOptions list of ModelOptions (see getDefaultModelOpts for what options can be set here). If NULL, default options are used.
#' @param TrainOptions list of TrainOptions (see getDefaultTrainOptions for what options can be set here). If NULL, default options are used.
#' @param outFile name of output file from Python MOFA
#' @param k number of latent factors to start with (default = 10)
#' @param MOFAdir directory of the MOFA Pyhton package installation
#' @details fill this
#' @return a untrained MOFA object with specified ModelOpts and TrainOpts 
#' @export

prepareMOFA <- function(object, DirOptions, DataOptions = NULL, ModelOptions = NULL, TrainOptions = NULL) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  # Create temporary folder to store data
  dir.create(DirOptions$dataDir, showWarnings = FALSE)
  
  # Store views as matrices in .txt files
  message(sprintf("Storing input views in tmp folder %s...", DirOptions$dataDir))
  for(view in viewNames(object)) {
    write.table(t(object@TrainData[[view]]), file=file.path(DirOptions$dataDir, paste0(view,".txt")),
                sep=" ", row.names=TRUE, col.names=TRUE, quote=F)
  }
  
  # Store covariates as a .txt file
  # if (!is.null(ModelOptions$covariates)) {
  # write.table(ModelOptions$covariates, file=file.path(DirOptions$dataDir, "covariates.txt"), 
  #             sep=" ", row.names=F, col.names=F, quote=F)
  # }
  
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
  
  return(object)
}



#' @title getDefaultTrainOpts: Get default training options
#' @name getDefaultTrainOpts
#' @description Function to obtain default training options
#' @return list with training options
#' @export
getDefaultTrainOpts <- function() {
  TrainOpts <- list(
    maxiter = 10000,              # Maximum number of iterations
    tolerance = 0.01,            # Convergence threshold based on change in the evidence lower bound
    DropFactorThreshold = 0.03   # Threshold on fraction of variance explained to drop a factor
  )
  return(TrainOpts)
}


#' @title getDefaultDataOpts: Get default data options
#' @name getDefaultDatasOpts
#' @description Function to obtain default data options
#' @return list with data options
#' @export
getDefaultDataOpts <- function() {
  DataOpts <- list(
    centerFeatures = F,   # Center features to zero mean (does not apply to binary or count views)
    scaleViews = T        # Scale views to unit variance (does not apply to binary or count views)
  )
  return(DataOpts)
}

#' @title getDefaultModelOpts: Get default model options
#' @name getDefaultModelOpts
#' @param object  untrained MOFA object to get model options for
#' @description Function to obtain default model options
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
    learnIntercept = TRUE,      # (bool) include a constant factor of 1s to learn the mean of features (intercept)? If not, you need to center the data
    likelihood = likelihood,    # (character vector) likelihood per view [gaussian/bernoulli/poisson]
    numFactors = 25,            # (numeric) initial number of latent factors
    sparsity=T                  # use feature-wise sparsity?
    # covariates = NULL
  )
  
  return(ModelOptions)
}
