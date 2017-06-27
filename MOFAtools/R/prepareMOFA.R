
#' @title prepareMOFA: Prepare an untrained MOFA object for running MOFA
#' @name prepareMOFA
#' @description Function to set the training and model options, produces .txt files that are used for python as input and 
#' creates an .sh file for calling MOFA with the specified options from the command line. These files are all stored in the specified directory.
#' @param object an untrained MOFA object
#' @param dir directory to store .txt and .sh files in
#' @param ModelOptions list of ModelOptions (see getDefaultModelOpts for what options can be set here). If none specified, default options are used.
#' @param TrainOptions list of TrainOptions (see getDefaultTrainOptions for what options can be set here). If none specified, default options are used.
#' @param outFile name of output file from Python MOFA
#' @param k number of latent factors to start with (default = 10)
#' @param MOFAdir directory of the MOFA Pyhton package installation
#' @details fill this
#' @return a untrained MOFA object with specified ModelOpts and TrainOpts 
#' @export

prepareMOFA <- function(object, DirOptions, ModelOptions = NULL, TrainOptions = NULL) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  # Create temporary folder to store the input matrices
  DirOptions$tmpDir <- tempdir()
  dir.create(DirOptions$tmpDir, showWarnings = FALSE)
  
  # Store views as matrices in .txt files
  message(sprintf("Storing input views in %s...", DirOptions$tmpDir))
  for(view in viewNames(object)) {
    write.table(t(object@TrainData[[view]]), file=file.path(DirOptions$tmpDir, paste0(view,".txt")),
                sep=" ", row.names=TRUE, col.names=TRUE, quote=F)
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
#' @details fill this
#' @return list with training options
#' @export

getDefaultTrainOpts <- function(silent=T) {
  if(!silent) message("Using default training options...")
  TrainOpts <- list(
    trials = 1,          # Number of trials
    
    maxiter = 2000,      # Maximum number of iterations
    tolerance = 0.01,    # Convergence threshold based on change in the evidence lower bound
    forceiter = 0,       # Do not stop when convergence criterion is met, force model to complete all iterations
    
    elbofreq = 1,        # Frequency of evidence lower bound calculation
    
    startdrop = 5,       # Initial iteration to start dropping factors
    freqdrop = 1,        # Frequency of dropping latent factors
    drop_by_norm = 0.00,  # Option to drop latent factors: minimum norm threshold
    # drop_by_cor = NA,   # (NOT IMPLEMENTED) Option to drop latent factors: maximum correlation between two factors
    drop_by_r2 = 0.00,   # Option to drop latent factors: minimum coefficient of determination (percentage of variance explained in at least one view)
    
    verbosity = 2,       # Verbosity (TO BE DEFINED)
    cores = 1            # Number of cores (usually it is one core per trial)
  )
  return(TrainOpts)
}

#' @title getDefaultModelOpts: Get default model options
#' @name getDefaultModelOpts
#' @param object  untrained MOFA object to get model options for
#' @param silent  boolean whether to print warnings
#' @description Function to obtain default model options
#' @details fill this
#' @return  list with training options
#' @export
#' 
getDefaultModelOpts <- function(object, silent=T) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  if(!.hasSlot(object,"Dimensions") | length(object@Dimensions) == 0)
    stop("Dimensions of object need to be defined before getting ModelOpts")
  if(!.hasSlot(object,"InputData"))
    stop("Input data needs to be specified before getting ModelOpts")
  
  # Guess likelihood type
  likelihood = rep("gaussian", object@Dimensions[["M"]]); names(likelihood) <- viewNames(object)
  for (view in viewNames(object)) {
    data <- MultiAssayExperiment::experiments(object@InputData)[[view]]
    if (all(data %in% c(0,1,NA))) {
      likelihood[view] <- "bernoulli"
    } else if (all(data%%1==0)) {
      likelihood[view] <- "poisson"
    }
  }
  
  # Define default model options
  ModelOptions <- list(
    learnMean = TRUE,
    learnTheta = rep(1,object@Dimensions[["M"]]),
    initTheta = rep(0.5,object@Dimensions[["M"]]),
    likelihood = likelihood,
    initialK = 10,
    schedule = c("Y","SW","Z","AlphaW","Theta","Tau"),
    covariatesFile = NULL,
    scale_covariates = NULL
  )
  
  return(ModelOptions)
}
