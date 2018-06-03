
###########################
## Functions to run MOFA ##
###########################

#' @title train a MOFA model
#' @name runMOFA
#' @description Function to train an untrained \code{\link{MOFAmodel}} object.
#' @details In this step the R package is calling the \code{mofa} Python package, where the the training is performed. \cr
#' The interface with Python is done with the \code{\link{reticulate}} package. 
#' If you have several versions of Python installed and Rstudio is not detecting the correct one, you can change it using
#' \code{reticulate::use_python}. \cr
#' This module is in beta testing so please, read our FAQ for troubleshooting and report any problems.
#' @param object an untrained \code{\link{MOFAmodel}} object
#' @return a trained \code{\link{MOFAmodel}} object
#' @import reticulate
#' @export

#' @examples
#' data("CLL_data")
#' #create and prepare the MOFaobject before using runMOFA
#' MOFAobject <- createMOFAobject(CLL_data)
#' MOFAobject <- prepareMOFA(MOFAobject)
#' # runMOFA takes some time for model fitting
#' \dontrun{
#' # MOFAobject <- runMOFA(MOFAobject)}

runMOFA <- function(object) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) 
    stop("'object' has to be an instance of MOFAmodel")
    
  if (object@Status=="trained") 
    stop("The model is already trained! If you want to retrain, create a new untrained MOFAmodel")
  
  if (object@TrainOptions$learnFactors == FALSE) {
    object@TrainOptions$DropFactorThreshold <- 0.00
  } else {
    if (object@TrainOptions$DropFactorThreshold==0) { 
      print("Warning: learnFactors is set to TRUE but dropFactorThreshold is 0, this is contradictory.")
      print("Please read the documentation in prepareMOFA about learning the number of factors.")
    }
  }
  
  if (file.exists(object@DirOptions$outFile))
    message("Warning: Output file already exists, it will be replaced")
  
  # Initiate reticulate
  mofa <- import("mofa")
  mofa_entrypoint <- mofa$core.init_asd2$entry_point()
  
  # Pass data options
  mofa_entrypoint$set_data_options(
    inFiles     = paste0(object@DirOptions$dataDir, "/", viewNames(object), ".txt"), 
    outFile     = object@DirOptions$outFile, 
    views       = viewNames(object), 
    delimiter   = "\t", 
    header_cols = TRUE, 
    header_rows = TRUE
  )
  
  # Pass training options  
  mofa_entrypoint$set_train_options(
    iter       = object@TrainOptions$maxiter,
    tolerance  = object@TrainOptions$tolerance,
    dropR2     = object@TrainOptions$DropFactorThreshold,
    seed       = object@TrainOptions$seed, 
    verbose    = object@TrainOptions$verbose
  )
  
  # Pass model options 
  mofa_entrypoint$set_model_options(
    factors        = object@ModelOptions$numFactors,
    likelihoods    = unname(object@ModelOptions$likelihood),
    learnIntercept = object@ModelOptions$learnIntercept,
    sparsity       = object@ModelOptions$sparsity
  )

  # Pass data processing options
  mofa_entrypoint$set_dataprocessing_options(
    center_features         = object@DataOptions$centerFeatures,
    scale_views             = object@DataOptions$scaleViews,
    RemoveIncompleteSamples = object@DataOptions$removeIncompleteSamples
    
  )

  # Load data
  mofa_entrypoint$load_data()
  
  # Define Priors
  mofa_entrypoint$define_priors()
  
  # Initialise variational distributions
  mofa_entrypoint$initialise_variational()
  
  # Parse the intercept factor
  mofa_entrypoint$parse_intercept()
  
  # Train the model
  mofa_entrypoint$train_model()
  
  # Load the trained model
  object <- loadModel(object@DirOptions$outFile, object)
  
  return(object)
}
