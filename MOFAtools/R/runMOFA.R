
###########################
## Functions to run MOFA ##
###########################

#' @title train a MOFA model
#' @name runMOFA
#' @description Function to train an untrained \code{\link{MOFAmodel}} object.
#' @details In this step the R package is calling the \code{mofa} Python package, where the the training is performed. \cr
#' The interface with Python is done with the \code{\link{reticulate}} package. 
#' If you have several versions of Python installed and Rstudio is not detecting the correct one, you can change it using
#' \code{reticulate::use_python}.
#' @param object an untrained \code{\link{MOFAmodel}} object
#' @param outfile output .hdf5 file
#' @return a trained \code{\link{MOFAmodel}} object
#' @import reticulate
#' @export
#' @examples
#' data("CLL_data")
#' # create and prepare the MOFAmodel
#' MOFAobject <- createMOFAobject(CLL_data)
#' MOFAobject <- prepareMOFA(MOFAobject)
#' # fit the model (takes some time)
#' \dontrun{
#' # MOFAobject <- runMOFA(MOFAobject)
#' }

runMOFA <- function(object, outfile) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) 
    stop("'object' has to be an instance of MOFAmodel")
    
  if (object@Status=="trained") 
    stop("The model is already trained! If you want to retrain, create a new untrained MOFAmodel")
  
  # Initiate reticulate
  mofa <- import("mofa")
  mofa_entrypoint <- mofa$core.entry_point$entry_point()
  
  # Pass data
  mofa_entrypoint$set_data(data=lapply(object@TrainData, function(x) r_to_py(t(x))))
  
  # Pass model options 
  mofa_entrypoint$set_model_options(
    factors        = object@ModelOptions$numFactors,
    likelihoods    = unname(object@ModelOptions$likelihood),
    learnIntercept = object@ModelOptions$learnIntercept,
    sparsity       = object@ModelOptions$sparsity
  )
  
  # Pass data processing options
  mofa_entrypoint$set_data_options(
    view_names              = viewNames(object), 
    center_features         = object@DataOptions$centerFeatures,
    scale_views             = object@DataOptions$scaleViews,
    RemoveIncompleteSamples = object@DataOptions$removeIncompleteSamples
  )
  
  # Parse data
  mofa_entrypoint$parse_data()
  
  # Pass training options  
  mofa_entrypoint$set_train_options(
    iter       = object@TrainOptions$maxiter,
    tolerance  = object@TrainOptions$tolerance,
    dropR2     = object@TrainOptions$DropFactorThreshold,
    seed       = object@TrainOptions$seed, 
    verbose    = object@TrainOptions$verbose
  )
  

  # Define priors
  mofa_entrypoint$define_priors()
  
  # Define initialisations
  mofa_entrypoint$define_init()
  
  # Parse the intercept factor
  mofa_entrypoint$parse_intercept()
  
  # Train the model
  mofa_entrypoint$train_model()
  
  # Save the model
  sample_names <- colnames(object@TrainData[[1]])
  feature_names <- unname(lapply(object@TrainData,rownames))
  mofa_entrypoint$save_model(outfile, sample_names=sample_names, feature_names=feature_names)
  
  # Load the model back into R
  object <- loadModel(outfile, object)
  
  return(object)
}
