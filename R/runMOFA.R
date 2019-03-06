
###########################
## Functions to run MOFA ##
###########################

#' @title train a MOFA model
#' @name runMOFA
#' @description Function to train an untrained \code{\link{MOFAmodel}} object.
#' @details In this step the R package is calling the \code{mofapy} Python package,
#'  where the the training is performed. \cr
#' The interface with Python is done with the \code{\link{reticulate}} package. 
#' If you have several versions of Python installed and Rstudio is not detecting
#'  the correct one, you can change it using \code{reticulate::use_python}.
#' @param object an untrained \code{\link{MOFAmodel}} object
#' @param outfile output .hdf5 file
#' @return a trained \code{\link{MOFAmodel}} object
#' @import reticulate
#' @export
#' @examples
#' data <- makeExampleData()
#' # create and prepare the MOFAmodel
#' MOFAobject <- createMOFAobject(data)
#' MOFAobject <- prepareMOFA(MOFAobject)
#' # fit the model (takes some time)
#' \dontrun{
#' # MOFAobject <- runMOFA(MOFAobject)
#' # MOFAobject
#' }

runMOFA <- function(object, outfile=NULL) {
  
  # Sanity checks on the model
  if (!is(object, "MOFAmodel")) 
    stop("'object' has to be an instance of MOFAmodel")
  if (length(object@ModelOptions)==0)
    stop("Can't find the model options, did you run prepareMOFA?")
  if (length(object@TrainOptions)==0)
    stop("Can't find the training options, did you run prepareMOFA?")
  if (length(object@DataOptions)==0)
    stop("Can't find the data options, did you run prepareMOFA?")
  
  # Sanity checks on the output file
  if (is.null(outfile)) {
    print("No output file provided, using a temporary file...")
    outfile <- tempfile()
  } else {
    if (!dir.exists(dirname(outfile))) {
      print("Output directory not found, creating it...")
      dir.create(dirname(outfile), recursive = TRUE, showWarnings = TRUE)
    }
  }
    
  if (object@Status=="trained") 
    stop("The model is already trained! If you want to retrain, create a new untrained MOFAmodel")
  
  # Initiate reticulate
  mofa <- tryCatch({
    import("mofapy")
  }, error = function(err) {
    import("mofa") }
  )
  mofa_entrypoint <- mofa$core.entry_point$entry_point()
  
  # Pass data
  mofa_entrypoint$set_data(data=unname(lapply(object@TrainData, function(x) r_to_py(t(x)))))
  
  # Pass model options 
  mofa_entrypoint$set_model_options(
    factors        = object@ModelOptions$numFactors,
    likelihoods    = unname(object@ModelOptions$likelihood),
    learnIntercept = FALSE,
    sparsity       = object@ModelOptions$sparsity
  )
  
  #TO BE FIXED: This should also be saved along with ModelOptions.
  numFactors <- object@ModelOptions$numFactors
  
  # Pass data processing options
  mofa_entrypoint$set_data_options(
    view_names              = viewNames(object), 
    center_features         = TRUE,
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
  mofa_entrypoint$save_model(
    outfile, sample_names=sample_names, feature_names=feature_names
  )
  
  # Load the model back into R
  object <- loadModel(outfile, object)
  
  #TO BE FIXED: Workaround as numFactors is not saved in ModelOpts
  object@ModelOptions$numFactors <- numFactors
  
  return(object)
}
