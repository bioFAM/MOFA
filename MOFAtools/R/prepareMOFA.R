
####################################################
## Functions to prepare a MOFA model for training ##
####################################################

#' @title prepare a MOFAobject for training
#' @name prepareMOFA
#' @description Function to prepare a \code{\link{MOFAmodel}} object for training.
#' Here, data, input/output option are specified and data, model and training options can be set.
#' @param object an untrained \code{\link{MOFAmodel}}
#' @param DataOptions list of DataOptions (see \code{\link{getDefaultDataOptions}} details). 
#' If NULL, default data options are used.
#' @param ModelOptions list of ModelOptions (see \code{\link{getDefaultModelOptions}} for details). 
#' If NULL, default model options are used.
#' @param TrainOptions list of TrainOptions (see \code{\link{getDefaultTrainOptions}} for details). 
#' If NULL, default training options are used.
#' @return Returns an untrained \code{\link{MOFAmodel}} with specified data, model and training options.
#' Next step is to train the model with \code{\link{runMOFA}}
#' @export
#' @examples 
#' # load data
#' data("CLL_Data")
#' #create a MOFAmodel object
#' MOFAobject <- createMOFAobject(CLL_data)
#' # set options
#' TrainOptions <- getDefaultTrainOptions()
#' ModelOptions <- getDefaultModelOptions(MOFAobject)
#' DataOptions <- getDefaultDataOptions()
#' # prepare MOFAmodel object for training
#' MOFAobject <- prepareMOFA(MOFAobject, 
#' DataOptions = DataOptions,
#' ModelOptions = ModelOptions,
#' TrainOptions = TrainOptions
#' )
#' # as a next step runMOFA can be used to train the model (can take some time)
#' # MOFAobject <- runMOFA(MOFAobject, outfile=tempfile())

prepareMOFA <- function(object, DataOptions = NULL, ModelOptions = NULL, TrainOptions = NULL) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) 
    stop("'object' has to be an instance of MOFAmodel")
  
  # Get data options
  message("Checking data options...")
  if (is.null(DataOptions)) {
    message("No data options specified, using default...")
    object@DataOptions <- getDefaultDataOptions()
  } else {
    if (!is(TrainOptions,"list") & !all(names(TrainOptions) == names(getDefaultTrainOptions())))
      stop("DataOptions are incorrectly specified, please read the documentation in getDefaultDataOptions")
    object@DataOptions <- DataOptions
  }
  
  # Get training options
  message("Checking training options...")
  if (is.null(TrainOptions)) {
    message("No training options specified, using default...")
    object@TrainOptions <- getDefaultTrainOptions()
  } else {
    if(!is(TrainOptions,"list") & !all(names(TrainOptions) == names(getDefaultTrainOptions())))
      stop("TrainOptions are incorrectly specified, please read the documentation in getDefaultTrainOptions")
    object@TrainOptions <- TrainOptions
  }
  
  # Get model options
  message("Checking model options...")
  if(is.null(ModelOptions)) {
    message("No model options specified, using default...")
    object@ModelOptions <- getDefaultModelOptions(object)
  } else {
    # (To-do) Check that ModelOptions is correct
    if(!is(ModelOptions,"list") & !all(names(ModelOptions) == names(getDefaultModelOptions(object))))
      stop("ModelOptions are incorrectly specified, please read the documentation in getDefaultModelOptions")
    object@ModelOptions <- ModelOptions
  }
  
  # Convert binary data to numeric
  idx <- names(which(object@ModelOptions$likelihood == "bernoulli"))
  if (length(idx)>0) {
    for (i in idx) {
      foo <- object@TrainData[[i]]
      object@TrainData[[i]] <- as.numeric(object@TrainData[[i]])
      dim(object@TrainData[[i]]) <- dim(foo)
      rownames(object@TrainData[[i]]) <- rownames(foo)
      colnames(object@TrainData[[i]]) <- colnames(foo)
    }
  }
  
  # Store feature-wise means
  object@FeatureMeans <- lapply(object@TrainData,rowMeans,na.rm=T)
  
  return(object)
}



#' @title Get default training options
#' @name getDefaultTrainOptions
#' @description Function to obtain the default training options.
#' @details The training options are the following: \cr
#' \itemize{
#'  \item{\strong{maxiter}:}{ numeric value indicating the maximum number of iterations. 
#'  Default is 5000, but we recommend using the 'tolerance' as convergence criteria.}
#'  \item{\strong{tolerance}:}{ numeric value indicating the convergence threshold based
#'   on the change in Evidence Lower Bound (deltaELBO). 
#'  For quick exploration we recommend this to be around 1.0,
#'   and for a thorough training we recommend a value of 0.01. Default is 0.1}
#'  \item{\strong{DropFactorThreshold}:}{ numeric hyperparamter to automatically learn the number of factors.
#'  It indicates the threshold on fraction of variance explained to consider a factor inactive and 
#'  automatically drop it from the model during training. 
#'  For example, a value of 0.01 implies that factors explaining less
#'  than 1\% of variance (in each view) will be dropped.
#'  Default is 0, which implies that only factors that explain no variance at all will be removed
#'  }
#'  \item{\strong{verbose}:}{ logical indicating whether to generate a verbose output.}
#'  \item{\strong{seed}:}{ random seed for reproducibility (default is NULL, which samples a random seed).}
#' }
#' @return Returns a list with default training options, which have to be passed
#'as an argument to \code{\link{prepareMOFA}}
#' @export
#' @examples 
#' TrainOptions <- getDefaultTrainOptions()
#' TrainOptions

getDefaultTrainOptions <- function() {
  TrainOptions <- list(
    maxiter = 5000,               # (numeric) Maximum number of iterations
    tolerance = 0.1,              # (numeric) Convergence threshold based on change in the evidence lower bound
    DropFactorThreshold = 0.00,   # (numeric) Threshold on fraction of variance explained to drop a factor
    verbose = FALSE,              # (logical) verbosity?
    seed = NULL                   # (numeric or NULL) random seed
  )
  return(TrainOptions)
}


#' @title Get default data options
#' @name getDefaultDataOptions
#' @description Function to obtain the default data options.
#' @details The data options are the following: \cr
#' \itemize{
#'  \item{\strong{centerFeatures}:}{ logical indicating whether to center the features to zero mean. 
#'  This only works for gaussian data. Default is TRUE.}
#'  \item{\strong{scaleViews}:}{ logical indicating whether to scale views to have the same unit variance. 
#'  As long as the scale differences between the data sets is not too high, this is not required.
#'   Default is FALSE.}
#'  \item{\strong{removeIncompleteSamples}:}{ logical indicating whether to remove samples that
#'   are not profiled in all omics. We recommend this only for testing,
#'    as the model can cope with samples having missing assays. Default is FALSE.}
#' }
#' @return Returns a list with the default data options, which have to be passed as
#'  an argument to \code{\link{prepareMOFA}}
#' @export
#' @examples 
#' DataOptions <- getDefaultDataOptions()
#' DataOptions

getDefaultDataOptions <- function() {
  DataOptions <- list(
    centerFeatures = TRUE,           # Center features to zero mean (does not apply to binary or count views)
    scaleViews = FALSE,              # Scale views to unit variance (does not apply to binary or count views)
    removeIncompleteSamples = FALSE  # Remove incomplete samples that are not profiled in all omics?
  )
  return(DataOptions)
}

#' @title Get default model options
#' @name getDefaultModelOptions
#' @param object an untrained \code{\link{MOFAmodel}} object
#' @description Function to obtain the default model options.
#' @details The model options are the following: \cr
#' \itemize{
#'  \item{\strong{likelihood}:}{ character vector with data likelihoods per view: 
#'  'gaussian' for continuous data, 'bernoulli' for binary data and 'poisson' for count data.
#'  By default, they are guessed internally.}
#'  \item{\strong{numFactors}:}{ numeric value indicating the initial number of factors. 
#'  If you want to learn the number of factors automatically we recommend
#'   setting this to a large value, default is 25.}
#'  \item{\strong{sparsity}:}{ logical indicating whether to use sparsity.
#'  This is always recommended, as it will make the loadings more interpretable. Default is TRUE.}
#'  \item{\strong{learnIntercept}:}{ logical indicating whether to learn an intercept
#'   term to capture differences in feature means. This is currently not functional. Therefore, we ask you to always center
#'    the data if it is gaussian (do not worry about the non-gaussian views) and keep learnIntercept to FALSE.}
#' }
#' @return Returns a list with the default model options, which have to be passed as
#'  an argument to \code{\link{prepareMOFA}}
#' @export
#' @examples 
#' # load data
#' data("CLL_Data")
#' #create a MOFAmodel object
#' MOFAobject <- createMOFAobject(CLL_data)
#' # set model options
#' ModelOptions <- getDefaultModelOptions(MOFAobject)
#' ModelOptions

getDefaultModelOptions <- function(object) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  if (!.hasSlot(object,"Dimensions") | length(object@Dimensions) == 0) 
    stop("Dimensions of object need to be defined before getting ModelOptions")
  if (!.hasSlot(object,"InputData")) 
    stop("InputData slot needs to be specified before getting ModelOptions")
  if (!.hasSlot(object,"TrainData")) 
    stop("TrainData slot needs to be specified before getting ModelOptions")
  
  # Guess likelihood type
  likelihood <- .inferLikelihoods(object)
  
  nsamples <- getDimensions(object)[["N"]]
  if (nsamples<10) print("Warning: too few samples for MOFA to be useful...")
  
  # Define default model options
  ModelOptions <- list(
    likelihood = likelihood,          # (character vector) likelihood per view [gaussian/bernoulli/poisson]
    learnIntercept = FALSE,           # (bool) include a constant factor of 1s to learn the mean of features (intercept)?
    numFactors = ceiling(nsamples/2), # (numeric) initial number of latent factors
    sparsity = TRUE                   # (logical) use feature-wise sparsity?
  )
  
  return(ModelOptions)
}
