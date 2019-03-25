
####################################################
## Functions to prepare a MOFA model for training ##
####################################################

#' @title regress out a covariate from the training data
#' @name regressCovariates
#' @description Function to regress out a covariate from the training data.\cr
#' If you have technical sources of variability (i.e. batch effects) that you do not want to be captured by factors in the model, 
#' you should regress them out before fitting MOFA. This function performs a simple linear regression model, extracts the residuals,
#' and replaces the original data in the TrainingData slot. \cr
#' Why is this important? If big technical factors exist, the model will "focus" on capturing the variability driven by these factors, and smaller sources of variability could be missed. \cr
#' But... can we not simply add those covariates to the model? Technically yes, but we extensively tested this functionality and it was not yielding good results. \cr 
#' The reason is that covariates are usually discrete labels that do not reflect the underlying molecular biology. 
#' For example, if you introduce age as a covariate, but the actual age is different from the "molecular age", 
#' the model will simply learn a new factor that corresponds to this "latent" molecular age, and it will drop the covariate from the model.\cr
#' We recommend factors to be learnt in a completely unsupervised manner and subsequently relate them to the covariates via visualisation or via a simple correlation analysis (see our vignettes for more details).
#' @param object an untrained \code{\link{MOFAmodel}}
#' @param views the view(s) to regress out the covariates.
#' @param covariates a vector (one covariate) or a data.frame (for multiple covariates) where each row corresponds to one sample, sorted in the same order as in the input data matrices. 
#' You can check the order by doing sampleNames(MOFAobject). If required, fill missing values with \code{NA}, which will be ignored when fitting the linear model.
#' @param min_observations number of non-missing observations required
#' @return Returns an untrained \code{\link{MOFAmodel}} where the specified covariates have been regressed out in the training data.
#' @importFrom stats lm
#' @export
#' @examples 
#' data("CLL_data", package = "MOFAdata")
#' data("CLL_covariates", package = "MOFAdata")
#' library(MultiAssayExperiment)
#' mae_CLL <- MultiAssayExperiment(
#' experiments = CLL_data, 
#' colData = CLL_covariates
#' )
#' MOFAobject <- createMOFAobject(mae_CLL)
#' MOFAobject <- prepareMOFA(MOFAobject)
#' MOFAobject_reg <- regressCovariates(
#' object = MOFAobject,
#' views = c("Drugs","Methylation","mRNA"),
#' covariates = MOFAobject@InputData$Gender
#' )
#' # MOFA object with training data after regressing out the specified covariate
#' MOFAobject_reg 

regressCovariates <- function(object, views, covariates, min_observations = 5) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) 
    stop("'object' has to be an instance of MOFAmodel")
  if (Status(object)=="trained")
    stop("Status(object) is 'trained'. regressCovariates has to be done before training the model")
  if (length(ModelOptions(object)$likelihood)==0) 
    stop("Run prepareMOFA before regressing out covariates") 
  if (any(ModelOptions(object)$likelihood[views]!="gaussian")) 
    stop("Some of the specified views contains discrete data. \nRegressing out covariates only works in views with continuous (gaussian) data")
  
  # Fetch data
  Y <- getTrainData(object, views=views)
  all_samples <- sampleNames(object)
  
  # Prepare data.frame with covariates
  if (!is(covariates,"data.frame"))
    covariates <- data.frame(x=covariates)
  stopifnot(nrow(covariates)==Dimensions(object)$N)
  
  Y_regressed <- list()
  for (m in views) {
    if (any(rowSums(!is.na(Y[[m]]))<min_observations) ) stop(sprintf("Some features do not have enough observations (N=%s) to fit the linear model",min_observations))
    Y_regressed[[m]] <- t( apply(Y[[m]], 1, function(y) {
      
      # Fit linear model
      df <- cbind(y,covariates)
      lm.out <- lm(y~., data=df)
      residuals <- lm.out[["residuals"]]+lm.out[["coefficients"]][1]
      
      # Fill missing values
      missing_samples <- all_samples[!all_samples %in% names(residuals)]
      residuals[missing_samples] <- NA
      residuals[all_samples]
    }))
  }
  TrainData(object)[views] <- Y_regressed
  
  return(object)
}


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
#' data("CLL_data", package = "MOFAdata")
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
#' MOFAobject


prepareMOFA <- function(object, DataOptions = NULL, ModelOptions = NULL, TrainOptions = NULL) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) 
    stop("'object' has to be an instance of MOFAmodel")
  if (getDimensions(object)[["N"]] < 15) {
      warning("This model is not appropriate for data sets with less than ~15 samples")
  }
  if (getDimensions(object)[["N"]] < getDimensions(object)[["K"]]) {
      warning("There are less samples than factors, likely to generate numerical errors")
  }
  if (min(getDimensions(object)[["D"]]) < getDimensions(object)[["K"]]) {
      warning("There are less factors than features, likely to generate numerical errors")
  }
  
  # Get data options
  message("Checking data options...")
  if (is.null(DataOptions)) {
    message("No data options specified, using default...")
    DataOptions(object) <- getDefaultDataOptions()
  } else {
    if (!is(TrainOptions,"list") & !all(names(TrainOptions) == names(getDefaultTrainOptions())))
      stop("DataOptions are incorrectly specified, please read the documentation in getDefaultDataOptions")
    DataOptions(object) <- DataOptions
  }
  if (any(nchar(sampleNames(object))>50))
    warning("Due to string size limitations in the HDF5 format, sample names will be trimmed to less than 50 characters")
  
  # Get training options
  message("Checking training options...")
  if (is.null(TrainOptions)) {
    message("No training options specified, using default...")
    TrainOptions(object) <- getDefaultTrainOptions()
  } else {
    if(!is(TrainOptions,"list") & !all(names(TrainOptions) == names(getDefaultTrainOptions())))
      stop("TrainOptions are incorrectly specified, please read the documentation in getDefaultTrainOptions")
    TrainOptions(object) <- TrainOptions
  }
  
  # Get model options
  message("Checking model options...")
  if(is.null(ModelOptions)) {
    message("No model options specified, using default...")
    ModelOptions(object) <- getDefaultModelOptions(object)
  } else {
    # (To-do) Check that ModelOptions is correct
    if(!is(ModelOptions,"list") & !all(names(ModelOptions) == names(getDefaultModelOptions(object))))
      stop("ModelOptions are incorrectly specified, please read the documentation in getDefaultModelOptions")
    ModelOptions(object) <- ModelOptions
  }
  
  # Convert binary data to numeric
  idx <- names(which(ModelOptions(object)$likelihood == "bernoulli"))
  if (length(idx)>0) {
    for (i in idx) {
      foo <- TrainData(object)[[i]]
      TrainData(object)[[i]] <- as.numeric(TrainData(object)[[i]])
      dim(TrainData(object)[[i]]) <- dim(foo)
      rownames(TrainData(object)[[i]]) <- rownames(foo)
      colnames(TrainData(object)[[i]]) <- colnames(foo)
    }
  }
  
  # Make sure that there are no features with zero variance
  for (m in seq_along(TrainData(object))) {
      if (!all(apply(TrainData(object)[[m]],1,var,na.rm=TRUE)>0, na.rm=TRUE))
        sprintf("Error: there are features with zero variance in view '%s', please remove them and create a new MOFAobject",viewNames(object)[m])
  }
  
  # Store feature-wise means
  # FeatureIntercepts(object) <- lapply(TrainData(object),rowMeans,na.rm=TRUE)
  
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
#' }
#' @return Returns a list with the default model options, which have to be passed as
#'  an argument to \code{\link{prepareMOFA}}
#' @export
#' @examples 
#' # load data
#' data("CLL_data", package = "MOFAdata")
#' #create a MOFAmodel object
#' MOFAobject <- createMOFAobject(CLL_data)
#' # set model options
#' ModelOptions <- getDefaultModelOptions(MOFAobject)
#' ModelOptions

getDefaultModelOptions <- function(object) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  if (!.hasSlot(object,"Dimensions") | length(Dimensions(object)) == 0) 
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
    numFactors = ceiling(nsamples/2), # (numeric) initial number of latent factors
    sparsity = TRUE                   # (logical) use feature-wise sparsity?
  )
  
  return(ModelOptions)
}
