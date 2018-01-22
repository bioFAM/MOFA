
############################################
## Functions to load a trained MOFA model ##
############################################


#' @title loading a trained MOFA model
#' @name loadModel
#' @description Method to load a trained MOFA model. \cr
#' The training of MOFA is done using a Python framework, and the model output is saved as an .hdf5 file, which has to be loaded in the R package.
#' @param file an hdf5 file saved by the MOFA python framework.
#' @param object either NULL (default) or an an existing untrained MOFA object. If NULL, the \code{\link{MOFAmodel}} object is created from the scratch.
#' @param sortFactors boolean indicating whether factors should be sorted by variance explained (default is TRUE)
#' @return a \code{\link{MOFAmodel}} model.
#' @importFrom rhdf5 h5read
#' @export

loadModel <- function(file, object = NULL, sortFactors = T) {
  
  # message(paste0("Loading the following MOFA model: ", file))
  
  if (is.null(object)) object <- new("MOFAmodel")
  
  # if(.hasSlot(object,"Status") & length(object@Status) !=0)
  #   if (object@Status == "trained") warning("The specified object is already trained, over-writing training output with new results!")
  
  # Load expectations
  object@Expectations <- h5read(file,"expectations")
  object@Status <- "trained"
  
  # Load training statistics
  tryCatch( {
    object@TrainStats <- h5read(file, 'training_stats',read.attributes=T);
    colnames(object@TrainStats$elbo_terms) <- attr(h5read(file,"training_stats/elbo_terms", read.attributes=T),"colnames")
  }, error = function(x) { print("Training stats not found, not loading it...") })

  
  # Load training options
  if (length(object@TrainOpts) == 0) {
    tryCatch(object@TrainOpts <- as.list(h5read(file, 'training_opts',read.attributes=T)), error = function(x) { print("Training opts not found, not loading it...") })
  }
    
  # Load model options
  if (length(object@ModelOpts) == 0) {
    tryCatch(object@ModelOpts <- as.list(h5read(file, 'model_opts',read.attributes=T)), error = function(x) { print("Model opts not found, not loading it...") })
  }
  
  # TO REMOVE....
  if ("learnMean" %in% names(object@ModelOpts)) {
    tmp <- names(object@ModelOpts)
    tmp[tmp=="learnMean"] <- "learnIntercept"
    names(object@ModelOpts) <- tmp
  }
  object@ModelOpts$learnIntercept <- as.logical(object@ModelOpts$learnIntercept)
  
  # Load training data
  tryCatch( {
    TrainData <- h5read(file,"data")
    featureData <- h5read(file,"features")
    sampleData <- h5read(file,"samples")
    for (m in names(TrainData)) {
      rownames(TrainData[[m]]) <- sampleData
      colnames(TrainData[[m]]) <- featureData[[m]]
      TrainData[[m]][is.nan(TrainData[[m]])] <- NA
    }
    TrainData <- lapply(TrainData, t)
    object@TrainData <- TrainData
    }, error = function(x) { print("Error loading the training data...") })
  
  
  # Load dimensions
  object@Dimensions[["M"]] <- length(object@TrainData)
  object@Dimensions[["N"]] <- ncol(object@TrainData[[1]])
  object@Dimensions[["D"]] <- sapply(object@TrainData,nrow)
  # K=tail(training_stats$activeK[!is.nan(training_stats$activeK)],n=1)
  object@Dimensions[["K"]] <- ncol(object@Expectations$Z$E)
  
  # Set view, sample, feature and factor names
  viewNames(object) <- names(object@TrainData)
  sampleNames(object) <- colnames(object@TrainData[[1]])
  featureNames(object) <- lapply(object@TrainData,rownames)
  factorNames(object) <- as.character(1:object@Dimensions[["K"]])
  
  #
  names(object@ModelOpts$likelihood) <- viewNames(object)
  
  # Rename covariates, including intercept
  # if (object@ModelOpts$learnIntercept == TRUE) factorNames(object) <- c("intercept",as.character(1:(object@Dimensions[["K"]]-1)))
  # if (!is.null(object@ModelOpts$covariates)) {
  #   if (object@ModelOpts$learnIntercept == TRUE) {
  #     factorNames(object) <- c("intercept", colnames(object@ModelOpts$covariates), as.character((ncol(object@ModelOpts$covariates)+1:(object@Dimensions[["K"]]-1-ncol(object@ModelOpts$covariates)))))
  #   } else {
  #     factorNames(object) <- c(colnames(object@ModelOpts$covariates), as.character((ncol(object@ModelOpts$covariates)+1:(object@Dimensions[["K"]]-1))))
  #   }
  # }
  
  # Rename factors
  if (object@ModelOpts$learnIntercept == TRUE) {
    intercept_idx <- names(which(sapply(apply(object@Expectations$Z$E,2,unique),length)==1))
    factornames <- as.character(1:(object@Dimensions[["K"]]))
    factornames[factornames==intercept_idx] <- "intercept"
    factorNames(object) <- factornames
    # object@Dimensions[["K"]] <- object@Dimensions[["K"]] - 1
  }
  # if (!is.null(object@ModelOpts$covariates)) {
  #   stop("Not working")
  # }
  
  # Parse factors
  if ((object@Dimensions$K-as.numeric(object@ModelOpts$learnIntercept))>0) {
    
    # (TO-DO) Mask passenger factors
    # object <- detectPassengers(object)
  
    # Order factors in order of variance explained
    if (sortFactors == T) {
      r2 <- rowSums(calculateVarianceExplained(object,totalVar=T,plotit=F)$R2PerFactor)
      order_factors <- c(names(r2)[order(r2, decreasing = T)])
      if (object@ModelOpts$learnIntercept==T) { order_factors <- c("intercept",order_factors) }
      object <- subsetFactors(object,order_factors)
      if (object@ModelOpts$learnIntercept==T) { 
        factorNames(object) <- c("intercept",1:(object@Dimensions$K-1))
      } else {
        factorNames(object) <- c(1:object@Dimensions$K) 
      }
    }
    return(object)
  } else {
    stop("The model has no active factors")
  }
  
  # Check for intercept factors
  # findInterceptFactors(object)
}

