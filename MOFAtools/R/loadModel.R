
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
#' @param minR2 minimum R2 threshold to call 'active' factors (default is 0.01).
#' @return a \code{\link{MOFAmodel}} model.
#' @importFrom rhdf5 h5read
#' @export
loadModel <- function(file, object = NULL, sortFactors = TRUE, minR2 = 0.01) {
  
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
  if (length(object@TrainOptions) == 0) {
    tryCatch(object@TrainOptions <- as.list(h5read(file, 'training_opts',read.attributes=T)), error = function(x) { print("Training opts not found, not loading it...") })
  }
  # different names in R and python package (drop_by_r2 in python corresponds to DropFactorThreshold in R) - needs to be adapted
  if("drop_by_r2" %in% names(object@TrainOptions)) {
    object@TrainOptions$DropFactorThreshold <- object@TrainOptions$drop_by_r2
    object@TrainOptions$drop_by_r2 <- NULL
  }
    
  # Load model options
  tryCatch(object@ModelOptions <- as.list(h5read(file, 'model_opts',read.attributes=T)), error = function(x) { print("Model opts not found, not loading it...") })
  object@ModelOptions$sparsity <- as.logical(object@ModelOptions$sparsity)
  
  
  # Load training data
  tryCatch( {
    TrainData <- h5read(file,"data")
    featureData <- h5read(file,"features")
    sampleData <- h5read(file,"samples")
    for (m in names(TrainData)) {
      rownames(TrainData[[m]]) <- sampleData
      colnames(TrainData[[m]]) <- featureData[[m]]
    }
    TrainData <- lapply(TrainData, t)
    object@TrainData <- TrainData
    }, error = function(x) { print("Error loading the training data...") })
  
  # Replace NaN by NA in the training data
  for (m in names(TrainData)) {
    TrainData[[m]][is.nan(TrainData[[m]])] <- NA
  }
  
  # Sanity check on the order of the likelihoods
  if (!is.null(attr(TrainData,"likelihood"))) {
    lik <- attr(TrainData,"likelihood")
    if (!all(object@ModelOptions$likelihood == lik)) {
      object@ModelOptions$likelihood <- lik
      names(object@ModelOptions$likelihood) <- names(TrainData)
    }
  }
  
  
  # Update old models
  object <- .updateOldModel(object)
  
  # Load dimensions
  object@Dimensions[["M"]] <- length(object@TrainData)
  object@Dimensions[["N"]] <- ncol(object@TrainData[[1]])
  object@Dimensions[["D"]] <- sapply(object@TrainData,nrow)
  object@Dimensions[["K"]] <- ncol(object@Expectations$Z)
  
  
  # Set view, sample, feature and factor names
  viewNames(object) <- as.character(names(object@TrainData))
  sampleNames(object) <- as.character(colnames(object@TrainData[[1]]))
  featureNames(object) <- lapply(object@TrainData, function(x) as.character(rownames(x)))
  factorNames(object) <- paste0("LF",as.character(1:object@Dimensions[["K"]]))
  
  # Add names to likelihood vector
  names(object@ModelOptions$likelihood) <- viewNames(object)
  
  # Rename covariates, including intercept
  # if (object@ModelOptions$learnIntercept == TRUE) factorNames(object) <- c("intercept",as.character(1:(object@Dimensions[["K"]]-1)))
  # if (!is.null(object@ModelOptions$covariates)) {
  #   if (object@ModelOptions$learnIntercept == TRUE) {
  #     factorNames(object) <- c("intercept", colnames(object@ModelOptions$covariates), as.character((ncol(object@ModelOptions$covariates)+1:(object@Dimensions[["K"]]-1-ncol(object@ModelOptions$covariates)))))
  #   } else {
  #     factorNames(object) <- c(colnames(object@ModelOptions$covariates), as.character((ncol(object@ModelOptions$covariates)+1:(object@Dimensions[["K"]]-1))))
  #   }
  # }
  
  
  # Rename factors if intercept is included
  if (object@ModelOptions$learnIntercept) {
    intercept_idx <- apply(object@Expectations$Z==1,2,all)
    nonconst_idx <- which(!intercept_idx)
    factornames <- factorNames(object)
    factornames[intercept_idx] <- "intercept"
    factornames[nonconst_idx] <- paste0("LF",as.character(1:length(nonconst_idx)))
    factorNames(object) <- factornames
  }
  
  # Parse factors: Mask passenger samples
  if(is.null(minR2)) minR2 <- object@TrainOptions$DropFactorThreshold
  object <- .detectPassengers(object, r2_threshold=minR2)

  # Parse factors: order factors in order of variance explained
  if (sortFactors == T) {
    r2 <- rowSums(calculateVarianceExplained(object)$R2PerFactor)
    order_factors <- order(r2, decreasing = T)
    object <- subsetFactors(object,order_factors)
    if (object@ModelOptions$learnIntercept==T) { 
      factorNames(object) <- c("intercept",paste0("LF",1:(object@Dimensions$K-1)))
    } else {
      factorNames(object) <- paste0("LF",c(1:object@Dimensions$K) )
    }
  }
  
  # Do quality control on the model
  qualityControl(object)
  
  return(object)
}

