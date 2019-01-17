
############################################
## Functions to load a trained MOFA model ##
############################################


#' @title loading a trained MOFA model
#' @name loadModel
#' @description Method to load a trained MOFA model. \cr
#' The training of MOFA is done using a Python framework, and the model output
#'  is saved as an .hdf5 file, which has to be loaded in the R package.
#' @param file an hdf5 file saved by the MOFA python framework.
#' @param object either NULL (default) or an an existing untrained MOFA object.
#'  If NULL, the \code{\link{MOFAmodel}} object is created from the scratch.
#' @param sortFactors boolean indicating whether factors should be sorted
#'  by variance explained (default is TRUE)
#' @param minR2 minimum R2 threshold to call 'active' factors (default is 0.01).
#' @return a \code{\link{MOFAmodel}} model.
#' @importFrom rhdf5 h5read
#' @export
#' @examples
#' # path to the hdf5 file
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' # load the model
#' MOFAobject <- loadModel(filepath)
loadModel <- function(file, object = NULL, sortFactors = TRUE, minR2 = 0.01) {
    
    # message(paste0("Loading the following MOFA model: ", file))
    
    if (is.null(object)) object <- new("MOFAmodel")
    
    if(.hasSlot(object,"Status") & length(object@Status) !=0) {
      if (object@Status == "trained") {
          warning("The specified object is already trained, over-writing training output with new results!")
          object@Expectations <- list()
          object@FeatureIntercepts <- list()
          object@TrainStats <- list()
          object@ImputedData <- list()
          object@Dimensions <- list()
      }
      
    }
    
    # remove options to make sure the loaded model contains the options it was trained with
    object@DataOptions <- list()
    object@TrainOptions <- list()
    object@ModelOptions <- list()
    
    
    # Load expectations
    object@Expectations <- h5read(file,"expectations")
    
    # Change model status from "untrained" to "trained"
    object@Status <- "trained"
    
    # Load training statistics
    tryCatch( {
        object@TrainStats <- h5read(file, 'training_stats',read.attributes=TRUE);
        colnames(object@TrainStats$elbo_terms) <- attr(h5read(file,"training_stats/elbo_terms", read.attributes=TRUE),
        "colnames")
    }, error = function(x) { print("Training stats not found, not loading it...") })
    
    
    # Load training options
    tryCatch(object@TrainOptions <- as.list(h5read(file, 'training_opts',read.attributes=TRUE)),
    error = function(x) { print("Training opts not found, not loading it...") })
    
    # different names in R and python package (drop_by_r2 in python corresponds to DropFactorThreshold in R)
    if("drop_by_r2" %in% names(object@TrainOptions)) {
        object@TrainOptions$DropFactorThreshold <- object@TrainOptions$drop_by_r2
        object@TrainOptions$drop_by_r2 <- NULL
    }
    
    # Load model options
    tryCatch(object@ModelOptions <- as.list(h5read(file, 'model_opts',read.attributes=TRUE)),
    error = function(x) { print("Model opts not found, not loading it...") })
    object@ModelOptions$sparsity <- as.logical(object@ModelOptions$sparsity)
    
    # Load training data
    tryCatch( {
        TrainData <- h5read(file,"data")
        features <- h5read(file,"features")
        samples <- h5read(file,"samples")
        
        # Add feature names
        if (length(features)>0) {
            for (m in names(TrainData)) colnames(TrainData[[m]]) <- features[[m]]
        } else {
            for (m in names(TrainData)) colnames(TrainData[[m]]) <- paste0("feature_",1:ncol(TrainData[[m]]),"_",m )
        }
        
        # Add sample names
        if (length(samples)>0) {
            if(class(samples) != "list") samples <- list(samples=samples) # for compatibility with older verisons
            for (m in names(TrainData))
            rownames(TrainData[[m]]) <- samples[[1]]
        } else {
            for (m in names(TrainData))
            rownames(TrainData[[m]]) <- paste0("sample_",1:nrow(TrainData[[m]]))
        }
        
        # Transpose the data so that rows are features and columns are samples
        TrainData <- lapply(TrainData, t)
        
        # Replace NaN by NA in the training data
        for (m in names(TrainData)) {
            TrainData[[m]][is.nan(TrainData[[m]])] <- NA
        }
        
        # Store training data in the corresponding slot
        object@TrainData <- TrainData
        
    }, error = function(x) { print("Error loading the training data...") })
    
    # Sanity check on the order of the likelihoods
    if (!is.null(attr(TrainData,"likelihood"))) {
        lik <- attr(TrainData,"likelihood")
        if (!all(object@ModelOptions$likelihood == lik)) {
            object@ModelOptions$likelihood <- lik
            names(object@ModelOptions$likelihood) <- names(TrainData)
        }
    }
    
    # Collect feature-wise intercepts (for gaussian data stored during model preparation)
    object@FeatureIntercepts <- tryCatch( {
        h5read(file,"intercept")[names(object@TrainData)]
    }, error = function(x) {
        print("Could not load the feature intercepts. Your features will be centered to zero mean")
        list()
    })
    
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
    
    
    # Remove zero factors
    nonzero_factors <- which(apply(object@Expectations$Z, 2,
    function(z) !all(z==0)))
    if (length(nonzero_factors) < object@Dimensions$K) {
        message("Removing ", object@Dimensions$K - length(nonzero_factors),
        " factors that are constant zero from the model...")
    }
    object <- subsetFactors(object, nonzero_factors)
    factorNames(object) <- paste0("LF",as.character(1:length(nonzero_factors)))
    
    
    # Parse factors: Mask passenger samples
    if(is.null(minR2)) minR2 <- object@TrainOptions$DropFactorThreshold
    object <- .detectPassengers(object, r2_threshold=minR2)
    
    # Parse factors: order factors in order of variance explained
    if (sortFactors == TRUE) {
        r2 <- rowSums(calculateVarianceExplained(object)$R2PerFactor)
        order_factors <- order(r2, decreasing = TRUE)
        object <- subsetFactors(object,order_factors)
        factorNames(object) <- paste0("LF",c(1:object@Dimensions$K) )
    }
    
    # Do quality control on the model
    qualityControl(object)
    
    return(object)
}

