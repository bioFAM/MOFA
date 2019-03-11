
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
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAdata")
#' # load the model
#' MOFAobject <- loadModel(filepath)
loadModel <- function(file, object = NULL, sortFactors = TRUE, minR2 = 0.01) {
    
    # message(paste0("Loading the following MOFA model: ", file))
    
    if (is.null(object)) object <- new("MOFAmodel")
    
    if(.hasSlot(object,"Status") & length(Status(object)) !=0) {
      if (Status(object) == "trained") {
          warning("The specified object is already trained, over-writing training output with new results!")
          Expectations(object) <- list()
          FeatureIntercepts(object) <- list()
          TrainStats(object) <- list()
          ImputedData(object) <- list()
          Dimensions(object) <- list()
      }
      
    }
    
    # remove options to make sure the loaded model contains the options it was trained with
    DataOptions(object) <- list()
    TrainOptions(object) <- list()
    ModelOptions(object) <- list()
    
    
    # Load expectations
    Expectations(object) <- h5read(file,"expectations")
    
    # Load training statistics
    tryCatch( {
        TrainStats(object) <- h5read(file, 'training_stats',read.attributes=TRUE);
        colnames(TrainStats(object)[["elbo_terms"]]) <- attr(h5read(file,"training_stats/elbo_terms", read.attributes=TRUE),
        "colnames")
    }, error = function(x) { print("Training stats not found, not loading it...") })
    
    
    # Load training options
    tryCatch(TrainOptions(object) <- as.list(h5read(file, 'training_opts',read.attributes=TRUE)),
    error = function(x) { print("Training opts not found, not loading it...") })
    
    # different names in R and python package (drop_by_r2 in python corresponds to DropFactorThreshold in R)
    if("drop_by_r2" %in% names(TrainOptions(object))) {
        TrainOptions(object)[["DropFactorThreshold"]] <- TrainOptions(object)[["drop_by_r2"]]
        TrainOptions(object)[["drop_by_r2"]] <- NULL
    }
    
    # Load model options
    tryCatch(ModelOptions(object) <- as.list(h5read(file, 'model_opts',read.attributes=TRUE)),
    error = function(x) { print("Model opts not found, not loading it...") })
    ModelOptions(object)[["sparsity"]] <- as.logical(ModelOptions(object)[["sparsity"]])
    
    # Load training data
    tryCatch( {
        TrainData <- h5read(file,"data")
        features <- h5read(file,"features")
        samples <- h5read(file,"samples")
        
        # Add feature names
        if (length(features)>0) {
            for (m in names(TrainData)) colnames(TrainData[[m]]) <- features[[m]]
        } else {
            for (m in names(TrainData)) colnames(TrainData[[m]]) <- paste0("feature_", seq_len(ncol(TrainData[[m]])),"_",m )
        }
        
        # Add sample names
        if (length(samples)>0) {
            if(!is(samples, "list")) samples <- list(samples=samples) # for compatibility with older verisons
            for (m in names(TrainData))
            rownames(TrainData[[m]]) <- samples[[1]]
        } else {
            for (m in names(TrainData))
            rownames(TrainData[[m]]) <- paste0("sample_",seq_len(nrow(TrainData[[m]])))
        }
        
        # Transpose the data so that rows are features and columns are samples
        TrainData <- lapply(TrainData, t)
        
        # Replace NaN by NA in the training data
        for (m in names(TrainData)) {
            TrainData[[m]][is.nan(TrainData[[m]])] <- NA
        }
        
        # Store training data in the corresponding slot
        TrainData(object) <- TrainData
        
    }, error = function(x) { print("Error loading the training data...") })
    
    # Sanity check on the order of the likelihoods
    if (!is.null(attr(TrainData,"likelihood"))) {
        lik <- attr(TrainData,"likelihood")
        if (!all(ModelOptions(options)[["likelihood"]] == lik)) {
            ModelOptions(options)[["likelihood"]] <- lik
            names(ModelOptions(options)[["likelihood"]]) <- names(TrainData)
        }
    }
    
    # Collect feature-wise intercepts (for gaussian data stored during model preparation)
    FeatureIntercepts(object) <- tryCatch( {
        h5read(file,"intercept")[names(TrainData(object))]
    }, error = function(x) {
        print("Could not load the feature intercepts. Your features will be centered to zero mean")
        list()
    })
    
    # Update old models
    object <- .updateOldModel(object)
    
    # Load dimensions
    Dimensions(object) <- list()
    Dimensions(object)[["M"]] <- length(TrainData(object))
    Dimensions(object)[["N"]] <- ncol(TrainData(object)[[1]])
    Dimensions(object)[["D"]] <- vapply(TrainData(object), nrow, numeric(1))
    Dimensions(object)[["K"]] <- ncol(Expectations(object)[["Z"]])
    
    # Set view, sample, feature and factor names
    viewNames(object) <- as.character(names(TrainData(object)))
    sampleNames(object) <- as.character(colnames(TrainData(object)[[1]]))
    featureNames(object) <- lapply(TrainData(object), function(x) as.character(rownames(x)))
    factorNames(object) <- paste0("LF",as.character(seq_len(Dimensions(object)[["K"]])))
    
    # Add names to likelihood vector
    names(ModelOptions(object)[["likelihood"]]) <- viewNames(object)
    
    
    # Remove zero factors
    nonzero_factors <- which(apply(Expectations(object)[["Z"]], 2,
    function(z) !all(z==0)))
    if (length(nonzero_factors) < Dimensions(object)[["K"]]) {
        message("Removing ", Dimensions(object)[["K"]] - length(nonzero_factors),
        " factors that are constant zero from the model...")
    }
    object <- subsetFactors(object, nonzero_factors)
    factorNames(object) <- paste0("LF",as.character(seq_len(length(nonzero_factors))))
    
    
    # Parse factors: Mask passenger samples
    if(is.null(minR2)) minR2 <- TrainOptions(object)[["DropFactorThreshold"]]
    object <- .detectPassengers(object, r2_threshold=minR2)
    
    # Parse factors: order factors in order of variance explained
    if (sortFactors == TRUE) {
        r2 <- rowSums(calculateVarianceExplained(object)$R2PerFactor)
        order_factors <- order(r2, decreasing = TRUE)
        object <- subsetFactors(object, order_factors)
        factorNames(object) <- paste0("LF", seq_len(Dimensions(object)[["K"]]))
    }
    
    # Change model status from "untrained" to "trained"
    Status(object) <- "trained"
    
    # Do quality control on the model
    qualityControl(object)
    
    return(object)
}

