
#' @title Initialize a \code{\link{MOFAmodel}} object with multi-omic data
#' @name createMOFAobject
#' @description Method to initialize a MOFAmodel object
#' @param data either MultiAssayExperiment or list of matrices with features as rows and samples as columns
#' @return an untrained MOFA model object
#' @export
createMOFAobject <- function(data) {
  
  if (class(data) == "MultiAssayExperiment") {
    message("Creating MOFA object from a MultiAssayExperiment object...")
    object <- .createMOFAobjectFromMAE(data)
  } else if (class(data) == "list") {
    message("Creating MOFA object from list of matrices, make sure that samples are columns and features are rows...")
    object <- .createMOFAobjectFromList(data)
  } else {
    stop("Error: input data has to be provided either as a list of matrices or as a MultiAssayExperiment object")
  }
  
  # Set dimensionalities
  object@Dimensions[["M"]] <- length(object@TrainData)
  object@Dimensions[["N"]] <- ncol(object@TrainData[[1]])
  object@Dimensions[["D"]] <- sapply(object@TrainData,nrow)
  object@Dimensions[["K"]] <- 0
  
  # Set view names
  if(!is.null(names(data))) {
    viewNames(object) <- names(data) 
  } else { 
    viewNames(object) <- paste("view",1:length(object@TrainData), sep="_")
    warning(paste0("View names are not specified in data, renaming them to: ",paste("view",1:length(object@TrainData), collapse=" "), "\n"))
  }
  
  print(object)
  
  return(object)
}

#' @title .createMOFAobjectFromMAE: initialise a MOFAmodel using a MultiAssayExperiment object
#' @name .createMOFAobjectFromMAE
#' @param data a \link{MultiAssayExperiment} object
#' @import MultiAssayExperiment
.createMOFAobjectFromMAE <- function(data) {

  # Initialise MOFA object
  object <- new("MOFAmodel")
  object@Status <- "untrained"
  object@InputData <- data
  
  # Re-arrange data for training in MOFA to matrices, fill in NAs and store in TrainData slot
  object@TrainData <- lapply(names(data), function(m) {
    
    # Extract general sample names
    primary <- unique(sampleMap(data)[,"primary"])
    
    # Extract view
    subdata <- assays(data)[[m]]
    
    # Rename view-specific sample IDs with the general sample names
    stopifnot(colnames(subdata)==sampleMap(data)[sampleMap(data)[,"assay"]==m,"colname"])
    colnames(subdata) <- sampleMap(data)[sampleMap(data)[,"assay"]==m,"primary"]
    
    # Fill subdata with NAs
    subdata_filled <- .subset_augment(subdata,primary)
    return(subdata_filled)
  })
  return(object)
}


# (Hidden) function to initialise a MOFAmodel object using a list of matrices
.createMOFAobjectFromList <- function(data) {
  
  # Initialise MOFA object
  object <- new("MOFAmodel")
  object@Status <- "untrained"
  
  # Fetch or assign sample names
  samples <- Reduce(union, lapply(data, colnames))
  if (is.null(samples)) {
    N <- unique(sapply(data,ncol))
    if (length(N)>1) { 
      stop("If the matrices have no column (samples) names that can be used to match the different views, all matrices must have the same number of columns")
    }
    samples <- as.character(1:N)
    for (m in 1:length(data)) { colnames(data[[m]]) <- samples }
  }
  
  object@TrainData <- lapply(data, function(view) .subset_augment(view, samples))
  return(object)
}




.subset_augment<-function(mat, pats) {
  pats <- unique(pats)
  mat <- t(mat)
  aug_mat<-matrix(NA, ncol=ncol(mat), nrow=length(pats))
  aug_mat<-mat[match(pats,rownames(mat)),,drop=FALSE]
  rownames(aug_mat)<-pats
  colnames(aug_mat)<-colnames(mat)
  return(t(aug_mat))
}