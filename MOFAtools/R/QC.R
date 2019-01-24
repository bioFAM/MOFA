#' @title qualityControl
#' @name qualityControl
#' @description Function to do quality control on a \code{\link{MOFAmodel}} object.
#' @param object a trained \code{\link{MOFAmodel}} object.
#' @param verbose logical indicating whether to generate a verbose output.
#' @return none
#' @export
#' @examples 
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAdata")
#' MOFAobject <- loadModel(filepath)
#' # do quality control
#' qualityControl(MOFAobject, verbose=TRUE)

qualityControl <- function(object, verbose = FALSE) {
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  if (object@Status != "trained") stop("This function only works in a trained MOFAmodel")
  
  # Check that the model has view names
  if (verbose) message("Checking view names...")
  stopifnot(!is.null(viewNames(object)))
  
  # Check that the model has sample names
  if (verbose) message("Checking sample names...")
  stopifnot(!is.null(sampleNames(object)))
  
  # Check that the model has feature names
  if (verbose) message("Checking feature names...")
  stopifnot(!is.null(featureNames(object)))

  # Check that the model has the right node names
  if (verbose) message("Checking nodes...")
  stopifnot(identical(sort(c("W","Z","Theta","Tau","Alpha","Y")), sort(names(object@Expectations))))
  
  # Check that all expectations are the correct object
  if (verbose) message("Checking expectations...")
  stopifnot(is.matrix(object@Expectations$Z))
  stopifnot(is.list(object@Expectations$W))
  stopifnot(all(vapply(object@Expectations$W, is.matrix, logical(1))))
  stopifnot(is.list(object@Expectations$Y))
  stopifnot(all(vapply(object@Expectations$Y, is.matrix, logical(1))))
  # stopifnot(is.list(object@Expectations$Theta))
  # stopifnot(all(vapply(object@Expectations$Theta, is.matrix, logical(1))))
  stopifnot(is.list(object@Expectations$Tau))
  stopifnot(all(vapply(object@Expectations$Tau, is.numeric, logical(1))))
  stopifnot(is.list(object@Expectations$Alpha))
  stopifnot(all(vapply(object@Expectations$Alpha, is.numeric, logical(1))))
  
  # Check that the dimensionalities match
  # TO-DO...
  if (verbose) message("Checking dimensionalities...")
  
  # Check that there are no features with complete missing values
  if (verbose) message("Checking there are no features with complete missing values...")
  for (view in viewNames(object)) {
    # FIX THIS
    if (!all(apply(object@TrainData[[view]],1, function(x) mean(is.na(x))) < 1, na.rm=TRUE)) {
      print("Warning: you have features which only contain missing values, consider removing them...")
    }
  }
  
  # Check that there are no features with zero variance
  if (verbose) message("Checking there are no features with zero variance...")
  for (view in viewNames(object)) {
    if (!all(apply(object@TrainData[[view]],1,var,na.rm=TRUE) > 0, na.rm=TRUE)) {
      print("Warning: you have features with zero variance, consider removing them...")
    }
  }
  
  # Check that the likelihoods match the data distribution
  if (verbose) message("Checking likelihooods...")
  predicted_lik <- .inferLikelihoods(object)
  for (view in viewNames(object)) {
    lk <- object@ModelOptions$likelihood[view]
    if (lk != predicted_lik[view])
      message(sprintf("Warning, view %s should follow a %s distribution rather than %s ",
                      view, predicted_lik[view], lk))
  }
  
}
