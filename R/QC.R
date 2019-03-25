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
  if (Status(object) != "trained") stop("This function only works in a trained MOFAmodel")
  
  # Check that the model has view names
  if (verbose) message("Checking view names...")
  stopifnot(!is.null(viewNames(object)))
  
  # Check that the model has sample names
  if (verbose) message("Checking sample names...")
  stopifnot(!is.null(sampleNames(object)))
  
  # Check that the model has feature names
  if (verbose) message("Checking feature names...")
  stopifnot(!is.null(featureNames(object)))

  # NOW DONE IN CLASS VALIDITY CHECK
  # # Check that the model has the right node names
  # if (verbose) message("Checking nodes...")
  # stopifnot(identical(sort(c("W","Z","Theta","Tau","Alpha","Y")), sort(names(Expectations(object)))))
  # 
  # # Check that all expectations are the correct object
  # if (verbose) message("Checking expectations...")
  # stopifnot(is.matrix(Expectations(object)[["Z"]]))
  # stopifnot(is.list(Expectations(object)[["W"]]))
  # stopifnot(all(vapply(Expectations(object)[["W"]], is.matrix, logical(1))))
  # stopifnot(is.list(Expectations(object)[["Y"]]))
  # stopifnot(all(vapply(Expectations(object)[["Y"]], is.matrix, logical(1))))
  # stopifnot(is.list(Expectations(object)[["Tau"]]))
  # stopifnot(all(vapply(Expectations(object)[["Tau"]], is.numeric, logical(1))))
  # stopifnot(is.list(Expectations(object)[["Alpha"]]))
  # stopifnot(all(vapply(Expectations(object)[["Alpha"]], is.numeric, logical(1))))
  # 
  # Check that the dimensionalities match
  if (verbose) message("Checking dimensionalities...")
    stopifnot(length(Expectations(object)[["Alpha"]]) == getDimensions(object)[["M"]])
    stopifnot(length(Expectations(object)[["W"]]) == getDimensions(object)[["M"]])
    stopifnot(length(Expectations(object)[["Y"]]) == getDimensions(object)[["M"]])
    stopifnot(length(Expectations(object)[["Theta"]]) == getDimensions(object)[["M"]])
    stopifnot(all(vapply(Expectations(object)[["W"]], dim, numeric(2)) == rbind(getDimensions(object)[["D"]],
                                                                getDimensions(object)[["K"]])))
     stopifnot(all(vapply(Expectations(object)[["Y"]], dim, numeric(2)) == rbind(getDimensions(object)[["N"]],
                                                                getDimensions(object)[["D"]])))
     stopifnot(all(vapply(Expectations(object)[["Alpha"]], length, numeric(1)) == getDimensions(object)[["K"]]))
     stopifnot(all(vapply(Expectations(object)[["Theta"]], length, numeric(1)) == getDimensions(object)[["K"]]))
     stopifnot(ncol(Expectations(object)[["Z"]]) == getDimensions(object)[["K"]])
     stopifnot(nrow(Expectations(object)[["Z"]]) == getDimensions(object)[["N"]])


  # Check that there are no features with complete missing values
  if (verbose) message("Checking there are no features with complete missing values...")
  for (view in viewNames(object)) {
    # FIX THIS
    if (!all(apply(TrainData(object)[[view]],1, function(x) mean(is.na(x))) < 1, na.rm=TRUE)) {
      print("Warning: you have features which only contain missing values, consider removing them...")
    }
  }
  
  # Check that there are no features with zero variance
  if (verbose) message("Checking there are no features with zero variance...")
  for (view in viewNames(object)) {
    if (!all(apply(TrainData(object)[[view]],1,var,na.rm=TRUE) > 0, na.rm=TRUE)) {
      print("Warning: you have features with zero variance, consider removing them...")
    }
  }
  
  # Check that the likelihoods match the data distribution
  if (verbose) message("Checking likelihooods...")
  predicted_lik <- .inferLikelihoods(object)
  for (view in viewNames(object)) {
    lk <- ModelOptions(object)[["likelihood"]][view]
    if (lk != predicted_lik[view])
      message(sprintf("Warning, view %s should follow a %s distribution rather than %s ",
                      view, predicted_lik[view], lk))
  }
  
}
