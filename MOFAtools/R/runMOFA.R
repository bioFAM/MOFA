
##############################################
## Functions to run MOFA from the R package ##
##############################################

#' @title runMOFA:
#' @name runMOFA
#' @description train a \code{\link{MOFAmodel}}
#' @param object an untrained \code{\link{MOFAmodel}}
#' @param DirOptions list with I/O options, should contain at least 'dataDir' where the input matrices as stored as .txt files and 'outFile' where the model is going to be stored as a .hdf5 file
#' @return a trained \code{\link{MOFAmodel}}
#' @export
runMOFA <- function(object, DirOptions) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(all(c("dataDir","outFile") %in% names(DirOptions)))
  
  # Prepare command
  command <- paste(sep=" ",
  # "mofa",
  "--inFiles", paste(paste0(DirOptions$dataDir, "/", viewNames(object), ".txt"), collapse = " "),
  "--header_cols --header_rows",
  "--outFile", DirOptions$outFile,
  "--views", paste(viewNames(object), collapse=" "),
  "--likelihoods", paste(object@ModelOpts$likelihood, collapse=" "),
  "--factors", object@ModelOpts$numFactors,
  "--iter", object@TrainOpts$maxiter,
  "--dropR2", object@TrainOpts$DropFactorThreshold,
  "--tolerance", object@TrainOpts$tolerance
  )
  if (!is.null(object@ModelOpts$covariates)) {
    command <- paste(command, sep=" ",
                     "--covariatesFile", file.path(DirOptions$dataDir, "covariates.txt"),
                     "--scale_covariates", rep(1,ncol(object@ModelOpts$covariates)))
  }
  if (object@ModelOpts$learnIntercept == T) { command <- paste(command, "--learnIntercept", sep=" ") }
  if (object@ModelOpts$sparsity == F) { command <- paste(command, "--learnTheta 0", sep=" ") }
  
  if (object@DataOpts$centerFeatures == T) { command <- paste(command, "--center_features", sep=" ") }
  if (object@DataOpts$scaleViews == T) { command <- paste(command, "--scale_views", sep=" ") }
  
  # Run!
  # system(command, ignore.stdout = F, ignore.stderr = T, wait=F)
  # system2(command="mofa", args=command, wait=F)
  system2(command="mofa", args=command, wait=T)
  
  # Load trained model
  object <- loadModel(DirOptions$outFile, object)
  
  return(object)
}
