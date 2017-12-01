
##############################################
## Functions to run MOFA from the R package ##
##############################################

#' @title runMOFA:
#' @name runMOFA
#' @description train a \code{\link{MOFAmodel}}
#' @param object an untrained \code{\link{MOFAmodel}}
#' @param DirOptions list with I/O options, should contain at least 'dataDir' where the input matrices as stored as .txt files and 'outFile' where the model is going to be stored as a .hdf5 file
#' @param ... Extra options to add to the mofa command
#' @return a trained \code{\link{MOFAmodel}}
#' @export
runMOFA <- function(object, DirOptions, ..., mofaPath="mofa") {
  
  # Sanity checks
  if (! is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  stopifnot(all(c("dataDir","outFile") %in% names(DirOptions)))

  arglist <- list(
    inFiles = paste0(DirOptions$dataDir, "/", viewNames(object), ".txt"),
    header_cols = TRUE,
    header_rows = TRUE,
    delimiter = object@DataOpts$delimiter,
    outFile = DirOptions$outFile,
    views = viewNames(object),
    likelihoods = object@ModelOpts$likelihood,
    factors = object@ModelOpts$numFactors,
    iter = object@TrainOpts$maxiter,
    dropR2 =  object@TrainOpts$DropFactorThreshold,
    tolerance = object@TrainOpts$tolerance
  )

  # Setting the below arguments to NULL doesn't actually add them to
  # the argument list, but reserves that argument name to prevent
  # extra.arglist from using it.
  if (!is.null(object@ModelOpts$covariates)) {
    arglist$covariatesFile <- file.path(DirOptions$dataDir, "covariates.txt")
    arglist$scale_covariates <- rep(1,ncol(object@ModelOpts$covariates))
  } else {
    arglist$covariatesFile <- NULL
    arglist$scale_covariates <- NULL
  }
  arglist$learnIntercept <- as.logical(object@ModelOpts$learnIntercept)
  if (! object@ModelOpts$sparsity) {
    arglist$learnTheta <- 0
  } else {
    arglist$learnTheta <- NULL
  }

  arglist$center_features <- as.logical(object@DataOpts$centerFeatures)
  arglist$scale_views <- as.logical(object@DataOpts$scaleViews)
  if (object@DataOpts$removeIncompleteSamples == T) { command <- paste(command, "--RemoveIncompleteSamples", sep=" ") }
  arglist$verbose <- as.logical(object@TrainOpts$verbose)

  extra.arglist <- list(...)
  if (any(is.na(names(extra.arglist)) | names(extra.arglist) == "")) {
    stop("All extra options must be named")
  }

  # Remove leading "--" from extra arg names if present (it will be
  # added back later)
  names(arglist) <- sub("^--", "", names(arglist))
  conflicting.argnames <- intersect(names(extra.arglist), names(arglist))
  if (length(conflicting.argnames) > 0)
    stop(paste0("You cannot pass the following arguments as extra options to runMOFA: ",
      deparse(conflicting.argnames)))

  # No conflicting argument names,
  arglist <- c(arglist, extra.arglist)

  argv <- character(0)
  for (argname in names(arglist)) {
    argval <- arglist[[argname]]
    argname <- paste0("--", argname)

    if (is.null(argval)) {
      # Placeholder option; don't add it
    }
    if (is.logical(argval)) {
      # Flag option
      if (length(argval) != 1) {
        stop(paste("Invalid argument value:", deprase(argval)))
      } else if (argval == FALSE || is.na(argval)) {
        # Unset flag: don't add it
      } else if (argval == TRUE) {
        # Set flag: add it
        argv <- c(argv, argname)
      }
    } else {
      # Option with arguments: add the option followed by it args
      argv <- c(argv, argname, argval)
    }
  }
  argv <- unlist(argv)

  if (length(mofaPath) != 1) stop("Invalid mofaPath")

  # If output already exists, remove it
  if (file.exists(DirOptions$outFile)) {
    if (arglist$verbose) {
      message("Deleting old output file")
      }
    file.remove(DirOptions$outFile)
  }

  if (arglist$verbose) {
    message("Running MOFA command: ", paste(collapse=" ", shQuote(c(mofaPath, argv))))
  }
  # Run!
  exitcode <- system2(command=mofaPath, args=shQuote(argv), wait=T)
  if (exitcode != 0) {
    stop(paste("mofa command failed with exit code", exitcode))
  }
  
  # Load trained model
  object <- loadModel(DirOptions$outFile, object)
  
  return(object)
}
