
###########################
## Functions to run MOFA ##
###########################

#' @title prepareMOFARunFile: Write .sh files to run MOFA in Python 
#' @name prepareMOFARunFile
#' @description Function to produce .sh files to run MOFA in Python from the command line with specified options.
#' @param object an untrained MOFA object
#' @param dir directory to store .txt in
#' @param outFile name of output file from Python MOFA
#' @param k number of latent factors to start with (default = 10)
#' @param MOFAdir directory of the MOFA Pyhton package installation
#' @details fill this
#' @return  
#' @export
#' 

runMOFA <- function(object, DirOptions) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(all(c("tmpDir","mofaDir","outFile") %in% names(DirOptions)))
  
  # Prepare command
  command <- paste(sep=" ",
  "/Users/ricard/anaconda2/bin/python", paste0(DirOptions$mofaDir,"/run/template_run.py"),
  "--inFiles", paste(paste0(DirOptions$tmpDir, "/", viewNames(object), ".txt"), collapse = " "),
  "--outFile", DirOptions$outFile,
  "--views", paste(viewNames(object), collapse=" "),
  "--likelihoods", paste(object@ModelOpts$likelihood, collapse=" "),
  "--learnTheta", paste(object@ModelOpts$learnTheta, collapse=" "),
  "--initTheta", paste(object@ModelOpts$initTheta, collapse=" "),
  "--schedule", paste(object@ModelOpts$schedule, collapse=" "),
  "--ntrials", object@TrainOpts$trials,
  "--ncores", object@TrainOpts$cores,
  "--iter", object@TrainOpts$maxiter,
  "--elbofreq", object@TrainOpts$elbofreq,
  "--startDrop", object@TrainOpts$startdrop,
  "--freqDrop", object@TrainOpts$freqdrop,
  "--startSparsity", object@TrainOpts$startSparsity,
  "--dropNorm", object@TrainOpts$drop_by_norm,
  "--dropR2", object@TrainOpts$drop_by_r2,
  "--factors", object@ModelOpts$initialK,
  "--tolerance", object@TrainOpts$tolerance
  )
  if (!is.null(object@ModelOpts$covariates)) {
    command <- paste(command, sep=" ",
                     "--covariatesFile", file.path(DirOptions$tmpDir, "covariates.txt"),
                     "--scale_covariates", paste(object@ModelOpts$scale_covariates, collapse=" ")
                     )
  }
  if (object@TrainOpts$forceiter == T) { command <- paste(command, "--nostop", sep=" ") }
  if (object@ModelOpts$learnMean == T) { command <- paste(command, "--learnMean", sep=" ") }
  
  # Run motherfuckers!!!!
  system(command)
  
  # Load trained model
  object <- loadModel(DirOptions$outFile, object)
  # object <- loadModel(DirOptions$outFile)
  
  return(object)
}
