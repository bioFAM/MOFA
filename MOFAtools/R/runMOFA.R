
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
  
  # cat(sprintf("#!/bin/bash\ninFolder='%s'\ninFiles=( %s)\noutFile=( '%s')\nlikelihoods=( %s)\nviews=( %s)\nschedule=( %s)\ntolerance=%f\nnostop=0\nntrials=%d\nncores=%d\niter=%d\nelbofreq=%d\nfactors=%d\nstartDrop=%d\nfreqDrop=%d\ndropNorm=%f\ndropR2=%f\nardW='mk'\nardZ=0\nlearnTheta=1\nscriptdir='%s'
  #             \ncmd='python $scriptdir/template_run.py\n--inFiles ${inFiles[@]}\n--outFile $outFile\n--likelihoods ${likelihoods[@]}\n--views ${views[@]}\n--schedule ${schedule[@]}\n--ntrials $ntrials\n--ncores $ncores\n--iter $iter\n--elbofreq $elbofreq\n--startDrop $startDrop\n--freqDrop $freqDrop\n--tolerance $tolerance\n--factors $factors\n--ardW $ardW\n--dropNorm $dropNorm\n--dropR2 $dropR2\n--learnTheta%s'
  #             \neval $cmd",
  #             dir, paste(paste("'","$inFolder/",viewNames(object),".txt", "'", sep=""), collapse=" "),
  #             file.path(dir,outFile), 
  #             paste(object@ModelOpts$likelihood,collapse=" "), 
  #             paste(viewNames(object),collapse=" "),
  #             paste(object@ModelOpts$schedule,collapse=" "), 
  #             object@TrainOpts$tolerance, object@TrainOpts$trials,
  #             object@TrainOpts$cores,object@TrainOpts$maxiter,object@TrainOpts$elbofreq,
  #             k,  object@TrainOpts$startdrop, object@TrainOpts$freqdrop,object@TrainOpts$drop_by_norm,
  #             object@TrainOpts$drop_by_r2, MOFAdir, ifelse(object@ModelOpts$learnMean,"\n--learnMean","")), file = file.path(dir,"run.sh"))
  # 
  command <- paste(sep=" ",
  "python", paste0(DirOptions$mofaDir,"/run/template_run.py"),
  "--inFiles", paste(paste0(DirOptions$tmpDir, "/", viewNames(object), ".txt"), collapse = " "),
  "--outFile", DirOptions$outFile,
  "--views", paste(viewNames(object), collapse=" "),
  "--likelihoods", paste(object@ModelOpts$likelihood, collapse=" "),
  "--schedule", paste(object@ModelOpts$schedule, collapse=" "),
  "--ntrials", object@TrainOpts$trials,
  "--ncores", object@TrainOpts$cores,
  "--iter", object@TrainOpts$maxiter,
  "--elbofreq", object@TrainOpts$elbofreq,
  "--startDrop", object@TrainOpts$startdrop,
  "--freqDrop", object@TrainOpts$freqdrop,
  "--dropNorm", object@TrainOpts$drop_by_norm,
  "--dropR2", object@TrainOpts$drop_by_r2,
  "--factors", object@ModelOpts$initialK,
  "--tolerance", object@TrainOpts$tolerance
  )
  
  if (!is.null(object@ModelOpts$covariatesFile)) {
    command <- paste(command, "--covariatesFile", object@ModelOpts$covariatesFile, sep=" ")
  }
  if (object@ModelOpts$learnTheta == T) {
    command <- paste(command, "--learnTheta", sep=" ")
  }
  if (object@TrainOpts$forceiter == T) {
    command <- paste(command, "--nostop", sep=" ")
  }
  if (object@ModelOpts$learnMean == T) {
    command <- paste(command, "--learnMean", sep=" ")
  }
  
  system(command)
  
  # Load trained model
  object <- loadModel(DirOptions$outFile, object)
  
  return(object)
}
