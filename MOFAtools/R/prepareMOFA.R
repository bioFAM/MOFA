
#' @title prepareMOFA: Prepare an untrained MOFA object for running MOFA
#' @name prepareMOFA
#' @description Function to set the training and model options, produces .txt files that are used for python as input and 
#' creates an .sh file for calling MOFA with the specified options from the command line. These files are all stored in the specified directory.
#' @param object an untrained MOFA object
#' @param dir directory to store .txt and .sh files in
#' @param ModelOptions list of ModelOptions (see getDefaultModelOpts for what options can be set here). If none specified, default options are used.
#' @param TrainOptions list of TrainOptions (see getDefaultTrainOptions for what options can be set here). If none specified, default options are used.
#' @param outFile name of output file from Python MOFA
#' @param k number of latent factors to start with (default = 10)
#' @param MOFAdir directory of the MOFA Pyhton package installation
#' @details fill this
#' @return a untrained MOFA object with specified ModelOpts and TrainOpts 
#' @export

prepareMOFA<- function(object, dir, ModelOptions = NULL, TrainOptions = NULL, outFile="MOFAout", k=10, MOFAdir) {
  
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  if(!dir.exists(dir)){
    warning("Directory did not exist and was created")
    dir.create(dir)
  }
  
  message("Preparing input .txt files...")
  prepareMOFAInputFiles(object, dir)
  
  # STILL NEED SANITY CHECKS ON IINDIVIDUAL ARGUMENTS OF TrainOpts
  message("Setting training options...")
  if(is.null(TrainOptions)){
  object@TrainOpts <- getDefaultTrainOpts()
  message("Using default Training Options as none specified")
  } 
  else {
    if(!class(TrainOptions) == "list" & !all(names(TrainOptions) == names(getDefaultTrainOpts()))) 
      stop("'TrainOptions' are misspecified, use the list format provided by getDefaultTrainOpts()")
    object@TrainOpts <- TrainOptions
  }
  
  # STILL NEED SANITY CHECKS ON IINDIVIDUAL ARGUMENTS OF MODELOPTS
  message("Setting model options...")
  if(is.null(ModelOptions)){
    object@ModelOpts <- getDefaultModelOpts(object)
    message("Using default Model Options as none specified")
  } 
  else {
    if(!class(ModelOptions) == "list" & !all(names(ModelOptions) == names(getDefaultModelOpts(object, silent=T)))) 
      stop("'TrainOpts' are misspecified, use the list format provided by getDefaultModelOpts()")
    object@ModelOpts <- ModelOptions
  }
  
  message("Preparing run file...")
  prepareMOFARunFile(object, dir,outFile=outFile, k=k, MOFAdir =MOFAdir) 
  message("Use 'run.sh' in ", dir, " to run the MOFA model with specified training and model options")
  
  return(object)
  
}



#' @title prepareMOFAInputFiles: Write omics .txt files for input to Python
#' @name prepareMOFAInputFiles
#' @description Function to produce .txt files that are used for python as input. These files are all stored in the specified directory.
#' @param object an untrained MOFA object
#' @param dir directory to store .txt in
#' @details fill this
#' @return  
#' @export
#' 
prepareMOFAInputFiles <- function(object, dir) {

  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  if(!dir.exists(dir)){
    warning("Directory did not exist and was created")
    dir.create(dir)
  }

  for(objnm in names(object@TrainData)){
    write.table(t(object@TrainData[[objnm]]), file=file.path(dir, paste(objnm, ".txt", sep="")),
                row.names=TRUE, col.names=TRUE, quote=F)
    message("Writing", file.path(dir, paste(objnm, ".txt", sep="")))
    }

}

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

prepareMOFARunFile <- function(object, dir, outFile="MOFAout", k=10, MOFAdir){

cat(sprintf("#!/bin/bash\ninFolder='%s'\ninFiles=( %s)\noutFile=( '%s')\nlikelihoods=( %s)\nviews=( %s)\nschedule=( %s)\ntolerance=%f\nnostop=0\nntrials=%d\nncores=%d\niter=%d\nelbofreq=%d\nfactors=%d\nstartDrop=%d\nfreqDrop=%d\ndropNorm=%f\ndropR2=%f\nardW='mk'\nardZ=0\nlearnTheta=1\nscriptdir='%s'
\ncmd='python $scriptdir/template_run.py\n--inFiles ${inFiles[@]}\n--outFile $outFile\n--likelihoods ${likelihoods[@]}\n--views ${views[@]}\n--schedule ${schedule[@]}\n--ntrials $ntrials\n--ncores $ncores\n--iter $iter\n--elbofreq $elbofreq\n--startDrop $startDrop\n--freqDrop $freqDrop\n--tolerance $tolerance\n--factors $factors\n--ardW $ardW\n--dropNorm $dropNorm\n--dropR2 $dropR2\n--learnTheta%s'
            \neval $cmd",
          dir, paste(paste("'","$inFolder/",viewNames(object),".txt", "'", sep=""), collapse=" "),
          file.path(dir,outFile), 
          paste(object@ModelOpts$likelihood,collapse=" "), 
          paste(viewNames(object),collapse=" "),
          paste(object@ModelOpts$schedule,collapse=" "), 
          object@TrainOpts$tolerance, object@TrainOpts$trials,
          object@TrainOpts$cores,object@TrainOpts$maxiter,object@TrainOpts$elbofreq,
          k,  object@TrainOpts$startdrop, object@TrainOpts$freqdrop,object@TrainOpts$drop_by_norm,
          object@TrainOpts$drop_by_r2, MOFAdir, ifelse(object@ModelOpts$learnMean,"\n--learnMean","")), file = file.path(dir,"run.sh"))
 

}

#' @title getDefaultTrainOpts: Get default training options
#' @name getDefaultTrainOpts
#' @description Function to obtain default training options
#' @details fill this
#' @return  list with training options
#' @export

getDefaultTrainOpts <- function(){
  TrainOptions <- list(freqdrop = 999999,
                       drop_by_norm = NaN,
                       forceiter = 0,
                       trials = 1,
                       verbosity = 2,
                       elbofreq = 1,
                       drop_by_pvar = NaN,
                       savefreq = NaN,
                       maxiter = 2000,
                       cores = 1,
                       startdrop = 999999,
                       drop_by_cor = NaN,
                       tolerance = 0.01,
                       savefolder = NaN,
                       drop_by_r2 = 0.05)
  return(TrainOptions)
}

#' @title getDefaultModelOpts: Get default model options
#' @name getDefaultModelOpts
#' @param object  untrained MOFA object to get model options for
#' @param silent  boolean whether to print warnings
#' @description Function to obtain default model options
#' @details fill this
#' @return  list with training options
#' @export
#' 
getDefaultModelOpts <- function(object, silent=F){
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  if(!.hasSlot(object,"Dimensions") | length(object@Dimensions) == 0)
    stop("Dimensions of object need to be defined before getting ModelOpts")
  
  ModelOptions <- list(learnMean = TRUE,
                       likelihood = rep("gaussian", object@Dimensions[["M"]]),
                       schedule =c("Y", "SW","Z","AlphaW","Theta","Tau")
  )
  
  if(!silent) message("Using gaussian liklihoods for all views by default, change liklihood entry in ModelOpts to accomodate non-gaussian views")
  return(ModelOptions)
}
