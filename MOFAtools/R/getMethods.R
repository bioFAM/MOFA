
################################################
## Get functions to fetch data from the model ##
################################################

#' @title getDimensions 
#' @name getDimensions
#' @description Extract dimensionalities from the model. 
#' @details K indicates the number of factors, D indicates the number of features, 
#' N indicates the (total) number of samples and M indicates the number of views.
#' @param object a \code{\link{MOFAmodel}} object.
#' @return list with the relevant dimensionalities of the model.
#'  N for the number of samples, M for the number of views, 
#'  D for the number of features of each view and K for the number of infered latent factors.
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' MOFAobject <- loadModel(filepath)
#' # get dimensions 
#' getDimensions(MOFAobject)

getDimensions <- function(object) {
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")  
  return(object@Dimensions)
}


#' @title getFactors
#' @name getFactors
#' @description Extract the latent factors from the model.
#' @param object a trained \code{\link{MOFAmodel}} object.
#' @param factors character vector with the factor name(s), or numeric vector with the factor index(es).
#' Default is "all".
#' @param as.data.frame logical indicating whether to return a long data frame instead of a matrix.
#' Default is \code{FALSE}.
#' @return By default it returns the latent factor matrix of dimensionality (N,K),
#'  where N is number of samples and K is number of factors. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE},
#'  returns a long-formatted data frame with columns (sample,factor,value).
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' MOFAobject <- loadModel(filepath)
#' # get factors as matrix
#' getFactors(MOFAobject, factors = 1:3)
#' # get factors as data.frame
#' head(getFactors(MOFAobject, factors = 1:5, as.data.frame = TRUE))

getFactors <- function(object, factors = "all", as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  
  # Get factors
  if (paste0(factors,collapse="") == "all") { 
    factors <- factorNames(object) 
  } else if (is.numeric(factors)) {
     factors <- factorNames(object)[factors]
  } else { 
    stopifnot(all(factors %in% factorNames(object))) 
  }

  # Collect factors
  Z <- getExpectations(object,"Z",as.data.frame)
  if (as.data.frame==FALSE) {
    Z <- Z[,factors, drop=FALSE]
  } else {
    Z <- Z[Z$factor %in% factors,]
  }

  return(Z)
}


#' @title getWeights
#' @name getWeights
#' @description Extract the weights from the model.
#' @param object a trained \code{\link{MOFAmodel}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es). 
#' Default is "all".
#' @param factors character vector with the factor name(s) or numeric vector with the factor index(es). \cr
#' Default is "all".
#' @param as.data.frame logical indicating whether to return a long data frame instead of a list of matrices. 
#' Default is \code{FALSE}.
#' @return By default it returns a list where each element is a loading matrix with dimensionality (D,K), 
#' where D is the number of features and K is the number of factors. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE},
#'  returns a long-formatted data frame with columns (view,feature,factor,value).
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' MOFAobject <- loadModel(filepath)
#' # get weights as a list of matrices
#' weightList <- getWeights(MOFAobject, view = "all", factors = 1:4)
#' # get weights as a data.frame
#' head(getWeights(MOFAobject, view = "Mutations", as.data.frame = TRUE))

getWeights <- function(object, views = "all", factors = "all", as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  
  # Get views and factors
  if (paste0(views,collapse="") == "all") { views <- viewNames(object) 
  } else { stopifnot(all(views %in% viewNames(object))) }

  # Get factors
  if (paste0(factors,collapse="") == "all") { 
    factors <- factorNames(object) 
  } else if (is.numeric(factors)) {
        factors <- factorNames(object)[factors]
    } else { stopifnot(all(factors %in% factorNames(object))) }
        
  # Fetch weights
  weights <- getExpectations(object,"W",as.data.frame)
  if (as.data.frame==TRUE) {
    weights <- weights[weights$view%in%views & weights$factor%in%factors, ]
  } else {
    weights <- lapply(views, function(m) weights[[m]][,factors,drop=FALSE])
    names(weights) <-  views
    # if (length(views)==1) { weights <- weights[[1]] }
  }
  return(weights)
}


#' @title getTrainData
#' @name getTrainData
#' @description Fetch the training data
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es). 
#' Default is "all".
#' @param features list of character vectors with the feature names or
#'  list of numeric vectors with the feature indices. 
#' Default is "all"
#' @param as.data.frame logical indicating whether to return a long data frame instead of a list of matrices.
#' Default is \code{FALSE}.
#' @return A list with one numeric matrix per view, containing the parsed data used to fit the MOFA model.
#' @details By default this function returns a list where each element
#'  is a data matrix with dimensionality (D,N) 
#' where D is the number of features and N is the number of samples. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE}, the function
#'  returns a long-formatted data frame with columns (view,feature,sample,value).
#' @export
#' @examples 
#' data("scMT_data")
#' MOFAobject <- createMOFAobject(scMT_data)
#' getTrainData(MOFAobject)
#' 
#' data("CLL_data")
#' MOFAobject <- createMOFAobject(CLL_data)
#' getTrainData(MOFAobject)

getTrainData <- function(object, views = "all", features = "all", as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  
  # Get views
  if (paste0(views,collapse="") == "all") { views <- viewNames(object) 
  } else { stopifnot(all(views %in% viewNames(object))) }
  
  # Get features
  if (is(features, "list")) {
    stopifnot(all(vapply(seq_along(features),
                         function(i) all(features[[i]] %in% featureNames(object)[[views[i]]]), logical(1))))
  } else {
    if (paste0(features,collapse="") == "all") { 
      features <- featureNames(object)[views]
    } else {
      stop("features not recognised, please read the documentation")
    }
  }
  
  # Fetch data
  trainData <- object@TrainData[views]
  trainData <- lapply(seq_along(trainData), function(m) trainData[[m]][features[[m]],,drop=FALSE])
  names(trainData) <- views
  
  # Convert to long data frame
  if (as.data.frame==TRUE) {
    tmp <- lapply(views, function(m) { tmp <- reshape2::melt(trainData[[m]])
    colnames(tmp) <- c("feature","sample","value"); tmp <- cbind(view=m,tmp); return(tmp) })
    trainData <- do.call(rbind,tmp)
    trainData[,c("view","feature","sample")] <- vapply(trainData[,c("view","feature","sample")], as.character, character(nrow(trainData)))
  }# else if ((length(views)==1) && (as.data.frame==FALSE)) {
  #  trainData <- trainData[[views]]
  #}
  
  return(trainData)
}


#' @title getImputedData
#' @name getImputedData
#' @description Function to get the imputed data. It requires the previous use of the
#'  \code{\link{impute}} method.
#' @param object a trained \code{\link{MOFAmodel}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es). 
#' Default is "all".
#' @param features list of character vectors with the feature names or
#'  list of numeric vectors with the feature indices. 
#' Default is "all"
#' @param as.data.frame logical indicating whether to return a long-formatted data frame
#'  instead of a list of matrices. 
#' Default is \code{FALSE}.
#' @return By default returns a list where each element is a matrix with dimensionality (D,N), 
#' where D is the number of features in this view and N is the number of samples. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE},
#'  returns a long-formatted data frame with columns (view,feature,sample,value).
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' MOFAobject <- loadModel(filepath)
#' # impute missing values
#' MOFAobject <- impute(MOFAobject)
#' # get imputations for a single view
#' imputedDrugs <- getImputedData(MOFAobject,view="Drugs")
#' head(imputedDrugs)

getImputedData <- function(object, views = "all", features = "all", as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  if (length(object@ImputedData)==0) {
    stop("No imputed data found. Please run impute(MOFAobject) first")
  }
  
  # Get views
  if (paste0(views,collapse="") == "all") { 
    views <- viewNames(object) 
  } else { 
    stopifnot(all(views %in% viewNames(object))) 
  }
  
  # Get features
  if (is(features, "list")) {
    stopifnot(all(vapply(seq_along(features), function(i) all(features[[i]] %in% featureNames(object)[[views[i]]]), logical(1))))
  } else {
    if (paste0(features,collapse="") == "all") { 
      features <- featureNames(object)[views]
    } else {
      stop("features not recognised, please read the documentation")
    }
  }
  
  # Fetch imputed data
  ImputedData <- object@ImputedData[views]
  ImputedData <- lapply(seq_along(ImputedData),
                        function(m) ImputedData[[m]][features[[m]],,drop=FALSE]) 
  names(ImputedData) <- views
  
  # Convert to long data frame
  if (as.data.frame==TRUE) {
    tmp <- lapply(views, function(m) { tmp <- reshape2::melt(ImputedData[[m]]) 
    colnames(tmp) <- c("feature","sample","value"); tmp <- cbind(view=m,tmp); return(tmp) })
    ImputedData <- do.call(rbind,tmp)
    ImputedData[,c("view","feature","sample")] <- vapply(ImputedData[,c("view","feature","sample")], as.character, character(nrow(ImputedData)))
  } 
  # else if ((length(views)==1) && (as.data.frame==FALSE)) {
  #   ImputedData <- ImputedData[[views]]
  # }
  
  return(ImputedData)
}

#' @name getCovariates
#' @title getCovariates
#' @description This function extracts covariates from the \code{colData}
#'  in the input \code{MultiAssayExperiment} object. \cr
#' Note that if you did not use \code{MultiAssayExperiment} to create
#'  your \code{\link{createMOFAobject}}, this function will not work.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param covariate names of the covariate
#' @return a vector containing the covariate
#' @import MultiAssayExperiment
#' @importFrom Biobase phenoData
#' @export
#' @examples
#' # Example on the CLL data
#' library(MultiAssayExperiment)
#' data("CLL_data")
#' data("CLL_covariates")
#' # Create MultiAssayExperiment object 
#' mae_CLL <- MultiAssayExperiment(CLL_data, colData=CLL_covariates)
#' MOFAobject  <- createMOFAobject(mae_CLL)
#' # Extract covariates from the colData of a MultiAssayExperiment
#' gender <- getCovariates(MOFAobject, "Gender")
#' diagnosis <- getCovariates(MOFAobject, "Diagnosis")
#' # Example on the scMT data
#' data("scMT_data")
#' MOFAobject  <- createMOFAobject(scMT_data)
#' # Extract covariates from the colData of a MultiAssayExperiment
#' culture <- getCovariates(MOFAobject, "culture")
#' # Extract covariates from the phenoData of the RNA assay
#' cdr <- getCovariates(MOFAobject, "cellular_detection_rate")
getCovariates <- function(object, covariate) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  if(!is(object@InputData, "MultiAssayExperiment")) {
    stop("To work with covariates, InputData has to be specified in form of a MultiAssayExperiment")
  }
  
  # Check that samples from MOFAobject and MultiAssayExperiment are consistent.
  samples <- sampleNames(object)
  if (!all(samples %in% rownames(colData(object@InputData)))) {
    stop("There are samples in the model which are not detected in the MultiAssayExperiment object")
  } else {
    mae <- object@InputData[,samples]
  }
  
  # Extract the covariate from colData or in the phenoData of specific assays
  out <- list()
  if (covariate %in% colnames(colData(mae))) { 
    out[[covariate]] <- colData(mae)[,covariate]
    names(out[[covariate]]) <- sampleNames(object)
  } else {
    for (view in viewNames(object)) {
      phenodata <- phenoData(mae@ExperimentList[[view]])
      if (covariate %in% colnames(phenodata)) {
        out[[paste(view,covariate,sep="_")]] <- phenodata[[covariate]]
        names(out[[paste(view,covariate,sep="_")]]) <- rownames(phenodata)
      }
    }
  }
  
  # Final sanity check
  if (length(out) == 0) {
      stop("Covariate not found in the colData or phenoData of the MultiAssayExperiment object")    
  } else if (length(out)>1) {
      stop("Covariate ambiguously found in the phenoData of multiple assays. Please, extract it manually.")    
  } else {
    return(out[[1]])
  }
  
}

#' @title getExpectations
#' @name getExpectations
#' @description Function to extract the expectations from the (variational) posterior
#'  distributions of a trained \code{\link{MOFAmodel}} object.
#' @param object a trained \code{\link{MOFAmodel}} object.
#' @param variable variable name: 'Z' for factors, 'W' for weights, 'Tau' for noise,
#' 'Y' for pseudodata, 'Theta' for feature-wise spike-and-slab sparsity,
#'  'AlphaW' for view and factor-wise ARD sparsity
#' @param as.data.frame logical indicating whether to output the result as a long data frame,
#'  default is \code{FALSE}.
#' @details Technical note: MOFA is a Bayesian model where each variable has a prior distribution
#'  and a posterior distribution. In particular, to achieve scalability we used the 
#'  variational inference framework, thus true posterior distributions are replaced
#'   by approximated variational distributions.
#' This function extracts the expectations of the variational distributions, 
#' which can be used as final point estimates to analyse the results of the model. \cr 
#' The priors and variational distributions of each variable are extensively
#'  described in the supplementary methods of the original paper.
#' @return the output varies depending on the variable of interest: \cr
#' \itemize{
#'  \item{"Z"}{a matrix with dimensions (samples,factors). 
#'  If \code{as.data.frame} is \code{TRUE}, a long-formatted data frame with columns (sample,factor,value)}
#'  \item{"W"}{a list of length (views) where each element is a matrix with dimensions (features,factors).
#'   If \code{as.data.frame} is \code{TRUE}, a long-formatted data frame with columns (view,feature,factor,value)}
#'  \item{"Y"}{a list of length (views) where each element is a matrix with dimensions (features,samples).
#'   If \code{as.data.frame} is \code{TRUE}, a long-formatted data frame with columns (view,feature,sample,value)}
#'  \item{"Theta"}{}
#'  \item{"Tau"}{}
#' }
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "CLL_model.hdf5", package = "MOFAtools")
#' MOFAobject <- loadModel(filepath)
#' # get expectations of Alpha as matrix
#' getExpectations(MOFAobject, variable="Alpha")

getExpectations <- function(object, variable, as.data.frame = FALSE) {
  
  # Sanity checks
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")
  stopifnot(variable %in% names(object@Expectations))
  
  # Get expectations in single matrix or list of matrices (for multi-view nodes)
  exp <- object@Expectations[[variable]]
  # if (variable=="Z") {
  #   exp <- object@Expectations$Z
  # } else {
  #   exp <- lapply(object@Expectations[[variable]], function(x) x$E)
  # }
  
  # Convert to long data frame
  if (as.data.frame==TRUE) {
    if (variable=="Z") {
      tmp <- reshape2::melt(exp)
      colnames(tmp) <- c("sample","factor","value")
    }
    else if (variable=="W") {
      tmp <- lapply(names(exp), function(m) { 
        tmp <- reshape2::melt(exp[[m]])
        colnames(tmp) <- c("feature","factor","value");
        tmp$view <- m
        tmp[c("view","feature","factor")] <- vapply(tmp[c("view","feature","factor")], as.character, character(nrow(tmp)))
        return(tmp) 
      })
      tmp <- do.call(rbind.data.frame,tmp)
    }
    else if (variable=="Y") {
      tmp <- lapply(names(exp), function(m) { 
        tmp <- reshape2::melt(exp[[m]])
        colnames(tmp) <- c("sample","feature","value")
        tmp$view <- m
        tmp[c("view","feature","factor")] <- vapply(tmp[c("view","feature","factor")], as.character, character(nrow(tmp)))
        return(tmp) 
      })
      tmp <- do.call(rbind,tmp)
    }
    else if (variable=="Tau") {
      stop("Not implemented")
      # tmp <- lapply(names(exp), function(m) { 
      #   data.frame(view=m, feature=names(exp[[m]]), value=unname(exp[[m]]))
      #   tmp[c("view","feature","factor")] <- vapply(tmp[c("view","feature","factor")], as.character, character(nrow(tmp)))
      #   return(tmp) 
      # })
      # tmp <- do.call(rbind,tmp)
    }
    else if (variable=="AlphaW") {
      tmp <- lapply(names(exp), function(m) { 
        tmp <- data.frame(view=m, factor=names(exp[[m]]), value=unname(exp[[m]]))
        tmp[c("view","feature","factor")] <- vapply(tmp[c("view","feature","factor")], as.character, character(nrow(tmp)))
        return(tmp) 
      })
      tmp <- do.call(rbind,tmp)
    }
    else if (variable=="Theta") {
      stop("Not implemented")
      # tmp <- lapply(names(exp), function(m) { tmp <- reshape2::melt(exp[[m]])
      # colnames(tmp) <- c("sample","feature","value")
      # tmp$view <- m; tmp[c("view","feature","factor")] <- vapply(tmp[c("view","feature","factor")], as.character, character(nrow(tmp)))
      # return(tmp) })
      # tmp <- do.call(rbind,tmp)
    }
    exp <- tmp
  }
  return(exp)
}


#' @title getELBO
#' @name getELBO
#' @description Extract the value of the ELBO statistics after model training.
#'  This can be useful for model selection.
#' @param object a \code{\link{MOFAmodel}} object.
#' @return value of the ELBO statistic at end of training
#' @importFrom utils tail
#' @export
#' @examples
#' # load a trained MOFAmodel object
#' filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAtools")
#' MOFAobject <- loadModel(filepath)
#' # get ELBO statistic
#' getELBO(MOFAobject)

getELBO <- function(object) {
  if (!is(object, "MOFAmodel")) stop("'object' has to be an instance of MOFAmodel")  
  return(tail(object@TrainStats$elbo,1))
}

