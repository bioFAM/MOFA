
################################################
## Get functions to fetch data from the model ##
################################################

#' @name getDimensions
#' @title Extract dimensionalities from the model. 
#' @description K indicates the number of factors, D indicates the number of features, N indicates the (total) number of samples and M indicates the number of views.
#' @param object a \code{\link{MOFAmodel}} object.
#' @export
getDimensions <- function(object) {
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")  
  return(object@Dimensions)
}


#' @name getFactors
#' @title Extract the latent factors from the model
#' @param object a \code{\link{MOFAmodel}} object.
#' @param factors character vector with the factor name(s), or numeric vector with the factor index(es). Default is "all".
#' @param include_intercept boolean indicating where to include the intercept term of the model, if present. Default is TRUE.
#' @param as.data.frame boolean indicating whether to return a long data frame instead of a matrix, default is FALSE.
#' @return By default returns the latent factor matrix of dimensionality (N,K), where N is number of samples and K is number of factors. \cr
#' Alternatively, if as.data.frame is TRUE, returns a long-formatted data frame with columns (sample,factor,value).
#' @export
#' 
getFactors <- function(object, factors = "all", as.data.frame = FALSE, include_intercept = TRUE) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")  
  
  # Get factors
  if (paste0(factors,collapse="") == "all") { factors <- factorNames(object) } else { stopifnot(all(factors %in% factorNames(object))) }

  # Collect factors
  Z <- getExpectations(object,"Z",as.data.frame)
    if (as.data.frame==FALSE) {
      Z <- Z[,factors, drop=FALSE]
    } else {
      Z <- Z[Z$factor %in% factors,]
    }

  # Remove intercept
  if (include_intercept == FALSE) {
    if (as.data.frame==FALSE) {
      if ("intercept" %in% colnames(Z)) Z <- Z[,colnames(Z)!="intercept"]
    } else {
      if ("intercept" %in% unique(Z$factor)) Z <- Z[Z$factor!="intercept",]
    }
  }
  return(Z)
}


#' @rdname getWeights
#' @name getWeights
#' @title Extract the weights from the model
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es). Default is "all".
#' @param factors character vector with the factor name(s) or numeric vector with the factor index(es). Default is "all".
#' @param as.data.frame boolean indicating whether to return a long data frame instead of a list of matrices. Default is FALSE.
#' @return By default returns a list where each element is a loading matrix with dimensionality (D,K), where D is the number of features in this view and K is the number of factors. \cr
#' Alternatively, if as.data.frame is TRUE, returns a long-formatted data frame with columns (view,feature,factor,value).
#' @export
#' 
getWeights <- function(object, views = "all", factors = "all", as.data.frame = FALSE) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # Get views and factors
  if (paste0(views,collapse="") == "all") { views <- viewNames(object) } else { stopifnot(all(views %in% viewNames(object))) }
  if (paste0(factors,collapse="") == "all") { factors <- factorNames(object) } else { stopifnot(all(factors %in% factorNames(object))) }
  
  # Fetch weights
  weights <- getExpectations(object,"SW",as.data.frame)
  if (as.data.frame==T) {
    weights <- weights[weights$view%in%views & weights$factor%in%factors, ]
  } else {
    weights <- lapply(views, function(m) weights[[m]][,factors,drop=F])
    # if (length(views)==1) { weights <- weights[[1]] }
    names(weights) <-  views
  }
  return(weights)
}


#' @rdname getTrainData
#' @name getTrainData
#' @title Fetch the training data
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es). Default is "all".
#' @param features list of character vectors with the feature names or list of numeric vectors with the feature indices. Default is "all"
#' @param as.data.frame boolean indicating whether to return a long data frame instead of a list of matrices. Default is FALSE.
#' @details By default this function returns a list where each element is a data matrix with dimensionality (D,N) where D is the number of features and N is the number of samples. \cr
#' Alternatively, if \code{as.data.frame} is \code{TRUE}, the function returns a long-formatted data frame with columns (view,feature,sample,value).
#' @export
getTrainData <- function(object, views = "all", features = "all", as.data.frame = F) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # Get views
  if (paste0(views,collapse="") == "all") { views <- viewNames(object) } else { stopifnot(all(views %in% viewNames(object))) }
  
  # Get features
  if (class(features)=="list") {
    stopifnot(all(sapply(1:length(features), function(i) all(features[[i]] %in% featureNames(object)[[views[i]]]))))
  } else {
    if (paste0(features,collapse="") == "all") { 
      features <- featureNames(object)[views]
    } else {
      stop("features not recognised, please read the documentation")
    }
  }
  
  # Fetch data
  trainData <- object@TrainData[views]
  trainData <- lapply(1:length(trainData), function(m) trainData[[m]][features[[m]],,drop=F]); names(trainData) <- views
  
  # Convert to long data frame
  if (as.data.frame==T) {
    tmp <- lapply(views, function(m) { tmp <- reshape2::melt(trainData[[m]]); colnames(tmp) <- c("feature","sample","value"); tmp <- cbind(view=m,tmp); return(tmp) })
    trainData <- do.call(rbind,tmp)
    trainData[,c("view","feature","sample")] <- sapply(trainData[,c("view","feature","sample")], as.character)
  }# else if ((length(views)==1) && (as.data.frame==F)) {
  #  trainData <- trainData[[views]]
  #}
  
  return(trainData)
}


#' @title Fetch the imputed data
#' @name getImputedData
#' @rdname getImputedData
#' @description Function to collect the imputed data. It requires the previous use of the \code{\link{imputeMissing}} method.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views character vector with the view name(s), or numeric vector with the view index(es). Default is "all".
#' @param features list of character vectors with the feature names or list of numeric vectors with the feature indices. Default is "all"
#' @param as.data.frame boolean indicating whether to return a long-formatted data frame instead of a list of matrices. Default is FALSE.
#' @return By default returns a list where each element is a matrix with dimensionality (D,N), where D is the number of features in this view and N is the number of samples. \cr
#' Alternatively, if as.data.frame is TRUE, returns a long-formatted data frame with columns (view,feature,sample,value).
#' @export
getImputedData <- function(object, views = "all", features = "all", as.data.frame = FALSE) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # Get views
  if (paste0(views,collapse="") == "all") { views <- viewNames(object) } else { stopifnot(all(views %in% viewNames(object))) }
  
  # Get features
  if (class(features)=="list") {
    stopifnot(all(sapply(1:length(features), function(i) all(features[[i]] %in% featureNames(object)[[views[i]]]))))
  } else {
    if (paste0(features,collapse="") == "all") { 
      features <- featureNames(object)[views]
    } else {
      stop("features not recognised, please read the documentation")
    }
  }
  
  # Fetch data
  ImputedData <- object@ImputedData[views]
  ImputedData <- lapply(1:length(ImputedData), function(m) ImputedData[[m]][features[[m]],,drop=F]); names(ImputedData) <- views
  
  # Convert to long data frame
  if (as.data.frame==T) {
    tmp <- lapply(views, function(m) { tmp <- reshape2::melt(ImputedData[[m]]); colnames(tmp) <- c("feature","sample","value"); tmp <- cbind(view=m,tmp); return(tmp) })
    ImputedData <- do.call(rbind,tmp)
    ImputedData[,c("view","feature","sample")] <- sapply(ImputedData[,c("view","feature","sample")], as.character)
  } else if ((length(views)==1) && (as.data.frame==F)) {
    ImputedData <- ImputedData[[views]]
  }
  
  return(ImputedData)
}

#' @rdname getCovariates
#' @name getCovariates
#' @title Extract covariates from the model
#' @description (NOT PROPERLY IMPLEMENTED YET) this function extracts the covariates defined by the user. 
#' Note that (for now) the model can be run with known covariates only if they are stored in the colData of the MultiAssayExperiment input object. See \code{\link{createMOFAobject}}.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param names names of the covariates
#' @export
#' 
getCovariates <- function(object, names) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")  
  if(class(object@InputData) != "MultiAssayExperiment") stop("To work with covariates, InputData has to be specified in form of a MultiAssayExperiment")  
  stopifnot(all(names %in% colnames(colData(object@InputData))))
  
  # Get covariates
  covariates <- colData(object@InputData)[,names]
  
  return(covariates)
}

#' @rdname getExpectations
#' @name getExpectations
#' @title Fetch expectations from the model
#' @description Function to extract the expectations from the (variational) posterior distributions of a trained MOFAmodel.
#' @param object a \code{\link{MOFAmodel}} object.
#' @param variable variable name, 'Z' for factors, 'SW' for weights, 'Tau' for noise, 'Y' for pseudodata, 'Theta' for feature-wise spike-and-slab sparsity, 'AlphaW' for view and factor-wise ARD sparsity
#' @param as.data.frame boolean indicating whether to output the result as a long data frame, default is FALSE.
#' @details Technical note: MOFA is a Bayesian model where each variable has a prior distribution and a posterior distribution. 
#' In particular, to gain speed we used the variational inference framework so true posterior distributions are replaced by approximated variational distributions.
#' The priors and variational distributions of each variable are extensively described in the supplementary methods of the paper.
#' @return the output varies depending on the variable of interest: \cr
#' Z: a matrix with dimensions (samples,factors). If as.data.frame is TRUE, a long-formatted data frame with columns (sample,factor,value). \cr
#' SW: a list of length (views) where each element is a matrix with dimensions (features,factors). If as.data.frame is TRUE, a long-formatted data frame with columns (view,feature,factor,value). \cr
#' Y: a list of length (views) where each element is a matrix with dimensions (features,samples). If as.data.frame is TRUE, a long-formatted data frame with columns (view,feature,sample,value). \cr
#' Theta: TO-FILL \cr
#' Tau: TO-FILL
#' @export
getExpectations <- function(object, variable, as.data.frame = FALSE) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(variable %in% names(object@Expectations))
  
  # Get expectations in single matrix or list of matrices (for multi-view nodes)
  if (variable=="Z") {
    exp <- object@Expectations$Z$E
  } else {
    exp <- lapply(object@Expectations[[variable]], function(x) x$E)
  }
  
  # Convert to long data frame
  if (as.data.frame==T) {
    if (variable=="Z") {
      tmp <- reshape2::melt(exp)
      colnames(tmp) <- c("sample","factor","value")
    }
    if (variable=="SW") {
      tmp <- lapply(names(exp), function(m) { tmp <- reshape2::melt(exp[[m]]); colnames(tmp) <- c("feature","factor","value"); tmp$view <- m;  tmp[c("view","feature","factor")] <- sapply(tmp[c("view","feature","factor")], as.character); return(tmp) })
      tmp <- do.call(rbind.data.frame,tmp)
    }
    if (variable=="Y") {
      tmp <- lapply(names(exp), function(m) { tmp <- reshape2::melt(exp[[m]]); colnames(tmp) <- c("sample","feature","value"); tmp$view <- m;  tmp[c("view","feature","factor")] <- sapply(tmp[c("view","feature","factor")], as.character); return(tmp) })
      tmp <- do.call(rbind,tmp)
    }
    if (variable=="Tau") {
      stop("Not implemented")
      tmp <- lapply(names(exp), function(m) { data.frame(view=m, feature=names(exp[[m]]), value=unname(exp[[m]])); tmp[c("view","feature","factor")] <- sapply(tmp[c("view","feature","factor")], as.character); return(tmp) })
      tmp <- do.call(rbind,tmp)
    }
    if (variable=="AlphaW") {
      tmp <- lapply(names(exp), function(m) { tmp <- data.frame(view=m, factor=names(exp[[m]]), value=unname(exp[[m]])); tmp[c("view","feature","factor")] <- sapply(tmp[c("view","feature","factor")], as.character); return(tmp) })
      tmp <- do.call(rbind,tmp)
    }
    if (variable=="Theta") {
      stop("Not implemented")
      tmp <- lapply(names(exp), function(m) { tmp <- reshape2::melt(exp[[m]]); colnames(tmp) <- c("sample","feature","value"); tmp$view <- m; tmp[c("view","feature","factor")] <- sapply(tmp[c("view","feature","factor")], as.character); return(tmp) })
      tmp <- do.call(rbind,tmp)
    }
    exp <- tmp
  }
  return(exp)
}


