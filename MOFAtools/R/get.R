
##########################################################
## User-friendly functions to fetch data from the model ##
##########################################################

#' @rdname getDimensions
#' @name getDimensions
#' @title wraper to extract dimensionalities from the MOFAmodel
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @export
#' 
getDimensions <- function(object) {
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")  
  return(object@Dimensions)
  
}

#' @rdname getCovariates
#' @name getCovariates
#' @title wraper to extract covariates from the MultiAssayExperiment stored in the InputData slot
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @param names: names of the covariates
#' @export
#' 
getCovariates <- function(object, names) {
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")  
  # TO-DO: CHECK THAT WE HAVE A MULTIASSAYEXPERIMENT SLOT
  stopifnot(all(names %in% colnames(colData(object@InputData))))
  return(colData(object@InputData)[,names])
}


#' @rdname getFactors
#' @name getFactors
#' @title wraper to extract the latent factors from the model
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @export
#' 
getFactors <- function(object, as.data.frame=F, include_intercept=T) {
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")  
  
  Z <- getExpectations(object,"Z","E",as.data.frame)
  
  if (!include_intercept) {
    if (as.data.frame==F) {
      Z <- Z[,colnames(Z)!="intercept"]
    } else {
      Z <- Z[Z$factor!="intercept",]
    }
  }
  
  return(Z)
}


#' @rdname getWeights
#' @name getWeights
#' @title wraper to extract the weights from the model
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @export
#' 
getWeights <- function(object, views="all", factors="all", as.data.frame=F) {
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  # Get views and factors
  if (paste0(views,collapse="") == "all") { views <- viewNames(object) } else { stopifnot(all(views %in% viewNames(object))) }
  if (paste0(factors,collapse="") == "all") { factors <- factorNames(object) } else { stopifnot(all(factors %in% factorNames(object))) }
  
  weights <- getExpectations(object,"SW","E",as.data.frame)
  if (as.data.frame) {
    weights <- weights[weights$view%in%views & weights$factor%in%factors, ]
  } else {
    weights <- lapply(views, function(m) weights[[m]][,factors,drop=F])
    if (length(views)==1) { weights <- weights[[1]] }
  }
  return(weights)
}


#' @rdname getTrainData
#' @name getTrainData
#' @title wraper to fetch training data from the model
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views: views (default is "all")
#' @param features: list with feature names, with the same order as views (default is "all")
#' @param as.data.frame: output the result as a long data frame?
#' @export
getTrainData <- function(object, views="all", features="all", as.data.frame=F) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # Get views
  if (paste0(views,collapse="") == "all") { 
    views <- viewNames(object) 
  } else {
    stopifnot(all(views %in% viewNames(object)))
  }
  
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
  
  # Get data
  trainData <- object@TrainData[views]
  trainData <- lapply(1:length(trainData), function(m) trainData[[m]][features[[m]],,drop=F]); names(trainData) <- views
  
  # Convert to long data frame
  if (as.data.frame) {
    tmp <- lapply(views, function(m) { tmp <- reshape2::melt(trainData[[m]]); colnames(tmp) <- c("feature","sample","value"); tmp <- cbind(view=m,tmp); return(tmp) })
    trainData <- do.call(rbind,tmp)
    trainData[,c("view","feature","sample")] <- sapply(trainData[,c("view","feature","sample")], as.character)
  } else if ((length(views)==1) && (as.data.frame==F)) {
    trainData <- trainData[[views]]
  }
  
  return(trainData)
}


#' @rdname getExpectations
#' @name getExpectations
#' @title wraper to fetch particular expectations from the model
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @param variable: variable name ('SW')
#' @param expectation: expectation name ('E')
#' @param as.data.frame: output the result as a long data frame?
#' @export
getExpectations <- function(object, variable, expectation, as.data.frame=F) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(variable %in% names(object@Expectations))
  
  # Get expectations in single matrix or list of matrices (for multi-view nodes)
  if (variable=="Z") {
    exp <- object@Expectations$Z[[expectation]]
  } else {
    exp <- lapply(object@Expectations[[variable]], function(x) x[[expectation]])
  }
  
  # Convert to long data frame
  if (as.data.frame) {
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
      tmp <- lapply(names(exp), function(m) { data.frame(view=m, feature=names(exp[[m]]), value=unname(exp[[m]])); tmp[c("view","feature","factor")] <- sapply(tmp[c("view","feature","factor")], as.character); return(tmp) })
      tmp <- do.call(rbind,tmp)
    }
    if (variable=="AlphaW") {
      tmp <- lapply(names(exp), function(m) { tmp <- data.frame(view=m, factor=names(exp[[m]]), value=unname(exp[[m]])); tmp[c("view","feature","factor")] <- sapply(tmp[c("view","feature","factor")], as.character); return(tmp) })
      tmp <- do.call(rbind,tmp)
    }
    if (variable=="Theta") {
      stop()
      tmp <- lapply(names(exp), function(m) { tmp <- reshape2::melt(exp[[m]]); colnames(tmp) <- c("sample","feature","value"); tmp$view <- m; tmp[c("view","feature","factor")] <- sapply(tmp[c("view","feature","factor")], as.character); return(tmp) })
      tmp <- do.call(rbind,tmp)
    }
    exp <- tmp
  }
  return(exp)
}


#' @rdname getParameters
#' @name getParameters
#' @title wraper to fetch particular parameters from the model
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @param variable: variable name ('AlphaW')
#' @param parameter: parameter name ('a')
#' @param as.data.frame: output the result as a long data frame?
#' @export
getParameters <- function(object, variable, parameter, as.data.frame=F) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  stopifnot(variable %in% names(object@Parameters))
  
  # Get expectations in single matrix or list of matrices (for multi-view nodes)
  if (variable=="Z") {
    par <- object@Parameters$Z[[parameter]]
  } else {
    par <- lapply(object@Parameters[[variable]], function(x) x[[parameter]])
  }
  
  # Convert to long data frame
  if (as.data.frame) {
    if (variable=="Z") {
      tmp <- reshape2::melt(par)
      colnames(tmp) <- c("sample","factor","value")
    }
    if (variable=="SW") {
      tmp <- lapply(names(par), function(m) { tmp <- reshape2::melt(par[[m]]); colnames(tmp) <- c("feature","factor","value"); tmp <- cbind(view=m,tmp); return(tmp) })
      tmp <- do.call(rbind,tmp)
    }
    if (variable=="Y") {
      tmp <- lapply(names(par), function(m) { tmp <- reshape2::melt(par[[m]]); colnames(tmp) <- c("sample","feature","value"); tmp <- cbind(view=m,tmp); return(tmp) })
      tmp <- do.call(rbind,tmp)
    }
    if (variable=="Tau") {
      tmp <- lapply(names(par), function(m) data.frame(view=m, feature=names(par[[m]]), value=unname(par[[m]])) )
      tmp <- do.call(rbind,tmp)
    }
    if (variable=="AlphaW") {
      tmp <- lapply(names(par), function(m) tmp <- data.frame(view=m, factor=names(par[[m]]), value=unname(par[[m]])) )
      tmp <- do.call(rbind,tmp)
    }
    if (variable=="Theta") {
      # PROBLEM: THETA CAN BE (D,K) OR (K) 
      stop("Not implemented")
    }
    # par <- cbind(variable = variable, tmp)
    par <- tmp
  }
  
}

