
##########################################################
## User-friendly functions to fetch data from the model ##
##########################################################

#' @rdname getDimensions
#' @name getDimensions
#' @title wraper to extract dimensionalities from the MOFAmodel. 
#' @description K indicates the number of factors, D indicates the number of features, N indicates the number of samples and M indicates the number of views.
#' @param object a \code{\link{MOFAmodel}} object.
#' @export
#' 
getDimensions <- function(object) {
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")  
  return(object@Dimensions)
}


#' @rdname getFactors
#' @name getFactors
#' @title wraper to extract the latent factors from the model
#' @param object a \code{\link{MOFAmodel}} object.
#' @param as.data.frame: boolean indicating whether to return factors as a long data frame with columns (sample,factor,value), default is FALSE.
#' @param include_intercept: boolean indicating where to include the intercept term of the model, default is TRUE.
#' @param factors: vector with the factors indices (numeric) or factor names (character) to fetch.
#' @return by default returns a matrix of dimensionality (N,K) where N is number of samples and k is number of factors. 
#' Alternatively, with as.data.frame=TRUE, it returns a data frame with columns (sample,factor,value)
#' @export
#' 
getFactors <- function(object, as.data.frame = FALSE, include_intercept = TRUE, factors = "all") {
  
  # Sanity check
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")  
  
  # Get factors
  if (paste0(factors,collapse="") == "all") { factors <- factorNames(object) } else { stopifnot(all(factors %in% factorNames(object))) }

  # Collect factors
  Z <- getExpectations(object,"Z","E",as.data.frame)
    if (as.data.frame==F) {
      Z <- Z[,factors, drop=FALSE]
    } else {
      Z <- Z[Z$factor %in% factors,]
    }

  # Remove intercept
  if (include_intercept == F & "intercept" %in% colnames(Z)) {
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
#' @param views vector containing the names of views (character) or index of views (numeric) to be fetched.
#' @param factors vector with the factors indices (numeric) or factor names (character) to fetch.
#' @param as.data.frame: boolean indicating whether to return factors as a long data frame with columns (view,feature,factor,value), default is FALSE.
#' @return by default returns a list of length M, where M is the number of views, of matrices with dimensionality (Dm,K) where Dm is number of features in view m and k is number of factors. 
#' Alternatively, with as.data.frame=TRUE, it returns a data frame with columns (view,feature,factor,value)
#' @export
#' 
getWeights <- function(object, views = "all", factors = "all", as.data.frame = F) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # Get views and factors
  if (paste0(views,collapse="") == "all") { views <- viewNames(object) } else { stopifnot(all(views %in% viewNames(object))) }
  if (paste0(factors,collapse="") == "all") { factors <- factorNames(object) } else { stopifnot(all(factors %in% factorNames(object))) }
  
  # Fetch weights
  weights <- getExpectations(object,"SW","E",as.data.frame)
  if (as.data.frame) {
    weights <- weights[weights$view%in%views & weights$factor%in%factors, ]
  } else {
    weights <- lapply(views, function(m) weights[[m]][,factors,drop=F])
    if (length(views)==1) { weights <- weights[[1]] }
    names(weights) <-  views
  }
  return(weights)
}


#' @rdname getTrainData
#' @name getTrainData
#' @title wraper to fetch training data from the model
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views vector containing the names of views (character) or index of views (numeric) to be fetched.
#' @param features: list with feature indices (numeric) or names (character), with the same order as views (default is "all")
#' @param as.data.frame: boolean indicating whether to return factors as a long data frame with columns (view,feature,sample,value), default is FALSE.
#' @return if as.data.frame == FALSE it return a list of length M, where M is the number of views, with matrices of dimensionality (N,Dm) where N is the number of samples and Dm is the numher of features in view m.
#' Alternatively, if as.data.frame == TRUE, it returns long data frame with columns (view,feature,sample,value).
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
  if (as.data.frame) {
    tmp <- lapply(views, function(m) { tmp <- reshape2::melt(trainData[[m]]); colnames(tmp) <- c("feature","sample","value"); tmp <- cbind(view=m,tmp); return(tmp) })
    trainData <- do.call(rbind,tmp)
    trainData[,c("view","feature","sample")] <- sapply(trainData[,c("view","feature","sample")], as.character)
  } else if ((length(views)==1) && (as.data.frame==F)) {
    trainData <- trainData[[views]]
  }
  
  return(trainData)
}

#' @rdname getImputedData
#' @name getImputedData
#' @title wraper to fetch imputed data from the model
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views vector containing the names of views (character) or index of views (numeric) to be fetched.
#' @param features: list with feature indices (numeric) or names (character), with the same order as views (default is "all")
#' @param as.data.frame: boolean indicating whether to return factors as a long data frame with columns (view,feature,sample,value), default is FALSE.
#' @return if as.data.frame == FALSE it return a list of length M, where M is the number of views, with matrices of dimensionality (N,Dm) where N is the number of samples and Dm is the numher of features in view m.
#' Alternatively, if as.data.frame == TRUE, it returns long data frame with columns (view,feature,sample,value).
#' @export
getImputedData <- function(object, views = "all", features = "all", as.data.frame = F) {
  
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
  if (as.data.frame) {
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
#' @title wraper to extract covariates from the MultiAssayExperiment stored in the InputData slot
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @param names: names of the covariates
#' @export
#' 
getCovariates <- function(object, names) {
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")  
  if(class(object@InputData) != "MultiAssayExperiment") stop("To work with covariates, InputData has to be specified in form of a MultiAssayExperiment")  
  stopifnot(all(names %in% colnames(colData(object@InputData))))
  return(colData(object@InputData)[,names])
}

#' @rdname getExpectations
#' @name getExpectations
#' @title wraper to fetch particular expectations from the model
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @param variable: variable name ('Z,'SW','Tau','Y','Theta' or 'Alpha')
#' @param expectation: expectation name ('E' for the first moment and 'E2' for the second moment (not implemented))
#' @param as.data.frame: output the result as a long data frame?
#' @return to-fill
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

