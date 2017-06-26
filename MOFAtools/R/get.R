
########################################################
## User-friendly functions to get data from the model ##
########################################################

#' @rdname getTrainData
#' @name getTrainData
#' @title wraper to fetch training data from the model
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @param views: views (default is "all")
#' @param as.data.frame: output the result as a long data frame?
#' @export
getTrainData <- function(object, views="all", as.data.frame=F) {
  
  # Sanity checks
  if (class(object) != "MOFAmodel") stop("'object' has to be an instance of MOFAmodel")
  
  # Get views
  if (paste0(views,collapse="") == "all") { 
    views <- viewNames(object) 
  } else {
    stopifnot(all(views %in% viewNames(object)))
  }
  
  # Get data
  trainData <- object@TrainData[views]
  
  # Convert to long data frame
  if (as.data.frame) {
    tmp <- lapply(views, function(m) { tmp <- reshape2::melt(trainData[[m]]); colnames(tmp) <- c("sample","feature","value"); tmp <- cbind(view=m,tmp); return(tmp) })
    trainData <- do.call(rbind,tmp)
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
      tmp <- lapply(names(exp), function(m) { tmp <- reshape2::melt(exp[[m]]); colnames(tmp) <- c("feature","factor","value"); tmp <- cbind(view=m,tmp); return(tmp) })
      tmp <- do.call(rbind,tmp)
    }
    if (variable=="Y") {
      tmp <- lapply(names(exp), function(m) { tmp <- reshape2::melt(exp[[m]]); colnames(tmp) <- c("sample","feature","value"); tmp <- cbind(view=m,tmp); return(tmp) })
      tmp <- do.call(rbind,tmp)
    }
    if (variable=="Tau") {
      tmp <- lapply(names(exp), function(m) data.frame(view=m, feature=names(exp[[m]]), value=unname(exp[[m]])) )
      tmp <- do.call(rbind,tmp)
    }
    if (variable=="AlphaW") {
      tmp <- lapply(names(exp), function(m) tmp <- data.frame(view=m, factor=names(exp[[m]]), value=unname(exp[[m]])) )
      tmp <- do.call(rbind,tmp)
    }
    if (variable=="Theta") {
      # PROBLEM: THETA CAN BE (D,K) OR (K) 
      stop("Not implemented")
    }
    # exp <- cbind(variable = variable, tmp)
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

