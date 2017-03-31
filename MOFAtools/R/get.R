
#' @rdname getExpectations
#' @name getExpectations
#' @title wraper to fetch particular expectations from the model
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @param variable: variable name ('SW')
#' @param expectation_name: expectation name ('E')
#' @export
getExpectations <- function(object, variable, expectation) {
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  if (variable=="Z") {
    object@Expectations$Z[[expectation]]
  } else {
    lapply(object@Expectations[[variable]], function(x) x[[expectation]])
  }
}

#' @rdname getParameters
#' @name getParameters
#' @title wraper to fetch particular parameters from the model
#' @description to-fill
#' @param object a \code{\link{MOFAmodel}} object.
#' @param variable: variable name ('SW')
#' @param parameter_name: parameter name ('E')
#' @export
getParameters <- function(object, variable, parameter) {
  if (class(object) != "MOFAmodel")
    stop("'object' has to be an instance of MOFAmodel")
  
  if (variable=="Z") {
    object@Parameters$Z[[parameter]]
  } else {
    lapply(object@Parameters[[variable]], function(x) x[[parameter]])
  }
}

