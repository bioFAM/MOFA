% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/AllGenerics.R, R/methods.R
\docType{methods}
\name{ModelOptions}
\alias{ModelOptions}
\alias{ModelOptions<-}
\alias{ModelOptions,MOFAmodel-method}
\alias{ModelOptions<-,MOFAmodel,list-method}
\title{ModelOptions: set and retrieve model options}
\usage{
ModelOptions(object)

ModelOptions(object) <- value

\S4method{ModelOptions}{MOFAmodel}(object)

\S4method{ModelOptions}{MOFAmodel,list}(object) <- value
}
\arguments{
\item{object}{a \code{\link{MOFAmodel}} object.}

\item{value}{a list with model options}
}
\value{
list of model options
}
\description{
Function to set and retrieve model options.
}
\examples{
# load a trained MOFAmodel object
filepath <- system.file("extdata", "scMT_model.hdf5", package = "MOFAdata")
MOFAobject <- loadModel(filepath)
ModelOptions(MOFAobject)
}
