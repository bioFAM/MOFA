
##################
## Factor Names ##
##################

#' @title factorNames: set and retrieve factor names
#' @name factorNames
#' @rdname factorNames
#' @export
setGeneric("factorNames", function(object) {standardGeneric("factorNames")})

#' @name factorNames<-
#' @rdname factorNames
#' @export
setGeneric("factorNames<-", function(object, value) {standardGeneric("factorNames<-")})


##################
## Sample Names ##
##################

#' @title sampleNames: set and retrieve sample names
#' @name sampleNames
#' @rdname sampleNames
#' @export
setGeneric("sampleNames", function(object) {standardGeneric("sampleNames")})

#' @name sampleNames<-
#' @rdname sampleNames
#' @export
setGeneric("sampleNames<-", function(object, value) {standardGeneric("sampleNames<-")})

###################
## Feature Names ##
###################

#' @title featureNames: set and retrieve feature names
#' @name featureNames
#' @rdname featureNames
#' @export
setGeneric("featureNames", function(object) {standardGeneric("featureNames")})

#' @name featureNames<-
#' @rdname featureNames
#' @export
setGeneric("featureNames<-", function(object, value) {standardGeneric("featureNames<-")})

################
## View Names ##
################

#' @title viewNames: set and retrieve view names
#' @name viewNames
#' @rdname viewNames
#' @export
setGeneric("viewNames", function(object) {standardGeneric("viewNames")})

#' @name viewNames<-
#' @rdname viewNames
#' @export
setGeneric("viewNames<-", function(object, value) {standardGeneric("viewNames<-")})


################
## Input Data ##
################

#' @title InputData: set and retrieve training data
#' @name InputData
#' @rdname InputData
#' @export
setGeneric("InputData", function(object) {standardGeneric("InputData")})

setGeneric(".InputData<-", function(object, value) {standardGeneric(".InputData<-")})


################
## Train Data ##
################

#' @title TrainData: set and retrieve training data
#' @name TrainData
#' @rdname TrainData
#' @export
setGeneric("TrainData", function(object) {standardGeneric("TrainData")})

setGeneric(".TrainData<-", function(object, value) {standardGeneric(".TrainData<-")})

###################
## Train Options ##
###################

#' @title TrainOpts: set and retrieve training opts
#' @name TrainOpts
#' @rdname TrainOpts
#' @export
setGeneric("TrainOpts", function(object) {standardGeneric("TrainOpts")})

setGeneric(".TrainOpts<-", function(object, value) {standardGeneric(".TrainOpts<-")})


###################
## Model Options ##
###################

#' @title ModelOpts: set and retrieve Model options
#' @name ModelOpts
#' @rdname ModelOpts
#' @export
setGeneric("ModelOpts", function(object) {standardGeneric("ModelOpts")})

setGeneric(".ModelOpts<-", function(object, value) {standardGeneric(".ModelOpts<-")})


######################
## Train Statistics ##
######################

#' @title TrainStats: set and retrieve training statistics
#' @name TrainStats
#' @rdname TrainStats
#' @export
setGeneric("TrainStats", function(object) {standardGeneric("TrainStats")})

setGeneric(".TrainStats<-", function(object, value) {standardGeneric(".TrainStats<-")})

##################
## Expectations ##
##################

#' @title Expectations: set and retrieve expectations
#' @name Expectations
#' @rdname Expectations
#' @export
setGeneric("Expectations", function(object) {standardGeneric("Expectations")})

setGeneric(".Expectations<-", function(object, value) {standardGeneric(".Expectations<-")})

################
## Parameters ##
################

#' @title Parameters: set and retrieve parameters
#' @name Parameters
#' @rdname Parameters
#' @export
setGeneric("Parameters", function(object) {standardGeneric("Parameters")})

setGeneric(".Parameters<-", function(object, value) {standardGeneric(".Parameters<-")})
