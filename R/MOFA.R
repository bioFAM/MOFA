#' Multi-Omics Factor Analysis (MOFA)
#' 
#' @description MOFA provides an unsupervised framework for the integration of multi-omics data sets.
#' Given several data matrices with measurements of multiple ‘omics data
#' types on the same or on overlapping sets of samples, MOFA infers an
#' interpretable low-dimensional data representation in terms of (hidden) factors.
#' These learnt factors represent the driving sources of variation across data modalities,
#' thus facilitating the identification of cellular states or disease subgroups.  
#' 
#' The package contains all function required for training MOFA on a multi-omics data set as well as
#' for different downstream analyes, such as visualisation of samples in factor space, annotation of factors 
#' to molecular markers or gene sets, outlier identification and imputation of missing values.
#' 
#' @details Please have a look at the vignette "MOFA" for a in-depth introduction to the package.
#' @docType package
#' @name MOFA
NULL