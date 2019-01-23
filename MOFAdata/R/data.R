#' CLL_data
#'
#' A list containing different omic measurements for CLL patient samples
#'
#' @format A list containing omic measurements for n=200 patient samples that are used as input data for MOFA
#' \itemize{
#'   \item{mRNA: }{normalized expression values of most variable genes, p=5000}
#'   \item{Methylation: }{methylation M-values in most varibale CpG sites, p=4248}
#'   \item{Drugs: }{viability values in response t different drugs adn concentrations, p=310}
#'   \item{Mutations: }{Mutation status for selected genes, p=69}
#' }
#' @name CLL_data
#' @usage data(CLL_data)
"CLL_data"


#' CLL_covariates
#'
#' Data frame containing additional information on the patient samples, i.e. diagnosis and gender.
#'
#' @format A  data frame diagnosis and gender for the n=200 patient samples in CLL_data
#' @name CLL_covariates
#' @usage data(CLL_covariates)
"CLL_covariates"

#' reactomeGS
#'
#' A matrix containing feature binary membership indictors for the Reactome Gene Sets.
#'
#' @format A matrix containing feature binary membership indictors for different genes (in columns)
#'  and Reactome Gene Sets (in rows).
#' @name reactomeGS
#' @usage data(reactomeGS)
 "reactomeGS"

#' scMT_data
#'
#' A MultiAssayExperiment containing data from a single cell multi-omics study on mESCs.
#'
#' @format A MultiAssayExperiment containing four Experiments:
#' \itemize{
#'   \item{RNA expression: }{ExpressionSet with normalized expression values of most variable genes, p=5000}
#'   \item{Met Enhancers: }{Methylation values at enhancers, p=5000}
#'   \item{Met CpG Islands: }{Methylation values at CpG Islands, p=5000}
#'   \item{Met Promoters: }{Methylation values at Promoters, p=5000}
#' } 
#' @name scMT_data
#' @usage data(scMT_data)
"scMT_data"