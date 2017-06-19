#' CLL multi-omic data
#'
#' A list containing differen omic measurements in CLL patient samples
#'
#' @format A list containing filtered and pre-processed omic measurements for 211 CLL patients from Dietrich, Oles, Sellner et al including partially missing patients.
#' \itemize{
#'   \item{lincRNA}{normalized expression of lincRNA, p=2363}
#'   \item{mRNA}{normalized expression of mRNA, p=5000}
#'   \item{miRNA}{normalized expression of miRNA, p=274}
#'   \item{meth}{methylation M-values of top variable sites excluding sex chromosomes, p=4248}
#'   \item{viab}{drug response data for 62 drugs in 5 concentrations, p=310}
#'   \item{mut}{somatic mutations, IGHV status and copy number aberrations, p=69}
#'   \item{covariates}{diganosis and sex, p=2}
#' }
#' @name CLL_views
#' @usage data(CLL_views)

NULL

# created using the following:
# 
# inputDir <- "~/Documents/MOFA/CLL_MOFA_data/views/all_small_noXY_alldrugs2/"
# files <- list.files(inputDir)
# viewfiles <- files[grepl(".txt", files) ]
# CLL_views <- lapply(viewfiles, function(filenm) t(as.matrix(read.table(file.path(inputDir, filenm)))))
# names(CLL_views) <- sub(".txt", "",viewfiles)
# save(CLL_views, file = "~/Documents/MOFA/MOFApackage/scGFA/MOFAtools/data/CLL_views.RData")
