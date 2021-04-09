#' A proteomics dataset with random protein expression values - with missing values.
#'
#' A longitudinal proteomics dataset with five timepoints, two conditions and
#' three replicates for each sample at each timepoint in each condition. The expression
#' values for the random data have been generated using the \code{rnorm} function. The
#' pattern of missing values has been directly copied from a semi-simulated spike-in
#' CPTAC dataset. The missing values in the used semi-simulated CPTAC dataset in turn have
#' been copied from the corresponding samples in the original CPTAC dataset. For more
#' details about the creation of the semi-simulated datasets, see the RolDE original
#' publication.
#'
#' @format A matrix with 1247 rows and 30 variables.
#' @source \url{https://cptac-data-portal.georgetown.edu/cptac/study/showDetails/10424}
"data2"
