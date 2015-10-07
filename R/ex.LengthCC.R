#' @name ex.LengthCC
#' @title Example dataset of length for use in catch curve
#'
#'
#' @description bla bla
#'
#' \itemize{
#'   \item bla1. Dateframe of lon/lat coordinates corresponding to columns of sst$field (5 deg resolution)
#'   \item bla2. Vector of monthly date values corresponding to rows of sst$field
#'   \item bla3. Matrix of sea level pressure values by month and lon/lat position.
#' }
#'
#' @docType data
#' @format A list consisting of: 1. a dataframe for lon/lat positions, 2. a vector of
#' monthly date values, and 3. a matrix of sst anomaly values by month and lon/lat position
#' (1906 rows, 264 columns)
#' @source \url{http://www.esrl.noaa.gov/psd/data/gridded/data.kaplan_sst.html}
#' @usage data(ex.LengthCC)
#' @keywords datasets length
#' @examples
#'
#' ### Ex 1. Plot of single month
#' data(ex.LengthCC)
#' head(ex.LengthCC)
#'
#'
NULL
