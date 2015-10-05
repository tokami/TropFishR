#' @name ex.LengthCC
#' @title Exemplary length composition data
#'
#'
#' @description Exemplary length composition data for Upeneus vittatus from Manila Bay, Philippines (from Ziegeler,1979) for the application of the length converted catch curve. Corresponding Linf and K are, 23.1cm and 0.59 per year, respectively. The \code{ex.LengthCC} data set is bull
#'
#' \itemize{
#'   \item midLengths Midpoints of length classes
#'   \item catch Catch in Numbers
#' }
#'
#' @docType data
#' @format A dataframe consisting of: 1.midLengths, and 2.catch (18 rows, 2 columns)
#' @source \url{Ziegler or Sparre?}
#' @usage data(ex.LengthCC)
#' @keywords datasets
#' @examples
#'
#' ### Ex 1.
#' data(ex.LengthCC)
#' output <- LengthConCatchCurve(midLengths = ex.LengthCC[,1], catch = ex.LengthCC[,2], Linf = 23.1, K = 0.59)
#'  output
#'
#'
NULL
