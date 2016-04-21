#' @name trammelnet
#'
#' @title Trammel net data
#'
#' @description Data of an experiment with several trammel nets with different mesh sizes.
#'    Can be used for function \code{\link{select_Millar}}.
#'
#' @docType data
#'
#' @format A list consiting of:
#' \itemize{
#'   \item \code{$midLengths}  the midlengths of size classes,
#'   \item \code{$meshSizes}  the meshsizes,
#'   \item \code{$catchPerNet_mat}  a matrix with the numbers in catch of the
#'      corresponding mesh sizes (same order),
#' }
#'
#' @source Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#'  54(3), 471-477.
#'
#' @usage data(trammelnet)
#'
#' @keywords data dataset selectivity trammelnet
#'
#' @examples
#' data(trammelnet)
#' str(trammelnet)
#' summary(trammelnet)
#'
NULL
