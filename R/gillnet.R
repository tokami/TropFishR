#' @name gillnet
#'
#' @title Gillnet data
#'
#' @description Data of an experiment with several gillnets with different mesh sizes.
#'    Can be used for function \code{\link{MillarsGillnetSelect}}.
#'
#' @docType data
#'
#' @format A list consiting of:
#' \itemize{
#'   \item \code{$midLengths}  the midlengths of size classes,
#'   \item \code{$meshSizes}  the meshsizes,
#'   \item \code{$catchPerNet_mat}  a matrix with the numbers in catch of the corresponding mesh sizes (same order),
#' }
#'
#' @source Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. ICES Journal of Marine Science: Journal du Conseil, 54(3), 471-477.
#'
#'  Holt, S. J. 1963. A method for determining gear selectivity and its application.
#'  ICNAF Special Publication, 5: 106-115.
#'
#' @usage data(gillnet)
#' @keywords dataset selectivity gillnet
#' @examples
#' data(gillnet)
#' str(gillnet)
#'
#'
NULL
