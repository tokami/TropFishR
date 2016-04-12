#' @name haddock
#'
#' @title Haddock data
#'
#' @description Data of a covered codend experimental catch of the haddock
#'    (\emph{Melanogrammus aeglefinus}).
#'    Can be used for function \code{\link{select_Millar}}.
#'
#' @docType data
#'
#' @format A list consiting of:
#' \itemize{
#'   \item \code{midLengths}  the midlengths of size classes,
#'   \item \code{numCodend}  the number of fish retained in codend,
#'   \item \code{numCover}  the number of fish retained in cover,
#' }
#'
#' @source Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#'  54(3), 471-477.
#'
#' @usage data(haddock)
#'
#' @keywords data dataset selectivity trawl
#'
#' @examples
#' data(haddock)
#' str(haddock)
#' summary(haddock)
#'
NULL
