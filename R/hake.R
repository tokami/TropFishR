#' @name hake
#' @title Hake data
#'
#' @description This dataset contains length-frequency data and biological
#'  characteristics about hake (\emph{Merluccius merluccius})
#'  and its fisheries off Senegal. It can be used for \code{\link{VPA}} or
#'  \code{\link{predict_mod}}.
#'
#' @docType data
#' @format A list consisting of: 1. a vector with midlengths of size classes,
#'  2. a vector with catch in numbers, 3. K value, 4. Linf value, 5. M value,
#'  6. a value, 7. b value, 8. a vector with fishing mortalities, and
#'  9. a vector with average value of fish per kg
#'
#' @source Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#' @usage data(hake)
#' @keywords data dataset hake VPA
#' @examples
#'
#' data(hake)
#' str(hake)
#' summary(hake)
#'
#'
NULL
