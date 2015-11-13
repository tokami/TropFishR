#' @name data_GillnetSelect
#'
#' @title Example dataset for the use of \code{\link[FishTropR]{GillnetSelect}}
#'
#' @description A list with characteristics of experimental catches of \textit{Tilapia esculenta}
#'   with gillnets of two mesh sizes. The Results of the experiment are the midlengths of
#'   size classes, the number of fish caught with net 1 & 2, and the meshsizes of both nets.
#'
#' @docType data
#'
#' @format A list consiting of:
#' \itemize{
#'   \item \code{$midLengths}  the midlengths of size classes,
#'   \item \code{$numNet1}  the number of fish caught with net 1,
#'   \item \code{$numNet2}  the number of fish caught with net 2,
#'   \item \code{$msNet1}  the meshsize of net 1,
#'   \item \code{$msNet2}  and the meshsize of net 2.
#' }
#'
#' @source Garrod, D.J., 1961. The selection characteristics of nylon gill nets for \textit{Tilapia esculenta} Graham. J.Cons.CIEM, 26(2):191-203
#'
#' @usage data(data_GillnetSelect)
#'
#' @keywords dataset selectivity
#'
#' @examples
#' data(data_GillnetSelect)
#' summary(data_GillnetSelect)
#' str(data_GillnetSelect)
#'
NULL
