#' @name synCPUE
#'
#' @title Synthetical Catch per unit of effort (CPUE) dataset
#'
#' @description Synthetical CPUE dataset from Exercise 4.3 in Sparre & Venema (1999).
#'    Can be used to caluclate Z by means of the function \code{\link{Z_CPUE}}.
#'
#' @docType data
#'
#' @format A dataframe consisting of following columns:
#' \itemize{
#'   \item cohort Name of cohort, e.g. 1982 S, meaning summer cohort of 1982,
#'   \item age Age of cohorts,
#'   \item CPUE Catch per unit of effort of cohorts.
#' }
#'
#' @source Sparre, P., Venema, S.C., 1999. Introduction to tropical fish stock
#'    assessment. Part 2. Excercises. FAO Fisheries Technical Paper, (306.2, Rev. 2).
#'    94 p.
#'
#' @usage data(synCPUE)
#' @keywords data dataset CPUE
#' @examples
#' data(synCPUE)
#' head(synCPUE)
#'
#'
NULL
