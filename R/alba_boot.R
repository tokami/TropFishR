#' @name alba_boot
#'
#' @title Bootstrapped VBGF estimates for the alba lfq dataset
#'
#' @description Bootstrapped VBGF estimates for the \link{alba} lfq dataset as
#'   estimated by \code{\link{ELEFAN_GA_boot}} (200 resamples).
#'
#' @docType data
#'
#' @format A list of class `lfqBoot`, see \code{\link{ELEFAN_GA_boot}}.
#'
#' @source The data set was produced with following code (takes several
#'   minutes):
#' \preformatted{
#' data(alba)
#' alba_boot <- ELEFAN_GA_boot(lfq = alba, MA = 7, seasonalised = FALSE,
#'   low_par = list(Linf = 8, K = 0.1, t_anchor = 0),
#'   up_par = list(Linf = 20, K = 5, t_anchor = 1),
#'   popSize = 100, maxiter = 50, run = 10, pmutation = 0.2,
#'   nresamp = 200, seed = 20, parallel = TRUE)
#' }
#'
#' @usage data(alba_boot)
#' @keywords data dataset length-frequency bootstrap
#'
#' @examples
#' data(alba_boot)
#'
#' head(alba_boot$bootRaw)
#' alba_boot$maxDen
#' alba_boot$CI
#'
#' # plot
#' univariate_density(alba_boot)
#' LinfK_scatterhist(alba_boot, phi.contour = TRUE)
#'
NULL
