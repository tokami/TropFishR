#' @title Millar's selectivity plot
#
#' @description  This function plots the selectivity estimates of Millar's selectivity model
#'    (\code{\link{select_Millar}}).
#'
#' @param x A list of the class \code{"select_Millar"} containing the results of Millar's selectivity model
#' @param plotlens A vector with lengths which should be used for drawing the selection curves
#' @param standardise A parameter indicating if the retention should be realtive to the maximum value (Default: \code{standardise = TRUE}).
#' @param ... Default parameter options from plot function
#'
#' @examples
#' # Gillnet
#' # load data
#' data(gillnet)
#'
#' # run model
#' output <- select_Millar(gillnet, x0 = c(60,4), rel.power = rep(1,8),
#'    rtype = "norm.loc")
#'
#' # plot results
#' plot(output,plotlens=seq(40,90,0.1))
#'
#' @details This function draws a selectivity plot for the object class
#'    \code{"select_Millar"}, which is created by applying Millar's selectivity model
#'    \code{\link{select_Millar}}.
#'
#' @references
#'  Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. ICES Journal of Marine Science: Journal du Conseil, 54(3), 471-477.
#'
#' @export

plot.select_Millar <- function(x,
                               plotlens=NULL,
                               standardise=TRUE,
                               ...
                               ){

  res <- x
  r <- rtypes_Millar(res$rtype) #Get selection curve function
  if(is.null(plotlens)) plotlens <- res$midLengths

  plot.title <- switch(res$rtype,
                    "norm.loc" = "Normal (common spread)",
                    "norm.sca" = "Normal",
                    "lognorm" = "Lognormal",
                    "binorm.sca" = "Bi-normal",
                    "bilognorm" = "Bi-lognormal",
                    "tt.logistic" = "Control and logistic","")

  rmatrix <- outer(plotlens,res$meshSizes,r,res$par)
  rmatrix=t(t(rmatrix)*res$rel.power)

  if(standardise) rmatrix=rmatrix/max(rmatrix, na.rm = TRUE)

  # Plot
  matplot(plotlens, rmatrix, type="l", las=1, ylim=c(0,1),
          xlab="Length (cm)", ylab="Relative retention",
          main = plot.title,...)
  abline(h=seq(0,1,0.25),lty=3)
}




