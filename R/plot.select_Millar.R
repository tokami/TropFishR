#' @title Millar's selectivity plot
#
#' @description  This function plots the selectivity estimates of Millar's
#'    selectivity model (\code{\link{select_Millar}}).
#'
#' @param x a list of the class \code{"select_Millar"} containing the results of
#'    Millar's selectivity model
#' @param plotlens A vector with lengths which should be used for drawing the
#'    selection curves
#' @param standardise A parameter indicating if the retention should be realtive
#'    to the maximum value (Default: TRUE).
#' @param deviance_plot logical; indicating whether a plot of deviance residuals should
#'    be displayed
#' @param ... additional parameter options from plot function
#'
#' @examples
#' \dontrun{
#' data(gillnet)
#'
#' output <- select_Millar(gillnet, x0 = c(60,4), rel.power = rep(1,8),
#'    rtype = "norm.loc", plot = FALSE)
#'
#' plot(output, plotlens = seq(40,90,0.1), deviance_plot = FALSE)
#'}
#' @details This function draws a selectivity plot for the object class
#'    \code{"select_Millar"}, which is created by applying Millar's selectivity model
#'    \code{\link{select_Millar}}.
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics abline axis matplot par plot points
#'
#' @references
#'  Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#'  54(3), 471-477.
#'
#' @export


plot.select_Millar <- function(x,
                               plotlens = NULL,
                               standardise = TRUE,
                               deviance_plot = TRUE,
                               ...){
# Adapted R code from Russell Millar (https://www.stat.auckland.ac.nz/~millar/selectware/)

  res <- x
  r <- rtypes_Millar(res$rtype)
  if(is.null(plotlens)) plotlens <- res$midLengths
  classes <- res$midLengths
  nlens <- length(classes)
  Dev.resids <- res$Dev.resids
  meshSizes <- res$meshSizes
  nmeshes <- length(meshSizes)
  AreLensUnique <- (length(classes)==length(unique(classes)))

  plot.title <- switch(res$rtype,
                    "norm.loc" = "Normal (common spread)",
                    "norm.sca" = "Normal",
                    "lognorm" = "Lognormal",
                    "binorm.sca" = "Bi-normal",
                    "bilognorm" = "Bi-lognormal",
                    "tt.logistic" = "Control and logistic","")

  rmatrix <- outer(plotlens, res$meshSizes, r, res$par)
  rmatrix = t(t(rmatrix) * res$rel.power)

  if(standardise) rmatrix = rmatrix / max(rmatrix, na.rm = TRUE)


  #dev.new()
  #create plot
  if(deviance_plot & ((nmeshes > 2 & AreLensUnique) | nmeshes == 2)){
    op <- par(mfrow = c(2,1), xpd = FALSE,
              mar = c(4, 4, 3, 1) + 0.1,
              oma = c(1, 1, 0.5, 2) + 0.1)

    if(nmeshes > 2 & AreLensUnique) {
      plot(1,1,xlim=range(classes),xlab="",ylab="Mesh size [cm]",
           ylim=range(meshSizes)+(1/50)*c(-1,1)*(max(meshSizes)-min(meshSizes)), # (cex/50)
           yaxt="n",type="n",main="Deviance residuals")
      axis(2,meshSizes,meshSizes,las=1)
      for(i in 1:nlens)
        for(j in 1:nmeshes)
          points(classes[i],meshSizes[j],pch=ifelse(Dev.resids[i,j]>0,16,1),
                 cex=3*abs(Dev.resids[i,j])*1/(abs(max(Dev.resids))))   # cex / (abs...)
    }else
      if(nmeshes == 2) {
        Dev.resids.len=sign(Dev.resids[,2])*sqrt(apply(Dev.resids^2,1,sum))
        plot(classes,Dev.resids.len,type=ifelse(AreLensUnique,"h","p"),las=1,
             main="Deviance residuals",xlab="",ylab="Mesh size [cm]",cex=1)   # cex = cex
        abline(h=0)
      }

  }else op <- par(mfrow = c(1,1),
                  mar = c(4, 4, 3, 1),
                  oma = c(3, 1, 1, 2))

  matplot(plotlens, rmatrix, type = "l", las = 1, ylim = c(0,1),
          xlab = "Length [cm]", ylab = "Relative retention",
          main = plot.title, ...)
  abline(h = seq(0,1,0.25), lty = 3)

  par(op)
}
