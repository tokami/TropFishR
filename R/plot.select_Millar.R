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
#' @param deviance_plot logical (Default: TRUE); indicating whether a plot of deviance residuals should
#'    be displayed
#' @param selectivity_plot logical (Default: TRUE); indicating whether a plot of relative retention selectivities
#' should be displayed
#' @param xlab_dev character string. Label for x axis of deviance plot. Default: "Length [cm]"
#' @param xlab_sel character string. Label for x axis of selectivity plot. Default: "Length [cm]"
#' @param ylab_dev character string. Label for y axis of deviance plot. Default: "Mesh size [cm]"
#' @param ylab_sel character string. Label for y axis of selectivity plot. Default: "Relative retention".
#' @param title_dev character string. Label for main title of deviance plot. Default: "Deviance residuals".
#' @param title_sel character string. Label for main title of selectivity plot. Default is
#' taken from the results of the select_Millar (e.g. res$rtype).
#' @param ... additional parameter options from plot function
#'
#' @examples
#' data(gillnet)
#'
#' output <- select_Millar(gillnet, x0 = c(60,4), rel.power = rep(1,8),
#'    rtype = "norm.loc", plot = FALSE)
#'
#' plot(output, plotlens = seq(40,90,0.1), deviance_plot = FALSE)
#'
#' @details This function draws a selectivity plot for the object class
#'    \code{"select_Millar"}, which is created by applying Millar's selectivity model
#'    \code{\link{select_Millar}}.
#'
#' @importFrom graphics abline axis matplot par plot points
#'
#' @references
#'  Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#'  54(3), 471-477.
#'
#' @method plot select_Millar
#' @export


plot.select_Millar <- function(x,
  plotlens = NULL,
  standardise = TRUE,
  deviance_plot = TRUE,
  selectivity_plot = TRUE,
  xlab_dev = "Length [cm]",
  xlab_sel = "Length [cm]",
  ylab_dev = "Mesh size [cm]",
  ylab_sel = "Relative retention",
  title_dev = "Deviance residuals",
  title_sel = NULL,
  ...
){
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

  if(is.null(title_sel)){
    title_sel <- switch(res$rtype,
      "norm.loc" = "Normal (common spread)",
      "norm.sca" = "Normal",
      "lognorm" = "Lognormal",
      "binorm.sca" = "Bi-normal",
      "bilognorm" = "Bi-lognormal",
      "tt.logistic" = "Control and logistic",
      ""
    )
  }


  rmatrix <- outer(plotlens, res$meshSizes, r, res$par)
  rmatrix = t(t(rmatrix) * res$rel.power)

  if(standardise) rmatrix = rmatrix / max(rmatrix, na.rm = TRUE)

  # define number of panels required for plots
  npanels <- sum(deviance_plot & ((nmeshes > 2 & AreLensUnique) | nmeshes == 2)) + sum(selectivity_plot)
  op <- par(mfrow = c(npanels,1))

  # deviance plot
  if(deviance_plot & ((nmeshes > 2 & AreLensUnique) | nmeshes == 2)){

    if(nmeshes > 2 & AreLensUnique) {
      plot(1, 1, xlim=range(classes), xlab=xlab_dev, ylab=ylab_dev,
           ylim=range(meshSizes)+(1/50)*c(-1,1)*(max(meshSizes)-min(meshSizes)), # (cex/50)
           yaxt="n", type="n", main=title_dev)
      axis(2,meshSizes,meshSizes,las=1)
      for(i in 1:nlens)
        for(j in 1:nmeshes)
          points(classes[i],meshSizes[j],pch=ifelse(Dev.resids[i,j]>0,16,1),
                 cex=3*abs(Dev.resids[i,j])*1/(abs(max(Dev.resids))))   # cex / (abs...)
    }else{
      if(nmeshes == 2) {
        Dev.resids.len=sign(Dev.resids[,2])*sqrt(apply(Dev.resids^2,1,sum))
        plot(classes, Dev.resids.len, type=ifelse(AreLensUnique, "h", "p"), las=1,
             main=title_dev, xlab=xlab_dev, ylab=ylab_dev, cex=1)   # cex = cex
        abline(h=0)
      }
    }

  }

  # selectivity plot
  if(selectivity_plot){
    matplot(plotlens, rmatrix, type = "l", las = 1, ylim = c(0,1),
      xlab=xlab_sel, ylab=ylab_sel,
      main = title_sel, ...)
    abline(h = seq(0,1,0.25), lty = 3)
  }

  par(op)
}
