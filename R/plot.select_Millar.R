#' @title Millar's selectivity plot
#
#' @description  This function plots the selectivity estimates of the function
#'    \code{\link{GillnetSelect}}. First sentence. second sentence.
#'
#' @param x A list of the class "GillnetSelect" containing the results of the gillnet selectivity function.
#' @param ret_dev draw both plots? Default = c(T,T)
#' @param cex  thickness of points Default = 1
#' @param title Title of plot Default = "Deviance residuals"
#' @param ... normal parameter options from plot function
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
#' # investigate results
#' output
#'
#' # plot results
#' plot(output,plotlens=seq(40,90,0.1))
#'
#' @details A function to plot the results of the gillnet selectivity estimation.
#'
#' @references
#' Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. ICES Journal of Marine Science: Journal du Conseil, 54(3), 471-477.
#'
#' @export

plot.select_Millar <- function(x,
                               plotlens=NULL,
                               standardize=TRUE,
                               ...
                               ){

  res <- x
  r <- rtypes_Millar(res$rtype) #Get selection curve function
  if(is.null(plotlens)) plotlens <- res$midLengths
  #if(is.null(meshSizes)) meshSizes <- res$meshSizes

  plot.title <- switch(res$rtype,
                    "norm.loc"="Normal (common spread)",
                    "norm.sca"="Normal",
                    "lognorm"="Lognormal",
                    "binorm.sca"="Bi-normal",
                    "bilognorm"="Bi-lognormal",
                    "tt.logistic"="Control and logistic","")

  rmatrix <- outer(plotlens,res$meshSizes,r,res$par)
  rmatrix=t(t(rmatrix)*res$rel.power)
  if(standardize) rmatrix=rmatrix/max(rmatrix)
  matplot(plotlens,rmatrix,type="l",las=1,ylim=c(0,1),
          xlab="Length (cm)",ylab="Relative retention",...)
  #abline(h=seq(0,1,0.25),lty=3)
  lenrmatrix=cbind(plotlens,rmatrix)
  colnames(lenrmatrix)=c("Length",res$meshSizes)
  invisible(lenrmatrix)
}




