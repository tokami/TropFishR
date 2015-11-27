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
#' # load data
#' data(gillnet) #data_GillnetSelect3
#'
#' # run model
#' output <- MillarsGillnetSelect(gillnet, model = "normal_fixed",
#'    plotlens = NULL, rel = NULL)
#'
#' # investigate results
#' output
#'
#' #plot results
#' plot(output)
#'
#' @details A function to plot the results of the gillnet selectivity estimation.
#'
#' @references
#' Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. ICES Journal of Marine Science: Journal du Conseil, 54(3), 471-477.
#'
#' @export



x = fit # from Netfit


plot.select_Millar <- function(x,
                               Meshsize=NULL,
                               plotlens=NULL,
                               standardize=TRUE,
                               ...
                               ){

  r=selncurves(fit$rtype) #Get selection curve function
  if(is.null(plotlens)) plotlens=fit$Data[,1]
  if(is.null(Meshsize)) Meshsize=fit$Meshsize
  plot.title=switch(fit$rtype,
                    "norm.loc"="Normal (common spread)",
                    "norm.sca"="Normal",
                    "lognorm"="Lognormal",
                    "binorm.sca"="Bi-normal",
                    "bilognorm"="Bi-lognormal",
                    "tt.logistic"="Control and logistic","")
  rmatrix=outer(plotlens,Meshsize,r,fit$par)
  rmatrix=t(t(rmatrix)*fit$rel.power)
  if(standardize) rmatrix=rmatrix/max(rmatrix)
  matplot(plotlens,rmatrix,type="l",las=1,ylim=c(0,1),
          xlab="Length (cm)",ylab="Relative retention",...)
  #abline(h=seq(0,1,0.25),lty=3)
  lenrmatrix=cbind(plotlens,rmatrix)
  colnames(lenrmatrix)=c("Length",Meshsize)
  invisible(lenrmatrix)

}




