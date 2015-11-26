#' @title Millar's gillnet selectivity plot
#
#' @description  This function plots the selectivity estimates of the function \code{\link{GillnetSelect}}.
#'
#' @param estim A list of the class "GillnetSelect" containing the results of the gillnet selectivity function.
#' @param ret_dev draw both plots? Default = c(T,T)
#' @param cex  thickness of points Default = 1
#' @param title Title of plot Default = "Deviance residuals"
#'
#' @examples
#' # load data
#' data(data_GillnetSelect3)
#'
#' # run model
#' output <- MillarsGillnetSelect(data_GillnetSelect3, model = "normal_fixed",
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

plot.MillarsGillnet <- function(estim, ret_dev = c(T,T), cex = 1,
                                title = "Deviance residuals"){  ###,...

  res <- estim
  classes.num <- res$midLengths
  model <- res$model
  plotlens <- res$plotlens
  rselect <- res$rselect
  residuals <- res$devres
  meshSizes <- res$meshSizes


  if(ret_dev[1] == T){
    #Plot the relative selection curves
    plot.title <- switch(model,
                         "normal_fixed"="Normal (common spread)",
                         "normal"="Normal",
                         "gamma"="Gamma",
                         "lognormal"="Log-normal")
    plot.title <- paste(plot.title,"retention curve")

    matplot(plotlens,rselect,type="l",lty=1,las=1,ylim=c(0,1),
            xlab="Length (cm)",ylab="Relative retention",main=plot.title)
  }
  if(ret_dev[2] == T){
    #Plot the deviance residuals matrix
    #if(is.null(classes.num)) lens <- 1:nrow(residuals)   #necessary? output should always have lengths and  mesh sizes!!
    #if(missing(meshSizes)) meshSizes <- 1:ncol(residuals)

    plot(c(min(classes.num),max(classes.num)),range(meshSizes), xlab="Length (cm)", ylab="Mesh size",
         ylim=range(meshSizes)+(cex/25)*c(-1,1)*(max(meshSizes)-min(meshSizes)),
         type="n", main=title)

    for(i in 1:nrow(residuals)){
      for(j in 1:ncol(residuals)){
        points(classes.num[i],meshSizes[j],pch=ifelse(residuals[i,j] > 0,16,1),
               cex=3*abs(residuals[i,j])*cex/(abs(max(residuals))))
      }
    }
  }
}
