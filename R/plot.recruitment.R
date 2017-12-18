#' @title Plot of recruitment patterns
#'
#' @description  This function plots the recruitment patterns from the
#'      \link{recruitment} model.
#'
#' @param x list of the class \code{"recruitment"}
#' @param percent logical; should number of recruits be relative (percentage)?
#' @param col colour of bars (default is "darkgreen")
#' @param xtitle title of x axis (default "rel. months" or no title, respectively)
#' @param ytitle title of y axis (default "# Recruits" or "Recruits [\%]", respectively)
#' @param ... standard parameters of \code{\link{barplot}}
#'
#' @examples
#' dat <- list(midLengths = seq(2,98,4),
#'                catch = c(0.6,17.6,93,83.2,12.6,0.3,0,0,0,1,17.1,51.4,
#'                26.1,2.2,0.2,4.5,21.6,17.6,3.7,8.7,10.6,6.3,5.6,2.9,0.8),
#'                Linf = 100,
#'                K = 0.5)
#' output <- recruitment(param = dat, tsample = 0.25)
#' plot(output, percent = FALSE)
#'
#' @importFrom grDevices dev.new
#' @importFrom graphics axis barplot par
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @method plot recruitment
#' @export

plot.recruitment <- function(x, percent = TRUE, col = "darkgreen",
                             xtitle = "default", ytitle = "default", ...){
  pes <- x

  if("t0" %in% names(pes)){
    xtitle <- ""
  }else xtitle <- "rel. months"
  if(percent){
    xvari <- pes$per_recruits
    ytitle <- "Recruits [%]"
  }else{
    xvari <- pes$mean_recruits
    ytitle <- "# Recruits"
  }
  if("t0" %in% names(pes)){
    months_num <- pes$months
    months_abb <- month.abb[months_num]
  }
  if(xtitle != "default") xtitle = xtitle
  if(ytitle != "default") ytitle = ytitle

  #save x axis positions if t0 != 0
  ylimits <- c(0,ceiling(max(xvari,na.rm=TRUE)*1.05))
  if("t0" %in% names(pes)){
    par(new = FALSE)
    mids <- barplot(xvari, xlab="", ann=TRUE, plot = FALSE)
  }

  #create plot
  #dev.new()
  #op <- par(mar = c(7, 5, 4, 5) + 0.3)
  barplot(xvari, col = col, ylab = ytitle,
          xlab = xtitle, ylim = ylimits)
  if("t0" %in% names(pes)){
    axis(1, at=mids, labels=months_abb)
  }else writeLines("If no t0 is provided only relative recruitment pattern can be estimated.")
  #par(op)
}
