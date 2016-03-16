#' @title Bhattacharya plot
#
#' @description  This function plots the seperated frequency distributions
#'      and selected regression lines of \link{Bhattacharya} method.
#'
#' @param x a list of the class \code{"Bhattacharya"} containing the results of
#'      Bhattacharya's method.
#'
#' @examples
#' \donttest{
#'  data(synLFQ1)
#'  output <- Bhattacharya(param = synLFQ1)
#'  plot(output)
#' }
#'
#' @details A function to plot the results of the Bhattacharya method.
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @export

plot.Bhattacharya <- function(x){
  pes <- x

  s.l.df <- pes$Lmean_SD_list
  a.b.df <- pes$regression_lines
  bhat.results <- pes$bhat_results
  temp.list <- pes$distributions


  midLengths <- bhat.results$mean.length.classes
  N <- bhat.results$N1.plus
  frequis <- rep(midLengths, N)

  colour.xy <- c('blue','darkgreen','red','goldenrod2','purple','orange',
                 'lightgreen','skyblue','brown','darkblue','darkorange','darkred')

  h <- hist(frequis, breaks = 50, plot = FALSE)

  # histogram
  maxTemps <- rep(NA,length(temp.list))
  for(maxis in 1:length(temp.list)){
    maxTemps[maxis] <- max(temp.list[[maxis]]$yyfit, na.rm = TRUE)
  }

  maxyyfit <- max(maxTemps, na.rm = TRUE)
  maxlim <- ifelse(maxyyfit > max(h$density), maxyyfit, max(h$density))

  par(mfrow = c(1,1))
  hist(frequis, breaks = 50, main = "", xlab = "L",
       probability = TRUE, ylim = c(0,maxlim))

  for(cohorti in 1:length(temp.list)){
    lines(temp.list[[cohorti]]$xfit, temp.list[[cohorti]]$yyfit,
          col=colour.xy[cohorti], lwd=1.5)
  }
}
