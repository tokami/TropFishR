#' @title Bhattacharya plot
#
#' @description  This function plots the seperated frequency distributions
#'      and selected regression lines of \link[TropFishR]{Bhattacharya} method.
#'
#' @param x a list of the class \code{"Bhattacharya"} containing the results of
#'      Bhattacharya's method.
#' @param ... additional options of the \code{\link{plot}} function
#' @param analysisPlot logical; indicating wheter the anaylsis graph with the regression
#'      lines should be created
#'
#' @examples
#' \donttest{
#'  data(synLFQ1)
#'  output <- Bhattacharya(param = synLFQ1)
#'  plot(output)
#' }
#'
#' @details This function plots the results of the Bhattacharya method.
#'
#' @importFrom graphics abline hist layout lines par plot points title
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. \emph{FAO Fisheries Technical Paper}, (306.1, Rev. 2). 407 p.
#'
#' @method plot Bhattacharya
#' @export

plot.Bhattacharya <- function(x, analysisPlot = TRUE, ...){
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

  if(analysisPlot == FALSE){
    par(mfrow = c(1,1))

    hist(frequis, breaks = 50, main = "", xlab = "L",
         probability = TRUE, ylim = c(0,maxlim))

    for(cohorti in 1:length(temp.list)){
      lines(temp.list[[cohorti]]$xfit, temp.list[[cohorti]]$yyfit,
            col=colour.xy[cohorti], lwd=1.5)
    }
  }

  if(analysisPlot == TRUE){
    par(mfrow = c(2,1),
        oma = c(6.5, 0.5, 2, 1) + 0.1,
        mar = c(0, 4, 0.5, 0) + 0.1)
    layout(matrix(c(1,2),nrow=2), heights = c(1,2.5))

    hist(frequis, breaks = 50, main = "", xaxt = 'n',
         probability = TRUE, ylim = c(0,maxlim))

    for(cohorti in 1:length(temp.list)){
      lines(temp.list[[cohorti]]$xfit, temp.list[[cohorti]]$yyfit,
            col=colour.xy[cohorti], lwd=1.5)
    }

    plot(bhat.results$L, bhat.results$delta.log.N1.plus, pch = 4,
         col = 'black', xlab = "",...,
         ylab = expression(paste(delta," log N+")))
    abline(h = 0)
    seqi <- seq(4,length(names(bhat.results)),by = 7) # depending on order of bhat.results table

    for(coli in 1:length(seqi)){
      sta <- a.b.df[coli,4]
      end <- a.b.df[coli,5]
      points(bhat.results$L[sta:end],
             bhat.results[sta:end,seqi[coli]],
             col = colour.xy[coli], pch = 16)
    }
    for(coli in 1:length(seqi)){
      ai <- a.b.df[coli,2]
      bi <- a.b.df[coli,3]
      abline(a = ai, b=bi, col=colour.xy[coli], lwd = 1.5)
    }
    title(xlab = "L", outer = TRUE, line = 2.5)
  }
}
