#' @title Plotting prediction models yield per recruit
#'
#' @description This function plots objects of the class "ypr_mod".
#'
#' @param x a object of the class 'ypr_mod',
#' @param ... optional parameters of plot function
#'
#' @examples
#'
#'
#' @references
#' Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock assessment.
#' Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev. 2). 407 p.
#'
#' @export

plot.ypr_mod <- function(x,...){
  pes <- x
  FM <- pes$FM
  if("tc" %in% names(pes)) tc_Lc <- pes$tc
  if("Lc" %in% names(pes)) tc_Lc <- pes$Lc
  if("list_tc_runs" %in% names(pes)) list_tc_Lc_runs <- pes$list_tc_runs
  if("list_Lc_runs" %in% names(pes)) list_tc_Lc_runs <- pes$list_Lc_runs

  #Plot
  op <- par(mfrow=c(1,1),new=F, mar = c(5, 4, 4, 4) + 0.3)

  #plot Y/R & B/R
  label <- ifelse("tc" %in% names(pes), "tc", "Lc")
  tc_Lc_start <- which(tc_Lc == max(tc_Lc,na.rm=T))
  offset_text <- list_tc_Lc_runs[[tc_Lc_start]]$Y_R[length(list_tc_Lc_runs[[tc_Lc_start]]$Y_R)] *
    0.02
  plot(list_tc_Lc_runs[[tc_Lc_start]]$FM, list_tc_Lc_runs[[tc_Lc_start]]$Y_R, type = 'l',
       ylab = "Y/R", xlab = "F",lty=tc_Lc_start,...)
  text(x = list_tc_Lc_runs[[tc_Lc_start]]$FM[length(list_tc_Lc_runs[[tc_Lc_start]]$FM)],
       y = (list_tc_Lc_runs[[tc_Lc_start]]$Y_R[length(list_tc_Lc_runs[[tc_Lc_start]]$Y_R)] +
              offset_text),
       labels = bquote(.(label)[.(names(list_tc_Lc_runs)[tc_Lc_start])]))
  seq_tc_Lc <- 1:length(list_tc_Lc_runs)
  seq_tc_Lc <- seq_tc_Lc[-tc_Lc_start]
  for(j in seq_tc_Lc){
    lines(list_tc_Lc_runs[[j]]$FM, list_tc_Lc_runs[[j]]$Y_R, type = 'l',
          ylab = "Y/R", xlab = "F",lty = j)
    text(x = list_tc_Lc_runs[[j]]$FM[length(list_tc_Lc_runs[[j]]$FM)],
         y = (list_tc_Lc_runs[[j]]$Y_R[length(list_tc_Lc_runs[[j]]$Y_R)] +
                offset_text),
         labels = bquote(.(label)[.(names(list_tc_Lc_runs)[j])]))
  } # only works with plotting whole graph if values are not getting bigger right? because otherwise graphs is not insed plotting area
  par(new=T)
  plot(list_tc_Lc_runs[[1]]$FM, list_tc_Lc_runs[[1]]$B_R, type = 'l',
       axes = F, ylab = '', xlab ='', lty = tc_Lc_start,
       col = 'blue')
  axis(side = 4, at = pretty(range(list_tc_Lc_runs[[1]]$B_R)), col = 'blue',
       col.axis = 'blue')
  mtext(side = 4, text = "B/R", line = 3, col = 'blue')
  for(j in seq_tc_Lc){
    lines(list_tc_Lc_runs[[j]]$FM, list_tc_Lc_runs[[j]]$B_R, type = 'l',
          ylab = "Y/R", xlab = "F", col='blue', lty = j)
  } # only works with plotting whole graph if values are not getting bigger right? because otherwise graphs is not insed plotting area
  par(op)
}


