#' @title Calculate age-based VBGF parameters from time-based estimates
#'
#' @description Translates t_anchor to t0 and ts to age-adjusted ts. In addition to
#' VBGF parameters calculated via e.g. ELEFAN, the length at recruitment
#' (\code{Lrecr}) must also be provided.
#'
#' @param par list. VBGF parameters.
#' @param Lrecr numeric. Length at recruitment
#' (default: \code{Lrecr = NULL}).
#' @param plot logical. Plot graphical representation of both
#' time and age based growth curves (default: \code{plot = TRUE}).
#'
#' @return a list containging age-based VBGF parameters.
#' @export
#'
#' @examples
#' lfq <- synLFQ4
#' lfq$par <- list(
#'   Linf = 80, K = 0.5, t_anchor = 0.25, C = 0.75, ts = 0.5
#' )
#' t_anchor2t0(par = lfq$par, Lrecr = 10)
#'
t_anchor2t0 <- function(par = NULL, Lrecr = NULL, plot = TRUE){
  if(is.null(par) | is.null(Lrecr)){
    stop("Error: must provide both 'par' and 'Lrecr'")
  }
  if(Lrecr > par$Linf){stop("Error: 'Lrecr' must be lower than 'par$Linf'")}

  # add seasonalized pars if needed
  seasonalized <- !is.null(par$C)
  if(!seasonalized){par$ts <- 0; par$C <- 0}

  t <- seq(0, 3, 0.01)
  Lt <- VBGF(param = par, t = t)
  trecr <- t[min(which((Lt - Lrecr)>0))]
  # calc_trecr(par = par, Lrecr = Lrecr, Lt = tail(Lt,1), t = tail(t,1)) # same
  t0 <- par$t_anchor - trecr
  if(seasonalized){
    ts <- (par$ts - trecr)%%1
  } else {
    ts <- 0
  }
  par_age <- par
  par_age$t_anchor <- NULL
  par_age$t0 <- t0
  par_age$ts <- ts
  age <- seq(t0, 3, 0.01)
  Lage <- VBGF(param = par_age, t = age)

  if(plot){
    op <- par(mfcol = c(1,2), mar = c(4,4,1,1), mgp = c(2,0.5,0), cex = 1)
    # plot by time
    plot(t, Lt, t="l", xlab = "time", ylim = c(0, max(Lt)*1.1), yaxs="i", ylab = bquote(L[t]))
    ylim <- par()$usr[3:4]
    grid(); abline(v = trecr, lty=3); box()
    lines(t, Lt)
    points(trecr, Lrecr, pch = 1, col = 2)
    points(trecr + t0, 0, pch = 4, col = 4)
    tmp <- data.frame(var = c(names(par), "trecr", "Lrecr"), val = c(unlist(par), trecr, Lrecr))
    tmp$pch = NA; tmp$col <- NA
    is_t0 <- which(tmp$var=="t_anchor")
    tmp$pch[is_t0] <- 4
    tmp$col[is_t0] <- 4
    is_Lrecr <- which(tmp$var=="Lrecr")
    tmp$pch[is_Lrecr] <- 1
    tmp$col[is_Lrecr] <- 2
    legend("bottomright",
      legend = c(paste(tmp$var, "=", round(tmp$val,2))),
      bty = "n", cex = 0.7, pt.cex = 1,
      pch = tmp$pch, col = tmp$col
    )
    # plot by age
    plot(age, Lage, t="l", ylim = ylim, yaxs = "i", ylab = bquote(L[age]))
    grid(); abline(v = 0, lty=3); box()
    lines(age, Lage)
    points(0, Lrecr, pch = 1, col = 2)
    points(t0, 0, pch = 4, col = 4)
    tmp <- data.frame(var = c(names(par_age), "Lrecr"), val = c(unlist(par_age), Lrecr))
    tmp$pch = NA; tmp$col <- NA
    is_t0 <- which(tmp$var=="t0")
    tmp$pch[is_t0] <- 4
    tmp$col[is_t0] <- 4
    is_Lrecr <- which(tmp$var=="Lrecr")
    tmp$pch[is_Lrecr] <- 1
    tmp$col[is_Lrecr] <- 2
    legend("bottomright",
      legend = c(paste(tmp$var, "=", round(tmp$val,2))),
      bty = "n", cex = 0.7, pt.cex = 1,
      pch = tmp$pch, col = tmp$col
    )
    par(op)
  }
  return(par_age)
}
