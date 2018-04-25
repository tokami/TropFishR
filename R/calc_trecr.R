
#' Calculate time at recruitment
#'
#' @description The \code{calc_trecr} function calculates the time at recruitment based
#' on \code{\link{VBGF}} parameters, a known length at time \code{t}
#' (i.e. date or numeric class) (\code{Lt}), and a defined length at
#' recruitment \code{Lrecr}. The solution is based on a iterative hindcast of
#' the soVBGF. \code{calc_trecr} is called internally by
#' \code{\link{recruitment2}}
#'
#' @param par \code{\link{VBGF}} parameters
#' @param Lrecr Length at recruitment. This can be e.g. the length at recruitment to the
#' population or fishery (default: \code{Lrecr = 0}).
#' @param Lt known length at time \code{t}
#' @param t time (class is numeric or date)
#' @param tincr resolution to use for hindcasting time at recruitment
#' (default: 0.01). Smaller values give higher precision, but hindcast calculation
#' is slower.
#'
#' @return a numeric value (typically decimal year). When \code{Lrecr > par$Linf} 
#' or \code{Lt > par$Linf} or \code{Lt < Lrecr}, \code{NA} is returned.
#' @export
#'
#' @examples
#' lfq <- synLFQ4
#' lfq$par <- list(
#'   Linf = 80, K = 0.5, t_anchor = 0.25, C = 0.75, ts = 0.5
#' )
#' age <- 1
#' Lt <- VBGF(param = lfq$par, t = age)
#' trecr <- calc_trecr(par = lfq$par, Lrecr = 0, Lt = Lt, t = age)
#' trecr
#' lfq$par$t_anchor # should be equivalent
#'
#'
calc_trecr <- function(
  par = NULL,
  Lrecr = 0,
  Lt = NULL,
  t = NULL,
  tincr = 0.01
){
  # if(Lrecr > par$Linf){stop("Error: 'Lrecr' must be lower than 'par$Linf'")}
  # if(Lt > par$Linf){stop("Error: 'Lt' must be lower than 'par$Linf'")}
  if(Lrecr > par$Linf | Lt > par$Linf | Lt < Lrecr){
    trecr <- NA
  }else{
    Lt.i <- rep(NaN, 500/tincr)
    Lt.i[1] <- Lt
    if(class(t)=="Date"){
      t.i <- seq(date2yeardec(t), by = -tincr, length.out = length(Lt.i))
    }else{
      if(class(t)=="numeric"){
        t.i <- seq(t, by = -tincr, length.out = length(Lt.i))
      }else{
        stop("Error: t must be of class 'Date' or 'numeric'" )
      }
    }
    seasonalized <- !is.null(par$C)
    if(!seasonalized){par$ts <- 0; par$C <- 0}

    i <- 2
    while(i <= length(t.i) & Lt.i[i-1] > Lrecr){
      t2 <- t.i[i-1]
      t1 <- t.i[i]
      slope <- {1 - exp(-(
        par$K*(t2-t1)
        - (((par$C*par$K)/(2*pi))*sin(2*pi*(t1-par$ts)))
        + (((par$C*par$K)/(2*pi))*sin(2*pi*(t2-par$ts)))
      ))}
      Lt.i[i] <- (par$Linf*slope - Lt.i[i-1]) / (slope - 1)
      i <- i + 1
    }
    trecr <- t.i[max(which(Lt.i > Lrecr))]
  }
  # return result
  return(trecr)
}


