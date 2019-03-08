#' calculate time for new length
#'
#' @param par VBGF parameters. If soVBGF, then ts should relate
#' to time of calendar year rather than time of age year.
#' @param Ltnew New length for which time is calculated
#' @param Lt Known length at time t
#' @param t known time t
#' @param tincr time increment (year fraction) used for iterative
#' hind-/fore-casting. Use smaller values for increased precision
#' (Default: tincr = 0.01)
#'
#' @return numeric value representing year decimal
#' @export
#'
#' @examples
#'
#' data(synLFQ4)
#' synLFQ4$par <- list(Linf = 80, K = 0.5, t_anchor = 0.25, C = 0.75, ts = 0.5)
#'
#' # time at length = 0
#' calc_tnew(par = synLFQ4$par, Lt = 70, t = 2000, Ltnew = 0)
#'
#' # time at length = 75
#' calc_tnew(par = synLFQ4$par, Lt = 70, t = 2000, Ltnew = 75)
#'
calc_tnew <- function(
  par = NULL,
  Ltnew = NULL,
  Lt = NULL,
  t = NULL,
  tincr = 0.01
){
  if(Ltnew == Lt | Ltnew >= par$Linf | Lt >= par$Linf | Lt < 0){
    tnew <- NA
  }else{
    Lt.i <- Lt
    if(class(t)=="Date"){
      t.i <- date2yeardec(t)
    }else{
      if(class(t)=="numeric"){
        t.i <- t
      }else{
        stop("Error: t must be of class 'Date' or 'numeric'" )
      }
    }
    seasonalized <- !is.null(par$C)
    if(!seasonalized){par$ts <- 0; par$C <- 0}

    if(Ltnew < Lt.i){
      while(Lt.i > Ltnew){
        t2 <- t.i
        t1 <- t.i-tincr
        slope <- {1 - exp(-(
          par$K*(t2-t1)
          - (((par$C*par$K)/(2*pi))*sin(2*pi*(t1-par$ts)))
          + (((par$C*par$K)/(2*pi))*sin(2*pi*(t2-par$ts)))
        ))}
        Lt.i <- (par$Linf*slope - Lt.i) / (slope - 1)
        t.i <- t1
      }
      tnew <- t.i
    }

    if(Ltnew > Lt){
      while(Lt.i < Ltnew){
        t2 <- t.i+tincr
        t1 <- t.i
        slope <- {1 - exp(-(
          par$K*(t2-t1)
          - (((par$C*par$K)/(2*pi))*sin(2*pi*(t1-par$ts)))
          + (((par$C*par$K)/(2*pi))*sin(2*pi*(t2-par$ts)))
        ))}
        Lt.i <- (par$Linf*slope + Lt.i) / (slope + 1)
        t.i <- t2
      }
      tnew <- t.i
    }
  }

  # return result
  return(tnew)
}
