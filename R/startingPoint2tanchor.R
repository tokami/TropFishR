#' @title Convert FiSAT's starting point to t_anchor value
#'
#' @description Starting points returned or chosen within FiSAT are not supported
#'    in TropFishR. Instead \code{t_anchor} takes on the job of anchoring VBGF growth curves on
#'    a temporal axis. This function allows to convert FiSAT's starting points to \code{t_anchor} values
#'
#' @param param list with dates, midLengths, and catch
#' @param par list with growth parameters 'Linf' and 'K' of VBGF
#' @param startingLength starting length as returned by FiSAT, indicating the length within
#'    the starting sample cut by a growth curve
#' @param startingSample starting sample as returned by FiSAT, indicating the sample which is
#'    cut by a growth curve
#'
#' @keywords function lfq startingPoints t_anchor
#'
#' @examples
#' data(synLFQ5)
#' lfqNEW <- startingPoint2tanchor(synLFQ5, par = list(Linf = 92, K = 0.37),
#'    startingLength = 31, startingSample = 4)
#' lfqRest <- lfqRestructure(lfqNEW, MA = 11)
#' plot(lfqRest,par=list(Linf=lfqRest$Linf,K=lfqRest$K,t_anchor=lfqRest$t_anchor))
#'
#' @return list with input elements and estimated t_anchor value
#'
#' @export

startingPoint2tanchor <- function(param, par, startingLength, startingSample){

  res <- param
  if(("Linf" %in% names(par)) == FALSE) stop(noquote("Please add Linf to param!"))
  if(("K" %in% names(par)) == FALSE) stop(noquote("Please add K to param!"))
  if(("dates" %in% names(res)) == FALSE) stop(noquote("Please add a dates to param!"))

  if(length(res$dates) < startingSample) stop(noquote(paste("Only",
                                                            length(res$dates),
                                                            "sample dates in param$dates, but startingSample is",
                                                            startingSample)))

  tx <- res$dates[startingSample]
  tx <- date2yeardec(tx)
  agex <- VBGF(L = startingLength, param = list(Linf = par$Linf, K = par$K, t0 = 0))
  t_anchor <- (tx-agex) %% floor(tx-agex)

  ret <- c(res, list(Linf = par$Linf, K = par$K, t_anchor = t_anchor))
  return(ret)
}
