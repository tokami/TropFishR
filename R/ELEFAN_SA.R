#' @title ELEFAN_SA
#'
#' @description Electronic LEngth Frequency ANalysis with simulated annealing
#'    for estimating growth parameters.
#'
#' @param param a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row) and sampling date (column);
#' }
#' @param seasonalised logical; FALSE
#' @param init_par default: list(Linf = 50, K = 0.5, t_anchor = 0.5, C = 0, ts = 0)
#' @param low_par list, default NULL, e.g. list(Linf = 1, K = 0.01, t_anchor = 0, C = 0, ts = 0)
#' @param up_par list, default NULL, e.g. list(Linf = 1000, K = 10, t_anchor = 1, C = 1, ts = 1)
#' @param SA_time default = 60 * 1
#' @param SA_temp default = 100000
#' @param verbose default   = TRUE
#' @param MA number indicating over how many length classes the moving average should be performed (defalut: 5, for
#'    more information see \link{lfqRestructure}).
#' @param addl.sqrt Passed to \link{lfqRestructure}. Applied an additional square-root transformation of positive values according to Brey et al. (1988).
#'    (default: FALSE, for more information see \link{lfqRestructure}).
#' @param plot logical; indicating if plot with restrucutred frequencies and growth curves should
#'    be displayed
#'
#' @examples
#' data(trout)
#' output <- ELEFAN_SA(trout, SA_time = 60*1)
#'
#' # view fit
#' plot(output, ylim = c(0,15))
#' Lt <- calcLt(output, col=1, par=output$par, draw=TRUE)
#'
#'
#' @details This is cool function
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{samplingPeriod}: length of sampling period in years,
#'   \item \strong{samplingDays}: time of sampling times in relation to first sampling time,
#'   \item \strong{delta_t}: array with time differences between relative sampling time set to zero and
#'      other sampling times,
#'   \item \strong{rcounts}: restructured frequencies,
#'   \item \strong{peaks_mat}: matrix with positive peaks with distinct values,
#'   \item \strong{ASP}: available sum of peaks, sum of posititve peaks which could be potential be hit by
#'      growth curves,
#'   \item \strong{C}: amplitude of growth oscillation,
#'   \item \strong{ts}: winter point winter point (WP = ts + 0.5);
#' when the K-Scan method is applied (fixed Linf) in addition following parameters:
#'   \item \strong{Rn_max}: highest score value,
#' }
#'
#' @importFrom graphics par plot title lines
#' @importFrom stats loess
#'
#' @references
#' Pauly, D. and N. David, 1981. ELEFAN I, a BASIC program for the objective extraction of
#' growth parameters from length-frequency data. \emph{Meeresforschung}, 28(4):205-211
#'
#' @export

ELEFAN_SA <- function(param,
                      seasonalised = FALSE,
                      init_par = list(Linf = 50, K = 0.5, t_anchor = 0.5, C = 0, ts = 0),
                      low_par = NULL, #list(Linf = 1, K = 0.01, t_anchor = 0, C = 0, ts = 0),
                      up_par = NULL, #list(Linf = 1000, K = 10, t_anchor = 1, C = 1, ts = 1),
                      SA_time = 60 * 3,
                      SA_temp = 100000,
                      verbose = TRUE,
                      MA = 5, addl.sqrt = FALSE,
                      plot = FALSE){

  res <- param
  classes <- res$midLengths
  n_classes <- length(classes)
  Linf_est <- classes[n_classes]


  # initial parameters
  init_par_ALL <- list(Linf = Linf_est,
                       K = 0.5, t_anchor = 0.5, C = 0, ts = 0)
  init_Linf <- ifelse("Linf" %in% names(init_par),
                      get("Linf", init_par),
                      get("Linf", init_par_ALL))
  init_K <- ifelse("K" %in% names(init_par),
                   get("K", init_par),
                   get("K", init_par_ALL))
  init_tanc <- ifelse("t_anchor" %in% names(init_par),
                      get("t_anchor", init_par),
                      get("t_anchor", init_par_ALL))
  init_C <- ifelse("C" %in% names(init_par),
                   get("C", init_par),
                   get("C", init_par_ALL))
  init_ts <- ifelse("ts" %in% names(init_par),
                    get("ts", init_par),
                    get("ts", init_par_ALL))

  # lower parameter bounds
  low_par_ALL <- list(Linf = Linf_est * 0.5,
                      K = 0.01,
                      t_anchor = 0,
                      C = 0,
                      ts = 0)
  low_Linf <- ifelse("Linf" %in% names(low_par),
                     get("Linf", low_par),
                     get("Linf", low_par_ALL))
  low_K <- ifelse("K" %in% names(low_par),
                  get("K", low_par),
                  get("K", low_par_ALL))
  low_tanc <- ifelse("t_anchor" %in% names(low_par),
                     get("t_anchor", low_par),
                     get("t_anchor", low_par_ALL))
  low_C <- ifelse("C" %in% names(low_par),
                  get("C", low_par),
                  get("C", low_par_ALL))
  low_ts <- ifelse("ts" %in% names(low_par),
                   get("ts", low_par),
                   get("ts", low_par_ALL))

  # upper parameter bounds
  up_par_ALL <- list(Linf = Linf_est * 1.5,
                     K = 1,
                     t_anchor = 1,
                     C = 1,
                     ts = 1)
  up_Linf <- ifelse("Linf" %in% names(up_par),
                    get("Linf", up_par),
                    get("Linf", up_par_ALL))
  up_K <- ifelse("K" %in% names(up_par),
                 get("K", up_par),
                 get("K", up_par_ALL))
  up_tanc <- ifelse("t_anchor" %in% names(up_par),
                    get("t_anchor", up_par),
                    get("t_anchor", up_par_ALL))
  up_C <- ifelse("C" %in% names(up_par),
                 get("C", up_par),
                 get("C", up_par_ALL))
  up_ts <- ifelse("ts" %in% names(up_par),
                  get("ts", up_par),
                  get("ts", up_par_ALL))


  # ELEFAN 0
  res <- lfqRestructure(res, MA = MA, addl.sqrt = addl.sqrt)
  catch_aAF_F <- res$rcounts
  peaks_mat <- res$peaks_mat
  ASP <- res$ASP

  # seasonalised cost function
  soSAfun <- function(par=c(init_Linf, init_K, init_tanc, init_C, init_ts), lfq){
    Lt <- calcLt(lfq, par=list(Linf=par[1], K=par[2],
                               t_anchor=par[3], C=par[4], ts=par[5]))
    return(-Lt$fESP)
  }
  # cost function
  SAfun <- function(par=c(init_Linf, init_K, init_tanc), lfq){
    Lt <- calcLt(lfq, par=list(Linf=par[1], K=par[2], t_anchor=par[3], C = 0, ts = 0))
    return(-Lt$fESP)
  }

  if(seasonalised){
    # Simulated annealing
    SAfit <- GenSA::GenSA(par = c(init_Linf, init_K, init_tanc, init_C, init_ts),
                          fn = soSAfun,
                          lower = c(low_Linf, low_K, low_tanc, low_C, low_ts),
                          upper = c(up_Linf, up_K, up_tanc, up_C, up_ts),
                          control = list(max.time = SA_time,
                                         temperature = SA_temp,
                                         verbose = verbose),
                          lfq = res)

    pars <- as.list(SAfit$par)
    names(pars) <- c("Linf","K","t_anchor","ts","C")
  }else{
    # Simulated annealing
    SAfit <- GenSA::GenSA(par = c(init_Linf, init_K, init_tanc),
                          fn = SAfun,
                          lower = c(low_Linf, low_K, low_tanc),
                          upper = c(up_Linf, up_K, up_tanc),
                          control = list(max.time = SA_time,
                                         temperature = SA_temp,
                                         verbose = verbose),
                          lfq = res)

    pars <- as.list(SAfit$par)
    names(pars) <- c("Linf","K","t_anchor")
  }

  # Score graph
  fitti <- as.data.frame(SAfit$trace.mat)
  fitti$score <- abs(fitti$current.minimum)
  fitti$jumpscore <- abs(fitti$function.value)
  lo <- stats::loess(jumpscore ~ nb.steps, data=fitti, span = 0.075)
  j <- order(fitti$nb.steps)

  op <- par(mfrow=c(1,1),mar=c(5,5,3,5))
  plot(fitti$nb.steps,fitti$jumpscore, type = "l",  lwd =3,
       xlab = "Time [steps]", ylab = "Score",col="grey")
  title("Score function")
  par(new = TRUE)
  plot(fitti$nb.steps,fitti$temperature, type = "l", col="blue",
       axes = FALSE, xlab="",ylab="", lwd=3)
  axis(4,col="blue", col.ticks = "blue", col.axis="blue")
  mtext(4,line = 2.5,text = "Temperature", col= "blue")
  par(new = TRUE)
  plot(fitti$nb.steps,fitti$jumpscore, type = "n",  lwd = 3,
       xlab = "", ylab = "",col="grey",axes=FALSE)
  lines(fitti$nb.steps[j], lo$fitted[j], col="black", lwd = 3)
  lines(fitti$nb.steps,fitti$score, lwd =3, col ="red")
  par(op)


  beepr::beep(10); beepr::beep(1)  ##beepr::beep(2) #

  ret <- c(res, list(par = pars,
                     Rn_max = abs(SAfit$value)))

  class(ret) <- "lfq"
  if(plot){
    plot(ret, Fname = "rcounts")
    Lt <- calcLt(ret, par = pars)
  }
  return(ret)
}
