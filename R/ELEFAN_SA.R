#' @title ELEFAN_SA
#'
#' @description Electronic LEngth Frequency ANalysis with simulated annealing
#'    for estimating growth parameters.
#'
#' @param lfq a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row)
#'      and sampling date (column);
#' }
#' @param seasonalised logical; indicating if the seasonalised von Bertalanffy
#'    growth function should be applied (default: FALSE).
#' @param init_par a list providing the Initial values for the components to be
#' optimized. When set to NULL the following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is the maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 0.5),
#'   \item \strong{t_anchor} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0.5),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param low_par a list providing the lower bounds for components. When set to
#' NULL the following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is calculated from maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 0.01),
#'   \item \strong{t_anchor} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param up_par a list providing the upper bounds for components. When set to
#' NULL the following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is calculated from maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 0.01),
#'   \item \strong{t_anchor} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param SA_time numeric; Maximum running time in seconds (default : 60 * 1).
#' @param maxit Integer. Maximum number of iterations of the
#'              algorithm. Default is NULL.
#' @param nb.stop.improvement Integer. The program will stop when
#'              there is no any improvement in 'nb.stop.improvement'
#'              steps. Default is NULL
#' @param SA_temp numeric; Initial value for temperature (default : 1e5).
#' @param verbose logical; TRUE means that messages from the algorithm
#'    are shown (default : TRUE).
#' @param MA number indicating over how many length classes the moving average
#'  should be performed (defalut: 5, for
#'    more information see \link{lfqRestructure}).
#' @param addl.sqrt Passed to \link{lfqRestructure}. Applied an additional
#'    square-root transformation of positive values according to Brey et al. (1988).
#'    (default: FALSE, for more information see \link{lfqRestructure}).
#' @param agemax maximum age of species; default NULL, then estimated from Linf
#' @param flagging.out logical; passed to \link{lfqFitCurves}. Default is TRUE
#' @param plot logical; Plot restructured counts with fitted lines using
#' \code{\link{plot.lfq}} and \code{\link{lfqFitCurves}} (default : FALSE).
#' @param plot.score logical; Plot simulated annealing score progression.
#'    (Default: plot.score=TRUE)
#'
#' @examples
#' \donttest{
#' ## synthetic lfq data example
#' data(synLFQ4)
#' plot(synLFQ4, Fname="catch")
#'
#' # ELEFAN_SA (takes approximately 2 minutes)
#' output <- ELEFAN_SA(synLFQ4, SA_time = 60*2, seasonalised = TRUE, MA = 11,
#'   init_par = list(Linf = 75, K = 0.5, t_anchor = 0.5, C = 0.5, ts = 0.5),
#'   low_par = list(Linf = 70, K = 0.3, t_anchor = 0, C = 0, ts = 0),
#'   up_par = list(Linf = 90, K = 0.7, t_anchor = 1, C = 1, ts = 1))
#' output$par
#' output$Rn_max
#'
#' # view fit
#' plot(output)
#'
#' # or
#' plot(output, draw = FALSE)
#' lfqFitCurves(output, col=1, par=output$par, draw=TRUE)$ESP
#'
#' # compare to original parameters
#' tmp <- lfqFitCurves(output, col=4, lty=1,
#'    par=list(Linf=80, K=0.5, t_anchor=0.25, C=0.75, ts=0.5), draw=TRUE)
#' tmp$fESP
#' output$Rn_max
#' }
#'
#' @details A more detailed description of the simulated annealing (SA) can be found in
#'    Xiang et al. (2013). The score value \code{cost_value} is not comparable with
#'    the score value of the other ELEFAN functions (\code{\link{ELEFAN}} or
#'    \code{\link{ELEFAN_GA}}).
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{rcounts}: restructured frequencies,
#'   \item \strong{peaks_mat}: matrix with positive peaks with distinct values,
#'   \item \strong{ASP}: available sum of peaks, sum of posititve peaks which
#'      could be potential be hit by
#'      growth curves,
#'   \item \strong{ncohort}: maximum age of species,
#'   \item \strong{agemax}: maximum age of species,
#'   \item \strong{par}: a list with the parameters of the von Bertalanffy growth
#'      function:
#'      \itemize{
#'        \item \strong{Linf}: length infinity in cm,
#'        \item \strong{K}: curving coefficient;
#'        \item \strong{t_anchor}: time point anchoring growth curves in year-length
#'          coordinate system, corrsponds to peak spawning month,
#'        \item \strong{C}: amplitude of growth oscillation
#'          (if \code{seasonalised} = TRUE),
#'        \item \strong{ts}: summer point of oscillation (ts = WP - 0.5)
#'          (if \code{seasonalised} = TRUE),
#'        \item \strong{phiL}: growth performance index defined as
#'          phiL = log10(K) + 2 * log10(Linf);
#'      }
#'   \item \strong{Rn_max}:  highest score value (absolute value of cost function,
#'   comparable with ELEFAN and ELEFAN_GA).
#' }
#'
#' @importFrom graphics par plot title lines grid
#' @importFrom grDevices adjustcolor
#' @importFrom stats median
#' @importFrom GenSA GenSA
#' @importFrom utils flush.console
#'
#' @references
#' Brey, T., Soriano, M., and Pauly, D. 1988. Electronic length frequency analysis: a revised and
#' expanded
#' user's guide to ELEFAN 0, 1 and 2.
#'
#' Pauly, D. and N. David, 1981. ELEFAN I, a BASIC program for the objective extraction of
#' growth parameters from length-frequency data. \emph{Meeresforschung}, 28(4):205-211
#'
#' Xiang, Y., Gubian, S., Suomela, B., & Hoeng, J. (2013). Generalized simulated
#' annealing for global optimization: the GenSA Package. R Journal, 5(1), 13-28.
#' @export

ELEFAN_SA <- function(lfq,
                      seasonalised = FALSE,
                      init_par = list(Linf = 50, K = 0.5, t_anchor = 0.5, C = 0, ts = 0),
                      low_par = NULL,
                      up_par = NULL,
                      SA_time = 60 * 1,
                      maxit = NULL,
                      nb.stop.improvement = NULL,
                      SA_temp = 1e5,
                      verbose = TRUE,
                      MA = 5, addl.sqrt = FALSE,
                      agemax = NULL,
                      flagging.out = TRUE,
                      plot = FALSE,
                      plot.score = TRUE){

    res <- lfq
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
    soSAfun <- function(lfq, par=c(init_Linf, init_K, init_tanc, init_C, init_ts),
                        agemax, flagging.out){
        Lt <- lfqFitCurves(lfq,
                           par=list(Linf=par[1], K=par[2], t_anchor=par[3], C=par[4], ts=par[5]),
                           flagging.out = flagging.out, agemax = agemax)
        return(-Lt$fESP)
    }
    # cost function
    SAfun <- function(lfq, par=c(init_Linf, init_K, init_tanc),
                      agemax, flagging.out){
        Lt <- lfqFitCurves(lfq,
                           par=list(Linf=par[1], K=par[2], t_anchor=par[3], C = 0, ts = 0),
                           flagging.out = flagging.out, agemax = agemax)
        return(-Lt$fESP)
    }


    ## control list
    control <- list(temperature = SA_temp,
                    verbose = verbose)
    if(!is.null(SA_time)) control$max.time = SA_time
    if(!is.null(maxit)) control$maxit = maxit
    if(!is.null(nb.stop.improvement)) control$nb.stop.improvement


    if(seasonalised){
        # Simulated annealing with seasonalised VBGF
        writeLines(paste(
            "Simulated annealing is running. \nThis will take approximately",
            round(SA_time/60,digits=2),
            "minutes.",
            sep=" "
        ))
        flush.console()
        SAfit <- GenSA::GenSA(
                            par = c(init_Linf, init_K, init_tanc, init_C, init_ts),
                            fn = soSAfun,
                            lower = c(low_Linf, low_K, low_tanc, low_C, low_ts),
                            upper = c(up_Linf, up_K, up_tanc, up_C, up_ts),
                            agemax = agemax,
                            flagging.out = flagging.out,
                            control = control,
                            lfq = res
                        )

        pars <- as.list(SAfit$par)
        names(pars) <- c("Linf","K","t_anchor", "C", "ts")
    }else{
        # Simulated annealing
        writeLines(paste(
            "Simulated annealing is running. \nThis will take approximately",
            round(SA_time/60,digits=2),
            "minutes.",
            sep=" "
        ))
        flush.console()
        SAfit <- GenSA::GenSA(par = c(init_Linf, init_K, init_tanc),
                              fn = SAfun,
                              lower = c(low_Linf, low_K, low_tanc),
                              upper = c(up_Linf, up_K, up_tanc),
                              agemax = agemax,
                              flagging.out = flagging.out,
                              control = control,
                              lfq = res
                              )

        pars <- as.list(SAfit$par)
        names(pars) <- c("Linf","K","t_anchor")
    }

    # Score graph in GA style
    tmp <- as.data.frame(SAfit$trace.mat)
    meani <- aggregate(tmp$function.value, list(step = tmp$nb.steps),mean, na.rm = TRUE)
    exe <- aggregate(tmp$current.minimum, list(step = tmp$nb.steps),mean, na.rm = TRUE)
    medi <- aggregate(tmp$function.value, list(step = tmp$nb.steps),median, na.rm = TRUE)
    ylim <- c(min(range(exe$x,na.rm = TRUE, finite = TRUE)),
              max(range(meani$x, na.rm = TRUE, finite = TRUE)))
    if(plot.score){
        op <- par(mar=c(5.1, 4.1, 1, 4.1))
        plot(tmp$nb.steps, tmp$function.value, type = "n", ylim = ylim, xlab = "Iteration",
             ylab = "Cost value")
        graphics::grid(equilogs = FALSE)
        points(tmp$nb.steps, tmp$current.minimum, type = "o", pch = 16, lty = 1,
               col = "green3", cex = 0.7)
        points(meani$step, meani$x, type = "o", pch = 1, lty = 2,
               col = "dodgerblue3", cex = 0.7)
        polygon(c(meani$step, rev(meani$step)),
                c(exe$x, rev(medi$x)),
                border = FALSE, col = adjustcolor("green3", alpha.f = 0.1))
        par(new=TRUE)
        plot(tmp$nb.steps, tmp$temperature, t="l", col=2, lty=2, log="y", axes = FALSE, xlab = "", ylab = "")
        axis(4, col=2, col.axis=2); mtext(text = "Temperature", side = 4, line = par()$mgp[1], col=2)
        legend("topright", legend = c("Best", "Mean", "Median", "Temperature"),
               col = c("green3", "dodgerblue3", adjustcolor("green3", alpha.f = 0.1), 2),
               pch = c(16, 1, NA, NA), lty = c(1,2,1,2),
               lwd = c(1, 1, 10, 1), pt.cex = c(rep(0.7,2), 2, NA),
               inset = 0.02)
        par(op)
    }

    final_res <- lfqFitCurves(lfq = res,par=pars,flagging.out = flagging.out,
                              agemax = agemax)

    # growth performance index
    phiL <- log10(pars$K) + 2 * log10(pars$Linf)
    pars$phiL <- phiL

    # Results
    res$ncohort <- final_res$ncohort
    res$agemax <- final_res$agemax
    res$par <- pars
    res$fESP <- abs(SAfit$value)
    res$Rn_max <- abs(SAfit$value)

    if(plot){
        plot(res, Fname = "rcounts")
        Lt <- lfqFitCurves(res, par = res$pars, draw=TRUE)
    }
    return(res)
}
