#' @title ELEFAN_GA
#'
#' @description Electronic LEngth Frequency ANalysis with genetic algorithm
#' used for estimating growth parameters.
#'
#' @param x a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row)
#'      and sampling date (column);
#' }
#' @param seasonalised logical; indicating if the seasonalised von Bertalanffy
#'    growth function should be applied (default: FALSE).
#' @param low_par a list providing the minimum of the search space in case
#' of real-valued or permutation encoded optimizations. When set to NULL the
#' following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is calculated from maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 0.01),
#'   \item \strong{t_anchor} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 0),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 0),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 0);
#' }
#' @param up_par a list providing the maximum of the search space in case of
#' real-valued or permutation encoded optimizations. When set to NULL the
#' following default values are used:
#'  \itemize{
#'   \item \strong{Linf} length infinity in cm (default is calculated from maximum
#'   length class in the data),
#'   \item \strong{K} curving coefficient (default: 1),
#'   \item \strong{t_anchor} time point anchoring growth curves in year-length
#'   coordinate system, corrsponds to peak spawning month (range: 0 to 1, default: 1),
#'   \item \strong{C} amplitude of growth oscillation (range: 0 to 1, default: 1),
#'   \item \strong{ts} summer point (ts = WP - 0.5) (range: 0 to 1, default: 1);
#' }
#' @param popSize the population size. Default: 50
#' @param maxiter the maximum number of iterations to run before the
#' GA search is halted. default:100
#' @param run the number of consecutive generations without any improvement
#' in the best fitness value before the GA is stopped. Default: equals maxiter
#' @param parallel a logical argument specifying if parallel computing
#' should be used (TRUE) or not (FALSE, default) for evaluating the
#' fitness function. See \code{\link[GA]{ga}} for details. Default:FALSE, but
#' setting to TRUE may substantially improve required calculation time. Use of
#' this functionality requires the following packages: parallel, doParallel.
#' @param pmutation the probability of mutation in a parent chromosome.
#' Usually mutation occurs with a small probability, and by default is set to 0.1.
#' @param pcrossover the probability of crossover between pairs of chromosomes.
#' Typically this is a large value and by default is set to 0.8.
#' @param elitism the number of best fitness individuals to survive at each generation.
#' By default the top 5\% individuals will survive at each iteration.
#' @param MA number indicating over how many length classes the moving average
#' should be performed (default: 5, for
#'    more information see \link{lfqRestructure})
#' @param addl.sqrt additional squareroot transformation of positive values
#' according to Brey et al. (1988) (default: FALSE, for
#'    more information see \link{lfqRestructure})
#' @param agemax maximum age of species; default NULL, then estimated from Linf
#' @param flagging.out logical; should positive peaks be flagged out? Original setting of
#' ELEFAN in TRUE. Default:TRUE
#' @param seed an integer value containing the random number generator state. This
#' argument can be used to replicate the results of a GA search. Note that
#' if parallel computing is required, the doRNG package must be installed.
<<<<<<< HEAD
#' (Default: `seed = NULL`)
=======
>>>>>>> no default printing of score values - potential crash fix
#' @param monitor a logical or an R function which takes as input the current
#'                state of the ‘ga-class’ object and show the evolution of the
#'                search. By default, ‘monitor = FALSE’ so any
#'                output is suppressed. Possible also, the functions
#'                ‘gaMonitor’ or ‘gaMonitor2’ (depending on whether or not is
#'                an RStudio session) which print the average and best fitness
#'                values at each iteration. If set to ‘plot’ these information
#'                are plotted on a graphical device. Other functions can be
#'                written by the user and supplied as argument.
#' @param plot logical; Plot restructured counts with fitted lines using
#' \code{\link{plot.lfq}} and \code{\link{lfqFitCurves}} (default : FALSE).
#' @param plot.score logical; Plot genetic algorithm fitness progression.
#'    (Default: plot.score=TRUE).
#' @param beep logical; Should termination of function result with an audible
#' notifucation sound (Default: FALSE).
#' @param ... additional parameters to pass to \code{\link[GA]{ga}}
#'
#'
#' @details A more detailed description of the generic algorithm (GA) can be found in
#'    Scrucca (2013). The score value \code{fitnessValue} is not comparable with
#'    the score value of the other ELEFAN functions (\code{\link{ELEFAN}} or
#'    \code{\link{ELEFAN_SA}}).
#'
#' @return A list with the input parameters and following list objects:
#' \itemize{
#'   \item \strong{samplingPeriod}: length of sampling period in years,
#'   \item \strong{samplingDays}: time of sampling times in relation to first sampling time,
#'   \item \strong{delta_t}: array with time differences between relative sampling time set to zero and
#'      other sampling times,
#'   \item \strong{rcounts}: restructured frequencies,
#'   \item \strong{peaks_mat}: matrix with positive peaks with distinct values,
#'   \item \strong{ASP}: available sum of peaks, sum of posititve peaks which could be
#'      potential be hit by
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
#'   \item \strong{Rn_max}: highest value of fitness function, (comparable with ELEFAN and ELEFAN_SA).
#' }
#'
#' @examples
#' \donttest{
#' # load data and view catch length frequencies
#' data(synLFQ4)
#' plot(synLFQ4, Fname="catch")
#'
#' # Genetic algorithm
#' # (if using a multicore processor,
#' #   consider adding the argument 'parallel=TRUE'
#' #   to reduce computation time)
#' output <- ELEFAN_GA(synLFQ4, seasonalised = TRUE,
#'    low_par = list(Linf = 70, K = 0.25, t_anchor = 0, C = 0, ts= 0),
#'    up_par = list(Linf = 90, K = 0.7, t_anchor = 1, C = 1, ts = 1),
#'    popSize = 40, maxiter = 50, run = 20,
#'    MA = 11, plot = TRUE, seed = 1111)
#' output$par
#' output$ASP
#' output$Rn_max
#'
#' # compare fitness score (fESP) to
#' # that calculated with "true" growth parameter values
#' plot(output, draw = FALSE)
#' lfqFitCurves(output, par=list(Linf=80, K=0.5, t_anchor=0.25, C=0.75, ts=0.5),
#'        draw = TRUE, col=1, flagging.out = FALSE)$fESP
#' lfqFitCurves(output, par=output$par, draw = TRUE, col=2, flagging.out = FALSE)$fESP
#' legend("top", legend=c("orig.", "GA"), lty=2, col=1:2, ncol=2)
#'}
#'
#' @import parallel
#' @import doParallel
#' @importFrom GA ga
#' @importFrom utils flush.console
#'
#' @references
#' Brey, T., Soriano, M., and Pauly, D. 1988. Electronic length frequency analysis: a
#'    revised and expanded
#' user's guide to ELEFAN 0, 1 and 2.
#'
#' Pauly, D. and N. David, 1981. ELEFAN I, a BASIC program for the objective extraction of
#' growth parameters from length-frequency data. \emph{Meeresforschung}, 28(4):205-211
#'
#' Scrucca, L. (2013). GA: a package for genetic algorithms in R. Journal of
#' Statistical Software, 53(4), 1-37.
#'
#' @export

ELEFAN_GA <- function(
  x,
  seasonalised = FALSE,
  low_par = NULL,
  up_par = NULL,
  popSize = 50,
  maxiter = 100,
  run = maxiter,
  parallel = FALSE,
  pmutation = 0.1,
  pcrossover = 0.8,
  elitism = base::max(1, round(popSize*0.05)),
  MA = 5,
  addl.sqrt = FALSE,
  agemax = NULL,
  flagging.out = TRUE,
  seed = NULL,
  monitor = FALSE,
  plot = FALSE,
  plot.score = TRUE,
  beep = FALSE,
  ...
){

  lfq <- x
  classes <- lfq$midLengths
  n_classes <- length(classes)
  Linf_est <- classes[n_classes]

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
  lfq <- lfqRestructure(lfq, MA = MA, addl.sqrt = addl.sqrt)

  # seasonalised fitness function
  sofun <- function(lfq, par, agemax, flagging.out){
    Lt <- lfqFitCurves(lfq,
                       par=list(Linf=par[1], K=par[2], t_anchor=par[3], C=par[4], ts=par[5]),
                       agemax = agemax, flagging.out = flagging.out)
    return(Lt$fESP)
  }
  # non-seasonalised fitness function
  fun <- function(lfq, par, agemax, flagging.out){
    Lt <- lfqFitCurves(lfq,
                       par=list(Linf=par[1], K=par[2], t_anchor=par[3], C = 0, ts = 0),
                       agemax = agemax, flagging.out = flagging.out)
    return(Lt$fESP)
  }


  # Genetic algorithm
  if(seasonalised){
    min = c(low_Linf, low_K, low_tanc, low_C, low_ts)
    max = c(up_Linf, up_K, up_tanc, up_C, up_ts)
    writeLines("Genetic algorithm is running. This might take some time. \nA beep tone will alert completion.")
    flush.console()
    fit <- GA::ga(
      type = "real-valued",
      fitness = sofun, lfq=lfq,
      min = min,
      max = max,
      agemax = agemax,
      flagging.out = flagging.out,
      popSize = popSize, maxiter = maxiter, run = run, parallel = parallel,
      pmutation = pmutation, pcrossover = pcrossover, elitism = elitism,
      seed = seed, monitor = FALSE,
      ...
    )
    pars <- as.list(fit@solution[1,])
    names(pars) <- c("Linf", "K", "t_anchor", "C", "ts")
  }else{
    min = c(low_Linf, low_K, low_tanc)
    max = c(up_Linf, up_K, up_tanc)
    writeLines("Genetic algorithm is running. This might take some time. \nA beep tone will alert completion.")
    flush.console()
    fit <- GA::ga(
      type = "real-valued",
      fitness = fun,
      lfq=lfq,
      min = min,
      max = max,
      agemax = agemax,
      flagging.out = flagging.out,
      popSize = popSize, maxiter = maxiter, run = run, parallel = parallel,
      pmutation = pmutation, pcrossover = pcrossover, elitism = elitism,
      seed = seed,
      monitor = FALSE,
      ...
    )
    pars <- as.list(fit@solution[1,])
    names(pars) <- c("Linf", "K", "t_anchor")
  }

  # Fitness graph
  if(plot.score){
    GA::plot(fit)
  }

  # notify completion
  if(beep) {beepr::beep(10); beepr::beep(1)}

  final_res <- lfqFitCurves(lfq = lfq, par=pars,
                            flagging.out = flagging.out,
                            agemax = agemax)

  # growth performance index
  phiL <- log10(pars$K) + 2 * log10(pars$Linf)
  pars$phiL <- phiL

  # Results
  ret <- c(lfq, list(ncohort = final_res$ncohort,
                     agemax = final_res$agemax,
                     par = pars,
                     #fitnessValue = fit@fitnessValue,
                     Rn_max = fit@fitnessValue))

  class(ret) <- "lfq"
  if(plot){
    plot(ret, Fname = "rcounts")
    Lt <- lfqFitCurves(ret, par = pars, draw=TRUE)
  }
  return(ret)
}
