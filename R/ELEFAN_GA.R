#' @title ELEFAN_GA
#'
#' @description Electronic LEngth Frequency ANalysis with genetic algorithm
#' used for estimating growth parameters.
#'
#' @param lfq a list consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with catches/counts per length class (row) and sampling date (column);
#' }
#' @param seasonalised logical. Default:FALSE
#' @param min a vector of length equal to the decision variables
#' providing the minimum of the search space in case of real-valued
#' or permutation encoded optimizations.
#' @param max a vector of length equal to the decision variables
#' providing the maximum of the search space in case of real-valued
#' or permutation encoded optimizations.
#' @param popSize the population size.
#' @param maxiter the maximum number of iterations to run before the
#' GA search is halted. default:100
#' @param run the number of consecutive generations without any improvement
#' in the best fitness value before the GA is stopped. Default:20
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
#' @param k number indicating over how many length classes the moving average
#' should be performed (default: 5)
#' @param addl.sqrt additional squareroot transformation of positive values
#' according to Brey et al. (1988) (default: FALSE).
#' @param flagging.out logical; should positive peaks be flagged out? Original setting of
#' ELEFAN in TRUE. Default:TRUE
#' @param plot logical. Plot restructured counts with fitted lines using
#' \code{\link{plot.lfq}} and \code{\link{calcLt}}
#' @param ... additional parameters to pass to \code{\link[GA]{ga}}
#'
#' @return lfq with restructured counts, fitted parameters, and fitness score
#' @export
#'
#' @examples
#' \dontrun{
#'
#' # load data and view catch length frequencies
#' data(synLFQ4)
#' plot(synLFQ4, Fname="catch")
#'
#' # Genetic algorithm
#' res <- ELEFAN_GA(
#'   synLFQ4,
#'   seasonalised = TRUE,
#'   min = c(70, 0.25, 0, 0, 0),
#'   max = c(90, 0.7, 1, 1, 1),
#'   parallel = TRUE,
#'   popSize = 40,
#'   maxiter = 50,
#'   run = 20,
#'   k = 11, addl.sqrt = FALSE,
#'   flagging.out = FALSE,
#'   plot = TRUE,
#'   seed = 1111
#' )
#' res$par
#' res$ASP
#' res$fitnessValue
#'
#' # compare fitness score (ESP) to
#' # that calculated with "true" growth parameter values
#' plot(res)
#' calcLt(res, par=list(Linf=80, K=0.5, t_anchor=0.25, C=0.75, ts=0),
#'        draw = TRUE, col=1, flagging.out = FALSE)$ESP
#' calcLt(res, par=res$par, draw = TRUE, col=2, flagging.out = FALSE)$ESP
#' legend("top", legend=c("orig.", "GA"), lty=2, col=1:2, ncol=2)
#'
#' }
#'


ELEFAN_GA <- function(
  lfq,
  seasonalised = FALSE,
  min = NULL,
  max = NULL,
  popSize = 100, maxiter = 100, run = 20, parallel = FALSE,
  pmutation = 0.3, pcrossover = 0.8, elitism = 10,
  k = 5, addl.sqrt = FALSE,
  flagging.out = TRUE,
  plot = FALSE,
  ...
){

  # ELEFAN 0
  lfq <- lfqRestructure2(lfq, k = k, addl.sqrt = addl.sqrt)

  # seasonalised fitness function
  sofun <- function(par, lfq){
    Lt <- calcLt(
      lfq,
      par=list(Linf=par[1], K=par[2], t_anchor=par[3], C=par[4], ts=par[5]),
      flagging.out = flagging.out
    )
    return(Lt$ESP)
  }
  # non-seasonalised fitness function
  fun <- function(par, lfq){
    Lt <- calcLt(
      lfq,
      par=list(Linf=par[1], K=par[2], t_anchor=par[3], C = 0, ts = 0),
      flagging.out = flagging.out
    )
    return(Lt$ESP)
  }


  # Genetic algorithm
  if(seasonalised){
    fit <- GA::ga(
      type = "real-valued",
      fitness = sofun, lfq=lfq,
      min = min, max = max,
      popSize = popSize, maxiter = maxiter, run = run, parallel = parallel,
      pmutation = pmutation, pcrossover = pcrossover, elitism = elitism,
      ...
    )
    pars <- as.list(fit@solution[1,])
    names(pars) <- c("Linf", "K", "t_anchor", "C", "ts")
  }else{
    fit <- GA::ga(
      type = "real-valued",
      fitness = fun, lfq=lfq,
      min = min, max = max,
      popSize = popSize, maxiter = maxiter, run = run, parallel = parallel,
      pmutation = pmutation, pcrossover = pcrossover, elitism = elitism,
      ...
    )
    pars <- as.list(fit@solution[1,])
    names(pars) <- c("Linf", "K", "t_anchor")
  }

  # Fitness graph
  GA::plot(fit)

  # notify completion
  beepr::beep(10); beepr::beep(1)

  # Results
  ret <- c(
    lfq,
    list(par = pars, fitnessValue = fit@fitnessValue)
  )
  class(ret) <- "lfq"
  if(plot){
    plot(ret, Fname = "rcounts")
    Lt <- calcLt(ret, par = pars, draw=TRUE)
  }
  return(ret)
}
