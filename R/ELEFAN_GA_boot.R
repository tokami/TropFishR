#' @title Bootstrapped ELEFAN_GA
#'
#' @description `ELEFAN_GA_boot` performs a bootstrapped fitting of the
#'   von Bertalanffy growth function (VBGF) via the \link{ELEFAN_GA} function.
#'   The length-frequency data is resampled `nresamp` times with
#'   \link{lfqResample} and \link{ELEFAN_GA} is fitted to each resampled data
#'   set. Most of the arguments are simply passed to \link{ELEFAN_GA}.
#'
#' @param lfq a length frequency object of the class `lfq`
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
#' @param parallel logical; should the resamples be fitted in parallel? This
#'    differs from the `parallel` argument in \link{ELEFAN_GA} in that it is
#'    not used within the `ga` function for calculation at the population
#'    level, but rather for the resamples. TropFishR has to be installed for
#'    the parallel version, because each worker loads it. Depending on the
#'    operating system, the argument `clusterType` can be adjusted.
#'    (Default: `parallel = FALSE`)
#' @param nresamp numeric; the number of resamples (Default: `nresamp = 200`)
#' @param no_cores numeric; number of cores used if `parallel = TRUE`
#'    (Default: one less than the number of available cores)
#' @param clusterType type of cluster, see \code{\link[parallel]{makeCluster}}
#'    (Default: `clusterType = "PSOCK"`)
#' @param outfile character; optional name of a text file, which records the
#'    progress of the resample completions. The file is overwritten.
#'    (Default: `outfile = NULL`, no file is written)
#' @param popSize the population size. (Default: `popSize = 60`)
#' @param maxiter the maximum number of iterations to run before the GA search is halted.
#'   (Default: `maxiter = 50`)
#' @param run the number of consecutive generations without any improvement
#'   in the best fitness value before the GA is stopped. (Default: `run = 10`)
#' @param pmutation numeric. A small fraction of 1.0. The probability of mutation in a
#'   parent chromosome. Usually mutation occurs with a small probability.
#'   (Default: `pmutation = 0.2`)
#' @param pcrossover the probability of crossover between pairs of chromosomes.
#' Typically this is a large value.  (Default: `pcrossover = 0.8`)
#' @param elitism the number of best fitness individuals to survive at each generation.
#' By default the top 5\% individuals will survive at each iteration.
#' @param MA number indicating over how many length classes the moving average
#' should be performed (default: 5, for more information see \link{lfqRestructure})
#' @param addl.sqrt additional squareroot transformation of positive values
#' according to Brey et al. (1988) (default: FALSE, for
#'    more information see \link{lfqRestructure})
#' @param agemax maximum age of species; default NULL, then estimated from Linf
#' @param flagging.out logical; should positive peaks be flagged out? Original setting of
#' ELEFAN in TRUE. Default:TRUE
#' @param seed seed value for random number reproducibility. Resample `x` uses
#'    the seed `seed + x`. (Default: NULL, a random seed is drawn)
#' @param CI percentage for confidence intervals (Default: 95)
#' @param bin_size optional; bin size passed to \link{lfqModify} to re-bin
#'    each resampled data set before fitting (Default: NA, no re-binning)
#'
#' @return an object of class `lfqBoot`, a list with following objects:
#'    \itemize{
#'    \item \strong{bootRaw} a data.frame of fitted VBGF parameters (columns,
#'       including `phiL` and `agemax`) by resample (rows),
#'    \item \strong{seed} the seed values used for the resamples,
#'    \item \strong{Rn_max} the scores of the fits,
#'    \item \strong{maxDen} the parameters with the maximum multivariate
#'       kernel density,
#'    \item \strong{median} the median parameters,
#'    \item \strong{CI} univariate confidence intervals of the parameters,
#'    \item \strong{multiCI} ranges of the parameters within the multivariate
#'       confidence region,
#'    \item \strong{misc} an empty list for additional information.
#'    }
#'    `phiL` and `agemax` of `maxDen`, `median`, `CI` and `multiCI` are derived
#'    from `Linf` and `K` of the respective estimates.
#'
#' @references
#' Schwamborn, R., Mildenberger, T. K., Taylor, M. H., 2019. Assessing sources
#' of uncertainty in length-based estimates of body growth in populations of
#' fishes and macroinvertebrates with bootstrapped ELEFAN.
#' \emph{Ecological Modelling}, 393, 37-51
#'
#' @importFrom parallel detectCores
#'
#' @examples
#' \donttest{
#' # load data
#' data(alba)
#'
#' # settings
#' MA <- 7
#' low_par <- list(Linf = 8, K = 0.1, t_anchor = 0)
#' up_par <- list(Linf = 15, K = 5, t_anchor = 1)
#' popSize <- 30
#' maxiter <- 20
#' run <- 10
#' pmutation <- 0.2
#' nresamp <- 20
#'
#' res <- ELEFAN_GA_boot(lfq = alba, MA = MA, seasonalised = FALSE,
#'   up_par = up_par, low_par = low_par,
#'   popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
#'   nresamp = nresamp, seed = 1
#' )
#' res$maxDen
#' res$CI
#'
#' # parallel version (TropFishR has to be installed)
#' res <- ELEFAN_GA_boot(lfq = alba, MA = MA, seasonalised = FALSE,
#'   up_par = up_par, low_par = low_par,
#'   popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
#'   nresamp = nresamp, seed = 1, parallel = TRUE, no_cores = 2
#' )
#'
#' # plot resulting distributions
#' univariate_density(res)
#' LinfK_scatterhist(res)
#' }
#'
#' @export
#'
ELEFAN_GA_boot <- function(
                           lfq, seasonalised = FALSE, low_par = NULL, up_par = NULL,
                           parallel = FALSE, nresamp = 200,
                           no_cores = max(1, parallel::detectCores() - 1, na.rm = TRUE),
                           clusterType = "PSOCK",
                           outfile = NULL,
                           popSize = 60, maxiter = 50, run = 10,
                           pmutation = 0.2, pcrossover = 0.8,
                           elitism = base::max(1, round(popSize * 0.05)),
                           MA = 5, addl.sqrt = FALSE, agemax = NULL,
                           flagging.out = TRUE, seed = NULL, CI = 95,
                           bin_size = NA
                           ){

    ## evaluate all arguments: fitOne is sent to the workers if parallel = TRUE,
    ## where unevaluated arguments (e.g. lfq = alba) cannot be found
    invisible(mget(names(formals()), envir = environment()))

    if(is.null(seed)){ seed <- round(runif(1, min = 0, max = 1e6)) }

    fitOne <- function(x){
        set.seed(seed + x)

        ## resample data
        lfqb <- lfqResample(lfq)
        if(!is.na(bin_size) && is.numeric(bin_size)){
            lfqb <- lfqModify(lfqb, bin_size = bin_size)
        }

        fitboot <- ELEFAN_GA(
            lfqb, seasonalised = seasonalised,
            low_par = low_par,
            up_par = up_par,
            popSize = popSize, maxiter = maxiter, run = run,
            pmutation = pmutation, pcrossover = pcrossover, elitism = elitism,
            MA = MA, parallel = FALSE, addl.sqrt = addl.sqrt,
            agemax = agemax, flagging.out = flagging.out,
            plot.score = FALSE, seed = NULL
        )

        c(unlist(fitboot$par), agemax = fitboot$agemax, seed + x, fitboot$Rn_max)
    }

    res <- bootFits(fitOne, nresamp = nresamp, parallel = parallel,
                    no_cores = no_cores, clusterType = clusterType,
                    outfile = outfile)

    bootSummary(res, CI = CI, agemax = agemax)
}
