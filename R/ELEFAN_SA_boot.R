#' @title Bootstrapped ELEFAN_SA
#'
#' @description `ELEFAN_SA_boot` performs a bootstrapped fitting of the
#'   von Bertalanffy growth function (VBGF) via the \link{ELEFAN_SA} function.
#'   The length-frequency data is resampled `nresamp` times with
#'   \link{lfqResample} and \link{ELEFAN_SA} is fitted to each resampled data
#'   set. Most of the arguments are simply passed to \link{ELEFAN_SA}.
#'
#' @param lfq a length frequency object of the class `lfq`
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
#' @param parallel logical; should the resamples be fitted in parallel?
#'    TropFishR has to be installed for the parallel version, because each
#'    worker loads it. Depending on the operating system, the argument
#'    `clusterType` can be adjusted. (Default: `parallel = FALSE`)
#' @param nresamp numeric; the number of resamples (Default: `nresamp = 200`)
#' @param no_cores numeric; number of cores used if `parallel = TRUE`
#'    (Default: one less than the number of available cores)
#' @param clusterType type of cluster, see \code{\link[parallel]{makeCluster}}
#'    (Default: `clusterType = "PSOCK"`)
#' @param outfile character; optional name of a text file, which records the
#'    progress of the resample completions. The file is overwritten.
#'    (Default: `outfile = NULL`, no file is written)
#' @param SA_time numeric; Maximum running time in seconds (default : 60 * 1).
#' @param SA_temp numeric; Initial value for temperature (default : 1e5).
#' @param maxit numeric; default: maxit = NULL.
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
#'
#' @return an object of class `lfqBoot`, see \link{ELEFAN_GA_boot}.
#'
#' @references
#' Schwamborn, R., Mildenberger, T. K., Taylor, M. H., 2019. Assessing sources
#' of uncertainty in length-based estimates of body growth in populations of
#' fishes and macroinvertebrates with bootstrapped ELEFAN.
#' \emph{Ecological Modelling}, 393, 37-51
#'
#' @importFrom parallel detectCores
#' @importFrom stats runif
#'
#' @examples
#' \donttest{
#' # load data
#' data(alba)
#'
#' # settings
#' MA <- 7
#' init_par <- list(Linf = 10, K = 1, t_anchor = 0.5)
#' low_par <- list(Linf = 8, K = 0.1, t_anchor = 0)
#' up_par <- list(Linf = 15, K = 5, t_anchor = 1)
#' SA_time <- 2
#' SA_temp <- 1e5
#' nresamp <- 10
#'
#' res <- ELEFAN_SA_boot(lfq = alba, MA = MA, seasonalised = FALSE,
#'   init_par = init_par, up_par = up_par, low_par = low_par,
#'   SA_time = SA_time, SA_temp = SA_temp, seed = 1,
#'   nresamp = nresamp
#' )
#' res$maxDen
#' res$CI
#'
#' # plot resulting distributions
#' univariate_density(res)
#' }
#'
#' @export
#'
ELEFAN_SA_boot <- function(
                           lfq, seasonalised = FALSE,
                           init_par = list(Linf = 50, K = 0.5, t_anchor = 0.5, C = 0, ts = 0),
                           low_par = NULL, up_par = NULL,
                           parallel = FALSE, nresamp = 200,
                           no_cores = max(1, parallel::detectCores() - 1, na.rm = TRUE),
                           clusterType = "PSOCK",
                           outfile = NULL,
                           SA_time = 60, SA_temp = 1e5, maxit = NULL,
                           MA = 5, addl.sqrt = FALSE, agemax = NULL,
                           flagging.out = TRUE, seed = NULL, CI = 95
                           ){

    ## evaluate all arguments: fitOne is sent to the workers if parallel = TRUE,
    ## where unevaluated arguments (e.g. lfq = alba) cannot be found
    invisible(mget(names(formals()), envir = environment()))

    if(is.null(seed)){ seed <- round(runif(1, min = 0, max = 1e6)) }

    fitOne <- function(x){
        set.seed(seed + x)

        ## resample data
        lfqb <- lfqResample(lfq)

        fitboot <- ELEFAN_SA(
            lfqb, seasonalised = seasonalised,
            init_par = init_par,
            low_par = low_par,
            up_par = up_par,
            SA_time = SA_time, SA_temp = SA_temp, maxit = maxit,
            MA = MA, addl.sqrt = addl.sqrt,
            agemax = agemax, flagging.out = flagging.out,
            plot.score = FALSE, verbose = FALSE
        )

        c(unlist(fitboot$par), agemax = fitboot$agemax, seed + x, fitboot$Rn_max)
    }

    res <- bootFits(fitOne, nresamp = nresamp, parallel = parallel,
                    no_cores = no_cores, clusterType = clusterType,
                    outfile = outfile)

    bootSummary(res, CI = CI, agemax = agemax)
}
