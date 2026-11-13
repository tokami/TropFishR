#' @title Empirical M formula by Then et al. (2015) with uncertainty
#'
#' @description Estimates the natural mortality rate (M) for each resample of
#'   a bootstrapped growth analysis with the growth-based estimator of Then et
#'   al. (2015) ("Then_growth" in \code{\link{M_empirical}}), including the
#'   uncertainty of the estimator itself: for each resample, the model of Then
#'   et al. (2015) with coefficients from a bootstrap of their data is used
#'   for each resample.
#'
#' @param boot an object of class `lfqBoot` with bootstrapped growth parameters
#'    (e.g. from \code{\link{ELEFAN_GA_boot}})
#' @param resid logical; should the residual error of the model of Then et al.
#'    (2015) be included (Default: `FALSE`)? See details.
#' @param seed optional; seed for the random numbers. Resample `x` uses
#'    `set.seed(seed + x)`. Default (`NULL`): the seeds of the resamples in
#'    `boot$seed`.
#' @param CI percentage for confidence intervals (Default: 95)
#'
#' @details The model \eqn{M = a K^b L_\infty^c} of Then et al. (2015) was
#'    refitted (with \code{nls}) to 2000 bootstrap samples of their data of
#'    230 stocks. The data themselves are not included in TropFishR, only the
#'    resulting 2000 sets of coefficients and the residuals of the fit to all
#'    data. For each resample, one of the coefficient sets is selected at
#'    random and M is predicted with the growth parameters of the resample.
#'    This accounts for the uncertainty of the estimated coefficients (i.e. of
#'    the mean M of species with the given growth parameters). With
#'    `resid = TRUE`, the scatter of individual species around the model is
#'    included as well: the prediction is multiplied with a randomly selected
#'    residual of the fit to all data (on the log scale,
#'    \eqn{M_{obs}/M_{fit}}), which results in much wider confidence
#'    intervals. The results are reproducible, because each resample uses its
#'    own seed (see `seed`).
#'
#' @return The \code{lfqBoot} object with the additional column \code{M_Then}
#'    (in \code{bootRaw}, \code{maxDen}, \code{median} and \code{CI}).
#'
#' @references
#' Then, A. Y., J. M. Hoenig, N. G. Hall, D. A. Hewitt. 2015. Evaluating the predictive
#' performance of empirical estimators of natural mortality rate using information on over
#' 200 fish species. ICES J. Mar. Sci. 72: 82-92.
#'
#' @examples
#' data(alba_boot)
#' res <- MstochThen(alba_boot)
#' res$maxDen
#' res$CI
#'
#' ## including the residual error of the model
#' res <- MstochThen(alba_boot, resid = TRUE)
#' res$CI
#'
#' @export
#'
MstochThen <- function(boot, resid = FALSE, seed = NULL, CI = 95){

    checkBoot(boot, c("Linf", "K"))

    ## bootstrapped coefficients of the model of Then et al. (2015) (internal
    ## data, see make_data/makedata_thenBoot.R)
    coefs <- thenBoot$coef
    logres <- thenBoot$logres

    bootRaw <- boot$bootRaw
    n <- nrow(bootRaw)
    seeds <- if(is.null(seed)) boot$seed else seed + seq_len(n)
    if(length(seeds) != n) stop("The boot object has to contain one seed per resample (boot$seed).")

    Ms <- vapply(seq_len(n), function(x){
        set.seed(seeds[x])
        co <- coefs[sample(nrow(coefs), 1), ]
        Mx <- co[["a"]] * bootRaw$K[x]^co[["b"]] * bootRaw$Linf[x]^co[["c"]]
        if(resid) Mx <- Mx * exp(sample(logres, 1))
        Mx
    }, numeric(1))

    bootAddPars(boot, list(M_Then = Ms), CI = CI)
}
