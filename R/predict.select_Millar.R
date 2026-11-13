#' @title Predict the retention of Millar's selectivity model
#'
#' @description Relative retention of a fitted selectivity model
#'    (\code{\link{select_Millar}}) for any lengths and mesh sizes.
#'
#' @param object a list of the class \code{"select_Millar"} with the results
#'    of \code{\link{select_Millar}}.
#' @param lengths lengths for which the retention is predicted (default: the
#'    length classes of the fitted data).
#' @param meshSizes mesh sizes for which the retention is predicted (default:
#'    the mesh sizes of the fitted data). The parameters of the model refer to
#'    the first (smallest) mesh size of the fitted data, which stays the
#'    reference mesh size for new mesh sizes.
#' @param rel.power optional; relative fishing power of the mesh sizes (one
#'    value per mesh size). By default (\code{NULL}), the relative power of the
#'    fit is used if \code{meshSizes} are the fitted mesh sizes, otherwise 1.
#' @param standardise logical; should the retention be divided by its maximum
#'    (Default: FALSE)?
#' @param ... unused.
#'
#' @details The retention curves have a maximum of 1 (except bi-modal
#'    curves and \code{"tt.logistic"}), multiplied with the relative power of
#'    the mesh size. For \code{"tt.logistic"}, the retention of the control
#'    gear (mesh size equal to the reference mesh size) is \code{1 - p} and
#'    of the experimental gear \code{p} times the logistic curve.
#'
#' @return A matrix with the relative retention with one row per length and
#'    one column per mesh size.
#'
#' @examples
#' data(gillnet)
#' output <- select_Millar(gillnet, x0 = NULL, rtype = "lognorm", plot = FALSE)
#'
#' ## retention at the fitted mesh sizes
#' predict(output, lengths = seq(40, 90, 10))
#'
#' ## retention of a new mesh size
#' predict(output, lengths = seq(40, 90, 10), meshSizes = 16)
#'
#' @seealso \code{\link{select_Millar}}, \code{\link{select_Millar_curve}}
#'
#' @importFrom stats predict
#' @method predict select_Millar
#' @export

predict.select_Millar <- function(object, lengths = object$midLengths,
                                  meshSizes = object$meshSizes,
                                  rel.power = NULL, standardise = FALSE, ...){

    if(is.null(rel.power)){
        rel.power <- if(length(meshSizes) == length(object$meshSizes) &&
                        all(meshSizes == object$meshSizes) &&
                        !is.null(object$rel.power))
                         object$rel.power else rep(1, length(meshSizes))
    }
    if(length(rel.power) != length(meshSizes))
        stop("'rel.power' has to have one value per mesh size.")

    ret <- millarRetention(lengths, meshSizes, object$rtype, object$par,
                           refMesh = object$meshSizes[1])
    ret <- t(t(ret) * rel.power)
    if(standardise) ret <- ret / max(ret, na.rm = TRUE)
    return(ret)
}



#' @title Selection curves of Millar's selectivity models
#'
#' @description Relative retention for given parameters of the selection
#'    curves of Millar's selectivity model (\code{\link{select_Millar}}), e.g.
#'    published gillnet selectivity parameters.
#'
#' @param lengths lengths for which the retention is calculated.
#' @param meshSizes mesh sizes for which the retention is calculated.
#' @param rtype type of the selection curve:
#'    \code{"norm.loc"} for a normal curve with common spread,
#'    \code{"norm.sca"} for a normal curve with variable spread,
#'    \code{"lognorm"} for a lognormal curve,
#'    \code{"binorm.sca"} for a bi-normal curve,
#'    \code{"bilognorm"} for a bi-lognormal curve,
#'    \code{"tt.logistic"} for a logistic curve,
#'    \code{"gamma"} for a gamma curve.
#' @param par parameters of the selection curve for the reference mesh size
#'    (\code{refMesh}), as reported by \code{\link{select_Millar}} in
#'    \code{estimates} (see Details). Either a numeric vector or the matrix
#'    \code{estimates} of \code{\link{select_Millar}}.
#' @param refMesh the mesh size to which the parameters refer (default: the
#'    first value of \code{meshSizes}).
#'
#' @details The parameters (\code{par}) are, in this order:
#' \itemize{
#'   \item \code{"norm.loc"}, \code{"norm.sca"}, \code{"lognorm"},
#'      \code{"gamma"}: mode and standard deviation of the selection curve of
#'      the reference mesh size;
#'   \item \code{"binorm.sca"}, \code{"bilognorm"}: mode and standard deviation
#'      of the first curve, mode and standard deviation of the second curve,
#'      proportion of the first curve (p);
#'   \item \code{"tt.logistic"}: length at 50\% retention (L50) and selection
#'      range (SR = L75 - L25); a third value (p, the split between control
#'      and experimental gear) is ignored.
#' }
#'    The selection curves of other mesh sizes follow from the principle of
#'    geometric similarity: the mode and standard deviation (\code{"norm.sca"},
#'    \code{"lognorm"}, \code{"gamma"}, \code{"binorm.sca"},
#'    \code{"bilognorm"}) or only the mode (\code{"norm.loc"}) are
#'    proportional to the mesh size. The \code{"tt.logistic"} curve is the same
#'    for all mesh sizes.
#'
#'    The curves are scaled to a maximum retention of 1 (except bi-modal
#'    curves). Note that the parameters in \code{par} of the result of
#'    \code{\link{select_Millar}} are on an internal scale (e.g. on the log
#'    scale for \code{"lognorm"}); use \code{estimates} or
#'    \code{\link{predict.select_Millar}} for a fitted model.
#'
#' @return A matrix with the relative retention with one row per length and
#'    one column per mesh size.
#'
#' @examples
#' ## lognormal selection curve with mode 55 cm and SD 4.5 cm for a mesh
#' ## size of 13.5 cm, for several mesh sizes
#' sel <- select_Millar_curve(lengths = 40:90, meshSizes = c(13.5, 15, 17),
#'                            rtype = "lognorm", par = c(55, 4.5))
#' matplot(40:90, sel, type = "l", lty = 1, xlab = "Length",
#'         ylab = "Relative retention")
#'
#' ## the same with the estimates of a fitted model
#' data(gillnet)
#' output <- select_Millar(gillnet, x0 = NULL, rtype = "lognorm", plot = FALSE)
#' select_Millar_curve(lengths = c(50, 60, 70), meshSizes = 16,
#'                     rtype = "lognorm", par = output$estimates,
#'                     refMesh = gillnet$meshSizes[1])
#'
#' @seealso \code{\link{select_Millar}}, \code{\link{predict.select_Millar}}
#'
#' @references
#'  Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#'  54(3), 471-477.
#'
#' @importFrom stats uniroot
#' @export

select_Millar_curve <- function(lengths, meshSizes, rtype, par,
                                refMesh = meshSizes[1]){

    if(is.matrix(par)) par <- if("par" %in% colnames(par)) par[, "par"] else par[, 1]
    par <- as.numeric(par)
    if(!is.numeric(refMesh) || length(refMesh) != 1 || !is.finite(refMesh) || refMesh <= 0)
        stop("'refMesh' has to be a single mesh size > 0.")
    th <- millarPar(rtype, par, refMesh)

    if(rtype == "tt.logistic"){
        wk <- exp(th[1] + th[2] * lengths)
        ret <- matrix(wk / (1 + wk), nrow = length(lengths), ncol = length(meshSizes))
        dimnames(ret) <- list(lengths, meshSizes)
        return(ret)
    }
    millarRetention(lengths, meshSizes, rtype, th, refMesh = refMesh)
}



## Retention matrix (lengths x mesh sizes) of a selection curve with the
## internal parameters th, relative to the reference mesh size refMesh
## (the curve functions of rtypes_Millar use the first mesh size as
## reference)
millarRetention <- function(lengths, meshSizes, rtype, th, refMesh){
    r <- rtypes_Millar(rtype)
    ret <- outer(lengths, c(refMesh, meshSizes), r, th)[, -1, drop = FALSE]
    dimnames(ret) <- list(lengths, meshSizes)
    ret
}


## Convert the parameters of a selection curve as reported in the
## estimates of select_Millar (mode, SD, p; L50, SR) to the internal
## parameters of rtypes_Millar
millarPar <- function(rtype, par, refMesh){
    nPar <- switch(rtype,
                   "norm.loc" = 2, "norm.sca" = 2, "lognorm" = 2, "gamma" = 2,
                   "binorm.sca" = 5, "bilognorm" = 5, "tt.logistic" = c(2, 3),
                   stop(paste0(rtype, " not recognised, possible curve types are ",
                               "\"norm.loc\", \"norm.sca\", \"lognorm\", \"binorm.sca\", ",
                               "\"bilognorm\", \"tt.logistic\" and \"gamma\".")))
    if(!(length(par) %in% nPar))
        stop(paste0("'par' has to have ", paste(nPar, collapse = " or "),
                    " values for rtype = \"", rtype, "\", but has ", length(par), "."))
    if(any(!is.finite(par)))
        stop("'par' has to contain finite values.")

    ## lognormal: mode = exp(mu - s^2), sd = sqrt(exp(2 mu + s^2) (exp(s^2) - 1))
    lognormPar <- function(mode, sd){
        if(mode <= 0 || sd <= 0) stop("Mode and standard deviation have to be > 0.")
        R2 <- (sd / mode)^2
        ## v = exp(s^2) solves v^3 (v - 1) = (sd / mode)^2
        v <- uniroot(function(v) v^3 * (v - 1) - R2, c(1, 2 + sqrt(R2)),
                     tol = 1e-14)$root
        s2 <- log(v)
        c(log(mode) + s2, sqrt(s2))
    }
    logit <- function(p){
        if(p <= 0 || p >= 1) stop("The proportion p has to be > 0 and < 1.")
        log(p / (1 - p))
    }

    switch(rtype,
           "norm.loc" = par,
           "norm.sca" = par,
           "lognorm" = lognormPar(par[1], par[2]),
           "binorm.sca" = c(par[1:4], logit(par[5])),
           "bilognorm" = c(lognormPar(par[1], par[2]), lognormPar(par[3], par[4]),
                           logit(par[5])),
           "tt.logistic" = {
               if(par[2] <= 0) stop("The selection range (SR) has to be > 0.")
               b <- 2 * log(3) / par[2]
               c(-par[1] * b, b, if(length(par) == 3) logit(par[3]) else 0)
           },
           "gamma" = {
               ## mode = (k - 1) alpha m, sd = sqrt(k) alpha m
               if(par[1] <= 0 || par[2] <= 0)
                   stop("Mode and standard deviation have to be > 0.")
               r <- par[2] / par[1]
               k <- ((1 + sqrt(1 + 4 * r^2)) / (2 * r))^2
               c(k, par[1] / ((k - 1) * refMesh))
           })
}
