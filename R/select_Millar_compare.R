#' @title Compare Millar's selectivity models
#'
#' @description Compares several selection curves fitted with
#'    \code{\link{select_Millar}} to the same data by the (quasi-) Akaike
#'    information criterion corrected for small sample sizes (AICc and
#'    QAICc).
#'
#' @param ... results of \code{\link{select_Millar}} (class
#'    \code{"select_Millar"}) fitted to the same data, or one list of such
#'    results. The names of the arguments (or of the list) are used as model
#'    names, otherwise the curve types (\code{rtype}).
#' @param c_hat optional; overdispersion factor. By default (\code{NULL}), it
#'    is estimated from the model with the most parameters (see Details).
#' @param c_hat_method method to estimate \code{c_hat}: \code{"deviance"}
#'    (default; deviance / degrees of freedom) or \code{"pearson"} (Pearson
#'    chi-square statistic / degrees of freedom).
#'
#' @details Catches of gillnet and hook experiments are often overdispersed,
#'    i.e. more variable than expected from the multinomial model of
#'    \code{\link{select_Millar}}. The QAICc accounts for this by dividing
#'    the log-likelihood by the overdispersion factor \code{c_hat}
#'    (Burnham and Anderson 2002):
#'
#'    QAICc = -2 logL / c_hat + 2 K + 2 K (K + 1) / (n - K - 1),
#'
#'    where logL is the log-likelihood of the model, K the number of
#'    parameters (plus 1 for \code{c_hat} if \code{c_hat > 1}), and n the
#'    number of observations, here the number of non-zero cells of the
#'    catch matrix (length classes x mesh sizes). The models can only be
#'    compared with the same \code{c_hat}, which is therefore estimated from
#'    the most general model (with the most parameters; if several models
#'    have the most parameters, the one with the lowest deviance) and set to 1
#'    if it is smaller than 1. The AICc is the same criterion without overdispersion
#'    (\code{c_hat = 1}, K without the additional parameter).
#'
#'    The log-likelihood is the multinomial log-likelihood of the catches
#'    conditional on the total catch per length class, without constant
#'    terms. Therefore, only differences between models fitted to the same
#'    data are meaningful. The model with the lowest QAICc (\code{delta = 0})
#'    is the best supported model; models with \code{delta} < 2 have similar
#'    support (Burnham and Anderson 2002).
#'
#' @return A data frame with one row per model, sorted by the QAICc, with the
#'    columns \code{model}, \code{rtype}, number of parameters (\code{npar}),
#'    \code{deviance}, degrees of freedom (\code{dof}), \code{c_hat},
#'    \code{AICc}, \code{QAICc}, and the difference to the lowest QAICc
#'    (\code{delta}).
#'
#' @examples
#' data(gillnet)
#' fits <- list(
#'   norm.loc = select_Millar(gillnet, x0 = NULL, rtype = "norm.loc", plot = FALSE),
#'   norm.sca = select_Millar(gillnet, x0 = NULL, rtype = "norm.sca", plot = FALSE),
#'   lognorm = select_Millar(gillnet, x0 = NULL, rtype = "lognorm", plot = FALSE),
#'   binorm.sca = select_Millar(gillnet, x0 = c(55, 4, 65, 4, 3),
#'                              rtype = "binorm.sca", plot = FALSE),
#'   bilognorm = select_Millar(gillnet, x0 = c(4, 0.2, 4.2, 0.1, 2),
#'                             rtype = "bilognorm", plot = FALSE))
#' select_Millar_compare(fits)
#'
#' ## AIC of a single fit (without overdispersion)
#' AIC(fits$lognorm)
#'
#' @seealso \code{\link{select_Millar}}
#'
#' @references
#'  Burnham, K. P., Anderson, D. R., 2002. Model selection and multimodel
#'  inference: a practical information-theoretic approach. 2nd edition.
#'  Springer, New York.
#'
#'  Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook selectivity
#'  using log-linear models. \emph{ICES Journal of Marine Science: Journal du Conseil},
#'  54(3), 471-477.
#'
#' @export

select_Millar_compare <- function(..., c_hat = NULL, c_hat_method = "deviance"){

    fits <- list(...)
    if(length(fits) == 1 && !inherits(fits[[1]], "select_Millar") && is.list(fits[[1]]))
        fits <- fits[[1]]
    if(length(fits) < 1 || !all(vapply(fits, inherits, logical(1), "select_Millar")))
        stop("Please provide results of select_Millar (class 'select_Millar').")
    c_hat_method <- match.arg(c_hat_method, c("deviance", "pearson"))

    ## all fits have to be based on the same data
    same <- function(a, b) isTRUE(all.equal(a, b, check.attributes = FALSE))
    ref <- fits[[1]]
    for(f in fits[-1]){
        if(!same(f$midLengths, ref$midLengths) || !same(f$meshSizes, ref$meshSizes) ||
           !same(f$CatchPerNet_mat, ref$CatchPerNet_mat) ||
           !same(f$rel.power, ref$rel.power))
            stop("The models have to be fitted to the same data (midLengths, meshSizes, catches and rel.power).")
    }

    mnames <- names(fits)
    rtypes <- vapply(fits, function(f) f$rtype, character(1))
    if(is.null(mnames)) mnames <- rep("", length(fits))
    mnames[mnames == ""] <- rtypes[mnames == ""]

    npar <- vapply(fits, function(f) length(f$par), numeric(1))
    nll <- vapply(fits, function(f) f$value, numeric(1))
    stat <- function(f, nm) f$out[nm, 1]
    dev <- vapply(fits, stat, numeric(1), nm = "Deviance")
    pchi <- vapply(fits, stat, numeric(1), nm = "Pearson.chisq")
    dof <- vapply(fits, stat, numeric(1), nm = "d.o.f.")
    n <- millarNobs(ref)

    ## common overdispersion factor from the most general model
    if(is.null(c_hat)){
        ## most parameters; with ties the best fitting of these models
        gen <- which(npar == max(npar))
        gen <- gen[which.min(dev[gen])]
        c_hat <- (if(c_hat_method == "deviance") dev[gen] else pchi[gen]) / dof[gen]
        c_hat <- max(1, unname(c_hat))
    }else if(!is.numeric(c_hat) || length(c_hat) != 1 || !is.finite(c_hat) || c_hat <= 0){
        stop("'c_hat' has to be a single number > 0.")
    }

    aicc <- function(m2ll, K){
        if(any(n - K - 1 <= 0))
            stop(paste0("The number of non-zero cells (", n,
                        ") is too small for the number of parameters."))
        m2ll + 2 * K + 2 * K * (K + 1) / (n - K - 1)
    }
    Kq <- npar + as.numeric(c_hat > 1)
    tab <- data.frame(model = unname(mnames), rtype = unname(rtypes),
                      npar = unname(npar), deviance = unname(dev),
                      dof = unname(dof), c_hat = c_hat,
                      AICc = unname(aicc(2 * nll, npar)),
                      QAICc = unname(aicc(2 * nll / c_hat, Kq)),
                      stringsAsFactors = FALSE)
    tab$delta <- tab$QAICc - min(tab$QAICc)
    tab <- tab[order(tab$QAICc), ]
    rownames(tab) <- NULL
    return(tab)
}


#' @title Log-likelihood of Millar's selectivity model
#'
#' @description Log-likelihood of a selection curve fitted with
#'    \code{\link{select_Millar}}, so that \code{AIC} and \code{BIC} can be
#'    used. The log-likelihood is the multinomial log-likelihood conditional
#'    on the total catch per length class, without constant terms; the
#'    number of observations is the number of non-zero cells of the catch
#'    matrix. For overdispersed data use \code{\link{select_Millar_compare}}.
#'
#' @param object a list of the class \code{"select_Millar"}.
#' @param ... unused.
#'
#' @return An object of class \code{"logLik"}.
#'
#' @importFrom stats logLik
#' @method logLik select_Millar
#' @export

logLik.select_Millar <- function(object, ...){
    structure(-object$value, df = length(object$par),
              nobs = millarNobs(object), class = "logLik")
}


## number of observations of a select_Millar fit: non-zero cells of the
## catch matrix
millarNobs <- function(fit){
    sum(fit$CatchPerNet_mat > 0, na.rm = TRUE)
}
