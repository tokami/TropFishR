#' @title VBGF plot and CI from bootstrapping results
#'
#' @description Plots the growth curves of all resamples of a bootstrapped
#'   ELEFAN fit (\link{ELEFAN_GA_boot}, \link{ELEFAN_SA_boot}) over relative
#'   time together with the confidence intervals and the growth curve with the
#'   maximum multivariate kernel density.
#'
#' @param res list of class `lfqBoot`
#' @param CI vector. Confidence interval levels to plot
#' @param agemax numeric. Maximum number of years to project.
#' @param add_legend logical. Should CI and max. density legend be added
#'   (Default: `add_legend = TRUE`).
#' @param add_max_dens_legend logical. Should maximum density line be added
#'   (Default: `add_max_dens_legend = TRUE`).
#' @param xlab Label for x-axis
#' @param ylab Label for y-axis
#' @param perm.col Color for each resample estimate line. See `?par`.
#' @param perm.lwd Line width for each resample estimate line. See `?par`.
#' @param ci.col Color for CI line. See `?par`.
#' @param ci.lty Line type for CI line. See `?par`.
#' @param ci.lwd Line width for CI line. See `?par`.
#' @param maxd.col Color for maximum density line. See `?par`.
#' @param maxd.lty Line type maximum density line. See `?par`.
#' @param maxd.lwd Line width maximum density line. See `?par`.
#'
#' @return plot and list containing: `limCI` - data.frame with CI limits by time,
#'   `inCi` -  data.frame with logical values defining whether bootstrapping samples
#'   are within each of the defined CIs, `density` - the multivariate kernel density
#'   estimates for each sample, and `max_dens` is a list with the VBGF parameter
#'   combination having the maximum density estimate.
#'
#' @importFrom grDevices adjustcolor
#' @importFrom graphics plot lines box legend
#' @importFrom stats cov quantile
#'
#' @examples
#' data(alba_boot)
#' CIinfo <- vbgfCI_time(
#'   res = alba_boot,
#'   agemax = 2, CI = 95,
#'   perm.col = adjustcolor("grey50",0.2)
#' )
#'
#' # plot more CI levels
#' CIinfo <- vbgfCI_time(
#'   res = alba_boot,
#'   agemax = 2, CI = c(95, 50),
#'   ci.lty = 1, ci.lwd = 2, ci.col = c("red", "orange"),
#'   perm.col = adjustcolor("grey50",0.2)
#' )
#'
#' # using output in lfq plot (see ?plot.lfq)
#' data(alba)
#' alba <- lfqRestructure(alba, MA = 7)
#' plot(alba, Fname = "rcounts")
#' for(i in seq(nrow(alba_boot$bootRaw))){
#'   x <- as.list(alba_boot$bootRaw[i,])
#'   tmp <- lfqFitCurves(
#'     lfq = alba, par = x,
#'     col = adjustcolor("grey50",0.1), draw = TRUE, lty=1
#'   )
#' }
#' tmp <- lfqFitCurves(lfq = alba, par = CIinfo$max_dens,
#'   col = 1, draw = TRUE, lty=1, lwd=2
#' )
#'
#' @export
#'
vbgfCI_time <- function(res, CI = 95, agemax = NULL,
                        add_legend = TRUE, add_max_dens_legend = TRUE,
                        xlab = "Relative time", ylab = "Length",
                        perm.col = adjustcolor("grey50",0.1), perm.lwd = 1,
                        ci.col = 1, ci.lty = 2, ci.lwd = 1,
                        maxd.col = 1, maxd.lty = 1, maxd.lwd = 2
                        ){

    res <- res$bootRaw

    if(is.null(agemax)){agemax <- max(ceiling((1/-res$K)*log(1-((res$Linf*0.95)/res$Linf))))}

    ## expand ci line attributes to length of CI (if needed, values are recycled)
    ci.col <- rep_len(ci.col, length(CI))
    ci.lty <- rep_len(ci.lty, length(CI))
    ci.lwd <- rep_len(ci.lwd, length(CI))

    ## remove derived parameters (phi', agemax) if included in res
    x <- res[,!(names(res) %in% c("phiL", "agemax"))]

    ## First a fitting round to look for shifts in t_anchor
    age = seq(0, agemax, 0.01)
    Lt0 <- matrix(NaN, nrow=length(age), ncol=nrow(x))
    Lt_minus <- matrix(NaN, nrow=length(age), ncol=nrow(x))
    Lt_plus <- matrix(NaN, nrow=length(age), ncol=nrow(x))
    for(i in seq(ncol(Lt0))){
        par0 <- par_minus <- par_plus <- as.list(x[i,])
        par_minus$t_anchor <- par0$t_anchor-1
        par_plus$t_anchor <- par0$t_anchor+1
        Lt0[,i] <- VBGF(param = par0, t = age)
        Lt_minus[,i] <- VBGF(param = par_minus, t = age)
        Lt_plus[,i] <- VBGF(param = par_plus, t = age)
    }

    ## replace negative lengths with NA
    Lt0 <- replace(Lt0, Lt0<0, NaN)
    Lt_minus <- replace(Lt_minus, Lt_minus<0, NaN)
    Lt_plus <- replace(Lt_plus, Lt_plus<0, NaN)

    ## determine if a positive or negative shift improves overall covariance for each permutation
    cov0 <- cov(Lt0, use = "pair")
    shift <- 0 * seq(ncol(Lt0))
    for(i in seq(ncol(Lt0))){
        cov_minus <- cov(x = Lt_minus[,i], y = Lt0, use = "pair")
        cov_plus <- cov(x = Lt_plus[,i], y = Lt0, use = "pair")
        shift[i] <- c(0,-1,1)[which.max(c(sum(cov0[i,]), sum(cov_minus), sum(cov_plus)))]
    }

    ## fixed predictions
    agenew <- seq(0+min(shift), agemax+max(shift), 0.01)
    Lt <- matrix(NaN, nrow=length(agenew), ncol=nrow(x))
    for(i in seq(ncol(Lt))){
        par0 <- as.list(x[i,])
        par0$t_anchor <- par0$t_anchor + shift[i]
        Lt[,i] <- VBGF(param = par0, t = agenew)
    }

    for(i in seq(ncol(Lt))){
        if(i == 1) {
            plot(agenew, Lt[,i], t="n",
                 xlim = c(min(x$t_anchor+shift), max(agenew))+c(-0.1,0),
                 ylim = c(0,max(Lt, na.rm = TRUE)*1.05),
                 xlab = xlab, ylab = ylab,
                 xaxs = "i", yaxs = "i"
                 )
        }
        lines(agenew, Lt[,i], col = perm.col, lwd = perm.lwd)
    }
    box()

    ## multivariate kernel density estimate for max. density
    fhat <- bootKde(x[,c("Linf","K","t_anchor")], what = "plotted")

    ## maximum density
    max_dens <- fhat$eval.points[which.max(fhat$estimate),]

    ## predict density estimate of original data
    x$estimate <- fhat$estimate

    ## determine which resamples are in the CI
    limCI <- vector(mode = "list", length(CI))
    inCI <- vector(mode = "list", length(CI))
    names(inCI) <- names(limCI) <- paste0("CI", CI)
    for(j in seq(limCI)){
        inCI[[j]] <- x$estimate > quantile(x$estimate, probs = (100-CI[j])/100 )
        limCI[[j]] <- data.frame(t = agenew, min = NaN, max = NaN)
        for(i in seq(agenew)){
            limCI[[j]]$min[i] <- min(Lt[i, which(inCI[[j]]) ], na.rm = TRUE)
            limCI[[j]]$max[i] <- max(Lt[i, which(inCI[[j]]) ], na.rm = TRUE)
        }
        lines(max ~ t, limCI[[j]], col = ci.col[j], lwd = ci.lwd[j], lty = ci.lty[j])
        lines(min ~ t, limCI[[j]], col = ci.col[j], lwd = ci.lwd[j], lty = ci.lty[j])
    }
    inCI <- as.data.frame(inCI)

    lines(
        agenew, VBGF(param = as.list(max_dens), t = agenew),
        col = maxd.col, lwd = maxd.lwd, lty = maxd.lty
    )

    if(add_legend){
        legend("bottomright", legend = c(paste0("CI = ", CI, "%"), "Max. Dens."),
               bty = "n", col = c(ci.col, maxd.col), lty = c(ci.lty, maxd.lty),
               lwd = c(ci.lwd, maxd.lwd)
               )
    }

    if(add_max_dens_legend){
        legend("topleft",
               legend = c(paste0(names(max_dens), " = ", sprintf("%.2f", round(unlist(max_dens),2) ))),
               bty = "n", title = "Max. Dens. \nparameters:",
               inset = c(0,0.1)
               )
    }

    RES <- list(limCI = limCI, inCI = inCI, density = x$estimate, max_dens = as.list(max_dens))
    return(RES)
}
