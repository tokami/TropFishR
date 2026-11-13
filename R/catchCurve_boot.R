## Bootstrapped catch curve (internal, called by catchCurve with boot)


## Automatic selection of the points for the regression line of a catch
## curve (used for bootstrapped catch curves): from the point after the
## maximum of y to the last point with y > 0.5, optionally limited to the
## x range `agerange`. With robustReg = TRUE, the last point is removed if the
## slope of the last three points (with y > 0.2) is much flatter than the
## slope of the whole descending limb (e.g. due to selectivity of large fish).
## Returns the indices of the first and last point in x and y (or NA).
ccAutoRegInt <- function(x, y, agerange = NULL, robustReg = FALSE){

    fin <- which(is.finite(x) & is.finite(y))
    if(length(fin) < 3) return(c(NA, NA))
    xv <- x[fin]
    yv <- y[fin]

    if(robustReg){
        keep <- which(yv > 0.2)
        xk <- xv[keep]
        yk <- yv[keep]
        last <- which.max(xk)
        first <- which.max(yk) + 1
        if(last >= 3 && last - first >= 2){
            ind3 <- (last-2):last
            sx <- abs(coef(lm(yk[ind3] ~ xk[ind3]))[2])
            indAll <- first:last
            sxx <- abs(coef(lm(yk[indAll] ~ xk[indAll]))[2])
            if(is.finite(sx) && is.finite(sxx) && sx/sxx < 0.8){
                drop <- keep[last]
                xv <- xv[-drop]
                yv <- yv[-drop]
                fin <- fin[-drop]
            }
        }
    }

    maxY <- which.max(yv)
    after <- which(xv > xv[maxY] & yv > 0.5)
    if(length(after) < 1){
        last <- which.max(xv)
    }else{
        last <- after[which.max(xv[after])]
    }
    cutter <- c(maxY + 1, last)

    ## limit regression line to a given interval of x
    if(!is.null(agerange)){
        agecut <- c(which.min(abs(xv - agerange[1])),
                    which.min(abs(xv - agerange[2])))
        cutter <- c(max(agecut[1], cutter[1]), min(agecut[2], cutter[2]))
    }

    if(any(is.na(cutter)) || cutter[1] > length(fin) || cutter[2] - cutter[1] < 2)
        return(c(NA, NA))
    return(fin[cutter])
}


catchCurve_boot <- function(param, boot, cumulative = FALSE, calc_ogive = FALSE,
                            reg_int = NULL, agerange = NULL, natMort = NULL,
                            robustReg = FALSE, yearSel = NA, yearCombine = FALSE,
                            binSize = NA, CI = 95){

    checkBoot(boot, c("Linf", "K", if(!is.null(natMort)) natMort))
    if(calc_ogive && cumulative)
        stop("It is not possible to estimate the selection ogive for the cumulative catch curve.")
    if(!("midLengths" %in% names(param)))
        stop("The bootstrapped catch curve requires length-frequency data (midLengths).")

    n <- nrow(boot$bootRaw)
    Z <- L50 <- L75 <- intercepts <- rep(NA_real_, n)
    datcc <- vector("list", n)
    skip <- rep(NA_character_, n)
    noOgive <- 0
    knifeEdge <- rep(FALSE, n)
    xname <- ifelse(cumulative, "ln_Linf_L", "t_midL")
    yname <- ifelse(cumulative, "lnC", "lnC_dt")

    for(x in seq_len(n)){

        lfqx <- bootLfq(param, boot, x, yearSel = yearSel,
                        yearCombine = yearCombine, binSize = binSize)
        ## growth parameters of the resample (top level: catchCurve uses
        ## the seasonal VBGF if C and ts are elements of the lfq list)
        lfqx <- c(lfqx, bootGrowthPar(boot, x))

        ## x and y of the catch curve
        cc0 <- try(suppressWarnings(
            catchCurve(lfqx, cumulative = cumulative, auto = TRUE, plot = FALSE)),
            silent = TRUE)
        if(inherits(cc0, "try-error")){
            skip[x] <- "catch curve failed"
            next
        }

        ## points for the regression line
        if(is.null(reg_int)){
            ri <- ccAutoRegInt(cc0[[xname]], cc0[[yname]], agerange = agerange,
                               robustReg = robustReg)
        }else{
            ri <- reg_int
        }
        if(any(is.na(ri))){
            skip[x] <- "no descending limb found"
            next
        }

        cc <- try(suppressWarnings(
            catchCurve(lfqx, cumulative = cumulative, calc_ogive = calc_ogive,
                       reg_int = ri, plot = FALSE)), silent = TRUE)
        if(inherits(cc, "try-error") && calc_ogive){
            ## selection ogive could not be estimated
            cc <- try(suppressWarnings(
                catchCurve(lfqx, cumulative = cumulative, calc_ogive = FALSE,
                           reg_int = ri, plot = FALSE)), silent = TRUE)
        }
        if(inherits(cc, "try-error")){
            skip[x] <- "regression failed"
            next
        }

        slope <- coef(cc$linear_mod)[2]
        if(!is.finite(slope) || slope > 0){
            skip[x] <- "positive slope"
            next
        }

        Z[x] <- cc$Z
        intercepts[x] <- coef(cc$linear_mod)[1]
        sel <- ri[1]:ri[2]
        datcc[[x]] <- data.frame(x = cc[[xname]][sel], y = cc[[yname]][sel])
        datcc[[x]] <- datcc[[x]][is.finite(datcc[[x]]$x) & is.finite(datcc[[x]]$y),]
        if(calc_ogive && (is.null(cc$L50) || !is.finite(cc$L50))){
            ## fewer than two length classes with catch below the regression
            ## range (e.g. the descending limb starts at the first length
            ## class): knife-edge selection at the lower limit of the first
            ## length class of the regression
            before <- seq_len(ri[1] - 1)
            if(sum(lfqx$catch[before] > 0, na.rm = TRUE) < 2){
                interval <- lfqx$midLengths[2] - lfqx$midLengths[1]
                L50[x] <- L75[x] <- lfqx$midLengths[ri[1]] - interval / 2
                knifeEdge[x] <- TRUE
            }else{
                noOgive <- noOgive + 1
            }
        }else if(calc_ogive){
            L50[x] <- cc$L50
            L75[x] <- cc$L75
        }
    }

    bootSkipWarning(skip, "catchCurve")
    if(calc_ogive && any(knifeEdge))
        message(sprintf("catchCurve: knife-edge selection assumed for %d of %d resamples (fewer than two length classes with catch below the regression line; see misc$knifeEdge).",
                        sum(knifeEdge), n))
    if(calc_ogive && noOgive > 0)
        warning(sprintf("catchCurve: the selection ogive could not be estimated for %d of %d resamples (L50 and L75 are NA).",
                        noOgive, n), call. = FALSE)

    new <- list(Z = Z)
    if(!is.null(natMort)){
        new$FM <- Z - boot$bootRaw[[natMort]]
        nneg <- sum(is.finite(new$FM) & new$FM <= 0)
        if(nneg > 0)
            warning(sprintf("catchCurve: %d of %d resamples have FM <= 0 (natural mortality '%s' >= Z). These resamples are skipped by VPA and predict_mod.",
                            nneg, n, natMort), call. = FALSE)
    }
    if(calc_ogive){
        new$L50 <- L50
        new$L75 <- L75
    }

    regx <- unlist(lapply(datcc, function(d) if(is.null(d)) NULL else range(d$x)))
    misc <- list(binSize = binSize, yearSel = yearSel, yearCombine = yearCombine,
                 cumulative = cumulative, intsCC = intercepts,
                 knifeEdge = if(calc_ogive) knifeEdge,
                 regint = if(length(regx) > 0) range(regx) else c(NA, NA),
                 datcc = datcc)

    bootAddPars(boot, new, CI = CI, misc = misc)
}
