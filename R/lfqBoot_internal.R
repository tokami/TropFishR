## Internal helpers for the bootstrapped assessment: ELEFAN_GA_boot,
## ELEFAN_SA_boot, and the boot argument of M_empirical, catchCurve, VPA and
## predict_mod

## Run `fitOne(x)` for x = 1, ..., nresamp, either serially or on a cluster.
## `fitOne` is a closure defined in the calling function; on the workers it
## resolves TropFishR functions from the TropFishR namespace, so the workers
## have to have TropFishR installed.
bootFits <- function(fitOne, nresamp, parallel = FALSE, no_cores = 1,
                     clusterType = "PSOCK", outfile = NULL){

    if(!is.null(outfile)) unlink(outfile) # delete old outfile

    run <- function(x){
        res <- fitOne(x)
        ## progress (for checking progress in outfile)
        if(!is.null(outfile)){
            cat(paste("resamp", x, "completed @", Sys.time()), "\n",
                file = outfile, append = TRUE)
        }
        return(res)
    }

    if(parallel){
        cl <- parallel::makeCluster(no_cores, type = clusterType)
        on.exit(parallel::stopCluster(cl))
        res <- parallel::parLapply(cl, seq_len(nresamp), run)
    }else{
        res <- lapply(seq_len(nresamp), run)
    }
    return(res)
}


## Summarise bootstrapped VBGF estimates.
## `res` is a list with one vector per resample: c(par, agemax, seed, Rn_max)
bootSummary <- function(res, CI = 95, agemax = NULL){

    tmp0 <- as.data.frame(do.call("rbind", res))
    nc <- ncol(tmp0)
    bootRaw <- tmp0[,-c(nc-1, nc)]
    rownames(bootRaw) <- NULL
    seeds <- as.numeric(tmp0[,nc-1])
    scores <- as.numeric(tmp0[,nc])

    ## multivariate kernel smoothing (without derived parameters)
    x <- bootRaw[,!(names(bootRaw) %in% c("phiL", "agemax"))]
    H <- ks::Hpi(x, nstage = 1)
    fhat <- ks::kde(x = x, H = H, eval.points = x)

    ## maximum density
    maxDens <- unlist(fhat$eval.points[which.max(fhat$estimate),])
    medians <- apply(fhat$eval.points, 2, median, na.rm = TRUE)

    ## confidence intervals (univariate)
    tmp <- (100-CI)/2/100
    limCIuni <- apply(fhat$eval.points, 2, quantile, na.rm = TRUE,
                      probs = c(tmp, 1-tmp))
    rownames(limCIuni) <- c("lo","up")

    ## confidence intervals (multivariate)
    inCI <- fhat$estimate > quantile(fhat$estimate, probs = (100-CI)/100)
    limCI <- apply(x[inCI,], 2, range)
    rownames(limCI) <- c("lo","up")

    ## derived parameters: growth performance index and maximum age
    phiL <- function(p) log10(p[["K"]]) + 2 * log10(p[["Linf"]])
    tmax <- function(p){
        if(is.null(agemax)){
            ceiling((1/-p[["K"]]) * log(1-((p[["Linf"]]*0.95)/p[["Linf"]])))
        }else agemax
    }
    maxDens <- c(maxDens, phiL = phiL(maxDens), agemax = tmax(maxDens))
    medians <- c(medians, phiL = phiL(medians), agemax = tmax(medians))
    ## limits of the derived parameters from the limits of Linf and K;
    ## sorted, because agemax decreases with K
    addCols <- function(m){
        lo <- as.list(m["lo",]); up <- as.list(m["up",])
        cbind(m, phiL = sort(c(phiL(lo), phiL(up))),
              agemax = sort(c(tmax(lo), tmax(up))))
    }
    limCIuni <- addCols(limCIuni)
    limCI <- addCols(limCI)

    ret <- list()
    ret$bootRaw <- bootRaw
    ret$seed <- seeds
    ret$Rn_max <- round(scores, 3)
    ret$maxDen <- maxDens
    ret$median <- medians
    ret$CI <- limCIuni
    ret$multiCI <- limCI
    ret$misc <- list()
    class(ret) <- "lfqBoot"
    return(ret)
}


## Check the boot argument
checkBoot <- function(boot, required = c("Linf", "K")){
    if(!inherits(boot, "lfqBoot"))
        stop("You provided an object for boot, but it does not have class 'lfqBoot'. Please check.")
    missing <- setdiff(required, names(boot$bootRaw))
    if(length(missing) > 0)
        stop(paste0("The boot object requires the column(s) ",
                    paste(missing, collapse = ", "),
                    " (in boot$bootRaw). ",
                    "Please run the corresponding step of the assessment first."))
    invisible(TRUE)
}


## Re-create the resampled lfq data set of resample x of a boot object (same
## seed and resampling as in ELEFAN_GA_boot/ELEFAN_SA_boot) and prepare it for
## the analysis of one year: optionally subset years, correct the catch,
## re-bin, and sum the catches over the samples.
bootLfq <- function(lfq, boot, x, yearSel = NA, yearCombine = FALSE,
                    binSize = NA, catch_corFac = NA, raise = NULL){

    set.seed(boot$seed[x])
    lfqb <- lfqResample(lfq)

    ## raise the resampled numbers to the total catch (the weight of the
    ## samples is re-estimated for each resample unless measured)
    if(!is.null(raise)){
        lfqb <- lfqRaise(lfqb, landings = raise$landings,
                         sampleWeight = raise$sampleWeight, by = raise$by,
                         a = raise$a, b = raise$b, LW_unit = raise$LW_unit)
    }

    ## subset years (after resampling, so that the resamples match)
    if(!is.na(yearSel[1])){
        sel <- format(lfqb$dates, "%Y") %in% as.character(yearSel)
        if(!any(sel)) stop(paste0("The selected year(s) ",
                                  paste(yearSel, collapse = ", "),
                                  " are not in the lfq data!"))
        lfqb$catch <- lfqb$catch[, sel, drop = FALSE]
        lfqb$dates <- lfqb$dates[sel]
    }

    ## correct catch with raising factor(s) (one per sample)
    if(!is.na(catch_corFac[1])){
        if(length(catch_corFac) != 1 && length(catch_corFac) != ncol(lfqb$catch))
            stop("catch_corFac has to be a single value or one value per sample.")
        lfqb$catch <- sweep(lfqb$catch, 2,
                            rep_len(catch_corFac, ncol(lfqb$catch)), "*")
    }

    ## sum catches per year (and re-bin)
    if(!is.na(binSize)){
        lfqb <- lfqModify(lfqb, vectorise_catch = TRUE, bin_size = binSize)
    }else{
        lfqb <- lfqModify(lfqb, vectorise_catch = TRUE)
    }
    if(is.matrix(lfqb$catch)){
        if(ncol(lfqb$catch) == 1){
            lfqb$catch <- as.numeric(lfqb$catch)
        }else if(yearCombine){
            lfqb$catch <- rowSums(lfqb$catch)
        }else{
            stop("The lfq data spans several years, please subset for one year (yearSel) at a time or use yearCombine = TRUE!")
        }
    }
    lfqb$par <- NULL
    return(lfqb)
}


## Growth parameters of resample x as list (t0 = 0, as in ELEFAN)
bootGrowthPar <- function(boot, x){
    br <- boot$bootRaw
    par <- list(Linf = br$Linf[x], K = br$K[x], t0 = 0)
    if(all(c("C", "ts") %in% names(br))){
        par$C <- br$C[x]
        par$ts <- br$ts[x]
    }
    return(par)
}


## Univariate summary of bootstrapped estimates: median, confidence interval
## (quantiles) and value with maximum density (mode of the kernel density)
bootUniSummary <- function(x, CI = 95){
    x <- as.numeric(x)
    xx <- x[is.finite(x)]
    citmp <- (100-CI)/2/100
    if(length(xx) == 0){
        return(list(maxDen = NA_real_, median = NA_real_, CI = c(NA_real_, NA_real_)))
    }
    med <- median(xx)
    ci <- unname(quantile(xx, probs = c(citmp, 1-citmp)))
    if(length(unique(xx)) == 1){
        maxDen <- xx[1]
    }else{
        k <- try(ks::kde(xx), silent = TRUE)
        if(inherits(k, "try-error")){
            maxDen <- NA_real_
        }else{
            maxDen <- k$eval.points[which.max(k$estimate)]
        }
    }
    list(maxDen = maxDen, median = med, CI = ci)
}


## Add bootstrapped estimates (data.frame `new`, one row per resample) to a
## boot object: columns of bootRaw and univariate summaries in maxDen, median
## and CI. Existing columns with the same names are replaced. Elements of
## `misc` are added to (or replace those in) boot$misc.
bootAddPars <- function(boot, new, CI = 95, misc = list()){
    new <- as.data.frame(new)
    if(nrow(new) != nrow(boot$bootRaw)) stop("Number of rows does not match the boot object.")
    ret <- boot
    if(is.null(ret$CI)){
        ret$CI <- matrix(nrow = 2, ncol = 0, dimnames = list(c("lo","up"), NULL))
    }
    for(nm in names(new)){
        s <- bootUniSummary(new[[nm]], CI = CI)
        ret$bootRaw[[nm]] <- new[[nm]]
        ret$maxDen[nm] <- s$maxDen
        ret$median[nm] <- s$median
        if(nm %in% colnames(ret$CI)){
            ret$CI[, nm] <- s$CI
        }else{
            ret$CI <- cbind(ret$CI, matrix(s$CI, ncol = 1, dimnames = list(c("lo","up"), nm)))
        }
    }
    if(is.null(ret$misc)) ret$misc <- list()
    for(nm in names(misc)) ret$misc[[nm]] <- misc[[nm]]
    class(ret) <- "lfqBoot"
    return(ret)
}


## Warn about resamples that were skipped. `reason` is a character vector with
## one element per resample: NA if the resample was analysed, otherwise the
## reason why it was skipped.
bootSkipWarning <- function(reason, fun){
    r <- reason[!is.na(reason)]
    if(length(r) > 0){
        tab <- table(r)
        warning(sprintf("%s: %d of %d resamples were skipped (estimates are NA): %s.",
                        fun, length(r), length(reason),
                        paste0(names(tab), " (", tab, ")", collapse = ", ")),
                call. = FALSE)
    }
    invisible(NULL)
}
