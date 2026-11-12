## Internal helpers for ELEFAN_GA_boot and ELEFAN_SA_boot

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
