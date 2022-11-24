#' Resampling of length-frequency data
#'
#' @param lfq a length frequency object of the class `lfq`.
#' @param boot optional; an object of class 'lfqBoot' for reproducing
#'    a lfq data by re-using seed values
#'
#' @description The function resamples the `lfq` data by sample date.
#'   Sampling is done in a non-parametric manner that follows the relative
#'   frequencies of the original data, allowing for individual counts to be
#'   selected more than once (i.e. `replace = TRUE` in \link[base]{sample}),
#'   and resulting in total counts (by sample) equal to the original data.
#'
#' @return a resampled version of the `lfq` class dataset.
#'
#' @examples
#' # load data
#' data(alba)
#'
#' # permutate data
#' alba_p <- lfqResample(alba)
#'
#' # side-by-side plot
#' op <- par( mfcol = c(2,1), mar=c(4,4,2,1) )
#' plot(lfqRestructure(alba)); mtext("original", side=3, line=0.25)
#' plot(lfqRestructure(alba_p)); mtext("resampled", side=3, line=0.25)
#' box()
#' par(op)
#'
#' # relative difference
#' alba_diff <- alba
#' alba_diff$delta <- (alba$catch - alba_p$catch) / alba$catch
#' alba_diff$delta[is.na(alba_diff$delta)] <- 0
#' with(alba_diff, image(x=dates, y=midLengths, z=t(delta), col=rev(cm.colors(21))) )
#' with(alba_diff, contour(x=dates, y=midLengths, z=t(delta), add=TRUE))
#'
#' @export
#'
lfqResample <- function(lfq, sampSize = NA, boot = NULL){

    ## reproduce used lfq data sets if boot object provided
    nI <- ifelse(!is.null(boot) & "seed" %in% names(boot), length(boot$seed),1)

    lfqList <- vector("list", nI)
    for(I in 1:nI){

        lfqb <- lfq

        if(!is.null(boot)) set.seed(boot$seed[I])

        for(i in seq(length(lfq$dates))){
            if(is.na(sampSize) | is.nan(sampSize)){
                sampSizi <- sum(lfq$catch[,i])
            }else{
                sampSizi <- sampSize
            }
            inds <- sample(x = lfq$midLengths, size = round(sampSizi),
                           prob = lfq$catch[,i], replace = TRUE)
            bin.width <- diff(lfq$midLengths) # bin width (should allow for uneven bin sizes)
            bin.lower <- lfq$midLengths - (c(bin.width[1], bin.width)/2) # upper bin limit
            bin.upper <- lfq$midLengths + (c(bin.width, bin.width[length(bin.width)])/2) # lower bin limit
            breaks <- unique(c(bin.lower, bin.upper))
            h <- hist(inds, breaks = breaks, plot = FALSE)
            lfqb$catch[,i] <- h$counts
        }
        lfqList[[I]] <- lfqb
    }

    if(is.null(boot)) lfqList <- lfqList[[1]]

    return(lfqList)
}





#' Bootstraped ELEFAN_SA
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
#' @param parallel logical; should parallelized computing be used. This differs from the
#'    `parallel` argument in `ELEFAN_GA` in that it is not used within the `ga` function for
#'    calculation at the population level, but rather for permutations. Depending on platform
#'    operating system, the argument `clusterType` can be adjusted (see argument description for
#'    details). (Default: `parallel = TRUE`)
#' @param nresamp numeric; the number of permutations to run (Default: `nresamp = 200`)
#' @param no_cores numeric (Default: `no_cores = detectCores() - 1`)
#' @param clusterType (Default: `clusterType = "PSOCK"`)
#' @param outfile character; text file name (Default: `outfile = "output.txt"`) which will
#'    records the progress of the permutation completions.
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
#' @param seed seed value for random number reproducibility (Default: NULL)
#' @param CI percentage for confidence intervals (Default: 95)
#'
#' @description `ELEFAN_SA_boot` performs a bootstrapped fitting of
#'   von Bertalanffy growth function (VBGF) via the \link{ELEFAN_SA} function.
#'   Most of the arguments are simply passed to the function within many
#'   permutations (resampling) of the original lfq data.
#'
#' @return a data.frame of fitted VBGF parameters (columns) by permutation (rows).
#'
#' @importFrom ks Hpi
#' @importFrom ks kde
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel parLapply
#'
#' @examples
#' # load data
#' data(alba)
#'
#' # settings
#' MA <- 7
#' init_par <- list(Linf = 10, K = 1, t_anchor = 0.5)
#' low_par <- list(Linf = 8, K = 0.1, t_anchor = 0)
#' up_par <- list(Linf = 15, K = 5, t_anchor = 1)
#' SA_time <- 5
#' SA_temp <- 1e5
#' nresamp <- 10
#'
#'
#' # parallel version
#' library(parallel)
#' t1 <- Sys.time()
#' res <- ELEFAN_SA_boot(lfq=alba, MA = MA, seasonalised = FALSE,
#'   init_par = init_par, up_par = up_par, low_par = low_par,
#'   SA_time = SA_time, SA_temp = SA_temp, seed = 1,
#'   nresamp = nresamp, parallel = TRUE, no_cores = detectCores()-1
#' )
#' t2 <- Sys.time()
#' t2 - t1
#' res
#'
#'
#' # non-parallel version
#' t1 <- Sys.time()
#' res <- ELEFAN_SA_boot(lfq=alba, seasonalised = FALSE,
#'   init_par = init_par, up_par = up_par, low_par = low_par,
#'   SA_time = SA_time, SA_temp = SA_temp, seed = 1,
#'   nresamp = nresamp, MA=MA, parallel = FALSE
#' )
#' t2 <- Sys.time()
#' t2 - t1
#' res
#'
#' # plot resulting distributions
#' op <- par(
#'   mfcol = c(1, ncol(res$bootRaw)-1),
#'   mar=c(1,2,2,1),
#'   mgp = c(2,0.5,0), tcl = -0.25, cex=1
#' )
#' for(i in seq(ncol(res$bootRaw)-1)){
#'   boxplot(res$bootRaw[[i]], boxwex = 0.25, col = 8)
#'   text(1, median(res$bootRaw[[i]]),
#'   labels = round(median(res$bootRaw[[i]]),2),
#'     pos = 4, col = 4)
#'   mtext(names(res$bootRaw)[i], line=0.25, side=3)
#' }
#' par(op)
#'
#' @export
#'
ELEFAN_SA_boot <- function(
                           lfq, seasonalised = FALSE,
                           init_par = list(Linf = 50, K = 0.5, t_anchor = 0.5, C = 0, ts = 0),
                           low_par = NULL, up_par = NULL,
                           parallel = TRUE, nresamp = 200, no_cores = detectCores() - 1,
                           clusterType = "PSOCK",
                           outfile = "output.txt",
                           SA_time = 60, SA_temp = 1e5, maxit = NULL,
                           MA = 5, addl.sqrt = FALSE, agemax = NULL,
                           flagging.out = TRUE, seed = NULL, CI = 95
                           ){

    if(!is.null(outfile)){unlink(outfile)} # delete old outfile

    if(is.null(seed)){ seed <- round(runif(1, min = 0, max = 1e6)) }

    if(parallel){ # Parallel version
        ARGS <- list(
            "lfqResample",
            "lfq", "seasonalised",
            "init_par", "low_par", "up_par",
            "parallel", "nresamp", "no_cores",
            "SA_temp", "SA_time", "maxit", "seed",
            "MA", "addl.sqrt", "agemax", "flagging.out",
            "outfile"
        )

        parFun <- function(x){
                                        # load packages and functions for each cluster
            library(TropFishR)

            set.seed(seed + x)

                                        # resample data
            lfqb <- lfqResample(lfq)

                                        # call ELEFAN_SA
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

                                        # print output (for checking progress in output.txt)
            if(!is.null(outfile)){
                sink(file=outfile, append = TRUE)
                print(paste("resamp", x, "completed @", Sys.time()))
                sink()
            }

                                        # return result
            return( c(unlist(fitboot$par), seed + x, fitboot$Rn_max) )
        }

        cl <- parallel::makeCluster(no_cores, type=clusterType)
        nn <- split(1:nresamp, 1:nresamp)
        parallel::clusterExport(cl, varlist = ARGS, envir=environment())
        res <- parallel::parLapply(cl, nn, parFun)
        parallel::stopCluster(cl)
    }

    if(!parallel){ # Non-parallel version
        res <- vector("list", nresamp) # empty results list
        for(x in seq(res)){
            set.seed(seed + x)

            ## resample data
            lfqb <- lfqResample(lfq)

                                        # call ELEFAN_SA
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

            if(!is.null(outfile)){
                sink(file=outfile, append = TRUE)
                print(paste("resamp", x, "completed @", Sys.time()))
                sink()
            }

            ## return result
            res[[x]] <- c(unlist(fitboot$par), seed + x, fitboot$Rn_max)
        }
    }

    tmp0 <- as.data.frame(do.call("rbind", res))
    tmp <- tmp0[,-c(ncol(tmp0)-1,ncol(tmp0))]

    seeds <- as.numeric(tmp0[,(ncol(tmp0)-1)])
    scores <- as.numeric(tmp0[,ncol(tmp0)])

    ## lfqboot object
    bootRaw <- tmp

    ## Conduct multivariate kernel smoothing
    phiLcol <- which(names(bootRaw) == "phiL")
    if(length(phiLcol > 0)){
        x <- bootRaw[,-phiLcol]
    } else {
        x <- bootRaw
    }

    H <- ks::Hpi(x, nstage = 1)
    fhat <- ks::kde(x = x, H = H, eval.points = x)

                                        # maximum density
    maxDens <- fhat$eval.points[which.max(fhat$estimate),]
    nami <- names(maxDens)
    maxDens <- as.numeric(maxDens)
    names(maxDens) <- nami

    medians <- apply(fhat$eval.points, 2, median, na.rm = TRUE)

    ## confidence intervals (univariate)
    tmp <- (100-CI)/2/100
    limCIuni <- apply(fhat$eval.points, 2, quantile, na.rm = TRUE, probs = c(tmp, 1-tmp))
    rownames(limCIuni) <- c("lo","up")

    ## confidence intervals (multivariate)
    inCI <- fhat$estimate > quantile(fhat$estimate, probs = (100-CI)/100)
    x_inCI <- x[inCI,]
    limCI <- apply(x_inCI, 2, range)
    rownames(limCI) <- c("lo","up")

    ## adding phiL to maxDen and CI assuming relationship
    maxDens[length(maxDens)+1] <- log10(maxDens[which(names(maxDens) == "K")]) +
        2 * log10(maxDens[which(names(maxDens) == "Linf")])
    names(maxDens) <- c(names(maxDens)[-length(maxDens)],"phiL")
    medians[length(medians)+1] <- log10(medians[which(names(medians) == "K")]) +
        2 * log10(medians[which(names(medians) == "Linf")])
    names(medians) <- c(names(medians)[-length(medians)],"phiL")
    limCIuni <- cbind(limCIuni,log10(limCIuni[,which(names(maxDens) == "K")]) +
                               2 * log10(limCIuni[,which(names(maxDens) == "Linf")]))
    colnames(limCIuni) <- c(colnames(limCIuni)[-ncol(limCIuni)],"phiL")
    limCI <- cbind(limCI,log10(limCI[,which(names(maxDens) == "K")]) +
                         2 * log10(limCI[,which(names(maxDens) == "Linf")]))
    colnames(limCI) <- c(colnames(limCI)[-ncol(limCI)],"phiL")

    ## save and return results
    ret <- list()
    ret$bootRaw <- bootRaw
    ret$seed <- seeds
    ret$Rn_max <- round(scores,3)
    ret$maxDen <- maxDens
    ret$median <- medians
    ret$CI <- limCIuni
    ret$multiCI <- limCI
    ret$misc <- list()
    class(ret) <- "lfqBoot"

    return(ret)
}





#' Bootstrapped ELEFAN_GA
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
#' @param parallel logical; should parallelized computing be used. This differs from the
#'    `parallel` argument in `ELEFAN_GA` in that it is not used within the `ga` function for
#'    calculation at the population level, but rather for permutations. Depending on platform
#'    operating system, the argument `clusterType` can be adjusted (see argument description for
#'    details). (Default: `parallel = TRUE`)
#' @param nresamp numeric; the number of resamples (Default: `nresamp = 200`)
#' @param no_cores numeric (Default: `no_cores = detectCores() - 1`)
#' @param clusterType (Default: `clusterType = "PSOCK"`)
#' @param outfile character; text file name (Default: `outfile = "output.txt"`) which will
#'    records the progress of the permutation completions.
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
#' @param seed seed value for random number reproducibility (Default: NULL)
#' @param CI percentage for confidence intervals (Default: 95)
#'
#' @description `ELEFAN_GA_boot` performs a bootstrapped fitting of
#'   von Bertalanffy growth function (VBGF) via the \link{ELEFAN_GA} function. Most of the arguments
#'   are simply passed to the function within many permutations (resampling) of the original
#'   lfq data.
#'
#' @return an object of class 'lfqBoot' with following objects
#'    \itemize{
#'    \item \strong{bootRaw} a data.frame of fitted VBGF parameters (columns) by
#'               permutation (rows),
#'    \item \strong{bootModes} a dataframe with the modes of the VBGF parameters,
#'    \item \strong{bootCIs} a dataframe with the lower and upper 95% confidence
#'               interval for the VBGF growth parameters.
#'    }
#'
#' @importFrom ks kde
#' @importFrom ks Hpi
#' @importFrom parallel makeCluster
#' @importFrom parallel stopCluster
#' @importFrom parallel parLapply
#'
#' @examples
#' # load data
#' data(alba)
#'
#' # settings
#' MA <- 7
#' low_par <- list(Linf = 8, K = 0.1, t_anchor = 0, C = 0, ts = 0)
#' up_par <- list(Linf = 15, K = 5, t_anchor = 1, C = 1, ts = 1)
#' popSize <- 60
#' maxiter <- 50
#' run <- 10
#' pmutation <- 0.2
#' nresamp <- 20
#'
#'
#' # parallel version
#' library(parallel)
#' t1 <- Sys.time()
#' res <- ELEFAN_GA_boot(lfq=alba, MA = MA, seasonalised = FALSE,
#'   up_par = up_par, low_par = low_par,
#'   popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
#'   nresamp = nresamp, parallel = TRUE, no_cores = detectCores()-1,
#'   seed = 1
#' )
#' t2 <- Sys.time()
#' t2 - t1
#' res
#'
#'
#' # non-parallel version
#' t1 <- Sys.time()
#' res <- ELEFAN_GA_boot(lfq=alba, seasonalised = FALSE, up_par = up_par, low_par = low_par,
#'   popSize = popSize, maxiter = maxiter, run = run,
#'   pmutation = pmutation, nresamp = nresamp, MA=MA, parallel = FALSE,
#'   seed = 1
#' )
#' t2 <- Sys.time()
#' t2 - t1
#' res
#'
#' # plot resulting distributions
#' op <- par(
#'   mfcol = c(1, ncol(res$bootRaw)-1),
#'   mar=c(1,2,2,1),
#'   mgp = c(2,0.5,0), tcl = -0.25, cex=1
#' )
#' for(i in seq(ncol(res$bootRaw)-1)){
#'   boxplot(res$bootRaw[[i]], boxwex = 0.25, col = 8)
#'   text(1, median(res$bootRaw[[i]]), labels = round(median(res$bootRaw[[i]]),2),
#'     pos = 4, col = 4)
#'   mtext(names(res$bootRaw)[i], line=0.25, side=3)
#' }
#' par(op)
#'
#' @export
#'
ELEFAN_GA_boot <- function(
                           lfq, seasonalised = FALSE, low_par = NULL, up_par = NULL,
                           parallel = TRUE, nresamp = 200, no_cores = detectCores() - 1,
                           clusterType = "PSOCK",
                           outfile = "output.txt",
                           popSize = 60, maxiter = 50, run = 10,
                           pmutation = 0.2, pcrossover = 0.8,
                           elitism = base::max(1, round(popSize * 0.05)),
                           MA = 5, addl.sqrt = FALSE, agemax = NULL,
                           flagging.out = TRUE, seed = NULL, CI = 95,
                           weight.by.sample.size = FALSE, bin_size = NA
                           ){

    if(!is.null(outfile)){unlink(outfile)} # delete old outfile

    if(is.null(seed)){ seed <- round(runif(1, min = 0, max = 1e6))}


    if(parallel){ # Parallel version
        ARGS <- list(
            "lfqResample",
            "lfq", "seasonalised", "low_par", "up_par",
            "parallel", "nresamp", "no_cores",
            "popSize", "maxiter", "run", "seed",
            "pmutation", "pcrossover", "elitism",
            "MA", "addl.sqrt", "agemax", "flagging.out",
            "outfile","bin_size"
        )


        parFun <- function(x){
                                        # load packages and functions for each cluster
            library(TropFishR)

            set.seed(seed + x)

                                        # permutate data
            lfqb <- lfqResample(lfq)

            if(!is.na(bin_size) && is.numeric(bin_size)){
                lfqb <- lfqModify(lfqb, bin_size = bin_size)
            }

                                        # call ELEFAN_GA
            fitboot <- ELEFAN_GA(
                lfqb, seasonalised = seasonalised,
                low_par = low_par,
                up_par = up_par,
                popSize = popSize, maxiter = maxiter, run = run,
                pmutation = pmutation, pcrossover = pcrossover, elitism = elitism,
                MA = MA, parallel = FALSE, addl.sqrt = addl.sqrt,
                agemax = agemax, flagging.out = flagging.out,
                plot.score = FALSE, seed = NULL,
                weight.by.sample.size = weight.by.sample.size
            )

                                        # print output (for checking progress in output.txt)
            if(!is.null(outfile)){
                sink(file=outfile, append = TRUE)
                print(paste("resamp", x, "completed @", Sys.time()))
                sink()
            }

                                        # return result
            return(c(unlist(fitboot$par), seed + x, fitboot$Rn_max))
        }

        cl <- parallel::makeCluster(no_cores, type=clusterType)
        nn <- split(seed + (0:(nresamp-1)), 1:nresamp)
        parallel::clusterExport(cl, varlist = ARGS, envir=environment())
        res <- parLapply(cl, nn, parFun)
        stopCluster(cl)
    }

    if(!parallel){ # Non-parallel version
        res <- vector("list", nresamp) # empty results list
        for(x in seq(res)){

            set.seed(seed + x)

            ## resample data
            lfqb <- lfqResample(lfq)

            if(!is.na(bin_size) && is.numeric(bin_size)){
                lfqb <- lfqModify(lfqb, bin_size = bin_size)
            }


                                        # call ELEFAN_GA
            fitboot <- ELEFAN_GA(
                lfqb, seasonalised = seasonalised,
                low_par = low_par,
                up_par = up_par,
                popSize = popSize, maxiter = maxiter, run = run,
                pmutation = pmutation, pcrossover = pcrossover, elitism = elitism,
                MA = MA, parallel = FALSE, addl.sqrt = addl.sqrt,
                agemax = agemax, flagging.out = flagging.out,
                plot.score = FALSE, seed = NULL,
                weight.by.sample.size = weight.by.sample.size
            )

            if(!is.null(outfile)){
                sink(file=outfile, append = TRUE)
                print(paste("resamp", x, "completed @", Sys.time()))
                sink()
            }

                                        # return result
            res[[x]] <- c(unlist(fitboot$par), seed + x, fitboot$Rn_max)
        }
    }
    tmp0 <- as.data.frame(do.call("rbind", res))
    tmp <- tmp0[,-c(ncol(tmp0)-1,ncol(tmp0))]

    seeds <- as.numeric(tmp0[,(ncol(tmp0)-1)])
    scores <- as.numeric(tmp0[,ncol(tmp0)])

    ## lfqboot object
    bootRaw <- tmp

    ## Conduct multivariate kernel smoothing
    phiLcol <- which(names(bootRaw) == "phiL")
    if(length(phiLcol) > 0){
        x <- bootRaw[,-phiLcol]
    } else {
        x <- bootRaw
    }

    H <- ks::Hpi(x, nstage = 1)
    fhat <- ks::kde(x = x, H = H, eval.points = x)

                                        # maximum density
    maxDens <- fhat$eval.points[which.max(fhat$estimate),]
    nami <- names(maxDens)
    maxDens <- as.numeric(maxDens)
    names(maxDens) <- nami

    medians <- apply(fhat$eval.points, 2, median, na.rm = TRUE)

    ## confidence intervals (univariate)
    tmp <- (100-CI)/2/100
    limCIuni <- apply(fhat$eval.points, 2, quantile, na.rm = TRUE, probs = c(tmp, 1-tmp))
    rownames(limCIuni) <- c("lo","up")

    ## confidence intervals (multivariate)
    inCI <- fhat$estimate > quantile(fhat$estimate, probs = (100-CI)/100)
    x_inCI <- x[inCI,]
    limCI <- apply(x_inCI, 2, range)
    rownames(limCI) <- c("lo","up")

    ## adding phiL to maxDen and CI assuming relationship
    maxDens[length(maxDens)+1] <- log10(maxDens[which(names(maxDens) == "K")]) +
        2 * log10(maxDens[which(names(maxDens) == "Linf")])
    names(maxDens) <- c(names(maxDens)[-length(maxDens)],"phiL")
    medians[length(medians)+1] <- log10(medians[which(names(medians) == "K")]) +
        2 * log10(medians[which(names(medians) == "Linf")])
    names(medians) <- c(names(medians)[-length(medians)],"phiL")
    limCIuni <- cbind(limCIuni,log10(limCIuni[,which(names(maxDens) == "K")]) +
                               2 * log10(limCIuni[,which(names(maxDens) == "Linf")]))
    colnames(limCIuni) <- c(colnames(limCIuni)[-ncol(limCIuni)],"phiL")
    limCI <- cbind(limCI,log10(limCI[,which(names(maxDens) == "K")]) +
                         2 * log10(limCI[,which(names(maxDens) == "Linf")]))
    colnames(limCI) <- c(colnames(limCI)[-ncol(limCI)],"phiL")

    ## save and return results
    ret <- list()
    ret$bootRaw <- bootRaw
    ret$seed <- seeds
    ret$Rn_max <- round(scores,3)
    ret$maxDen <- maxDens
    ret$median <- medians
    ret$CI <- limCIuni
    ret$multiCI <- limCI
    ret$misc <- list()
    class(ret) <- "lfqBoot"

    return(ret)
}




#' VBGF plot and CI
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
#' @importFrom ks Hpi
#' @importFrom ks kde
#'
#' @examples
#'
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
#' # using output in lfq plot (see ?TropFishR::plot.lfq)
#' library(TropFishR)
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

                                        # expand ci line attributes to length of CI (if needed, values are recycled)
    ci.col <- rep_len(ci.col, length(CI))
    ci.lty <- rep_len(ci.lty, length(CI))
    ci.lwd <- rep_len(ci.lwd, length(CI))

                                        # remove phi' if included in res
    phiLcol <- which(names(res) == "phiL")
    if(length(phiLcol > 0)){
        x <- res[,-phiLcol]
    } else {
        x <- res
    }
    d <- ncol(x)

                                        # First a fitting round to look for shifts in t_anchor
    age = seq(0, agemax, 0.01)
    Lt0 <- matrix(NaN, nrow=length(age), ncol=nrow(x))
    Lt_minus <- matrix(NaN, nrow=length(age), ncol=nrow(x))
    Lt_plus <- matrix(NaN, nrow=length(age), ncol=nrow(x))
    for(i in seq(ncol(Lt0))){
        par0 <- par_minus <- par_plus <- as.list(x[i,])
        par_minus$t_anchor <- par0$t_anchor-1
        par_plus$t_anchor <- par0$t_anchor+1
        Lt0[,i] <- TropFishR::VBGF(param = par0, t = age)
        Lt_minus[,i] <- TropFishR::VBGF(param = par_minus, t = age)
        Lt_plus[,i] <- TropFishR::VBGF(param = par_plus, t = age)
    }

                                        # replace negative lengths with NA
    Lt0 <- replace(Lt0, Lt0<0, NaN)
    Lt_minus <- replace(Lt_minus, Lt_minus<0, NaN)
    Lt_plus <- replace(Lt_plus, Lt_plus<0, NaN)

                                        # determine if a positive or negative shift improves overall covariance for each permutation
    cov0 <- cov(Lt0, use = "pair")
    shift <- 0 * seq(ncol(Lt0))
    for(i in seq(ncol(Lt0))){
        cov_minus <- cov(x = Lt_minus[,i], y = Lt0, use = "pair")
        cov_plus <- cov(x = Lt_plus[,i], y = Lt0, use = "pair")
        shift[i] <- c(0,-1,1)[which.max(c(sum(cov0[i,]), sum(cov_minus), sum(cov_plus)))]
    }

                                        # fixed predictions
    agenew <- seq(0+min(shift), agemax+max(shift), 0.01)
    Lt <- matrix(NaN, nrow=length(agenew), ncol=nrow(x))
    for(i in seq(ncol(Lt))){
        par0 <- as.list(x[i,])
        par0$t_anchor <- par0$t_anchor + shift[i]
        Lt[,i] <- TropFishR::VBGF(param = par0, t = agenew)
    }

    for(i in seq(ncol(Lt))){
        if(i == 1) {
            plot(agenew, Lt[,i], t="n",
                 xlim = c(min(x$t_anchor+shift), max(agenew))+c(-0.1,0),
                 ylim = c(0,max(Lt)*1.05),
                 xlab = xlab, ylab = ylab,
                 xaxs = "i", yaxs = "i"
                 )
        }
        lines(agenew, Lt[,i], col = perm.col, lwd = perm.lwd)
    }
    box()

                                        # multivariate kernel density estimate for max. density
    H <- ks::Hpi(x[,c("Linf","K","t_anchor")], nstage = 1)
    fhat <- ks::kde(x = x[,c("Linf","K","t_anchor")], H = H, eval.points = x[,c("Linf","K","t_anchor")])

                                        # maximum density
    max_dens <- fhat$eval.points[which.max(fhat$estimate),] # maximum density

                                        # predict density estimate of original data
    x$estimate <- fhat$estimate


                                        # determin which resamples are in the CI
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
        agenew, TropFishR::VBGF(param = as.list(max_dens), t = agenew),
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
               legend = c(paste0(names(max_dens), " = ", sprintf("%.2f", round(max_dens,2) ))),
               bty = "n", title = "Max. Dens. \nparameters:",
               inset = c(0,0.1)
               )
    }

    RES <- list(limCI = limCI, inCI = inCI, density = x$estimate, max_dens = as.list(max_dens))
    return(RES)
}



#' add phiprime contours to Linf_K scatterplot
#'
#' @param gridsize default = 20
#' @param ... Additional arguments passed to `par`.
#'
#' @return plot
#' @export
#'
add_phiprime <- function(gridsize = 20, ...){
    usr <- par()$usr
    Linf = seq(usr[1], usr[2], length.out = gridsize)
    K = seq(usr[3], usr[4], length.out = gridsize)
    grd <- expand.grid(
        Linf = Linf,
        K = K
    )
    grd$phiL <- log10(grd$K) + 2 * log10(grd$Linf)

    M <- list(x = Linf, y = K, z = matrix(grd$phiL, nrow = gridsize, ncol = gridsize))
    contour(x = M, add = TRUE, ...)
}




#' Univariate kernel density estimate plot of VBGF parameter from bootstrapping results
#'
#' @param res list of class `lfqBoot`
#' @param CI numeric. Confidence interval level (Default: 95)
#' @param use_hist logical. Plot histogram in addition to smoothed kernel density.
#' @param nbreaks numeric. Number of breaks in the histogram.
#' @param truePar if true parameters known lines can be added to graph
#' @param mar vector. Inner margins settings. See `?par`.
#' @param oma vector. Outer margins settings.See `?par`.
#' @param mgp vector. See `?par`.
#' @param tcl vector. See `?par`.
#' @param mfrow vector. See `?par` (default: NA).
#' @param cex numeric. See `?par`.
#' @param ... Additional arguments passed to `par`.
#'
#' @return plot
#' @export
#'
#' @examples
#'
#' data(alba_boot)
#' univariate_density(alba_boot)
#'
#'
univariate_density <- function(res, CI=95, use_hist = FALSE, nbreaks = 10,
                               truePar = NA, display_val=TRUE, display_median=FALSE,
                               display_legend = FALSE, display_rugs = FALSE,
                               mar = c(1.5,2,2,0), oma = c(1.5,0,0,0.5),
                               mgp = c(2,0.5,0), tcl = -0.25, cex = 1, mfrow=NA,
                               ylim=NA,...
                               ){
    res <- res$bootRaw
    na <- ifelse(display_legend,ncol(res)+1,ncol(res))
    nc <- min(c(6,na))
    nr <- ceiling(na/nc)
    if(is.na(mfrow[1])) mfrow <- c(nr,nc)
    op <- par(no.readonly = TRUE)
    par(
                                        # mfcol = c(floor(sqrt(ncol(res))), ceiling(sqrt(ncol(res)))),
        mfcol = c(1, ncol(res)),
        mar = mar, oma = oma, mfrow = mfrow,
        mgp = mgp, tcl = tcl, cex = cex, ...
    )

    VARS <- list(
        Linf = expression(bolditalic(L)[infinity]),
        K = expression(bolditalic(K)),
        t_anchor = expression(bolditalic(t)[anchor]),
        C = expression(bolditalic(C)),
        ts = expression(bolditalic(t)[s]),
        phiL = expression(bolditalic("phi'")),
        M_Then = expression(bolditalic(M)),
        M_Pauly = expression(bolditalic(M)),
        Z = expression(bolditalic(Z)),
        E = expression(bolditalic(E)),
        FM = expression(bolditalic(F)),
        L50 = expression(bolditalic("L50")),
        L75 = expression(bolditalic("L75")),
        N = expression(bolditalic(N)),
        B = expression(bolditalic(B)),
        F01 = expression(bolditalic("F01")),
        Fmax = expression(bolditalic(Fmax)),
        F05 = expression(bolditalic("F05")),
        FF01 = expression(bolditalic("FF01")),
        FFmax = expression(bolditalic(FFmax)),
        FF05 = expression(bolditalic("FF05")),
        ypr = expression(bolditalic(ypr)),
        yprrel = expression(bolditalic(yprrel)),
        bpr = expression(bolditalic(bpr)),
        bprrel = expression(bolditalic(bprrel))
    )


    ## univariate plots
    for(i in seq(ncol(res))){
        tmp <- as.numeric(na.omit(res[,i]))
        x <- try(kde(tmp),silent=TRUE)
        if(class(x) == "try-error"){
            plot(seq(0,1,0.2),
                 seq(min(tmp)*0.9,max(tmp)*1.1,length.out = 6),
                 t="n",
                 xlim = xlim,
                 xaxs = "i",
                 ylab="", xlab="", col=1, lty=1
                 )
        }else{
            h = hist(tmp, plot=FALSE, breaks = nbreaks)
            xlim <- c(0, max(x$estimate))
            if(use_hist){
                xlim <- range(c(xlim, max(h$density)))
            }
            xlim <- xlim * c(0,1.1)

            ylims <- range(x$eval.points)
            if(any(!is.null(ylim)) && any(!is.na(ylim))){
                if(!is.na(ylim[[i]][1])){
                    ylims[1] <- ylim[[i]][1]
                    ylims[2] <- ylim[[i]][2]
                }
            }


            plot(x$estimate, x$eval.points, t="n",
                 xlim = xlim,
                 ylim = ylims,
                 xaxs = "i",
                 ylab="", xlab="", col=1, lty=1
                 )
            usr <- par()$usr

            CItxt <- paste0(round(100-CI), "%")
            inCI <- rle( x$estimate > x$cont[CItxt] )
            start.idx <- c(1, cumsum(inCI$lengths[-length(inCI$lengths)])+1)
            end.idx <- cumsum(inCI$lengths)
            limCI <- range(x$eval.points[start.idx[min(which(inCI$values))]:end.idx[max(which(inCI$values))]])

            in1 <- which(x$estimate > x$cont["99%"])
            mean1 <- mean(x$eval.points[in1])

            if(use_hist){
                rect(
                    xleft = 0, ybottom = h$breaks[-length(h$breaks)],
                    xright = h$density, ytop = h$breaks[-1],
                    col = "grey90", border = "grey50"
                )
            }else{
                for(j in seq(inCI$lengths)){
                    if(inCI$values[j]){
                        polygon(
                            y = c(x$eval.points[start.idx[j]:end.idx[j]], rev(x$eval.points[start.idx[j]:end.idx[j]])),
                            x = c(x$estimate[start.idx[j]:end.idx[j]], rep(0, length(x$estimate[start.idx[j]:end.idx[j]]))),
                            col = "grey90", #col = rgb(0,1,0,0.2),
                            border = NA, lty = 3, lwd = 1
                        )
                    }
                }
            }
        }

        ## abline(v = x$cont[CItxt], lty=2, col="grey50")
        if(class(x) != "try-error"){
            lines(x$estimate, x$eval.points, lwd = 1, col = "grey50")

            ## rug
            if(display_rugs)
                segments(x0 = 0, x1 = par()$cxy[1]*0.3, y0 = x$x, y1 = x$x, col=rgb(0,0,0,0.5), lwd=0.3)

            ## true pars
            if(!is.na(truePar[1])){
                trupari <- truePar[match(names(res)[i], names(truePar))]
                abline(h=trupari, lty=2, col=2)
            }

            ## range of CI
            abline(h = limCI, lty = 2, lwd=1, col = 1)
            if(display_val) text(y =c(limCI), x = mean(usr[1:2]),
                                 labels = paste(sprintf("%.2f", round(c(limCI),2))),
                                 pos = c(1,3), offset = 0.25, col = 1)

            ## median
            if(display_median) abline(h=median(x$x), col=4, lty=2)

            ## maximum density
            abline(h = mean1, lty = 1, lwd=1, col = 1)
            if(display_val) text(y =  mean1, x = mean(usr[1:2]),
                                 labels = sprintf("%.2f", round(mean1,2)),
                                 pos = 3,
                                 offset = 0.25, col = 1
                                 )


        }else{
            abline(h=unique(tmp), lty = 3)
        }
        varlab <- VARS[[match(names(res)[i], names(VARS))]]
        mtext(varlab, line=0.25, side=3, font = 2)
        box(lwd=1.2)
    }
    mtext("Density", side = 1, line = 0, outer = TRUE)

    if(display_legend){
        plot.new()
        legend("center",
               legend=c("true","median","max density","CI"),
               lty=c(2,2,1,2), col=c(2,4,1,1), bty="n")
    }

    par(op)
}







#' Link/K scatterplot of bootstrapping results
#'
#' @param res list of class `lfqBoot`
#' @param Linf.breaks vector. Breaks for Linf histogram.
#' @param K.breaks vector. Breaks for K histogram.
#' @param gridsize vector. 2 values for defining the resolution of the grid
#' @param H object from \code{\link[ks]{Hpi}} (Default: `ks::Hpi(res[,c("Linf", "K")])`)
#' @param shading logical. Should 2d field of density estimates be colored with colors
#'   specified by `shading.cols` argument (Default: `shading = TRUE`).
#' @param shading.cols vector. Colors for background shading of 2d field of density estimates
#'   (Default: `shading.cols = colorRampPalette(c("white", blues9))(50)`).
#' @param dens.contour logical. Should contour lines be added (Default: `dens.contour = TRUE`)
#' @param probs vector. Density probability cutoffs to be plotted by contours when
#'   `dens.contour = TRUE` (Default: `probs = c(25,50,75,95)`).
#' @param phi.contour logical. Should phi prime isolines be displayed
#'   (Default: `phi.contour = FALSE`)
#' @param phi.levels vector. Phi prime values to display when `phi.contour = TRUE`
#'   (Default: `phi.levels = NULL`). When not defined (`phi.levels = NULL`), the default
#'   levels are chosen automatically by the \code{\link[graphics]{contour}} function.
#' @param phi.contour.col vector. Color to use for phi prime contours.
#'   Passed to \code{\link[graphics]{contour}}.
#' @param phi.contour.lty vector. Line type to use for phi prime contours.
#'   Passed to \code{\link[graphics]{contour}}.
#' @param phi.contour.lwd vector. Line width to use for phi prime contours.
#'   Passed to \code{\link[graphics]{contour}}.
#' @param phi.contour.labcex vector. Labels for the contour lines.
#'   Passed to \code{\link[graphics]{contour}}.
#' @param pt.pch pch value to use for resampling points
#' @param pt.col color to use for resampling points
#' @param pt.cex size value to use for resampling points
#' @param pt.bg background color to use for resampling points
#' @param xlab label for x-axis
#' @param ylab lavel for y axis
#' @param xlim limits for x-axis
#' @param ylim limits for y-axis
#'
#' @return plot
#'
#' @importFrom ks Hpi
#'
#' @examples
#' data(alba_boot)
#' LinfK_scatterhist(alba_boot)
#'
#' @export
#'
LinfK_scatterhist <- function(res, Linf.breaks = "Sturges", K.breaks = "Sturges",
                              gridsize = 151, H = ks::Hpi(res[,c("Linf", "K")]),
                              shading = TRUE, shading.cols = colorRampPalette(c("white", blues9))(50),
                              dens.contour = TRUE, probs = c(25,50,75,95),
                              phi.contour = FALSE, phi.levels = NULL,
                              phi.contour.col = 8, phi.contour.lty = 2, phi.contour.lwd = 1,
                              phi.contour.labcex = 0.75,
                              pt.pch = 16, pt.col = rgb(0,0,0,0.25), pt.cex = 0.5, pt.bg = 4,
                              xlab=expression(italic("L")[infinity]), ylab=expression(italic("K")),
                              xlim = NULL, ylim = NULL
                              ){

    res <- res$bootRaw
    zones <- matrix(c(2,0,1,3), ncol=2, byrow=TRUE)
    op <- par(no.readonly = TRUE)
    nf <- layout(zones, widths=c(4/5,1/5), heights=c(1/5,4/5), respect = FALSE)
                                        # layout.show(nf)
    par(cex = 1)

                                        # histogram data
    xhist = hist(res[,"Linf"], plot=FALSE, breaks = Linf.breaks)
    yhist = hist(res[,"K"], plot=FALSE, breaks = K.breaks)
    top = max(c(xhist$counts, yhist$counts))

    ## density estimation
    binFlag <- ifelse(nrow(res) >= 500 | nrow(res) < 20, FALSE, TRUE)   ## necessary for large sample sizes
    par(mar=c(4,4,0,0), mgp = c(2,0.5,0), tcl = -0.25)
    kk <- try(ks::kde(
                      x = res[,c("Linf", "K")], gridsize = gridsize, H = H,
                      binned = binFlag
                  ),silent = TRUE)
    if(class(kk) == "try-error"){
        kk <- ks::kde(
                      x = res[,c("Linf", "K")], gridsize = gridsize, H = H,
                      binned = FALSE
                  )
    }



                                        # plot limits (defaults are eval.points of kde)
    if(is.null(xlim)){xlim <- range(kk$eval.points[[1]])}
    if(is.null(ylim)){ylim <- range(kk$eval.points[[2]])}

                                        # 2d density plot
    image(
        x = kk$eval.points[[1]], y = kk$eval.points[[2]], z = kk$estimate,
        col = if(shading){shading.cols}else{NA},
        xlab = xlab, ylab = ylab, xlim=xlim, ylim=ylim
    )
    if(phi.contour){
        if(is.null(phi.levels)){
            add_phiprime(
                col = phi.contour.col, lty = phi.contour.lty, lwd = phi.contour.lwd,
                labcex = phi.contour.labcex
            )
        }else{
            add_phiprime(
                levels = phi.levels,
                col = phi.contour.col, lty = phi.contour.lty, lwd = phi.contour.lwd,
                labcex = phi.contour.labcex
            )
        }
    }
    points(kk$x[,1], kk$x[,2], pch = pt.pch, cex = pt.cex, col = pt.col, bg = pt.bg)
    if(dens.contour){
        plot(kk, type = "slice", add = TRUE, cont = probs)
    }
    box()

                                        # x histogram
    par(mar=c(0,4,1,0))
    plot(xhist$mids, xhist$counts, axes=FALSE, xlab="", ylab="", t="n", ylim=c(0, top), xlim = xlim, yaxs="i", xaxs="i")
    rect(
        xleft = xhist$breaks[-length(xhist$breaks)], ybottom = 0,
        xright = xhist$breaks[-1], ytop = xhist$counts,
        col = 8, border = 1
    )

                                        # y histogram
    par(mar=c(4,0,0,1))
    plot(yhist$counts, yhist$mids, axes=FALSE, xlab="", ylab="", t="n", xlim=c(0, top), ylim = ylim, yaxs="i", xaxs="i")
    rect(
        xleft = 0, ybottom = yhist$breaks[-length(yhist$breaks)],
        xright = yhist$counts, ytop = yhist$breaks[-1],
        col = 8, border = 1
    )

    par(op)
}




## x = lfq
## boot = resBoot
## Fname = "catch"

## plotting bootstrapped lfq data
plotLFQboot <- function(x,    ## lfq object
                        boot, ## boot object
                        Fname = "rcounts",  # alternative : "catch"
                        par = NULL,
                        agemax = NULL,
                        rel = FALSE,
                        curve.col = 1,
                        hist.sc = 0.5,
                        hist.col = c("white", "black", "orange","darkgreen"),
                        image.col = NULL,
                        region.col = NULL,
                        zlim = NULL,
                        zlimtype = "balanced",   # alternative : "range"
                        date.axis = "traditional",  # alternative : "modern"
                        date.at = seq(as.Date("1500-01-01"), as.Date("2500-01-01"), by="months"),
                        date.format = "'%y-%b", xlab = "", ylab = "Length classes",
                        draw = TRUE,
                        ...){


    dates <- x$dates
    classes <- x$midLengths
    catch <- get(Fname, x)


    for(i in 1:ncol(catch)){

        catchlist <- lapply(tmp, function(z) z$catch[,i])
        catchmat <- do.call(cbind, catchlist)


        r <- rep(NA,nrow(catchmat))
        for(ii in 1:nrow(catchmat)){
            z <- catchmat[ii,]

            if(all(z == 0)){
                r[ii] <- 0
            }else{
                tmp <- ks::kde(z)
                r[ii] <- round(mean(tmp$eval.points[which(tmp$estimate > tmp$cont["99%"])]))
            }
        }



        barplot(catch[,1], col = rgb(1,0.1,0.4,0.2))
        barplot(r,add=TRUE, col = rgb(0.1,1,0.4,0.2))


    }



    ## display relative catches (relative to number of samples per month)
    if(rel){
        catchrel <- catch
        for(i in 1:ncol(catch)){
            catchrel[,i] <- catch[,i]/colSums(catch, na.rm = TRUE)[i]
        }
        catch <- catchrel
        catch[is.nan(catch)] <- 0

        ## combined lfq plot
        if(any(!is.na(y))){
            catchrelY <- catchY
            for(i in 1:ncol(catchY)){
                catchrelY[,i] <- catchY[,i]/colSums(catchY, na.rm = TRUE)[i]
            }
            catchY <- catchrelY
            catchY[is.nan(catchY)] <- 0
        }
    }


    ## bin height scaling
    sc <- unclass(min(diff(dates)) * hist.sc / max(abs(catch)))

    if(any(!is.na(y))){
        ## bin height scaling
        scY <- unclass(min(diff(dates)) * hist.sc / max(abs(catchY)))
    }

    bin.width <- diff(classes)
    bin.lower <- classes - c(bin.width[1], bin.width)/2
    bin.upper <- classes + c(bin.width, bin.width[length(bin.width)])/2

                                        # image colour
    if(is.null(image.col)){
        pal <- colorRampPalette(c(rgb(1,0.8,0.8), rgb(1,1,1), rgb(0.8,0.8,1)))
        image.col <- pal(21)
    }
    if(!is.null(region.col)){
        image.col <- NA
    }

                                        # zlim value + type
    if(is.null(zlim) & zlimtype == "balanced"){
        zlim = c(-1,1) * max(abs(catch), na.rm=TRUE)
    }
    if(is.null(zlim) & zlimtype == "range"){
        zlim = range(catch, na.rm = TRUE)
    }

    ## Initial plot
    if(any(!is.na(y))){
        catchComb <- catch + catchY
        image(
            x=mergi$dates, y=mergi2$classes, z=t(catchComb), col=image.col, zlim=zlim,
            xaxt="n", xlab = xlab, ylab = ylab, ...)
    }else{
        image(
            x=dates, y=classes, z=t(catch), col=image.col, zlim=zlim,
            xaxt="n", xlab = xlab, ylab = ylab, ...)
    }

    if(!is.null(region.col)){
        usr <- par()$usr
        if(par()$xlog) usr[1:2] <- 10^usr[1:2]
        if(par()$ylog) usr[3:4] <- 10^usr[3:4]
        rect(usr[1], usr[3], usr[2], usr[4], col=region.col)
    }

                                        # add time axis
    if(date.axis == "modern"){
        axis.Date(side = 1, x=dates, at=date.at, format = date.format)
    }else if(date.axis == "traditional"){
        axis.Date(side = 1, x = dates, at = date.at, format = "%b")
        year <- seq(min(as.numeric(format(dates, "%Y"))), max(as.numeric(format(dates, "%Y"))), 1)
        date_seq <- seq.Date(dates[1],dates[length(dates)], by = "month")
        date_label <- format(date_seq, "%m")
        year_pre <- which(date_label %in% "01")
        if(!(1 %in% year_pre)) year_pre <- c(1,which(date_label %in% "01"))
        dates_for_years <- as.Date(paste(format(date_seq,"%Y"),date_label,"01",sep="-"))
        year_ticks <- dates_for_years[year_pre]
        mtext(side = 1, at = year_ticks, text = year, line = 2.5)
    }

    ## Histograms
    if(any(!is.na(y))){
        bin.width <- diff(mergi2$classes)
        bin.lower <- mergi2$classes - c(bin.width[1], bin.width)/2
        bin.upper <- mergi2$classes + c(bin.width, bin.width[length(bin.width)])/2
        for(i in seq(length(mergi$dates))){
            score.sc <- catch[,i] * sc
            score.scY <- catchY[,i] * scY
            for(j in seq(mergi2$classes)){
                polygon(
                    x = c(mergi$dates[i], mergi$dates[i], mergi$dates[i]-score.sc[j], mergi$dates[i]-score.sc[j]),
                    y = c(bin.lower[j], bin.upper[j], bin.upper[j], bin.lower[j]),
                    col = hist.col[(score.sc[j]>0)+1],
                    border = "grey20", lwd = 1)
                polygon(
                    x = c(mergi$dates[i], mergi$dates[i], mergi$dates[i]-score.scY[j], mergi$dates[i]-score.scY[j]),
                    y = c(bin.lower[j], bin.upper[j], bin.upper[j], bin.lower[j]),
                    col = hist.col[(score.scY[j]>0)+3],
                    border = "grey20", lwd = 1)

            }
        }
    }else{
        for(i in seq(length(dates))){
            score.sc <- catch[,i] * sc
            for(j in seq(classes)){
                                        # if(score.sc[j] != 0){
                polygon(
                    x = c(dates[i], dates[i], dates[i]-score.sc[j], dates[i]-score.sc[j]),
                    y = c(bin.lower[j], bin.upper[j], bin.upper[j], bin.lower[j]),
                    col = hist.col[(score.sc[j]>0)+1],
                    border = "grey20", lwd = 1)
                                        # }
            }
        }
    }

                                        # optional addition of cohort growth curves
    if("par" %in% names(x) & is.null(par) & draw){
        Lt <- lfqFitCurves(lfq = x, par = x$par,
                           agemax = x$agemax, draw = TRUE, col=curve.col
                           )
    }
    if(!is.null(par) & draw){
        Lt <- lfqFitCurves(x, par = par,
                           agemax = agemax, draw = TRUE, col=curve.col
                           )
    }

                                        # frame
    box()
}





#' @name alba_boot
#'
#' @title Bootstrapped VBGF estimates for the alba lfq dataset
#'
#' @description Bootstrapped VBGF estimates for the alba lfq dataset as estimated by
#'   \code{\link[TropFishR]{ELEFAN_GA}}.
#'
#'
#' @docType data
#'
#' @format A list of class `lfqBoot`
#'
#'
#' @usage data(alba_boot)
#' @keywords data dataset length-frequency bootstrap
#'
#' @examples
#'
#' \donttest{
#' #### For documentation only - How data was produced ####
#' library(parallel)
#' library(TropFishR)
#' library(ks)
#'
#' data(alba)
#'
#' # ELEFAN_GA_boot settings
#' MA <- 7
#' low_par <- list(Linf = 8, K = 0.1, t_anchor = 0, C = 0, ts = 0)
#' up_par <- list(Linf = 15, K = 5, t_anchor = 1, C = 1, ts = 1)
#' seasonalised <- FALSE
#' popSize <- 100
#' maxiter <- 50
#' run <- 10
#' pmutation <- 0.2
#' nresamp <- 200
#'
#' # Bootstrapped estimates
#' alba_boot <- ELEFAN_GA_boot(lfq=alba, MA = MA, seasonalised = seasonalised,
#'   up_par = up_par, low_par = low_par,
#'   popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
#'   nresamp = nresamp, parallel = TRUE, no_cores = detectCores(),
#'   seed = 20
#' )
#' }
#'
#'
#' data(alba_boot)
#'
#' head(alba_boot$bootRaw)
#' head(alba_boot$seed)
#'
#' # plot
#' univariate_density(alba_boot)
#' LinfK_scatterhist(alba_boot, phi.contour = TRUE)
#'
#'
NULL




#' @name plotBoot
#'
#' @title Bootstrapped VBGF estimates for the alba lfq dataset
#'
#' @description Bootstrapped VBGF estimates for the alba lfq dataset as estimated by
#'   \code{\link[TropFishR]{ELEFAN_GA}}.
#'
#' @examples
#' \donttest{
#' #### For documentation only - How data was produced ####
#' library(parallel)
#' library(TropFishR)
#' library(ks)
#'
#' data(alba)
#'
#' # ELEFAN_GA_boot settings
#' MA <- 7
#' low_par <- list(Linf = 8, K = 0.1, t_anchor = 0, C = 0, ts = 0)
#' up_par <- list(Linf = 15, K = 5, t_anchor = 1, C = 1, ts = 1)
#' seasonalised <- FALSE
#' popSize <- 100
#' maxiter <- 50
#' run <- 10
#' pmutation <- 0.2
#' nresamp <- 200
#'
#' # Bootstrapped estimates
#' alba_boot <- ELEFAN_GA_boot(lfq=alba, MA = MA, seasonalised = seasonalised,
#'   up_par = up_par, low_par = low_par,
#'   popSize = popSize, maxiter = maxiter, run = run, pmutation = pmutation,
#'   nresamp = nresamp, parallel = TRUE, no_cores = detectCores(),
#'   seed = 20
#' )
#' }
#'
#' plotBoot(lfq=alba, boot=alba_boot, plot="catchCurve", display="CI")
#' plotBoot(lfq=alba, boot=alba_boot, plot="catchCurve", display="bootRaw")
#'
#' @export
#'




plotBoot <- function(lfq,    ## lfq object
                     boot, ## boot object
                     pType="lfq",   ## alternative: catchCurve
                     display="CI",  ## alternative: bootRaw
                     best="maxDen",  ## alternative: "median"
                     Fname = "catch",  # alternative : "catch"
                     pCCori = FALSE,
                     pRef="fmax",
                     par = NULL,
                     agemax = NULL,
                     rel = FALSE,
                     curve.col = 1,
                     hist.sc = 0.5,
                     hist.col = c("white", "black", "orange","darkgreen"),
                     image.col = NULL,
                     region.col = NULL,
                     zlim = NULL,
                     zlimtype = "balanced",   # alternative : "range"
                     date.axis = "traditional",  # alternative : "modern"
                     date.at = seq(as.Date("1500-01-01"), as.Date("2500-01-01"), by="months"),
                     date.format = "'%y-%b", xlab = "default", ylab = "default",
                     draw = TRUE, ylim=NULL, xlim=NULL,
                     CI = 95,
                     add_legend = TRUE, add_max_dens_legend = TRUE,
                     points.col = "black", points.pch=19,
                     perm.col = adjustcolor("grey50",0.1), perm.lwd = 1,
                     ci.col = 1, ci.lty = 2, ci.lwd = 1,
                     maxd.col = 1, maxd.lty = 1, maxd.lwd = 2,
                     ...){

    dates <- lfq$dates
    classes <- lfq$midLengths

    if(pType=="lfq"){

        catch <- get(Fname, lfq)


        xlabel=""
        ylabel="Length classes"

        ## use user defined labels if given
        if(xlab != "default") xlabel = xlab
        if(ylab != "default") ylabel = ylab

        ## display relative catches (relative to number of samples per month)
        if(rel){
            catchrel <- catch
            for(i in 1:ncol(catch)){
                catchrel[,i] <- catch[,i]/colSums(catch, na.rm = TRUE)[i]
            }
            catch <- catchrel
            catch[is.nan(catch)] <- 0
        }
        ## bin height scaling
        sc <- unclass(min(diff(dates)) * hist.sc / max(abs(catch)))

        bin.width <- diff(classes)
        bin.lower <- classes - c(bin.width[1], bin.width)/2
        bin.upper <- classes + c(bin.width, bin.width[length(bin.width)])/2

        ## image colour
        if(is.null(image.col)){
            pal <- colorRampPalette(c(rgb(1,0.8,0.8), rgb(1,1,1), rgb(0.8,0.8,1)))
            image.col <- pal(21)
        }
        if(!is.null(region.col)){
            image.col <- NA
        }

        ## zlim value + type
        if(is.null(zlim) & zlimtype == "balanced"){
            zlim = c(-1,1) * max(abs(catch), na.rm=TRUE)
        }
        if(is.null(zlim) & zlimtype == "range"){
            zlim = range(catch, na.rm = TRUE)
        }


        ## Initial plot
        image(
            x=dates, y=classes, z=t(catch), col=image.col, zlim=zlim,
            xaxt="n", xlab = xlab, ylab = ylab, ...)


        if(!is.null(region.col)){
            usr <- par()$usr
            if(par()$xlog) usr[1:2] <- 10^usr[1:2]
            if(par()$ylog) usr[3:4] <- 10^usr[3:4]
            rect(usr[1], usr[3], usr[2], usr[4], col=region.col)
        }
        ## add time axis
        if(date.axis == "modern"){
            axis.Date(side = 1, x=dates, at=date.at, format = date.format)
        }else if(date.axis == "traditional"){
            axis.Date(side = 1, x = dates, at = date.at, format = "%b")
            year <- seq(min(as.numeric(format(dates, "%Y"))), max(as.numeric(format(dates, "%Y"))), 1)
            date_seq <- seq.Date(dates[1],dates[length(dates)], by = "month")
            date_label <- format(date_seq, "%m")
            year_pre <- which(date_label %in% "01")
            if(!(1 %in% year_pre)) year_pre <- c(1,which(date_label %in% "01"))
            dates_for_years <- as.Date(paste(format(date_seq,"%Y"),date_label,"01",sep="-"))
            year_ticks <- dates_for_years[year_pre]
            mtext(side = 1, at = year_ticks, text = year, line = 2.5)
        }

        ## Histograms
        for(i in seq(length(dates))){
            score.sc <- catch[,i] * sc
            for(j in seq(classes)){
                                        # if(score.sc[j] != 0){
                polygon(
                    x = c(dates[i], dates[i], dates[i]-score.sc[j], dates[i]-score.sc[j]),
                    y = c(bin.lower[j], bin.upper[j], bin.upper[j], bin.lower[j]),
                    col = hist.col[(score.sc[j]>0)+1],
                    border = "grey20", lwd = 1)
                                        # }
            }
        }


        ## plotting cohort growth curves with uncertainty
        res <- boot$bootRaw
        if(is.null(agemax))
            agemax <- max(ceiling((1/-res$K)*log(1-((res$Linf*0.95)/res$Linf))))

        ## expand ci line attributes to length of CI (if needed, values are recycled)
        ci.col <- rep_len(ci.col, length(CI))
        ci.lty <- rep_len(ci.lty, length(CI))
        ci.lwd <- rep_len(ci.lwd, length(CI))

        ## remove phi' if included in res
        colind <- which(names(res) %in% c("Linf","K","t_anchor"))
        if("C" %in% names(res)) colind <- c(colind, which(names(res) %in% c("C","ts")))
        x <- res[,colind]
        d <- ncol(x)

        ## First a fitting round to look for shifts in t_anchor
        age = seq(0, agemax, 0.01)
        Lt0 <- matrix(NaN, nrow=length(age), ncol=nrow(x))
        Lt_minus <- matrix(NaN, nrow=length(age), ncol=nrow(x))
        Lt_plus <- matrix(NaN, nrow=length(age), ncol=nrow(x))
        for(i in seq(ncol(Lt0))){
            par0 <- par_minus <- par_plus <- as.list(x[i,])
            par_minus$t_anchor <- par0$t_anchor-1
            par_plus$t_anchor <- par0$t_anchor+1
            Lt0[,i] <- TropFishR::VBGF(param = par0, t = age)
            Lt_minus[,i] <- TropFishR::VBGF(param = par_minus, t = age)
            Lt_plus[,i] <- TropFishR::VBGF(param = par_plus, t = age)
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
            Lt[,i] <- TropFishR::VBGF(param = par0, t = agenew)
        }


        ## multivariate kernel density estimate for max. density
        H <- ks::Hpi(x, nstage = 1)
        fhat <- ks::kde(x = x, H = H, eval.points = x)

        ## maximum density
        max_dens <- fhat$eval.points[which.max(fhat$estimate),] ## maximum density

        ## predict density estimate of original data
        x$estimate <- fhat$estimate

        ## determine which resamples are in the CI
        inCI <- x$estimate > quantile(x$estimate, probs = (100-CI)/100 )
        limCI <- data.frame(t = agenew, min = NaN, max = NaN)
        for(i in seq(agenew)){
            limCI$min[i] <- min(Lt[i, which(inCI) ], na.rm = TRUE)
            limCI$max[i] <- max(Lt[i, which(inCI) ], na.rm = TRUE)
        }


        years <- as.numeric(as.character(unique(format(dates, "%Y"))))
        minyear <- min(years) - agemax
        maxyear <- max(years)
        years2 <- seq(minyear,maxyear,1)
        for(i in 1:length(years2)){
            tmp <- years2[i] + max_dens$t_anchor + agenew
            xvals <- TropFishR::yeardec2date(tmp)

            lines(xvals, limCI$max, col = ci.col, lwd = ci.lwd, lty = ci.lty)
            lines(xvals, limCI$min, col = ci.col, lwd = ci.lwd, lty = ci.lty)

            lines(xvals,
                  TropFishR::VBGF(param = as.list(max_dens), t = agenew),
                  col = maxd.col, lwd = maxd.lwd, lty = maxd.lty)
        }

        ## frame
        box(lwd=1.2)


        ## catch curve plot
    }else if(pType == "cc"){

        binSize <- boot$misc$binSize
        yearSel <- boot$misc$yearSel
        yearCombine <- boot$misc$yearCombine


        ints <- boot$misc$intsCC
        zs <- boot$bootRaw$Z
        intMed <- median(ints, na.rm = TRUE)

        ## confidence intervals
        citmp <- (100-95)/2/100
        intCI <- quantile(intMed,  probs = c(citmp, 1-citmp), na.rm = TRUE)
        zCI <- boot$CI[,which(colnames(boot$CI)=="Z")]

        ## max densities
        x <- try(ks::kde(as.numeric(na.omit(ints))), TRUE)
        if(class(x) != "try-error"){
            ind <- which(x$estimate > x$cont["99%"])
            intMaxDen <- mean(x$eval.points[ind])
        }else stop("Cannot estimate MaxDen of intercepts of the catch curve. Use median or single line display.")

        if(best=="maxDen"){
            grPars <- boot$maxDen
            intX <- intMed
            zX <- boot$maxDen[which(names(boot$maxDen)=="Z")]
        }else if(best=="median"){
            grPars <- boot$median
            intX <- intMaxDen
            zX <- boot$median[which(names(boot$median)=="Z")]
        }

        midLengths <- lfq$midLengths
        interval <- midLengths[2] - midLengths[1]
        ## L and t of lower length classes
        lowerLengths <- midLengths - (interval / 2)
        ## seasonalised
        if("C" %in% names(grPars)){
            t_midL <- VBGF(param = list(Linf = grPars[1], K = grPars[2],
                                        t0 = 0, C=grPars[4], ts=grPars[5]), L = midLengths)
            t_L1 <- VBGF(param = list(Linf = grPars[1], K = grPars[2],
                                      t0 = 0, C=grPars[4], ts=grPars[5]), L = lowerLengths)
        }else{
            t_midL <- VBGF(param = list(Linf = grPars[1], K = grPars[2],
                                        t0 = 0), L = midLengths)
            t_L1 <- VBGF(param = list(Linf = grPars[1], K = grPars[2],
                                      t0 = 0), L = lowerLengths)
        }

        ## delta t
        deti <- rep(NA,length(midLengths))
        for(x1 in 1:(length(deti)-1)){
            deti[x1] <- t_L1[x1+1] - t_L1[x1]
        }


        xplot <- t_midL
        xlabel <- "Relative age [years - t0]"
        xplotAGE <- t_midL
        ## yaxis
        if(!is.na(yearSel)){
            yearSel <- as.character(yearSel)
            ## warning if year not in dates
            dateYears <- format(lfq$dates, "%Y")
            if(all(yearSel %in% dateYears == FALSE)) stop(paste0("The selected year ", yearSel, " is not in the LFQ data!"))
            lfq$catch <- lfq$catch[,which(format(lfq$dates,"%Y") %in% yearSel)]
            lfq$dates <- lfq$dates[format(lfq$dates,"%Y") %in% yearSel]
        }
        ## vectorise
        if(!is.na(binSize)){
            lfq2 <- lfqModify(lfq, vectorise_catch = TRUE, bin_size = binSize)
        }else{
            lfq2 <- lfqModify(lfq, vectorise_catch = TRUE)
        }
        if(yearCombine && inherits(lfq2$catch, "matrix")) lfq2$catch <- rowSums(lfq2$catch)
        ## error if lfq data spans several years!
        if(class(lfq2$catch) == "matrix") stop("The lfq data spans several years, please subset for one year at a time!")
        catch <- lfq2$catch
        lnC_dt <- log(catch / deti)
        lnC_dt[which(lnC_dt == -Inf)] <- NA   ### OR zero???
        yplot <- lnC_dt
        ylabel <- "ln(N / dt)"


        ## remove all NAs and Infs
        temp <- cbind(xplot,yplot)
        temp <- as.matrix(na.exclude(temp))
        temp <- temp[(!(temp[,1] == Inf | temp[,1] == -Inf)),]
        temp <- temp[(!(temp[,2] == Inf | temp[,2] == -Inf)),]
        xplot <- temp[,1]
        yplot <- temp[,2]


        ##for final plot

        ## use user defined labels if given
        if(xlab != "default") xlabel = xlab
        if(ylab != "default") ylabel = ylab


        maxi1 <- which.min(abs(xplot-boot$misc$regint[1]))
        maxi2 <- which.min(abs(xplot-boot$misc$regint[2]))

##        ind <- 1:length(yplot)
##        ind <- ind[is.finite(yplot)]
##        maxi1 <- ind[which.max(yplot[is.finite(yplot)])]
##        maxi2 <- ind[which.max(xplot[is.finite(yplot)])]


        datcc <- boot$misc$datcc

        xrang <- range(unlist(lapply(datcc, function(x) ifelse(!is.na(x),x$x,NA))),na.rm=TRUE)
        yrang <- range(unlist(lapply(datcc, function(x) ifelse(!is.na(x),x$y,NA))),na.rm=TRUE)
        maxa1 <- xrang[1]
        maxa2 <- xrang[2]

        step <- 0.01

        if(is.null(ylim)){
##            minyplot <- ifelse(min(yplot,na.rm=TRUE) < 0, min(yplot,na.rm=TRUE),0)
##            maxyplot <- max(yplot[is.finite(yplot)],na.rm=TRUE) + 1
##            ylims <- c(minyplot,maxyplot)
            ylims <- c(0,max(yrang))
        }else ylims <- ylim


        if(is.null(xlim)){
##            xlims <- c(min(xplot[which(yplot > 0 & is.finite(yplot))], na.rm = TRUE),
##                       max(xplot[which(yplot > 0 & is.finite(yplot))], na.rm = TRUE))
            xlims <- xrang
        }else xlims <- xlim


        ## final plot
        plot(x = xrang, y = yrang, ylim = ylims, ty="n",
             xlab = xlabel, ylab = ylabel, xlim = xlims,
             cex = 1)

        if(pCCori) points(x = xplot, y = yplot)
        if("bootRaw" %in% display){
            for(i in 1:length(ints)){
                tmp <- range(datcc[[i]]$x)
                seqi <- seq(tmp[1]-step,tmp[2]+step,step)
                lines(seqi, ints[i]-seqi*zs[i], col="grey80",lwd=2)
            }
        }
        for(i in 1:length(datcc)){
            if(!all(is.na(datcc[[i]])))
                points(datcc[[i]]$x,datcc[[i]]$y, pch=points.pch, col=points.col)
        }
        if("CI" %in% display){
            seqi <- seq(maxa1-step,maxa2+step,step)
            polygon(x=c(seqi,rev(seqi)), y=c(intCI[1]-seqi*zCI[1],
                                             rev(intCI[2]-seqi*zCI[2])),
                    col=rgb(t(col2rgb("blue"))/255,alpha=0.3),border=NA)
            lines(seqi, intX-seqi*zX, col="blue",lwd=2.5)
        }
        if(pCCori) points(xplot[(maxi1):maxi2],yplot[(maxi1):maxi2],
                          pch=19,col="blue")
        box(lwd=1.2)

        ## yield per recruit plot
    }else if(pType == "ypr"){

        px <- boot$misc$FM_change
        xlabel <- "Fishing mortality"
        ylabel1 <- "Yield per recruit"
        ylabel2 <- "Biomass"
        ylabel3 <- "Value"

        py <- boot$misc$totY[[1]]

        ylims <- range(unlist(lapply(boot$misc$totY, range, na.rm=TRUE)))
        xlims <- range(px)


        ## fixed predictions
        totys <- matrix(NaN, nrow=length(px), ncol=nrow(boot$bootRaw))
        for(i in seq(ncol(totys))){
            totys[,i] <- boot$misc$totY[[i]]
        }


        toti <- boot$misc$totY


        if(pRef=="fmax"){
            nmax <- unlist(lapply(toti, max, na.rm=TRUE))
            fmaxR <- range(boot$bootRaw$Fmax, na.rm=TRUE)
            ## median
            medi <- median(nmax, na.rm = TRUE)

            ## confidence intervals
            citmp <- (100-CI)/2/100
            ciList <- quantile(nmax,  probs = c(citmp, 1-citmp), na.rm = TRUE)

            ## max densities
            x <- ks::kde(as.numeric(na.omit(nmax)))
            ind <- which(x$estimate > x$cont["99%"])
            maxden <- mean(x$eval.points[ind])

            ## determine which resamples are in the CI
            inCI <- x$estimate > quantile(x$estimate, probs = (100-CI)/100 )
            limCI <- data.frame(t = px, min = NaN, max = NaN)
            for(i in seq(px)){
                limCI$min[i] <- min(totys[i, which(inCI) ], na.rm = TRUE)
                limCI$max[i] <- max(totys[i, which(inCI) ], na.rm = TRUE)
            }

            ## final plot
            plot(x = px, y = py, ylim = ylims, ty="n",
                 xlab = xlabel, ylab = ylabel1, xlim = xlims,
                 cex = 1)

            rect(fmaxR[1],-1,fmaxR[2],ylims[2]*2, border=NA, col=rgb(t(col2rgb("blue"))/255,alpha=0.3))
            text(mean(fmaxR), y=ylims[2]*0.98, labels="Fmax", cex=1.2, font=2)

            nbreaks = 10

            tmp <- as.numeric(na.omit(boot$bootRaw$Fmax))
            x <- try(kde(tmp),silent=TRUE)
            if(class(x) != "try-error"){
                h = hist(tmp, plot=FALSE, breaks = nbreaks)
                tmp <- seq(fmaxR[1],fmaxR[2], length.out = length(h$breaks))
                tmp2 <- h$density
                tmp2 <- tmp2/max(tmp2)
                rect(xleft = tmp[-length(tmp)], ybottom = 0,
                     xright = tmp[-1], ytop = tmp2*(0.2*ylims[2]),
                     col = rgb(t(col2rgb("black"))/255,alpha=0.8),
                     border = "grey50")
            }

            for(i in 1:length(boot$misc$totY)){
                lines(px, boot$misc$totY[[i]], col = rgb(t(col2rgb("black"))/255,alpha=0.3))
            }

            ##        lines(px, limCI$max, col = ci.col, lwd = ci.lwd, lty = ci.lty)
            ##        lines(px, limCI$min, col = ci.col, lwd = ci.lwd, lty = ci.lty)

            lines(px, boot$misc$totY[[which.min(abs(nmax-maxden))]],
                  col = maxd.col, lwd = maxd.lwd, lty = maxd.lty)
            box(lwd=1.2)
        }else if(pRef=="f01"){
            n01 <- unlist(lapply(seq(length(boot$misc$nmax)),
                          function(x) boot$misc$totY[[x]][boot$misc$n01[x]]))

            f01R <- range(boot$bootRaw$F01, na.rm=TRUE)
            ## median
            medi <- median(n01, na.rm = TRUE)

            ## confidence intervals
            citmp <- (100-CI)/2/100
            ciList <- quantile(n01,  probs = c(citmp, 1-citmp), na.rm = TRUE)

            ## max densities
            x <- ks::kde(as.numeric(na.omit(n01)))
            ind <- which(x$estimate > x$cont["99%"])
            maxden <- mean(x$eval.points[ind])

            ## determine which resamples are in the CI
            inCI <- x$estimate > quantile(x$estimate, probs = (100-CI)/100 )
            limCI <- data.frame(t = px, min = NaN, max = NaN)
            for(i in seq(px)){
                limCI$min[i] <- min(totys[i, which(inCI) ], na.rm = TRUE)
                limCI$max[i] <- max(totys[i, which(inCI) ], na.rm = TRUE)
            }

            ## final plot
            plot(x = px, y = py, ylim = ylims, ty="n",
                 xlab = xlabel, ylab = ylabel1, xlim = xlims,
                 cex = 1)

            rect(f01R[1],-1,f01R[2],ylims[2]*2, border=NA, col=rgb(t(col2rgb("blue"))/255,alpha=0.3))
            text(mean(f01R), y=ylims[2]*0.98, labels="F01", cex=1.2, font=2)

            nbreaks = 10

            tmp <- as.numeric(na.omit(boot$bootRaw$Fmax))
            x <- try(kde(tmp),silent=TRUE)
            if(class(x) != "try-error"){
                h = hist(tmp, plot=FALSE, breaks = nbreaks)
                tmp <- seq(f01R[1],f01R[2], length.out = length(h$breaks))
                tmp2 <- h$density
                tmp2 <- tmp2/max(tmp2)
                rect(xleft = tmp[-length(tmp)], ybottom = 0,
                     xright = tmp[-1], ytop = tmp2*(0.2*ylims[2]),
                     col = rgb(t(col2rgb("black"))/255,alpha=0.8),
                     border = "grey50")
            }

            for(i in 1:length(boot$misc$totY)){
                lines(px, boot$misc$totY[[i]], col = rgb(t(col2rgb("black"))/255,alpha=0.3))
            }

            ##        lines(px, limCI$max, col = ci.col, lwd = ci.lwd, lty = ci.lty)
            ##        lines(px, limCI$min, col = ci.col, lwd = ci.lwd, lty = ci.lty)

            lines(px, boot$misc$totY[[which.min(abs(n01-maxden))]],
                  col = maxd.col, lwd = maxd.lwd, lty = maxd.lty)
            box(lwd=1.2)
        }

    }
}
