#' @title Resampling of length-frequency data
#'
#' @description The function resamples the `lfq` data by sample date.
#'   Sampling is done in a non-parametric manner that follows the relative
#'   frequencies of the original data, allowing for individual counts to be
#'   selected more than once (i.e. `replace = TRUE` in \link[base]{sample}),
#'   and resulting in total counts (by sample) equal to the original data
#'   (unless `sampSize` is provided).
#'
#' @param lfq a length frequency object of the class `lfq`.
#' @param sampSize optional; number of individuals drawn per sample. By default
#'    (`sampSize = NA`), the number of individuals of each sample in the
#'    original data is used.
#' @param boot optional; an object of class `lfqBoot` (e.g. from
#'    \code{\link{ELEFAN_GA_boot}}) for reproducing the resampled lfq data sets
#'    by re-using its seed values.
#'
#' @return a resampled version of the `lfq` class dataset, or, if `boot` is
#'    provided, a list with one resampled `lfq` dataset per seed value in `boot`.
#'
#' @examples
#' # load data
#' data(alba)
#'
#' # resample data
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
#' @importFrom graphics hist
#'
#' @export
#'
lfqResample <- function(lfq, sampSize = NA, boot = NULL){

    ## reproduce used lfq data sets if boot object provided
    nI <- ifelse(!is.null(boot) && "seed" %in% names(boot), length(boot$seed), 1)

    lfqList <- vector("list", nI)
    for(I in 1:nI){

        lfqb <- lfq

        if(!is.null(boot)) set.seed(boot$seed[I])

        for(i in seq(length(lfq$dates))){
            if(is.na(sampSize) || is.nan(sampSize)){
                sampSizi <- sum(lfq$catch[,i])
            }else{
                sampSizi <- sampSize
            }
            inds <- sample(x = lfq$midLengths, size = round(sampSizi),
                           prob = lfq$catch[,i], replace = TRUE)
            ## bin limits (allowing for uneven bin sizes)
            bin.width <- diff(lfq$midLengths)
            bin.lower <- lfq$midLengths - (c(bin.width[1], bin.width)/2)
            bin.upper <- lfq$midLengths + (c(bin.width, bin.width[length(bin.width)])/2)
            breaks <- unique(c(bin.lower, bin.upper))
            h <- hist(inds, breaks = breaks, plot = FALSE)
            lfqb$catch[,i] <- h$counts
        }
        lfqList[[I]] <- lfqb
    }

    if(is.null(boot)) lfqList <- lfqList[[1]]

    return(lfqList)
}
