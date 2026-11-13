#' @title Raise length-frequency samples to the total catch
#'
#' @description The function raises the sampled numbers of fish per length
#'    class (`lfq$catch`) to the total catch, using the total landings in
#'    weight and the weight of the samples. This is required for methods that
#'    need absolute numbers, such as the length-based virtual population
#'    analysis (\code{\link{VPA}}). Methods that only use the relative length
#'    composition (e.g. ELEFAN, \code{\link{catchCurve}}) and the resampling
#'    (\code{\link{lfqResample}}) work with the sampled numbers.
#'
#' @param lfq a list of the class "lfq" consisting of following parameters:
#' \itemize{
#'   \item \strong{midLengths} midpoints of the length classes,
#'   \item \strong{dates} dates of sampling times (class Date),
#'   \item \strong{catch} matrix with the sampled numbers of fish per length
#'      class (rows) and sampling date (columns),
#'   \item \strong{par} (optional) a list with the parameters of the
#'      length-weight relationship (\code{a} and \code{b}), required if
#'      \code{sampleWeight} is not provided.
#' }
#' @param landings total landings (or catch) in kg: one value per sample
#'    (\code{by = "sample"}) or per year (\code{by = "year"}). With
#'    \code{by = "year"}, the values are matched to the years by their names
#'    if named, otherwise they are used in the order of the years.
#' @param sampleWeight optional; the weight of the samples in kg (one value
#'    per sample or year, as \code{landings}). By default
#'    (\code{sampleWeight = NULL}), the weight is estimated from the sampled
#'    numbers and the length-weight relationship: sum of
#'    \code{catch * a * midLengths^b}. Use the measured weight of the samples
#'    if available.
#' @param by \code{"sample"} (default) to raise each sample (column of
#'    \code{lfq$catch}) separately, or \code{"year"} to raise all samples of
#'    a year with one factor (landings of the year / weight of all samples of
#'    the year).
#' @param a,b optional; parameters of the length-weight relationship
#'    (W = a * L^b). By default taken from \code{lfq$par}.
#' @param LW_unit unit of the weight in the length-weight relationship
#'    (\code{a}): \code{"kg"} (default) or \code{"g"}.
#'
#' @details The raising factor of each sample (or year) is the landings
#'    divided by the weight of the sample(s). Raising each sample separately
#'    (\code{by = "sample"}) gives each sample the weight of its landings, e.g.
#'    months with large landings get more weight than months with small
#'    landings, independent of the number of fish measured.
#'
#'    Keep the original lfq data with the sampled numbers for all other
#'    analyses. In particular, the bootstrap functions
#'    (\code{\link{lfqResample}}, \code{\link{ELEFAN_GA_boot}}, ...) require
#'    the sampled numbers, because they resample the measured fish. Instead of
#'    raising the data beforehand, the landings can be passed directly to
#'    \code{\link{VPA}} (argument \code{landings}), which raises the data (or
#'    each resample with \code{boot}) internally.
#'
#' @return The lfq list with the raised catch (numbers of fish) and the
#'    element \code{raising}, a data frame with the landings, the sample
#'    weight, whether the sample weight was estimated, and the raising factor
#'    per sample or year.
#'
#' @examples
#' data(synLFQ7)
#' ## length-weight relationship (W in g)
#' synLFQ7$par <- list(a = 0.015, b = 3)
#' ## landings in kg at the sampling dates (one value per sample)
#' landings <- seq(20000, 40000, length.out = length(synLFQ7$dates))
#' lfqR <- lfqRaise(synLFQ7, landings = landings, LW_unit = "g")
#' lfqR$raising
#'
#' ## one raising factor per year
#' lfqY <- lfqRaise(synLFQ7, landings = c("2016" = 240000), by = "year",
#'                  LW_unit = "g")
#' lfqY$raising
#'
#' @export
#'
lfqRaise <- function(lfq, landings, sampleWeight = NULL, by = "sample",
                     a = NULL, b = NULL, LW_unit = "kg"){

    by <- match.arg(by, c("sample", "year"))
    if(!is.null(lfq$raising))
        stop("The catch in lfq has already been raised (lfq$raising).")
    if(is.null(lfq$midLengths) || is.null(lfq$catch))
        stop("lfq has to contain midLengths and catch.")

    catch <- lfq$catch
    isVec <- !is.matrix(catch) && !is.data.frame(catch)
    catch <- as.matrix(catch)
    if(nrow(catch) != length(lfq$midLengths))
        stop("The number of rows of the catch matrix does not match the number of length classes.")

    ## groups of samples with one raising factor each
    if(by == "sample"){
        group <- seq_len(ncol(catch))
        groupNames <- if(!is.null(lfq$dates) && length(lfq$dates) == ncol(catch))
                          as.character(lfq$dates) else as.character(group)
    }else{
        if(is.null(lfq$dates) || length(lfq$dates) != ncol(catch))
            stop("by = 'year' requires the dates of the samples (lfq$dates).")
        years <- format(lfq$dates, "%Y")
        groupNames <- sort(unique(years))
        group <- match(years, groupNames)
    }
    nG <- length(groupNames)

    ## match landings (and sample weights) to groups
    matchGroups <- function(x, name){
        if(length(x) != nG)
            stop(paste0("'", name, "' has to have one value per ", by, " (",
                        nG, "), but has ", length(x), "."))
        if(by == "year" && !is.null(names(x))){
            if(!all(groupNames %in% names(x)))
                stop(paste0("The names of '", name, "' do not match the years of the lfq data (",
                            paste(groupNames, collapse = ", "), ")."))
            x <- x[groupNames]
        }
        x <- as.numeric(x)
        if(any(!is.finite(x) | x < 0))
            stop(paste0("'", name, "' has to be finite and >= 0."))
        x
    }
    landings <- matchGroups(landings, "landings")

    ## weight of the samples in kg
    estimated <- is.null(sampleWeight)
    if(estimated){
        if(is.null(a)) a <- if(!is.null(lfq$par$a)) lfq$par$a else lfq$a
        if(is.null(b)) b <- if(!is.null(lfq$par$b)) lfq$par$b else lfq$b
        if(is.null(a) || is.null(b))
            stop("Please provide the sample weight (sampleWeight) or the parameters of the length-weight relationship (a and b).")
        a <- lwUnitA(list(a = a), LW_unit)$a
        w <- a * lfq$midLengths^b
        sampleWeight <- as.numeric(tapply(colSums(catch * w, na.rm = TRUE),
                                          factor(group, levels = seq_len(nG)), sum))
    }else{
        sampleWeight <- matchGroups(sampleWeight, "sampleWeight")
    }

    noFish <- sampleWeight == 0 & landings > 0
    if(any(noFish))
        stop(paste0("The following samples have no fish (or weight) but positive landings and cannot be raised: ",
                    paste(groupNames[noFish], collapse = ", "), "."))
    fac <- ifelse(sampleWeight > 0, landings / sampleWeight, 0)

    catch <- sweep(catch, 2, fac[group], "*")
    lfq$catch <- if(isVec) as.numeric(catch) else catch
    lfq$raising <- data.frame(group = groupNames, landings = landings,
                              sampleWeight = sampleWeight,
                              estimated = estimated, factor = fac,
                              stringsAsFactors = FALSE)
    names(lfq$raising)[1] <- by
    return(lfq)
}


## Stop if lfq contains raised catches (resampling requires the sampled
## numbers of fish) and warn if the catch is not in whole numbers
checkSampledLfq <- function(lfq, fun){
    if(!is.null(lfq$raising))
        stop(paste0(fun, " requires the sampled numbers of fish, but the catch in lfq has been raised (lfq$raising). Use the original lfq data and pass the landings to VPA (argument 'landings') instead."))
    cc <- as.numeric(as.matrix(lfq$catch))
    if(any(abs(cc - round(cc)) > 1e-8, na.rm = TRUE))
        warning(paste0("The catch in lfq contains non-integer values. ", fun,
                       " resamples the measured fish and requires the sampled numbers of fish, not raised catches or weights."),
                call. = FALSE)
    invisible(TRUE)
}
