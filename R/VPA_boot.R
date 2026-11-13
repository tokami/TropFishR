## Bootstrapped length-based VPA / cohort analysis (internal, called by VPA
## with boot)

VPA_boot <- function(param, boot, catch_unit = NA, catch_corFac = NA,
                     LW_unit = "kg", landings = NULL, sampleWeight = NULL,
                     raiseBy = "sample",
                     analysis_type = "VPA", algorithm = "new",
                     natMort = NULL, yearSel = NA, yearCombine = FALSE,
                     binSize = NA, CI = 95){

    if(is.null(natMort))
        stop("Please provide the name of the column with natural mortality estimates in boot$bootRaw (natMort).")
    checkBoot(boot, c("Linf", "K", natMort, "FM"))
    if(!("midLengths" %in% names(param)))
        stop("The bootstrapped VPA requires length-frequency data (midLengths).")

    ## length-weight relationship
    a <- if(!is.null(param$par$a)) param$par$a else param$a
    b <- if(!is.null(param$par$b)) param$par$b else param$b
    if(is.null(a) || is.null(b))
        stop("VPA requires information about the length-weight relationship. Please provide 'a' and 'b' estimates in param.")
    checkSampledLfq(param, "VPA with boot")

    ## raising of each resample to the total catch
    raise <- NULL
    if(!is.null(landings)){
        if(!is.na(catch_unit))
            stop("With 'landings', the catch is raised to numbers of fish: please leave catch_unit = NA.")
        raise <- list(landings = landings, sampleWeight = sampleWeight,
                      by = raiseBy, a = a, b = b, LW_unit = LW_unit)
    }

    bootRaw <- boot$bootRaw
    n <- nrow(bootRaw)
    N <- B <- rep(NA_real_, n)
    FMvecVPA <- vector("list", n)
    skip <- rep(NA_character_, n)

    for(x in seq_len(n)){

        M <- bootRaw[[natMort]][x]
        terminalF <- bootRaw$FM[x]
        if(!is.finite(M)){
            skip[x] <- "M not available"
            next
        }
        if(!is.finite(terminalF) || terminalF <= 0){
            skip[x] <- "FM not available or <= 0"
            next
        }

        lfqx <- bootLfq(param, boot, x, yearSel = yearSel,
                        yearCombine = yearCombine, binSize = binSize,
                        catch_corFac = catch_corFac, raise = raise)
        growth <- bootGrowthPar(boot, x)

        ## plus group at Linf of the resample
        interval <- lfqx$midLengths[2] - lfqx$midLengths[1]
        upperLength <- lfqx$midLengths + interval / 2
        if(growth$Linf < max(upperLength)){
            pg <- lfqx$midLengths[which.min(abs(upperLength - floor(growth$Linf)))]
            lfqx <- lfqModify(lfqx, plus_group = pg)
            plus_group <- TRUE
        }else{
            plus_group <- FALSE
        }

        p <- list(midLengths = lfqx$midLengths, catch = lfqx$catch,
                  Linf = growth$Linf, K = growth$K, t0 = 0, M = M, a = a, b = b)
        ## catch_unit = NA: catch in numbers (no warning for every resample)
        v <- try(withCallingHandlers(
            VPA(p, terminalF = terminalF, catch_unit = catch_unit,
                analysis_type = analysis_type, algorithm = algorithm,
                plus_group = plus_group, LW_unit = LW_unit, plot = FALSE),
            warning = function(w){
                if(grepl("did not specify catch_unit", conditionMessage(w)))
                    invokeRestart("muffleWarning")
            }), silent = TRUE)
        if(inherits(v, "try-error")){
            skip[x] <- "VPA failed"
            next
        }

        N[x] <- sum(v$annualMeanNr, na.rm = TRUE)
        B[x] <- sum(v$meanBiomassTon, na.rm = TRUE)
        FMvecVPA[[x]] <- data.frame(midLengths = p$midLengths, FM = v$FM_calc)
    }

    bootSkipWarning(skip, "VPA")

    bootAddPars(boot, list(N = N, B = B), CI = CI,
                misc = list(FMvecVPA = FMvecVPA))
}
