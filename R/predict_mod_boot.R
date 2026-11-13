## Bootstrapped prediction models (internal, called by predict_mod with boot)

predict_mod_boot <- function(param, boot, type, FM_change = NA,
                             stock_size_1 = NA, age_unit = "year",
                             plus_group = NA, Lmin = NA, Lincr = NA,
                             natMort = NULL, bootSel = "catchCurve",
                             LW_unit = "kg",
                             binSize = NA, CI = 95){

    if(is.null(natMort))
        stop("Please provide the name of the column with natural mortality estimates in boot$bootRaw (natMort).")
    checkBoot(boot, c("Linf", "K", natMort))
    if(!(type %in% c("ThompBell", "ypr")))
        stop("type has to be either 'ThompBell' or 'ypr'.")
    if(!("midLengths" %in% names(param)))
        stop("The bootstrapped prediction models require length-frequency data (midLengths).")
    if(length(FM_change) == 1 && is.na(FM_change[1]))
        stop("Please provide FM_change (E_change is not available with boot, because the corresponding fishing mortalities would differ between resamples).")

    getPar <- function(name){
        if(!is.null(param$par[[name]])) param$par[[name]] else param[[name]]
    }
    a <- getPar("a"); b <- getPar("b")
    if(type == "ThompBell" && (is.null(a) || is.null(b)))
        stop("The Thompson and Bell model requires information about the length-weight relationship. Please provide 'a' and 'b' estimates in param.")
    Lmat <- getPar("Lmat"); wmat <- getPar("wmat")
    Lr <- getPar("Lr"); tr <- getPar("tr")
    if(type == "ypr" && is.null(Lr) && is.null(tr))
        stop("The yield per recruit model requires the length or age at recruitment ('Lr' or 'tr' in param).")

    bootRaw <- boot$bootRaw
    bootSel <- match.arg(bootSel, c("catchCurve", "VPA"))
    useSel <- bootSel == "catchCurve"
    if(useSel && !all(c("L50", "L75", "FM") %in% names(bootRaw)))
        stop("bootSel = 'catchCurve' requires selectivity and fishing mortality estimates (L50, L75, FM): run catchCurve with boot, calc_ogive = TRUE and natMort first.")
    if(!useSel){
        if(type == "ypr")
            stop("bootSel = 'VPA' is only available for type = 'ThompBell'.")
        if(is.null(boot$misc$FMvecVPA))
            stop("bootSel = 'VPA' requires the fishing mortality at length of the VPA: run VPA with boot first.")
    }

    ## length classes (same for all resamples unless F at length from VPA)
    midLengths <- param$midLengths
    if(!is.na(binSize)) midLengths <- lfqModify(param, bin_size = binSize)$midLengths

    n <- nrow(bootRaw)
    refs <- c("F01", "Fmax", "F05", "F04")
    res <- matrix(NA_real_, nrow = n, ncol = length(refs) + 1,
                  dimnames = list(NULL, c("SPR", refs)))
    totY <- meanB <- vector("list", n)
    FMcurr <- rep(NA_real_, n)
    FM_grid <- NULL
    skip <- rep(NA_character_, n)

    for(x in seq_len(n)){

        M <- bootRaw[[natMort]][x]
        FM <- if("FM" %in% names(bootRaw)) bootRaw$FM[x] else NA
        growth <- bootGrowthPar(boot, x)
        if(!is.finite(M)){
            skip[x] <- "M not available"
            next
        }

        p <- list(midLengths = midLengths, Linf = growth$Linf, K = growth$K,
                  t0 = 0, M = M, a = a, b = b)
        if(!is.null(Lmat)) p$Lmat <- Lmat
        if(!is.null(wmat)) p$wmat <- wmat
        s_list <- NA

        ## selectivity: selection ogive of the catch curve, or F at length of
        ## the VPA
        if(useSel){
            L50 <- bootRaw$L50[x]; L75 <- bootRaw$L75[x]
            if(!is.finite(FM) || FM <= 0){
                skip[x] <- "FM not available or <= 0"
                next
            }
            if(!is.finite(L50) || !is.finite(L75) || L50 < 0 ||
               L75 < L50 || L75 >= growth$Linf){
                skip[x] <- "selectivity not available or invalid"
                next
            }
            s_list <- if(L50 == L75){
                list(selecType = "knife_edge", L50 = L50)
            }else{
                list(selecType = "trawl_ogive", L50 = L50, L75 = L75)
            }
            p$FM <- FM
            FMcurr[x] <- FM
            if(type == "ypr") p$Lc <- L50
        }else{
            vpa <- boot$misc$FMvecVPA[[x]]
            if(is.null(vpa) || all(!is.finite(vpa$FM))){
                skip[x] <- "VPA estimates not available"
                next
            }
            ## (VPA with boot already uses a plus group at Linf)
            p$midLengths <- vpa$midLengths
            p$FM <- vpa$FM
            ## FM_change refers to the maximum fishing mortality at length
            FMcurr[x] <- max(vpa$FM[is.finite(vpa$FM)])
        }
        if(type == "ypr"){
            if(!is.null(Lr)) p$Lr <- Lr
            if(!is.null(tr)) p$tr <- tr
        }

        r <- NULL
        out <- utils::capture.output(
            r <- try(suppressWarnings(
                predict_mod(p, type = type, FM_change = FM_change,
                            s_list = s_list, stock_size_1 = stock_size_1,
                            age_unit = age_unit, plus_group = plus_group,
                            Lmin = Lmin, Lincr = Lincr, plot = FALSE,
                            LW_unit = LW_unit,
                            hide.progressbar = TRUE)), silent = TRUE))
        if(inherits(r, "try-error")){
            skip[x] <- "prediction model failed"
            next
        }

        Es <- r$df_Es
        for(ref in intersect(refs, names(Es))) res[x, ref] <- Es[[ref]][1]
        if(type == "ThompBell"){
            FM_grid <- r$FM_change
            totY[[x]] <- r$totY
            meanB[[x]] <- r$meanB
            ## SPR at the current fishing mortality
            if(!is.null(r$SPR) && any(is.finite(r$SPR))){
                res[x, "SPR"] <- stats::approx(r$FM_change, r$SPR,
                                               xout = FMcurr[x])$y
            }
        }else{
            run <- r$list_Lc_runs[[1]]
            FM_grid <- run$FM
            totY[[x]] <- run$Y_R
            meanB[[x]] <- run$B_R
        }
    }

    bootSkipWarning(skip, "predict_mod")

    new <- as.data.frame(res)
    ## current fishing mortality relative to the reference points
    for(ref in refs) new[[paste0("F", ref)]] <- FMcurr / new[[ref]]
    ## omit SPR-based estimates if they cannot be estimated (Beverton and Holt
    ## model, or without Lmat and wmat); columns of resamples that were
    ## skipped or failed are kept as NA
    drop <- character(0)
    if(type == "ypr" || is.null(Lmat) || is.null(wmat))
        drop <- c("SPR", "F04", "FF04")
    new <- new[, setdiff(names(new), drop), drop = FALSE]

    bootAddPars(boot, new, CI = CI,
                misc = list(FM_change = FM_grid, totY = totY, meanB = meanB,
                            predType = type))
}
