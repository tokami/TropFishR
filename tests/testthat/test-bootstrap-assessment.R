## Bootstrapped assessment: M_empirical, MstochThen, catchCurve, VPA and
## predict_mod with boot, based on the first resamples of alba_boot

subset_boot <- function(boot, i) {
    boot$bootRaw <- boot$bootRaw[i, ]
    boot$seed <- boot$seed[i]
    boot$Rn_max <- boot$Rn_max[i]
    boot
}
boot15 <- subset_boot(alba_boot, 1:15)
alba_lw <- alba
alba_lw$a <- 0.0001
alba_lw$b <- 3
alba_lw$Lmat <- 5
alba_lw$wmat <- 1

## assessment chain, computed once and re-used
chain_cache <- new.env()
boot_chain <- function() {
    if (is.null(chain_cache$cc)) {
        m <- M_empirical(method = "Then_growth", boot = boot15)
        ## Then's M is larger than Z for alba (FM < 0); use a user-defined
        ## column with lower M for the rest of the assessment
        m$bootRaw$M_low <- m$bootRaw$M_Then / 3
        chain_cache$m <- m
        chain_cache$cc <- catchCurve(alba_lw, boot = m, calc_ogive = TRUE,
                                     natMort = "M_low", yearCombine = TRUE)
    }
    chain_cache
}

## resample x prepared as in the bootstrapped functions
transform_lr <- function(lfq) { lfq$Lr <- 2; lfq }

resample_x <- function(boot, x) {
    set.seed(boot$seed[x])
    lfqb <- lfqResample(alba_lw)
    lfqb <- lfqModify(lfqb, vectorise_catch = TRUE)
    lfqb$catch <- rowSums(lfqb$catch)
    lfqb
}

## natural mortality --------------------------------------------------

test_that("M_empirical with boot estimates M per resample", {
    res <- M_empirical(method = c("Then_growth", "Pauly_Linf"), temp = 20, boot = boot15)
    expect_s3_class(res, "lfqBoot")
    for (x in c(1, 7)) {
        expect_equal(res$bootRaw$M_Then[x],
                     M_empirical(Linf = boot15$bootRaw$Linf[x], K_l = boot15$bootRaw$K[x],
                                 method = "Then_growth")[1, 1])
        expect_equal(res$bootRaw$M_Pauly[x],
                     M_empirical(Linf = boot15$bootRaw$Linf[x], K_l = boot15$bootRaw$K[x],
                                 temp = 20, method = "Pauly_Linf")[1, 1])
    }
    ## summaries for both methods, everything else kept
    expect_equal(names(res$maxDen), c(names(boot15$maxDen), "M_Then", "M_Pauly"))
    expect_equal(colnames(res$CI), names(res$maxDen))
    expect_equal(res$median[["M_Then"]], median(res$bootRaw$M_Then))
    expect_equal(res$seed, boot15$seed)
    expect_equal(res$Rn_max, boot15$Rn_max)
    expect_equal(res$multiCI, boot15$multiCI)
})

test_that("M_empirical with boot rejects unsupported input", {
    expect_error(M_empirical(method = "Hoenig", boot = boot15), "Then_growth")
    expect_error(M_empirical(method = "Pauly_Linf", boot = boot15), "temp")
    expect_error(M_empirical(method = "Then_growth", boot = list()), "lfqBoot")
})

test_that("re-running a step replaces its columns", {
    res <- M_empirical(method = "Then_growth", boot = boot15)
    res2 <- M_empirical(method = "Then_growth", boot = res)
    expect_equal(sum(names(res2$bootRaw) == "M_Then"), 1)
    expect_equal(sum(colnames(res2$CI) == "M_Then"), 1)
})

test_that("MstochThen is reproducible and close to the deterministic estimate", {
    ## seeds of the resamples by default, independent of the global seed
    set.seed(1); a <- MstochThen(boot15)
    set.seed(2); b <- MstochThen(boot15)
    expect_identical(a$bootRaw$M_Then, b$bootRaw$M_Then)
    expect_true(all(is.finite(a$bootRaw$M_Then)))
    c1 <- MstochThen(boot15, seed = 10)
    expect_identical(c1$bootRaw$M_Then, MstochThen(boot15, seed = 10)$bootRaw$M_Then)
    expect_false(identical(c1$bootRaw$M_Then, a$bootRaw$M_Then))
    det <- M_empirical(method = "Then_growth", boot = boot15)
    expect_equal(median(a$bootRaw$M_Then), median(det$bootRaw$M_Then), tolerance = 0.3)
})

test_that("MstochThen with resid includes the residual error", {
    a <- MstochThen(alba_boot)
    r <- MstochThen(alba_boot, resid = TRUE)
    expect_true(all(r$bootRaw$M_Then > 0))
    expect_gt(diff(r$CI[, "M_Then"]), 2 * diff(a$CI[, "M_Then"]))
})

test_that("the bootstrapped coefficients match Then et al. (2015)", {
    tb <- TropFishR:::thenBoot
    ## fit to all data: published coefficients
    expect_equal(unname(tb$coefFull), c(4.118, 0.73, -0.33), tolerance = 0.01)
    expect_equal(dim(tb$coef), c(2000, 3))
    expect_true(all(is.finite(tb$coef)))
    ## published coefficients within the bootstrap distribution
    med <- apply(tb$coef, 2, median)
    expect_equal(unname(med), c(4.118, 0.73, -0.33), tolerance = 0.1)
    expect_length(tb$logres, tb$nobs)
})

test_that("MstochThen uses one bootstrapped coefficient set per resample", {
    res <- MstochThen(boot15)
    tb <- TropFishR:::thenBoot
    x <- 4
    set.seed(boot15$seed[x])
    co <- tb$coef[sample(nrow(tb$coef), 1), ]
    expect_equal(res$bootRaw$M_Then[x],
                 unname(co["a"] * boot15$bootRaw$K[x]^co["b"] * boot15$bootRaw$Linf[x]^co["c"]))
})

## catch curve --------------------------------------------------------

test_that("ccAutoRegInt chooses the descending limb", {
    x <- 1:12
    y <- c(1, 4, 6, 5, 4.2, 3.3, 2.5, 1.6, 0.9, 0.6, 0.3, NA)
    ## from the point after the maximum to the last point with y > 0.5
    expect_equal(TropFishR:::ccAutoRegInt(x, y), c(4, 10))
    expect_equal(TropFishR:::ccAutoRegInt(x, y, agerange = c(5, 8)), c(5, 8))
    ## last point with much flatter slope is removed with robustReg
    y2 <- c(1, 4, 6, 5, 4, 3, 2, 1.95, 1.9, NA, NA, NA)
    expect_equal(TropFishR:::ccAutoRegInt(x, y2), c(4, 9))
    expect_equal(TropFishR:::ccAutoRegInt(x, y2, robustReg = TRUE), c(4, 8))
    ## not enough points
    expect_equal(TropFishR:::ccAutoRegInt(1:2, c(3, 2)), c(NA, NA))
})

test_that("catchCurve with boot equals the catch curve of each resample", {
    ch <- boot_chain()
    res <- ch$cc
    expect_true(all(c("Z", "FM", "L50", "L75") %in% names(res$bootRaw)))
    expect_equal(res$bootRaw$FM, res$bootRaw$Z - res$bootRaw$M_low)
    for (x in c(2, 9)) {
        lfqx <- c(resample_x(boot15, x),
                  list(Linf = boot15$bootRaw$Linf[x], K = boot15$bootRaw$K[x], t0 = 0))
        cc0 <- suppressWarnings(catchCurve(lfqx, auto = TRUE, plot = FALSE))
        ri <- TropFishR:::ccAutoRegInt(cc0$t_midL, cc0$lnC_dt)
        cc <- suppressWarnings(catchCurve(lfqx, reg_int = ri, calc_ogive = TRUE, plot = FALSE))
        expect_equal(res$bootRaw$Z[x], cc$Z)
        expect_equal(res$bootRaw$L50[x], cc$L50)
        expect_equal(res$misc$intsCC[x], unname(coef(cc$linear_mod)[1]))
    }
    expect_length(res$misc$datcc, 15)
    expect_true(all(res$bootRaw$L50 < res$bootRaw$L75, na.rm = TRUE))
})

test_that("catchCurve with boot assumes knife-edge selection without ascending limb", {
    ## regression line from the first length class: no length classes below
    expect_message(res <- suppressWarnings(
        catchCurve(alba_lw, boot = boot15, calc_ogive = TRUE, yearCombine = TRUE,
                   reg_int = c(1, 8))), "knife-edge")
    ke <- res$misc$knifeEdge
    expect_true(any(ke))
    lower1 <- alba_lw$midLengths[1] - diff(alba_lw$midLengths)[1] / 2
    expect_equal(res$bootRaw$L50[ke], rep(lower1, sum(ke)))
    expect_equal(res$bootRaw$L75[ke], res$bootRaw$L50[ke])
})

test_that("catchCurve with boot requires one year or yearCombine", {
    expect_error(catchCurve(alba_lw, boot = boot15), "several years")
    expect_error(catchCurve(alba_lw, boot = boot15, natMort = "M_Then", yearCombine = TRUE),
                 "M_Then")
})

test_that("catchCurve with binSize re-bins the length classes", {
    lfq <- list(midLengths = synLFQ7$midLengths, catch = synLFQ7$catch,
                Linf = 130, K = 0.2, t0 = 0)
    rb <- TropFishR:::rebinCatch(lfq$midLengths, lfq$catch, 4)
    ref <- lfqModify(synLFQ7, bin_size = 4)
    expect_equal(rb$midLengths, ref$midLengths)
    expect_equal(unname(rb$catch), unname(ref$catch))
    expect_equal(sum(rb$catch), sum(lfq$catch))
    res <- suppressWarnings(catchCurve(lfq, catch_columns = 1:3, binSize = 4,
                                       auto = TRUE, plot = FALSE))
    lfq2 <- lfq
    lfq2$midLengths <- rb$midLengths
    lfq2$catch <- rb$catch
    expect_equal(res$Z, suppressWarnings(catchCurve(lfq2, catch_columns = 1:3,
                                                    auto = TRUE, plot = FALSE))$Z)
    ## vector catch
    expect_equal(TropFishR:::rebinCatch(1:6 + 0.5, rep(1, 6), 2)$catch, c(1, 2, 2, 1))
    expect_error(TropFishR:::rebinCatch(seq(1, 11, 2), rep(1, 6), 1), "smaller")
})

test_that("catchCurve handles infinite values of ln(C/dt)", {
    ## seasonal growth: length classes above Linf used to get the same age
    ## (dt = 0, infinite ln(C/dt)), which made lm fail
    lfq <- list(midLengths = synLFQ7$midLengths, catch = rowSums(synLFQ7$catch),
                Linf = 110, K = 0.25, t0 = 0, C = 0.3, ts = 0)
    res <- suppressWarnings(catchCurve(lfq, auto = TRUE, plot = FALSE))
    expect_true(is.finite(res$Z))
    expect_true(all(is.na(res$lnC_dt[lfq$midLengths > 110])))
})

## VPA and prediction models -----------------------------------------

test_that("VPA with boot equals the VPA of each resample", {
    ch <- boot_chain()
    res <- VPA(alba_lw, boot = ch$cc, natMort = "M_low", yearCombine = TRUE)
    expect_true(all(c("N", "B") %in% names(res$bootRaw)))
    x <- which(is.finite(res$bootRaw$N))[1]
    lfqx <- resample_x(boot15, x)
    Linf <- boot15$bootRaw$Linf[x]
    upper <- lfqx$midLengths + diff(lfqx$midLengths)[1] / 2
    pg <- Linf < max(upper)
    if (pg) lfqx <- lfqModify(lfqx, plus_group = lfqx$midLengths[which.min(abs(upper - floor(Linf)))])
    v <- suppressWarnings(VPA(list(midLengths = lfqx$midLengths, catch = lfqx$catch,
                                   Linf = Linf, K = boot15$bootRaw$K[x], t0 = 0,
                                   M = ch$cc$bootRaw$M_low[x], a = 0.0001, b = 3),
                              terminalF = ch$cc$bootRaw$FM[x], plus_group = pg, plot = FALSE))
    expect_equal(res$bootRaw$N[x], sum(v$annualMeanNr))
    expect_equal(res$bootRaw$B[x], sum(v$meanBiomassTon))
    expect_equal(res$misc$FMvecVPA[[x]]$FM, v$FM_calc)
    expect_equal(res$misc$FMvecVPA[[x]]$midLengths, lfqx$midLengths)
    ## resamples with FM <= 0 or NA are skipped
    skipped <- !is.finite(ch$cc$bootRaw$FM) | ch$cc$bootRaw$FM <= 0
    expect_true(all(is.na(res$bootRaw$N[skipped])))
})

test_that("VPA with boot and landings raises each resample", {
    ch <- boot_chain()
    landings <- c(100, 200, 150, 50, 100, 250, 150)  # kg, one per sample
    res <- VPA(alba_lw, boot = ch$cc, natMort = "M_low", yearCombine = TRUE,
               landings = landings)
    x <- which(is.finite(res$bootRaw$N))[1]
    ## resample x raised with the weight of the resampled fish
    set.seed(boot15$seed[x])
    lfqx <- lfqRaise(lfqResample(alba_lw), landings = landings)
    lfqx <- lfqModify(lfqx, vectorise_catch = TRUE)
    lfqx$catch <- rowSums(lfqx$catch)
    ## the raised catch of the resample weighs as much as the landings
    expect_equal(sum(lfqx$catch * 0.0001 * lfqx$midLengths^3), sum(landings))
    Linf <- boot15$bootRaw$Linf[x]
    upper <- lfqx$midLengths + diff(lfqx$midLengths)[1] / 2
    pg <- Linf < max(upper)
    if (pg) lfqx <- lfqModify(lfqx, plus_group = lfqx$midLengths[which.min(abs(upper - floor(Linf)))])
    v <- suppressWarnings(VPA(list(midLengths = lfqx$midLengths, catch = lfqx$catch,
                                   Linf = Linf, K = boot15$bootRaw$K[x], t0 = 0,
                                   M = ch$cc$bootRaw$M_low[x], a = 0.0001, b = 3),
                              terminalF = ch$cc$bootRaw$FM[x], plus_group = pg, plot = FALSE))
    expect_equal(res$bootRaw$N[x], sum(v$annualMeanNr))
    expect_equal(res$misc$FMvecVPA[[x]]$FM, v$FM_calc)

    expect_error(VPA(alba_lw, boot = ch$cc, natMort = "M_low", yearCombine = TRUE,
                     landings = landings, catch_unit = "'000"), "catch_unit = NA")
    raised <- lfqRaise(alba_lw, landings = landings)
    expect_error(VPA(raised, boot = ch$cc, natMort = "M_low", yearCombine = TRUE),
                 "raised")
})

test_that("predict_mod with boot equals the model of each resample", {
    ch <- boot_chain()
    FM_change <- seq(0, 10, 0.2)
    res <- predict_mod(alba_lw, type = "ThompBell", FM_change = FM_change,
                       boot = ch$cc, natMort = "M_low")
    expect_true(all(c("F01", "Fmax", "F05", "F04", "SPR", "FF01", "FFmax") %in% names(res$bootRaw)))
    br <- ch$cc$bootRaw
    x <- which(is.finite(res$bootRaw$Fmax))[1]
    p <- list(midLengths = alba_lw$midLengths, Linf = br$Linf[x], K = br$K[x], t0 = 0,
              M = br$M_low[x], a = 0.0001, b = 3, Lmat = 5, wmat = 1, FM = br$FM[x])
    direct <- function(p){
        r <- NULL
        capture.output(r <- suppressWarnings(predict_mod(
            p, type = "ThompBell", FM_change = FM_change,
            s_list = list(selecType = "trawl_ogive", L50 = br$L50[x], L75 = br$L75[x]),
            hide.progressbar = TRUE)))
        r
    }
    r <- direct(p)
    expect_equal(res$bootRaw$Fmax[x], r$df_Es$Fmax)
    expect_equal(res$bootRaw$F01[x], r$df_Es$F01)
    expect_equal(res$bootRaw$SPR[x], approx(FM_change, r$SPR, xout = br$FM[x])$y)
    expect_equal(res$bootRaw$FFmax, br$FM / res$bootRaw$Fmax)
    expect_equal(res$misc$totY[[x]], r$totY)
    expect_equal(res$misc$FM_change, FM_change)
    ## LW_unit = "g": the same reference points, yield in kg instead of g
    res_g <- predict_mod(alba_lw, type = "ThompBell", FM_change = FM_change,
                         boot = ch$cc, natMort = "M_low", LW_unit = "g")
    expect_equal(res_g$bootRaw$Fmax, res$bootRaw$Fmax)
    expect_equal(res_g$misc$totY[[x]], res$misc$totY[[x]] / 1000)
})

test_that("predict_mod with boot uses F at length from VPA", {
    ch <- boot_chain()
    v <- VPA(alba_lw, boot = ch$cc, natMort = "M_low", yearCombine = TRUE)
    v$bootRaw$L50 <- NULL
    v$bootRaw$L75 <- NULL
    expect_error(predict_mod(alba_lw, type = "ThompBell", FM_change = seq(0, 10, 0.2),
                             boot = v, natMort = "M_low"), "bootSel")
    res <- predict_mod(alba_lw, type = "ThompBell", FM_change = seq(0, 10, 0.2),
                       boot = v, natMort = "M_low", bootSel = "VPA")
    ok <- is.finite(v$bootRaw$N)
    expect_true(all(is.finite(res$bootRaw$Fmax[ok])))
    ## current F: maximum F at length (FM_change refers to it)
    Fapical <- vapply(v$misc$FMvecVPA, function(d)
        if(is.null(d)) NA_real_ else max(d$FM[is.finite(d$FM)]), numeric(1))
    expect_equal(res$bootRaw$FFmax, Fapical / res$bootRaw$Fmax)
    expect_error(predict_mod(transform_lr(alba_lw), type = "ypr", FM_change = seq(0, 2, 0.1),
                             boot = v, natMort = "M_low", bootSel = "VPA"), "ThompBell")
})

test_that("predict_mod with boot: Beverton and Holt yield per recruit", {
    ch <- boot_chain()
    lfq <- alba_lw
    lfq$Lr <- 2
    res <- predict_mod(lfq, type = "ypr", FM_change = seq(0, 10, 0.2),
                       boot = ch$cc, natMort = "M_low")
    expect_true(all(c("F01", "Fmax", "F05") %in% names(res$bootRaw)))
    expect_false("F04" %in% names(res$bootRaw))
    expect_true(all(res$bootRaw$F01 <= res$bootRaw$Fmax, na.rm = TRUE))
})

test_that("predict_mod with boot checks its input", {
    ch <- boot_chain()
    expect_error(predict_mod(alba_lw, type = "ThompBell", E_change = seq(0, 0.9, 0.1),
                             boot = ch$cc, natMort = "M_low"), "FM_change")
    expect_error(predict_mod(alba_lw, type = "ThompBell", FM_change = seq(0, 2, 0.1),
                             boot = ch$m, natMort = "M_low"), "selectivity")
})

## skipped resamples ---------------------------------------------------

test_that("skipped resamples give a warning and NA estimates", {
    ch <- boot_chain()
    ## Then's M exceeds Z for most resamples of alba
    expect_warning(cc <- catchCurve(alba_lw, boot = ch$m, calc_ogive = TRUE,
                                    natMort = "M_Then", yearCombine = TRUE),
                   "FM <= 0")
    nok <- sum(cc$bootRaw$FM > 0, na.rm = TRUE)
    expect_lt(nok, 15)
    expect_warning(res <- predict_mod(alba_lw, type = "ThompBell", FM_change = seq(0, 10, 0.2),
                                      boot = cc, natMort = "M_Then"),
                   "resamples were skipped")
    ## columns are kept, also if all resamples are skipped
    expect_true(all(c("F01", "Fmax", "F05", "F04", "SPR", "FF01", "FFmax", "FF05", "FF04")
                    %in% names(res$bootRaw)))
    expect_equal(sum(is.finite(res$bootRaw$Fmax)), nok)
    cc$bootRaw$FM <- -abs(cc$bootRaw$FM)
    res0 <- suppressWarnings(predict_mod(alba_lw, type = "ThompBell", FM_change = seq(0, 10, 0.2),
                                         boot = cc, natMort = "M_Then"))
    expect_true(all(is.na(res0$bootRaw$Fmax)))
    expect_true(is.na(res0$maxDen[["Fmax"]]))
    ## SPR-based estimates are omitted without maturity information
    lfq <- alba_lw
    lfq$Lmat <- NULL
    lfq$wmat <- NULL
    res1 <- predict_mod(lfq, type = "ThompBell", FM_change = seq(0, 10, 0.2),
                        boot = ch$cc, natMort = "M_low")
    expect_false(any(c("SPR", "F04", "FF04") %in% names(res1$bootRaw)))
    expect_true("F05" %in% names(res1$bootRaw))
})

test_that("bootstrapped summaries use the mode of the kernel density", {
    set.seed(2)
    x <- c(rnorm(300, 1, 0.1), rnorm(100, 3, 0.1))
    s <- TropFishR:::bootUniSummary(x)
    k <- ks::kde(x)
    expect_equal(s$maxDen, k$eval.points[which.max(k$estimate)])
    expect_equal(s$maxDen, 1, tolerance = 0.05)
    expect_equal(TropFishR:::bootUniSummary(rep(2, 5))$maxDen, 2)
    expect_true(is.na(TropFishR:::bootUniSummary(c(NA, NA))$maxDen))
})

test_that("plotBoot shows catch curves and yield per recruit curves", {
    ch <- boot_chain()
    res <- predict_mod(alba_lw, type = "ThompBell", FM_change = seq(0, 10, 0.2),
                       boot = ch$cc, natMort = "M_low")
    expect_equal(count_plot_pages(plotBoot(boot = res, pType = "cc")), 1)
    expect_equal(count_plot_pages(plotBoot(boot = res, pType = "cc", best = "median",
                                           display = "bootRaw")), 1)
    expect_equal(count_plot_pages(plotBoot(boot = res, pType = "ypr", pRef = "F01")), 1)
    expect_error(plotBoot(boot = ch$m, pType = "cc"), "catchCurve")
})

test_that("univariate_density plots selected parameters", {
    expect_equal(count_plot_pages(univariate_density(alba_boot, pars = c("Linf", "K"))), 1)
    expect_error(univariate_density(alba_boot, pars = "Z"), "Z")
})
