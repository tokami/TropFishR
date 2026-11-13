## Beverton & Holt yield per recruit ---------------------------------

threadfin <- list(Winf = 286, K = 0.37, t0 = -0.2, M = 1.1, tr = 0.4)

test_that("Beverton-Holt YPR reference points are ordered and plausible", {
    res <- predict_mod(threadfin, FM_change = seq(0, 6, 0.1),
                       tc_change = seq(0.2, 1, 0.2), type = "ypr", plot = FALSE)
    ref <- res$df_Es
    expect_equal(ref$tc, seq(0.2, 1, 0.2))
    expect_true(all(ref$F01 <= ref$Fmax))
    expect_true(all(ref$E01 > 0 & ref$Emax < 1))
    expect_equal(ref$Emax, ref$Fmax / (ref$Fmax + threadfin$M), tolerance = 1e-6)
    ## later age at first capture allows higher F
    expect_true(all(diff(ref$Fmax) >= 0))
    ## regression values
    expect_equal(ref$F01, c(0.6, 0.7, 0.9, 1.2, 1.6))
    expect_equal(ref$Fmax, c(0.7, 0.9, 1.2, 1.6, 2.3))
})

test_that("Beverton-Holt yield per recruit is zero without fishing", {
    res <- predict_mod(threadfin, FM_change = seq(0, 3, 0.1), tc_change = 0.6,
                       type = "ypr", plot = FALSE)
    run <- res$list_Lc_runs[[1]]
    expect_equal(run$Y_R[1], 0)
    expect_true(all(run$Y_R >= 0))
    ## relative biomass per recruit starts at 1 (100 %) and decreases with F
    expect_equal(run$B_R.rel[1], 1)
    expect_true(all(diff(run$B_R) < 0))
    expect_equal(run$E, run$FM / (run$FM + threadfin$M))
})

## Thompson & Bell ----------------------------------------------------

test_that("stock_sim totals are the sums over age classes", {
    res <- stock_sim(shrimps, age_unit = "month")
    expect_equal(res$totals$totC, sum(res$C, na.rm = TRUE))
    expect_equal(res$totals$totY, sum(res$Y, na.rm = TRUE))
    expect_equal(res$totals$totV, sum(res$V, na.rm = TRUE))
    ## numbers can only decrease over age
    expect_true(all(diff(res$N[!is.na(res$N)]) <= 0))
})

## the plus group should behave like an infinite chain of age classes with
## the plus group's F, M, weight and value
plus_group_reference <- function(param, k, n = 600) {
    ext <- param
    ext$age <- c(param$age[1:k], param$age[k] + seq_len(n))
    for (v in c("meanWeight", "meanValue", "FM")) {
        ext[[v]] <- c(param[[v]][1:k], rep(param[[v]][k], n))
    }
    stock_sim(ext, age_unit = "month")
}

test_that("stock_sim plus group matches an infinite chain of age classes", {
    res <- stock_sim(shrimps, age_unit = "month", plus_group = 11)
    ref <- plus_group_reference(shrimps, k = 11)
    expect_equal(res$totals$totC, ref$totals$totC)
    expect_equal(res$totals$totY, ref$totals$totY)
    expect_equal(res$totals$totV, ref$totals$totV)
    ## number-years lived, weighted by mean weight
    expect_equal(sum(res$B * res$dt, na.rm = TRUE),
                 sum(ref$B * ref$dt, na.rm = TRUE))
})

test_that("stock_sim plus group ignores age classes beyond the next one", {
    ## dt of the plus group needs the following age; anything older is unused
    cut <- function(p, n) {
        for (v in c("age", "meanWeight", "meanValue", "FM")) p[[v]] <- p[[v]][1:n]
        p
    }
    full <- stock_sim(shrimps, age_unit = "month", plus_group = 11)
    short <- stock_sim(cut(shrimps, 12), age_unit = "month", plus_group = 11)
    expect_equal(full$totals, short$totals)
})

test_that("stock_sim yield is catch times mean weight in every class", {
    res <- stock_sim(shrimps, age_unit = "month", plus_group = 11)
    expect_equal(res$Y, res$C * shrimps$meanWeight[1:11])
    expect_equal(res$V, res$Y * shrimps$meanValue[1:11])
})

test_that("Thompson-Bell yield increases and biomass decreases with F", {
    res <- predict_mod(shrimps, FM_change = seq(0.1, 3, 0.1), type = "ThompBell",
                       age_unit = "month", plot = FALSE)
    expect_true(all(diff(res$totY) > 0))   # shrimps: no Fmax below F = 3
    expect_true(all(diff(res$meanB) < 0))
})

test_that("Thompson-Bell gives the unfished biomass at F = 0", {
    res <- predict_mod(shrimps, FM_change = seq(0, 3, 0.1), type = "ThompBell",
                       age_unit = "month", plot = FALSE)
    expect_equal(res$totY[1], 0)
    expect_gt(res$meanB[1], res$meanB[2])
    ## biomass at F = 0 from a direct simulation without fishing
    unfished <- shrimps
    unfished$FM <- shrimps$FM * 0
    expect_equal(res$meanB[1],
                 stock_sim(unfished, age_unit = "month")$totals$meanB)
})

test_that("stock_sim biomass is unchanged for F > 0 and finite for F = 0", {
    res <- stock_sim(shrimps, age_unit = "month")
    expect_equal(res$B[-length(res$B)],
                 (res$Y / (shrimps$FM * res$dt))[-length(res$B)])
    unfished <- shrimps
    unfished$FM <- shrimps$FM * 0
    B0 <- stock_sim(unfished, age_unit = "month")$B
    expect_true(all(is.finite(B0[-length(B0)])))
    expect_true(all(B0[-length(B0)] > 0))
    ## same with a plus group
    pg <- stock_sim(unfished, age_unit = "month", plus_group = 11)
    expect_true(all(is.finite(pg$B[1:11]) & pg$B[1:11] > 0))
})

test_that("Thompson-Bell F0.5 halves the biomass at F = 0", {
    res <- predict_mod(shrimps, FM_change = seq(0, 10, 0.05), type = "ThompBell",
                       age_unit = "month", plot = FALSE)
    B_F05 <- approx(res$FM_change, res$meanB, xout = res$df_Es$F05)$y
    expect_equal(B_F05 / res$meanB[1], 0.5, tolerance = 0.02)
})

test_that("Thompson-Bell reference points use E = F / (F + M)", {
    res <- predict_mod(shrimps, FM_change = seq(0.1, 5, 0.1), type = "ThompBell",
                       age_unit = "month", plot = FALSE)
    ref <- res$df_Es
    expect_true(all(c(ref$E01, ref$Emax, ref$E05) <= 1))
})

test_that("Thompson-Bell reference table has a BPR_F01 column", {
    res <- predict_mod(shrimps, FM_change = seq(0.1, 5, 0.1), type = "ThompBell",
                       age_unit = "month", plot = FALSE)
    expect_true("BPR_F01" %in% names(res$df_Es))
    expect_false(anyDuplicated(names(res$df_Es)) > 0)
})

test_that("Thompson-Bell reference points are searched over the FM_change range", {
    res <- predict_mod(shrimps, FM_change = seq(0.1, 20, 0.1), type = "ThompBell",
                       age_unit = "month", plot = FALSE)
    ## YPR peaks inside the range, so Fmax should equal the F with maximal yield
    expect_equal(res$df_Es$Fmax, res$FM_change[which.max(res$totY)], tolerance = 0.1)
})

## length-based stock_sim: plus group at Linf

lb_stock <- function(Linf = 60) {
    ml <- seq(5, 79, 2)
    list(midLengths = ml, Linf = Linf, K = 0.3, t0 = 0, M = 0.3,
         FM = rep(0.5, length(ml)), a = 0.01, b = 3)
}

test_that("length-based stock_sim catches all fish with a plus group at Linf", {
    ## constant F: all fish die eventually, the fraction F / Z by fishing
    p <- lb_stock()
    res <- stock_sim(p)
    expect_equal(res$totals$totC, 1000 * 0.5 / 0.8)
    ## the class including Linf is the plus group, classes above are NA
    lower <- p$midLengths - 1
    above <- lower >= p$Linf
    expect_true(all(is.na(res$C[above])))
    expect_true(all(is.finite(res$C[!above])))
    ## the same as with the length classes up to the one including Linf
    cut <- p
    cut$midLengths <- p$midLengths[!above]
    cut$FM <- p$FM[!above]
    expect_equal(stock_sim(cut)$totals, res$totals)
    ## Linf at a class limit
    expect_equal(stock_sim(lb_stock(Linf = 62))$totals$totC, 1000 * 0.5 / 0.8)
})

test_that("predict_mod LW_unit = 'g' gives yield in kg and the same reference points", {
    p <- lb_stock()
    p$FM <- 0.5
    sl <- list(selecType = "trawl_ogive", L50 = 20, L75 = 24)
    run <- function(...) {
        r <- NULL
        capture.output(r <- predict_mod(p, type = "ThompBell", FM_change = seq(0, 3, 0.1),
                                        s_list = sl, hide.progressbar = TRUE, ...))
        r
    }
    kg <- run()
    g <- run(LW_unit = "g")
    expect_equal(g$totY, kg$totY / 1000)
    expect_equal(g$df_Es$Fmax, kg$df_Es$Fmax)
    expect_error(run(LW_unit = "t"), "arg")
})
