## M_empirical: compare with the published equations ------------------

test_that("M_empirical Pauly (1980) length equation", {
    Linf <- 80; K <- 0.5; temp <- 25
    expected <- 10^(-0.0066 - 0.279 * log10(Linf) + 0.6543 * log10(K) +
                    0.4634 * log10(temp))
    res <- M_empirical(Linf = Linf, K_l = K, temp = temp, method = "Pauly_Linf")
    expect_equal(res[1, 1], round(expected, 3))
})

test_that("M_empirical Then (2015) tmax equation", {
    res <- M_empirical(tmax = 20, method = "Then_tmax")
    expect_equal(res[1, 1], round(4.899 * 20^-0.916, 3))
})

test_that("M_empirical complains about missing inputs", {
    expect_error(M_empirical(tmax = 10, method = "AlversonCarney"), "K_l")
})

## catchCurve: synthetic data with known Z ----------------------------

test_that("age-based catch curve recovers Z exactly", {
    age <- 0:10
    Z <- 0.8
    catch <- 1000 * exp(-Z * age) * (1 - exp(-Z)) * 0.5
    res <- suppressWarnings(  # perfect fit warning from summary.lm
        catchCurve(list(age = age, catch = catch), reg_int = c(2, 10), plot = FALSE))
    expect_equal(res$Z, Z)
})

## length classes of 2 cm with numbers dying in each class for known Z
make_length_catch <- function(Z = 0.8, Linf = 60, K = 0.4) {
    ml <- seq(5, 55, 2)
    t_low <- -log(1 - (ml - 1) / Linf) / K
    t_up <- -log(1 - (ml + 1) / Linf) / K
    list(midLengths = ml, Linf = Linf, K = K,
         catch = 1e5 * (exp(-Z * t_low) - exp(-Z * t_up)))
}

test_that("length-converted catch curve recovers Z", {
    res <- catchCurve(make_length_catch(), reg_int = c(3, 20), plot = FALSE)
    expect_equal(res$Z, 0.8, tolerance = 0.01)
    expect_true(res$confidenceInt[1] <= res$Z && res$Z <= res$confidenceInt[2])
})

test_that("cumulative catch curve approximately recovers Z", {
    res <- catchCurve(make_length_catch(), cumulative = TRUE, reg_int = c(3, 20),
                      plot = FALSE)
    ## catches above the largest class are missing, hence the small bias
    expect_equal(res$Z, 0.8, tolerance = 0.05)
})

test_that("catchCurve with selection ogive on synLFQ3 (regression values)", {
    res <- catchCurve(synLFQ3, reg_int = c(9, 21), calc_ogive = TRUE, plot = FALSE)
    expect_equal(res$Z, 0.9847515, tolerance = 1e-6)
    expect_equal(res$L50, 14.60952, tolerance = 1e-6)
    expect_equal(res$L75, 17.56116, tolerance = 1e-6)
    expect_lt(res$L50, res$L75)
})

test_that("catchCurve needs growth parameters for length data", {
    lfq <- make_length_catch()
    lfq$Linf <- NULL
    expect_error(catchCurve(lfq, reg_int = c(3, 20), plot = FALSE), "Linf and K")
    ## the error explains where to put them (issue #48)
    expect_error(catchCurve(lfq, reg_int = c(3, 20), plot = FALSE), "lfq\\$par")
    lfq$Linf <- NA
    expect_error(catchCurve(lfq, reg_int = c(3, 20), plot = FALSE), "Linf and K")
})

test_that("catchCurve checks reg_int (point numbers, issue #48)", {
    lfq <- make_length_catch()
    n <- length(lfq$midLengths)
    ## lengths instead of point numbers
    expect_error(catchCurve(lfq, reg_int = c(10.5, 14.5), plot = FALSE),
                 "not lengths or ages")
    expect_error(catchCurve(lfq, reg_int = c(3, n + 1), plot = FALSE),
                 paste0("between 1 and ", n))
    expect_error(catchCurve(lfq, reg_int = c(3, 3), plot = FALSE), "first and last")
    expect_error(catchCurve(lfq, reg_int = list(c(3, 10), c(0, 5)), reg_num = 2,
                            plot = FALSE), "first and last")
    ## the order of the two points does not matter
    expect_equal(catchCurve(lfq, reg_int = c(20, 3), plot = FALSE)$Z,
                 catchCurve(lfq, reg_int = c(3, 20), plot = FALSE)$Z)
    expect_equal(catchCurve(lfq, reg_int = c(20, 3), plot = FALSE)$reg_int, c(3, 20))
})

test_that("catchCurve checks catch_columns (issue #48)", {
    lfq <- synLFQ7
    lfq$par <- list(Linf = 123, K = 0.2, t0 = 0)
    lfq <- lfqModify(lfq, plus_group = lfq$midLengths[max(which(lfq$midLengths < 123))])
    vec <- lfqModify(lfq, vectorise_catch = TRUE)
    expect_error(catchCurve(vec, catch_columns = c(1, 2), reg_int = c(9, 21),
                            plot = FALSE), "catch is a vector")
    expect_error(catchCurve(lfq, catch_columns = 0, reg_int = c(9, 21),
                            plot = FALSE), "column numbers")
    expect_error(catchCurve(lfq, catch_columns = ncol(lfq$catch) + 1,
                            reg_int = c(9, 21), plot = FALSE), "column numbers")
    ## several columns are summed
    one <- lfq
    one$catch <- rowSums(lfq$catch[, 1:2])
    expect_equal(catchCurve(lfq, catch_columns = 1:2, reg_int = c(9, 21),
                            plot = FALSE)$Z,
                 catchCurve(one, reg_int = c(9, 21), plot = FALSE)$Z)
})

## Z_BevertonHolt -----------------------------------------------------

test_that("Z_BevertonHolt (length) follows Z = K (Linf - Lmean) / (Lmean - L')", {
    res <- Z_BevertonHolt(synLFQ2, catch_columns = 2, Lprime_tprime = 47.5)
    expect_equal(res$Z, synLFQ2$K * (synLFQ2$Linf - res$Lmean) /
                        (res$Lmean - res$Lprime))
    ## L' is the lower bound of the class whose mid length is closest to 47.5
    interval <- diff(synLFQ2$midLengths)[1] / 2
    expect_equal(res$Lprime, 47.5 - interval)
    expect_equal(res$Z, 0.8872215, tolerance = 1e-6)
})

test_that("Z_BevertonHolt needs catch columns for a catch matrix", {
    expect_error(Z_BevertonHolt(synLFQ2, Lprime_tprime = 47.5), "column")
})

## Z_CPUE: synthetic CPUE with known Z --------------------------------

test_that("Z_CPUE standard method recovers Z for every cohort pair", {
    age <- 1:6
    cpue <- list(cohort = paste0("c", age), age = age, CPUE = 500 * exp(-0.6 * age))
    zmat <- as.matrix(Z_CPUE(cpue)$Z_mat[, -1])
    expect_true(all(zmat[!is.na(zmat)] == 0.6))
})

test_that("Z_CPUE Heincke and Robson-Chapman follow their equations", {
    age <- 1:6
    N <- 500 * exp(-0.6 * age)
    cpue <- list(cohort = paste0("c", age), age = age, CPUE = N)
    expect_equal(Z_CPUE(cpue, method = "Heincke")$Z,
                 log(sum(N) / sum(N[-1])))
    Tsum <- sum((seq_along(N) - 1) * N)
    expect_equal(Z_CPUE(cpue, method = "RobsonChapman")$Z,
                 log((sum(N) + Tsum - 1) / Tsum))
})

## Powell-Wetherall ---------------------------------------------------

test_that("Powell-Wetherall estimates Linf of synLFQ3 (true Linf = 50)", {
    res <- quiet(powell_wetherall(synLFQ3, catch_columns = 1, reg_int = c(10, 28)))
    expect_equal(res$Linf_est, 50, tolerance = 0.05)
    expect_equal(res$Linf_est, 49.34458, tolerance = 1e-6)
    expect_true(res$confidenceInt_Linf[1] < res$Linf_est &&
                res$Linf_est < res$confidenceInt_Linf[2])
})

test_that("catchCurve uses seasonal growth parameters in par", {
    lfq <- list(midLengths = synLFQ7$midLengths, catch = rowSums(synLFQ7$catch))
    gp <- list(Linf = 130, K = 0.2, t0 = 0, C = 0.75, ts = 0.5)
    top <- suppressWarnings(catchCurve(c(lfq, gp), reg_int = c(20, 40), plot = FALSE))
    inpar <- suppressWarnings(catchCurve(c(lfq, list(par = gp)), reg_int = c(20, 40),
                                         plot = FALSE))
    nonseasonal <- suppressWarnings(catchCurve(c(lfq, list(par = gp[1:3])),
                                               reg_int = c(20, 40), plot = FALSE))
    expect_equal(inpar$Z, top$Z)
    expect_false(isTRUE(all.equal(inpar$Z, nonseasonal$Z)))
})

test_that("catchCurve auto ends at the last point with ln(C/dt) > 0.2", {
    ## empty length class in the descending limb: excluded by y > 0.2, which
    ## used to shift the end of the regression line by one point
    lfq <- make_length_catch()
    lfq$catch[15] <- 0
    res <- suppressWarnings(catchCurve(lfq, auto = TRUE, plot = FALSE))
    ok <- which(res$lnC_dt > 0.2)
    expect_equal(res$reg_int[2], max(ok))
    expect_equal(res$reg_int[1], which.max(res$lnC_dt) + 1)
})
