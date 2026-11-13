## ELEFAN runs are slow (10-30 s each), so they are skipped on CRAN.
## synLFQ5 was simulated with Linf = 80, K = 0.5, t_anchor = 0.25;
## synLFQ4 with the same parameters plus seasonality (C = 0.75, ts = 0.5).

test_that("ELEFAN K-scan picks the true parameters from a small grid", {
    skip_on_cran()
    res <- quiet(ELEFAN(synLFQ5, Linf_range = c(76, 80, 84),
                        K_range = c(0.3, 0.5, 0.7), MA = 11,
                        hide.progressbar = TRUE))
    expect_s3_class(res, "lfq")
    expect_equal(res$par$Linf, 80)
    expect_equal(res$par$K, 0.5)
    expect_true(res$par$t_anchor >= 0 && res$par$t_anchor <= 1)
    ## score matrix: K in rows, Linf in columns; Rn_max is its maximum
    expect_equal(dim(res$score_mat), c(3, 3))
    expect_equal(res$Rn_max, max(res$score_mat))
    expect_equal(res$par$phiL, log10(0.5) + 2 * log10(80))
})

test_that("ELEFAN with fixed Linf only scans K", {
    skip_on_cran()
    res <- quiet(ELEFAN(synLFQ5, Linf_fix = 80, K_range = seq(0.3, 0.7, 0.1),
                        MA = 11, hide.progressbar = TRUE))
    expect_equal(res$par$Linf, 80)
    expect_equal(res$par$K, 0.5)
})

test_that("ELEFAN graphs are controlled by plot.score and plot", {
    skip_on_cran()
    ## response surface analysis (Linf and K) and K-scan (fixed Linf)
    rsa <- function(...) ELEFAN(synLFQ5, Linf_range = c(78, 80), K_range = 0.5,
                                MA = 11, hide.progressbar = TRUE, ...)
    kscan <- function(...) ELEFAN(synLFQ5, Linf_fix = 80, K_range = c(0.4, 0.5),
                                  MA = 11, hide.progressbar = TRUE, ...)
    for (fit in list(rsa, kscan)) {
        ## default: score graph only
        expect_equal(count_plot_pages(fit()), 1)
        expect_equal(count_plot_pages(fit(plot.score = FALSE)), 0)
        ## score graph plus restructured data with growth curves
        expect_equal(count_plot_pages(fit(plot = TRUE)), 2)
        expect_equal(count_plot_pages(fit(plot = TRUE, plot.score = FALSE)), 1)
    }
})

test_that("score graph of ELEFAN can be annotated (as in the ELEFAN vignette)", {
    skip_on_cran()
    pages <- count_plot_pages({
        fit <- ELEFAN(synLFQ5, Linf_range = c(78, 80), K_range = c(0.4, 0.5),
                      MA = 11, hide.progressbar = TRUE)
        points(fit$par$Linf, fit$par$K, pch = "*")
    })
    expect_equal(pages, 1)
})

ga_bounds <- list(low = list(Linf = 70, K = 0.3, t_anchor = 0),
                  up = list(Linf = 90, K = 0.7, t_anchor = 1))

test_that("ELEFAN_GA returns parameters within bounds and is reproducible", {
    skip_on_cran()
    fit <- function() quiet(ELEFAN_GA(synLFQ5, low_par = ga_bounds$low,
                                      up_par = ga_bounds$up, popSize = 20,
                                      maxiter = 10, MA = 11, seed = 1,
                                      plot.score = FALSE))
    a <- fit()
    b <- fit()
    expect_identical(a$par, b$par)
    for (p in c("Linf", "K", "t_anchor")) {
        expect_gte(a$par[[p]], ga_bounds$low[[p]])
        expect_lte(a$par[[p]], ga_bounds$up[[p]])
    }
    expect_gt(a$Rn_max, 0)
    expect_equal(a$Rn_max, a$fESP)
})

test_that("seasonalised ELEFAN_GA also estimates C and ts", {
    skip_on_cran()
    res <- quiet(ELEFAN_GA(synLFQ4, seasonalised = TRUE,
                           low_par = c(ga_bounds$low, C = 0, ts = 0),
                           up_par = c(ga_bounds$up, C = 1, ts = 1),
                           popSize = 20, maxiter = 5, MA = 11, seed = 1,
                           plot.score = FALSE))
    expect_true(all(c("C", "ts") %in% names(res$par)))
    expect_true(res$par$C >= 0 && res$par$C <= 1)
    expect_true(res$par$ts >= 0 && res$par$ts <= 1)
})

test_that("ELEFAN_SA draws the score graph and passes nb.stop.improvement", {
    skip_on_cran()
    fit <- function(...) {
        set.seed(2)
        ELEFAN_SA(synLFQ5, init_par = list(Linf = 78, K = 0.45, t_anchor = 0.3),
                  low_par = ga_bounds$low, up_par = ga_bounds$up,
                  SA_time = 1, maxit = 50, MA = 11, verbose = FALSE, ...)
    }
    ## issue #54: the score graph was never drawn
    expect_equal(count_plot_pages(fit()), 1)
    expect_equal(count_plot_pages(fit(plot.score = FALSE)), 0)
    expect_equal(count_plot_pages(fit(plot = TRUE)), 2)
    ## nb.stop.improvement reaches GenSA
    ctrl <- NULL
    realGenSA <- GenSA::GenSA
    local_mocked_bindings(GenSA = function(...) {
        args <- list(...)
        ctrl <<- args$control
        realGenSA(...)
    }, .package = "GenSA")
    quiet(fit(plot.score = FALSE, nb.stop.improvement = 7))
    expect_equal(ctrl$nb.stop.improvement, 7)
})

test_that("ELEFAN_SA returns parameters within bounds", {
    skip_on_cran()
    set.seed(2)
    res <- quiet(ELEFAN_SA(synLFQ5,
                           init_par = list(Linf = 78, K = 0.45, t_anchor = 0.3),
                           low_par = ga_bounds$low, up_par = ga_bounds$up,
                           SA_time = 3, maxit = 100, MA = 11, verbose = FALSE,
                           plot.score = FALSE))
    for (p in c("Linf", "K", "t_anchor")) {
        expect_gte(res$par[[p]], ga_bounds$low[[p]])
        expect_lte(res$par[[p]], ga_bounds$up[[p]])
    }
    expect_gt(res$Rn_max, 0)
})
