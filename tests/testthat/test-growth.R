## length-at-age data generated from a known VBGF
true <- list(Linf = 60, K = 0.3, t0 = -0.5)
age <- rep(1:8, each = 3)
exact <- list(age = age, length = VBGF(true, t = age))

test_that("growth_length_age recovers VBGF parameters from exact data", {
    for (m in c("FordWalford", "Chapman", "BertalanffyPlot")) {
        res <- quiet(suppressWarnings(  # perfect-fit warnings from summary.lm
            growth_length_age(exact, method = m, Linf_est = 60)))
        expect_equal(res$Linf, 60, tolerance = 1e-6, label = m)
        expect_equal(res$K, 0.3, tolerance = 1e-6, label = m)
    }
    res <- quiet(suppressWarnings(
        growth_length_age(exact, method = "BertalanffyPlot", Linf_est = 60)))
    expect_equal(res$t0, -0.5, tolerance = 1e-6)
})

test_that("growth_length_age Gulland-Holt approximates K", {
    res <- quiet(suppressWarnings(
        growth_length_age(exact, method = "GullandHolt")))
    expect_equal(res$K, 0.3, tolerance = 0.01)
})

test_that("growth_length_age LSM fits noisy data", {
    set.seed(1)
    noisy <- list(age = age, length = exact$length + rnorm(length(age), 0, 0.5))
    res <- quiet(growth_length_age(noisy, method = "LSM", Linf_init = 50,
                                   K_init = 0.2, age_plot = NULL))
    expect_equal(res$Linf, 60, tolerance = 0.05)
    expect_equal(res$K, 0.3, tolerance = 0.1)
})

test_that("growth_length_age LSM confidence band (delta method, issue #40)", {
    set.seed(1)
    noisy <- list(age = age, length = exact$length + rnorm(length(age), 0, 0.5))
    ages <- seq(1, 7, 0.5)
    res <- quiet(growth_length_age(noisy, method = "LSM", Linf_init = 50,
                                   K_init = 0.2, CI = TRUE, age_plot = ages))
    band <- res$predCI
    expect_equal(band$age_plot, ages)
    expect_equal(band$fit, unname(predict(res$mod, newdata = data.frame(t = ages))))
    expect_true(all(band$lower < band$fit & band$fit < band$upper))
    ## standard error from a numerical gradient
    cf <- coef(res$mod)
    vb <- function(p) p[1] * (1 - exp(-p[2] * (ages - p[3])))
    G <- sapply(1:3, function(i) {
        h <- 1e-6 * max(1, abs(cf[i])); d <- replace(numeric(3), i, h)
        (vb(cf + d) - vb(cf - d)) / (2 * h)
    })
    se <- sqrt(rowSums((G %*% vcov(res$mod)) * G))
    q <- qt(0.975, df.residual(res$mod))
    expect_equal(band$upper - band$fit, q * se, tolerance = 1e-6)
    ## the band follows ci.level
    res90 <- quiet(growth_length_age(noisy, method = "LSM", Linf_init = 50,
                                     K_init = 0.2, CI = TRUE, ci.level = 0.9,
                                     age_plot = ages))
    expect_equal((res90$predCI$upper - res90$predCI$fit) / (band$upper - band$fit),
                 rep(qt(0.95, df.residual(res$mod)) / q, length(ages)))
    expect_warning(quiet(growth_length_age(noisy, method = "LSM", Linf_init = 50,
                                           K_init = 0.2, do.sim = TRUE)),
                   "deprecated")
})

test_that("growth_length_age LSM confidence band covers the true curve", {
    set.seed(2)
    ages <- c(2, 4, 6)
    covered <- replicate(300, {
        d <- list(age = age, length = exact$length + rnorm(length(age), 0, 1))
        b <- lsmBand(nls(Lt ~ Linf * (1 - exp(-K * (t - t0))),
                         data = list(Lt = d$length, t = d$age),
                         start = list(Linf = 50, K = 0.2, t0 = 0)),
                     ages, 0.95)
        Lt <- VBGF(true, t = ages)
        b$lower <= Lt & Lt <= b$upper
    })
    expect_true(all(abs(rowMeans(covered) - 0.95) < 0.04))
})

## tagging data with exact VBGF increments
tag <- local({
    L1 <- c(20, 30, 40, 50, 60, 70)
    dt <- c(0.5, 1, 0.3, 0.8, 1.2, 0.4)
    list(L1 = L1, L2 = 100 - (100 - L1) * exp(-0.5 * dt), delta_t = dt)
})

test_that("growth_tagging Gulland-Holt recovers Linf and K", {
    res <- quiet(growth_tagging(tag, "GullandHolt", time_unit = "year"))
    expect_equal(res$Linf, 100, tolerance = 0.01)
    expect_equal(res$K, 0.5, tolerance = 0.01)
})

test_that("growth_tagging Munro recovers Linf and K", {
    res <- quiet(growth_tagging(tag, "Munro", time_unit = "year",
                                Linf_range = c(80, 150)))
    expect_equal(res$Linf, 100, tolerance = 1e-3)
    expect_equal(res$indi_Ks, rep(0.5, 6), tolerance = 1e-3)
    expect_equal(res$K, 0.5, tolerance = 1e-3)
})

test_that("startingPoint2tanchor converts a cohort start to t_anchor", {
    lfq <- list(dates = as.Date(c("2020-03-01", "2020-06-01")))
    res <- startingPoint2tanchor(lfq, par = list(Linf = 80, K = 0.5),
                                 startingLength = 20, startingSample = 1)
    age_at_20 <- VBGF(list(Linf = 80, K = 0.5, t0 = 0), L = 20)
    expect_equal(res$t_anchor, (date2yeardec(lfq$dates[1]) - age_at_20) %% 1)
})
