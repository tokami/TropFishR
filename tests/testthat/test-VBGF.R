test_that("VBGF length-at-age matches the closed-form equation", {
    par <- list(Linf = 80, K = 0.5, t0 = -0.2)
    t <- c(0, 0.5, 1, 2, 5, 10)
    expect_equal(VBGF(par, t = t),
                 80 * (1 - exp(-0.5 * (t + 0.2))))
})

test_that("VBGF age-at-length is the inverse of length-at-age", {
    par <- list(Linf = 80, K = 0.5, t0 = -0.2)
    t <- c(0.3, 1, 2.5, 7)
    expect_equal(VBGF(par, L = VBGF(par, t = t)), t)
})

test_that("VBGF with C = 0 equals the non-seasonal curve", {
    t <- seq(0.1, 5, 0.3)
    expect_equal(VBGF(list(Linf = 60, K = 0.4, t0 = 0, C = 0, ts = 0.3), t = t),
                 VBGF(list(Linf = 60, K = 0.4, t0 = 0), t = t))
})

test_that("seasonal VBGF stays monotone for C <= 1 and approaches Linf", {
    t <- seq(0, 30, 0.01)
    L <- VBGF(list(Linf = 60, K = 0.4, t0 = 0, C = 1, ts = 0.5), t = t)
    expect_true(all(diff(L) >= -1e-10))
    expect_lt(abs(tail(L, 1) - 60), 0.01)
})

test_that("seasonal VBGF age-at-length is the inverse of length-at-age", {
    par <- list(Linf = 100, K = 0.3, t0 = 0, C = 0.5, ts = 0.2)
    t <- c(0.5, 2.3, 7.1)
    expect_equal(VBGF(par, L = VBGF(par, t = t)), t, tolerance = 1e-3)
})

test_that("age-at-length at and above Linf is the same with and without seasonality", {
    L <- c(100, 101)
    nonseasonal <- suppressWarnings(VBGF(list(Linf = 100, K = 0.3, t0 = 0), L = L))
    seasonal <- suppressWarnings(VBGF(list(Linf = 100, K = 0.3, t0 = 0, C = 0.5, ts = 0.2), L = L))
    expect_equal(nonseasonal, c(Inf, NaN))
    expect_equal(seasonal, c(Inf, NaN))
})

test_that("VBGF accepts t_anchor in place of t0", {
    L <- VBGF(list(Linf = 80, K = 0.5, t_anchor = 0.25), t = c(1, 2))
    expect_equal(L, 80 * (1 - exp(-0.5 * (c(1, 2) - 0.25))))
})

test_that("VBGF na.rm drops missing ages", {
    par <- list(Linf = 80, K = 0.5, t0 = 0)
    expect_equal(VBGF(par, t = c(1, NA, 2), na.rm = TRUE),
                 VBGF(par, t = c(1, 2)))
})
