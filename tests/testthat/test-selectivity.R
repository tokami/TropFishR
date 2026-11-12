test_that("trawl ogive passes through 25/50/75 % at L25, L50, L75", {
    s <- list(selecType = "trawl_ogive", L50 = 30, L75 = 35)
    expect_equal(select_ogive(s, Lt = c(25, 30, 35)), c(0.25, 0.5, 0.75))
    expect_true(all(diff(select_ogive(s, Lt = 10:60)) > 0))
})

test_that("trawl ogive shifts with a new Lc keeping the L75/L50 ratio", {
    s <- list(selecType = "trawl_ogive", L50 = 30, L75 = 35)
    L75_new <- 35 * 40 / 30
    expect_equal(select_ogive(s, Lt = c(40, L75_new), Lc = 40), c(0.5, 0.75))
})

test_that("knife-edge selection is 0 below and 1 from L50", {
    s <- list(selecType = "knife_edge", L50 = 30)
    expect_equal(select_ogive(s, Lt = c(29, 30, 31)), c(0, 1, 1))
    expect_equal(select_ogive(s, Lt = c(29, 30, 31), Lc = 31), c(0, 0, 1))
})

test_that("select_ogive rejects unknown curve types", {
    expect_error(select_ogive(list(selecType = "foo"), Lt = 1:3), "not recognised")
})

test_that("trawl selectivity from covered codend data (bream)", {
    res <- select(bream)
    expect_lt(res$L25, res$L50)
    expect_lt(res$L50, res$L75)
    ## meshSizes = c(cover, codend); SF refers to the codend mesh
    expect_equal(res$SF, res$L50 / bream$meshSizes[2], tolerance = 1e-6)
    expect_equal(c(res$L50, res$L75), c(13.168332, 14.701427), tolerance = 1e-6)
})

test_that("gillnet selectivity from two mesh sizes (tilapia)", {
    res <- select(tilapia)
    expect_equal(res$LmNet1 / res$LmNet2,
                 tilapia$meshSizes[1] / tilapia$meshSizes[2], tolerance = 1e-6)
    expect_equal(res$SF, 2.610694, tolerance = 1e-6)
})

## production models

test_that("Schaefer and Fox MSY follow from the regression coefficients", {
    res <- prod_mod(trawl_fishery_Java, plot = FALSE)
    a <- res$Schaefer_lm[1]; b <- res$Schaefer_lm[2]
    expect_equal(unname(res$Schaefer_MSY), unname(-a^2 / (4 * b)))
    expect_equal(unname(res$Schaefer_fMSY), unname(-a / (2 * b)))
    c <- res$Fox_lm[1]; d <- res$Fox_lm[2]
    expect_equal(unname(res$Fox_fMSY), unname(-1 / d))
    expect_equal(unname(res$Fox_MSY), unname(-(1 / d) * exp(c - 1)))
})

test_that("Fox virgin stock index is positive", {
    res <- prod_mod(trawl_fishery_Java, plot = FALSE)
    expect_gt(res$Fox_Bv, 0)
    expect_equal(unname(res$Fox_Bv), unname(exp(res$Fox_lm[1])))
})
