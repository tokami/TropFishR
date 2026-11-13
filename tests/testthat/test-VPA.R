test_that("age-based VPA satisfies Z = F + M and the catch equation", {
    res <- suppressWarnings(  # catch_unit not specified
        VPA(whiting, catch_columns = 1, terminalE = 0.5, analysis_type = "VPA"))
    expect_equal(res$Z, res$FM_calc + whiting$M)

    n <- length(res$Z)
    ## catch equation C = F/Z * N * (1 - exp(-Z)) for all but the plus group
    ## (VPA solves it iteratively, hence the tolerance)
    pred <- res$FM_calc / res$Z * res$survivors_L1 * (1 - exp(-res$Z))
    expect_equal(pred[-n], res$catch.cohort[-n], tolerance = 1e-3)
    ## survivors of one class are the starting numbers of the next
    expect_equal(res$survivors_L2[-n], res$survivors_L1[-1])
    ## terminal exploitation rate
    expect_equal(res$FM_calc[n] / res$Z[n], 0.5)
})

test_that("Pope's cohort analysis approximates VPA", {
    vpa <- suppressWarnings(VPA(whiting, catch_columns = 1, terminalE = 0.5,
                                analysis_type = "VPA"))
    ca <- suppressWarnings(VPA(whiting, catch_columns = 1, terminalE = 0.5,
                               analysis_type = "CA"))
    expect_equal(ca$FM_calc, vpa$FM_calc, tolerance = 0.05)
})

test_that("length-based cohort analysis on hake (regression values)", {
    res <- VPA(hake, terminalE = 0.5, analysis_type = "CA",
               catch_unit = "'000", plus_group = TRUE)
    expect_equal(round(res$FM_calc, 4),
                 c(0.0402, 0.3874, 1.0661, 0.6468, 0.4911, 0.5910, 0.6462,
                   0.3855, 0.2882, 0.3075, 0.4013, 0.3894, 0.1101, 0.2800))
    expect_equal(res$Z, res$FM_calc + hake$M)
})

test_that("VPA with LW_unit = 'g' equals a in kg", {
    hake_g <- hake
    hake_g$a <- hake$a * 1000
    kg <- VPA(hake, terminalE = 0.5, analysis_type = "CA", catch_unit = "'000",
              plus_group = TRUE)
    g <- VPA(hake_g, terminalE = 0.5, analysis_type = "CA", catch_unit = "'000",
             plus_group = TRUE, LW_unit = "g")
    expect_equal(g$meanBiomassTon, kg$meanBiomassTon)
    expect_equal(g$FM_calc, kg$FM_calc)
})

test_that("VPA checks terminalF and terminalE (issue #47)", {
    run <- function(...) VPA(hake, analysis_type = "CA", catch_unit = "'000", ...)
    ## several values, e.g. Z - M with M of two methods
    expect_error(run(terminalF = c(0.5, 0.6)), "single value, but has 2")
    expect_error(run(terminalE = c(0.5, 0.6)), "single value, but has 2")
    expect_error(run(terminalF = NULL), "either the terminal")
    expect_error(run(), "either the terminal")
    expect_error(run(terminalF = -0.2), "> 0")
    expect_error(run(terminalF = Inf), "finite")
    expect_error(run(terminalE = 1), "< 1")
    expect_error(run(terminalE = 0), "> 0")
    ## both provided: terminalF is used
    expect_warning(both <- run(terminalF = 0.5, terminalE = 0.2), "terminalE is ignored")
    expect_equal(both$FM_calc, run(terminalF = 0.5)$FM_calc)
    expect_equal(tail(both$FM_calc, 1), 0.5)
    ## terminalE and the equivalent terminalF give the same result
    M <- hake$M
    expect_equal(run(terminalE = 0.5)$FM_calc, run(terminalF = M)$FM_calc)
    ## NULL is treated as not provided; a named value is fine
    expect_equal(run(terminalF = NULL, terminalE = 0.5)$FM_calc,
                 run(terminalE = c(E = 0.5))$FM_calc)
})
