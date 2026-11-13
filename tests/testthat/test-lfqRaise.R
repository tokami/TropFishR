## Raising of length-frequency samples to the total catch (lfqRaise) and the
## landings argument of VPA

## two samples in two years, length classes 10 and 20 cm; W = 0.01 * L^3 kg
## gives 10 kg and 80 kg per fish
lfq_toy <- function() {
    structure(list(midLengths = c(10, 20),
                   dates = as.Date(c("2020-06-01", "2021-06-01")),
                   catch = matrix(c(2, 1,
                                    4, 0), nrow = 2),
                   par = list(a = 0.01, b = 3)),
              class = "lfq")
}

test_that("lfqRaise raises each sample to its landings", {
    ## sample weights: 2*10 + 1*80 = 100 kg and 4*10 = 40 kg
    res <- lfqRaise(lfq_toy(), landings = c(1000, 200))
    expect_equal(res$raising$sampleWeight, c(100, 40))
    expect_equal(res$raising$factor, c(10, 5))
    expect_true(all(res$raising$estimated))
    expect_equal(res$catch, matrix(c(20, 10, 20, 0), nrow = 2))
    ## raised catch in weight equals the landings
    expect_equal(colSums(res$catch * 0.01 * res$midLengths^3), c(1000, 200))
    expect_s3_class(res, "lfq")
})

test_that("lfqRaise by year matches named landings to the years", {
    lfq <- lfq_toy()
    lfq$dates <- as.Date(c("2021-03-01", "2020-06-01"))
    res <- lfqRaise(lfq, landings = c("2020" = 400, "2021" = 1000), by = "year")
    expect_equal(res$raising$year, c("2020", "2021"))
    expect_equal(res$raising$factor, c(10, 10))
    ## two samples in the same year share one factor
    lfq$dates <- as.Date(c("2020-03-01", "2020-06-01"))
    res <- lfqRaise(lfq, landings = 1400, by = "year")
    expect_equal(res$raising$sampleWeight, 140)
    expect_equal(res$catch, lfq$catch * 10)
    expect_error(lfqRaise(lfq, landings = c("2019" = 1), by = "year"), "years")
})

test_that("lfqRaise uses measured sample weights and the unit of a", {
    res <- lfqRaise(lfq_toy(), landings = c(1000, 200), sampleWeight = c(50, 40))
    expect_equal(res$raising$factor, c(20, 5))
    expect_false(any(res$raising$estimated))
    ## a in g: same result as a in kg
    lfq <- lfq_toy()
    lfq$par$a <- 10
    expect_equal(lfqRaise(lfq, landings = c(1000, 200), LW_unit = "g")$catch,
                 lfqRaise(lfq_toy(), landings = c(1000, 200))$catch)
    ## a and b as arguments
    lfq$par <- NULL
    expect_error(lfqRaise(lfq, landings = c(1000, 200)), "length-weight")
    expect_equal(lfqRaise(lfq, landings = c(1000, 200), a = 0.01, b = 3)$catch,
                 lfqRaise(lfq_toy(), landings = c(1000, 200))$catch)
})

test_that("lfqRaise checks its input", {
    lfq <- lfq_toy()
    expect_error(lfqRaise(lfq, landings = 1000), "one value per sample")
    expect_error(lfqRaise(lfq, landings = c(1000, NA)), "finite")
    expect_error(lfqRaise(lfqRaise(lfq, landings = c(1, 1)), landings = c(1, 1)),
                 "already been raised")
    lfq$catch[, 2] <- 0
    expect_error(lfqRaise(lfq, landings = c(1000, 200)), "2021-06-01")
    ## no landings for an empty sample are fine
    expect_equal(lfqRaise(lfq, landings = c(1000, 0))$catch[, 2], c(0, 0))
    ## catch vector (one sample)
    v <- list(midLengths = c(10, 20), catch = c(2, 1), a = 0.01, b = 3)
    expect_equal(lfqRaise(v, landings = 1000)$catch, c(20, 10))
})

test_that("raised catches are refused by the resampling", {
    raised <- lfqRaise(lfq_toy(), landings = c(1000, 200))
    expect_error(lfqResample(raised), "raised")
    expect_error(checkSampledLfq(raised, "ELEFAN_GA_boot"), "raised")
    noninteger <- lfq_toy()
    noninteger$catch <- noninteger$catch * 1.5
    expect_warning(checkSampledLfq(noninteger, "ELEFAN_GA_boot"), "non-integer")
    expect_silent(checkSampledLfq(lfq_toy(), "ELEFAN_GA_boot"))
})

test_that("VPA with landings equals VPA of the raised catch per year", {
    data(synLFQ7)
    lfq <- synLFQ7
    lfq$par <- list(Linf = 123, K = 0.2, t0 = 0, M = 0.2, a = 0.015, b = 3)
    pg <- lfq$midLengths[max(which(lfq$midLengths < lfq$par$Linf))]
    lfq <- lfqModify(lfq, plus_group = pg)
    landings <- seq(20000, 40000, length.out = length(lfq$dates))

    expect_no_warning(res <- VPA(lfq, landings = landings, terminalE = 0.5,
                                 LW_unit = "g"))
    raised <- lfqRaise(lfq, landings = landings, LW_unit = "g")
    ref <- VPA(list(midLengths = lfq$midLengths, catch = rowSums(raised$catch),
                    Linf = 123, K = 0.2, t0 = 0, M = 0.2, a = 0.015 / 1000, b = 3),
               terminalE = 0.5, catch_unit = "'000")
    expect_equal(res$FM_calc, ref$FM_calc)
    expect_equal(res$annualMeanNr, ref$annualMeanNr / 1000)
    expect_equal(res$raising, raised$raising)

    ## one raising factor: same F at length, numbers proportional to landings
    r1 <- VPA(lfq, landings = rep(30000, 8), terminalE = 0.5, LW_unit = "g")
    r2 <- VPA(lfq, landings = rep(60000, 8), terminalE = 0.5, LW_unit = "g")
    expect_equal(r1$FM_calc, r2$FM_calc)
    expect_equal(r2$annualMeanNr, 2 * r1$annualMeanNr)

    expect_error(VPA(lfq, landings = landings, terminalE = 0.5, catch_unit = "kg"),
                 "catch_unit = NA")
    data(whiting)
    expect_error(VPA(whiting, landings = 1, catch_columns = 1, terminalE = 0.5),
                 "length-based")
})
