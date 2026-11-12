## small raw data set: 8 fish over two sampling dates
raw <- data.frame(
    len  = c(1.2, 1.7, 2.5, 3.1, 3.3, 4.9, 5.5, 2.2),
    date = as.Date(c(rep("2020-01-10", 4), rep("2020-02-10", 4))),
    n    = c(1, 2, 3, 4, 5, 6, 7, 8)
)

## lfqCreate ----------------------------------------------------------

test_that("lfqCreate bins individual measurements correctly", {
    lfq <- lfqCreate(raw, Lname = "len", Dname = "date", bin_size = 1)
    expect_s3_class(lfq, "lfq")
    expect_equal(lfq$midLengths, seq(0.5, 6.5, 1))
    expect_equal(lfq$dates, as.Date(c("2020-01-10", "2020-02-10")))
    expect_equal(unname(lfq$catch[, 1]), c(0, 2, 1, 1, 0, 0, 0))
    expect_equal(unname(lfq$catch[, 2]), c(0, 0, 1, 1, 1, 1, 0))
})

test_that("lfqCreate weights measurements by the frequency column", {
    lfq <- lfqCreate(raw, Lname = "len", Dname = "date", Fname = "n", bin_size = 1)
    expect_equal(sum(lfq$catch), sum(raw$n))
    expect_equal(unname(lfq$catch[, 2]), c(0, 0, 8, 5, 6, 7, 0))
})

test_that("lfqCreate converts mm and m to cm", {
    mm <- transform(raw, len = len * 10)
    expect_equal(lfqCreate(mm, "len", "date", length_unit = "mm")$catch,
                 lfqCreate(raw, "len", "date")$catch)
    m <- transform(raw, len = len / 100)
    expect_equal(lfqCreate(m, "len", "date", length_unit = "m")$catch,
                 lfqCreate(raw, "len", "date")$catch)
})

test_that("lfqCreate aggregate_dates lumps samples to the 15th of each month", {
    lfq <- lfqCreate(raw, "len", "date", aggregate_dates = TRUE)
    expect_equal(lfq$dates, as.Date(c("2020-01-15", "2020-02-15")))
    expect_equal(sum(lfq$catch), nrow(raw))
})

test_that("lfqCreate requires Date class", {
    bad <- transform(raw, date = as.character(date))
    expect_error(lfqCreate(bad, "len", "date"), "Date")
})

test_that("lfqCreate plus group keeps the total number of fish", {
    lfq <- lfqCreate(raw, "len", "date", bin_size = 1, plus_group = c(TRUE, 3.5))
    expect_equal(max(lfq$midLengths), 3.5)
    expect_equal(sum(lfq$catch), nrow(raw))
    expect_equal(unname(lfq$catch[4, ]), c(1, 3))
})

## lfqModify ----------------------------------------------------------

test_that("lfqModify re-binning keeps total catch", {
    m <- lfqModify(synLFQ7, bin_size = 4)
    expect_equal(sum(m$catch), sum(synLFQ7$catch))
    expect_equal(unique(diff(m$midLengths)), 4)
    expect_equal(ncol(m$catch), ncol(synLFQ7$catch))
})

test_that("lfqModify vectorise_catch sums over samples", {
    v <- lfqModify(synLFQ7, vectorise_catch = TRUE)
    expect_null(dim(v$catch))
    expect_equal(length(v$catch), length(v$midLengths))
    expect_equal(sum(v$catch), sum(synLFQ7$catch))
})

test_that("lfqModify aggregates by quarter", {
    q <- lfqModify(synLFQ7, aggregate = "quarter")
    expect_equal(ncol(q$catch), 4)
    expect_equal(length(q$dates), 4)
    expect_equal(sum(q$catch), sum(synLFQ7$catch))
})

test_that("lfqModify filters by date and length", {
    d <- lfqModify(synLFQ7, minDate = as.Date("2016-03-01"),
                   maxDate = as.Date("2016-09-01"))
    expect_true(all(d$dates >= as.Date("2016-03-01") &
                    d$dates <= as.Date("2016-09-01")))
    expect_equal(ncol(d$catch), length(d$dates))

    l <- lfqModify(synLFQ7, Lmin = 20, Lmax = 100)
    expect_equal(range(l$midLengths), c(20, 100))
    expect_equal(nrow(l$catch), length(l$midLengths))
})

test_that("lfqModify numeric plus group collects all larger fish", {
    pg <- lfqModify(synLFQ7, plus_group = 100)
    expect_equal(max(pg$midLengths), 100)
    expect_equal(sum(pg$catch), sum(synLFQ7$catch))
    expect_equal(sum(pg$catch[nrow(pg$catch), ]),
                 sum(synLFQ7$catch[synLFQ7$midLengths >= 100, ]))
})

## lfqRestructure / lfqFitCurves --------------------------------------

test_that("lfqRestructure adds restructured counts of the right shape", {
    rs <- lfqRestructure(synLFQ5, MA = 11)
    expect_s3_class(rs, "lfq")
    expect_equal(dim(rs$rcounts), dim(synLFQ5$catch))
    expect_equal(dim(rs$peaks_mat), dim(synLFQ5$catch))
    expect_true(all(is.finite(rs$rcounts)))
    expect_true(rs$ASP > 0)
    expect_equal(rs$MA, 11)
})

test_that("lfqRestructure treats NA catches as zero", {
    lfq <- synLFQ5
    lfq$catch[lfq$catch == 0] <- NA
    expect_equal(lfqRestructure(lfq, MA = 11)$rcounts,
                 lfqRestructure(synLFQ5, MA = 11)$rcounts)
})

test_that("lfqFitCurves scores the true growth curve above a wrong one", {
    ## synLFQ5 was simulated with Linf = 80, K = 0.5, t_anchor = 0.25
    rs <- lfqRestructure(synLFQ5, MA = 11)
    true <- lfqFitCurves(rs, par = list(Linf = 80, K = 0.5, t_anchor = 0.25))
    wrong <- lfqFitCurves(rs, par = list(Linf = 60, K = 1, t_anchor = 0.7))
    expect_gt(true$fESP, wrong$fESP)
    expect_gt(true$fESP, 0)
    expect_lte(true$fESP, 1)
})
