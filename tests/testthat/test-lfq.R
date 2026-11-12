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

## Straightforward (slow) implementation of the ELEFAN score for the growth
## curve points `fit$Lt` returned by lfqFitCurves
reference_ESP <- function(rs, fit, flagging.out = TRUE) {
    mid <- rs$midLengths
    w <- diff(mid)
    lower <- mid - c(w[1], w) / 2
    upper <- mid + c(w, w[length(w)]) / 2
    yd <- date2yeardec(rs$dates)
    hits <- matrix(0, length(mid), length(yd))
    for (k in seq_len(nrow(fit$Lt))) {
        j <- which(yd == fit$Lt$t[k])
        i <- which(lower <= fit$Lt$Lt[k] & upper > fit$Lt$Lt[k])
        hits[i, j] <- hits[i, j] + 1
    }
    if (!flagging.out) hits <- (hits > 0) * 1
    if (flagging.out) {
        ## each positive peak (id != 0) counts at most once, at its maximum
        for (id in setdiff(unique(c(rs$peaks_mat)), 0)) {
            cells <- which(rs$peaks_mat == id)
            if (sum(hits[cells]) > 1) {
                hits[cells] <- 0
                hits[cells[which.max(rs$rcounts[cells])]] <- 1
            }
        }
    }
    sum(rs$rcounts * hits)
}

test_that("lfqFitCurves score equals a straightforward implementation", {
    set.seed(3)
    ## truncated data: the first cell belongs to a positive peak
    trunc <- lfqRestructure(lfqModify(synLFQ5, Lmin = 23.5), MA = 11)
    expect_true(trunc$peaks_mat[1, 1] != 0)
    for (rs in list(lfqRestructure(synLFQ5, MA = 11),
                    lfqRestructure(synLFQ4, MA = 11),
                    lfqRestructure(alba, MA = 7),
                    trunc)) {
        Lmax <- max(rs$midLengths)
        for (i in 1:5) {
            p <- list(Linf = runif(1, 0.8, 1.4) * Lmax, K = runif(1, 0.2, 1.5),
                      t_anchor = runif(1), C = runif(1), ts = runif(1))
            for (fo in c(TRUE, FALSE)) {
                fit <- lfqFitCurves(rs, par = p, flagging.out = fo)
                expect_equal(fit$ESP, reference_ESP(rs, fit, flagging.out = fo))
            }
        }
    }
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
