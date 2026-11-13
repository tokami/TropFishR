## Bhattacharya: the interactive selection (identify, readline) is replaced
## by a list of selections (first and last point of each regression line)

run_bhat <- function(lfq, sel) {
    msgs <- character(0)
    k <- 0
    e <- new.env(parent = asNamespace("TropFishR"))
    e$identify <- function(...) {
        k <<- k + 1
        if (k > length(sel)) return(list(ind = integer(0), pos = integer(0)))
        list(ind = sel[[k]], pos = c(1, 1))
    }
    e$readline <- function(...) "y"
    e$interactive <- function() TRUE
    e$dev.new <- function(...) grDevices::pdf(NULL)
    e$writeLines <- function(text, ...) msgs <<- c(msgs, text)
    f <- Bhattacharya
    environment(f) <- e
    on.exit(grDevices::graphics.off())
    res <- f(lfq)
    list(res = res, msgs = msgs, n_identify = k)
}

## length-frequency of normal distributions (not rounded)
lfq_norm <- function(midLengths, means, sds, n = 1e6) {
    catch <- rowSums(sapply(seq_along(means), function(i)
        n * dnorm(midLengths, means[i], sds[i])))
    list(midLengths = midLengths, catch = catch)
}

test_that("Bhattacharya recovers mean and SD of a normal distribution", {
    lfq <- lfq_norm(10:50, means = 30, sds = 3)
    ## points around the mean
    out <- run_bhat(lfq, list(c(18, 26)))
    ls <- out$res$Lmean_SD_list
    expect_equal(unname(ls[1, "Lmean"]), 30, tolerance = 1e-3)
    expect_equal(unname(ls[1, "SD"]), 3, tolerance = 1e-3)
})

test_that("Bhattacharya separates two distributions with class width 0.5", {
    lfq <- lfq_norm(seq(10, 50, 0.5), means = c(20, 35), sds = c(2, 2))
    L <- lfq$midLengths
    sel <- list(range(which(L >= 17 & L <= 23)), range(which(L >= 32 & L <= 38)))
    out <- run_bhat(lfq, sel)
    ls <- out$res$Lmean_SD_list
    expect_equal(unname(ls[, "Lmean"]), c(20, 35), tolerance = 1e-2)
    expect_equal(unname(ls[, "SD"]), c(2, 2), tolerance = 1e-2)
})

test_that("bhatMidChoice returns the middle row of the selection", {
    expect_equal(bhatMidChoice(c(2, 8)), 5)     # 7 points
    expect_equal(bhatMidChoice(c(12, 17)), 15)  # 6 points: one after the midpoint
    expect_equal(bhatMidChoice(c(3, 5)), 4)
})

test_that("Bhattacharya rejects selections without falling line (issue #52)", {
    data(synLFQ1)
    ## points 19 to 23 are the third cohort: selected first, the line rises
    out <- run_bhat(synLFQ1, list(c(19, 23)))
    expect_true(any(grepl("falling line", out$msgs)))
    expect_true(any(grepl("did not choose any points", out$msgs)))
    ## the rejected selection is followed by a new one
    out <- run_bhat(synLFQ1, list(c(19, 23), c(2, 8)))
    expect_equal(out$n_identify, 3)
    expect_equal(nrow(out$res$Lmean_SD_list), 1)
    expect_true(out$res$Lmean_SD_list[1, "SD"] > 0)
    ## skipping the second cohort: rejected instead of SD = NaN
    out <- run_bhat(synLFQ1, list(c(2, 8), c(19, 23)))
    expect_true(any(grepl("falling line", out$msgs)))
    expect_equal(nrow(out$res$Lmean_SD_list), 1)
    expect_false(any(is.nan(out$res$Lmean_SD_list[, "SD"])))
})

test_that("Bhattacharya example selection of the help page works", {
    data(synLFQ1)
    out <- run_bhat(synLFQ1, list(c(2, 8), c(12, 17), c(19, 23), c(26, 30)))
    ls <- out$res$Lmean_SD_list
    expect_equal(nrow(ls), 4)
    expect_true(all(is.finite(ls[, "SD"]) & ls[, "SD"] > 0))
    expect_true(all(diff(ls[, "Lmean"]) > 0))
    expect_false(any(grepl("not possible", out$msgs)))
})

test_that("Bhattacharya is reproducible and returns the input data", {
    data(synLFQ1)
    sel <- list(c(2, 8), c(12, 17), c(19, 23), c(26, 30))
    r1 <- run_bhat(synLFQ1, sel)$res
    r2 <- run_bhat(synLFQ1, sel)$res
    expect_identical(r1$Lmean_SD_list, r2$Lmean_SD_list)
    expect_identical(r1$bhat_results, r2$bhat_results)
    expect_equal(r1$midLengths, synLFQ1$midLengths)
    expect_equal(r1$catch, synLFQ1$catch)
    expect_s3_class(r1, "Bhattacharya")
    expect_true(all(nzchar(names(r1))))
})

test_that("bhatRange gives the quantiles of the expected extremes", {
    expect_equal(bhatRange(30, 3, 1000), 30 + 3 * qnorm(c(0.625, 999.625) / 1000.25))
    ## close to the mean of the extremes of n random draws (previous method)
    set.seed(1)
    ext <- replicate(500, range(rnorm(1000, 30, 3)))
    expect_equal(bhatRange(30, 3, 1000), rowMeans(ext), tolerance = 0.01)
})
