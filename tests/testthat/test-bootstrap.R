## lfqResample --------------------------------------------------------

test_that("lfqResample keeps the structure and the number of fish per sample", {
    set.seed(1)
    res <- lfqResample(alba)
    expect_s3_class(res, "lfq")
    expect_equal(res$midLengths, alba$midLengths)
    expect_equal(res$dates, alba$dates)
    expect_equal(dim(res$catch), dim(alba$catch))
    expect_equal(colSums(res$catch), colSums(alba$catch))
    ## only length classes observed in a sample can be drawn
    expect_true(all(res$catch[alba$catch == 0] == 0))
    expect_false(identical(res$catch, alba$catch))
})

test_that("lfqResample is reproducible and allows a fixed sample size", {
    set.seed(5); a <- lfqResample(alba)
    set.seed(5); b <- lfqResample(alba)
    expect_identical(a, b)
    s <- lfqResample(alba, sampSize = 50)
    expect_equal(unname(colSums(s$catch)), rep(50, ncol(alba$catch)))
})

## bootstrapped ELEFAN ------------------------------------------------

## small bootstrap runs (alba, few resamples, short optimisation); computed
## once and re-used by the tests below
boot_cache <- new.env()
alba_bounds <- list(low = list(Linf = 8, K = 0.1, t_anchor = 0),
                    up = list(Linf = 20, K = 5, t_anchor = 1))
ga_boot_args <- list(lfq = alba, MA = 7, low_par = alba_bounds$low,
                     up_par = alba_bounds$up, popSize = 20, maxiter = 5,
                     nresamp = 10, seed = 1)
small_ga_boot <- function() {
    if (is.null(boot_cache$ga)) boot_cache$ga <- quiet(do.call(ELEFAN_GA_boot, ga_boot_args))
    boot_cache$ga
}

test_that("ELEFAN_GA_boot returns one fit per resample", {
    skip_on_cran()
    res <- small_ga_boot()
    expect_s3_class(res, "lfqBoot")
    expect_equal(nrow(res$bootRaw), 10)
    expect_equal(names(res$bootRaw), c("Linf", "K", "t_anchor", "phiL", "agemax"))
    expect_equal(res$seed, 1 + 1:10)
    expect_length(res$Rn_max, 10)
    for (p in c("Linf", "K", "t_anchor")) {
        expect_true(all(res$bootRaw[[p]] >= alba_bounds$low[[p]] &
                        res$bootRaw[[p]] <= alba_bounds$up[[p]]))
    }
})

test_that("ELEFAN_GA_boot resamples are ELEFAN_GA fits to lfqResample(seed + x)", {
    skip_on_cran()
    res <- small_ga_boot()
    for (x in 1:2) {
        set.seed(res$seed[x])
        lfqb <- lfqResample(alba)
        fit <- quiet(ELEFAN_GA(lfqb, MA = 7, low_par = alba_bounds$low,
                               up_par = alba_bounds$up, popSize = 20, maxiter = 5,
                               run = 10, pmutation = 0.2, plot.score = FALSE))
        expect_equal(unlist(res$bootRaw[x, c("Linf", "K", "t_anchor", "phiL")]),
                     unlist(fit$par))
        expect_equal(res$Rn_max[x], round(fit$Rn_max, 3))
    }
    ## lfqResample with boot reproduces the resampled data sets
    lfqs <- lfqResample(alba, boot = res)
    expect_length(lfqs, 10)
    set.seed(res$seed[2])
    expect_equal(lfqs[[2]], lfqResample(alba))
})

test_that("ELEFAN_GA_boot summaries are consistent", {
    skip_on_cran()
    res <- small_ga_boot()
    pars <- c("Linf", "K", "t_anchor")
    ## maximum density and median are within the range of the resamples
    for (p in pars) {
        expect_true(res$maxDen[[p]] >= min(res$bootRaw[[p]]) &&
                    res$maxDen[[p]] <= max(res$bootRaw[[p]]))
        expect_equal(res$median[[p]], median(res$bootRaw[[p]]))
        expect_equal(unname(res$CI[, p]),
                     unname(quantile(res$bootRaw[[p]], c(0.025, 0.975))))
    }
    expect_true(all(res$CI["lo", ] <= res$CI["up", ]))
    expect_true(all(res$multiCI["lo", ] <= res$multiCI["up", ]))
    ## the maximum density is one of the resamples
    expect_true(any(apply(res$bootRaw[, pars], 1, function(r)
        isTRUE(all.equal(unname(r), unname(res$maxDen[pars]))))))
    ## derived parameters
    expect_equal(res$maxDen[["phiL"]],
                 log10(res$maxDen[["K"]]) + 2 * log10(res$maxDen[["Linf"]]))
    expect_equal(res$maxDen[["agemax"]], ceiling(-log(0.05) / res$maxDen[["K"]]))
    expect_equal(colnames(res$CI), c(pars, "phiL", "agemax"))
})

test_that("ELEFAN_GA_boot writes progress to outfile", {
    skip_on_cran()
    out <- withr::local_tempfile(fileext = ".txt")
    args <- modifyList(ga_boot_args, list(nresamp = 4, outfile = out))
    quiet(do.call(ELEFAN_GA_boot, args))
    expect_length(readLines(out), 4)
})

test_that("ELEFAN_GA_boot gives the same results in parallel", {
    skip_on_cran()
    skip_if_not(parallel::detectCores() >= 2, "needs 2 cores")
    ## the workers load the installed TropFishR, so install this version
    ## into a temporary library first
    lib <- withr::local_tempdir()
    pkg <- normalizePath(testthat::test_path("..", ".."))
    skip_if_not(file.exists(file.path(pkg, "DESCRIPTION")), "package sources not found")
    status <- system2(file.path(R.home("bin"), "R"),
                      c("CMD", "INSTALL", "--no-docs", "--no-multiarch",
                        "--no-test-load", "-l", shQuote(lib), shQuote(pkg)),
                      stdout = FALSE, stderr = FALSE)
    skip_if(status != 0, "could not install TropFishR into a temporary library")
    withr::local_envvar(R_LIBS = lib)

    par <- quiet(do.call(ELEFAN_GA_boot,
                         c(ga_boot_args, parallel = TRUE, no_cores = 2)))
    ser <- small_ga_boot()
    expect_equal(par$bootRaw, ser$bootRaw)
    expect_equal(par$seed, ser$seed)
    expect_equal(par$maxDen, ser$maxDen)
})

test_that("ELEFAN_SA_boot returns one fit per resample", {
    skip_on_cran()
    res <- quiet(ELEFAN_SA_boot(alba, MA = 7,
                                init_par = list(Linf = 12, K = 1, t_anchor = 0.5),
                                low_par = alba_bounds$low, up_par = alba_bounds$up,
                                SA_time = 1, maxit = 20, nresamp = 10, seed = 3))
    expect_s3_class(res, "lfqBoot")
    expect_equal(nrow(res$bootRaw), 10)
    expect_equal(names(res$bootRaw), c("Linf", "K", "t_anchor", "phiL", "agemax"))
    expect_equal(res$seed, 3 + 1:10)
    for (p in c("Linf", "K", "t_anchor")) {
        expect_true(all(res$bootRaw[[p]] >= alba_bounds$low[[p]] &
                        res$bootRaw[[p]] <= alba_bounds$up[[p]]))
    }
})

## example data and plots ---------------------------------------------

test_that("alba_boot has the structure returned by ELEFAN_GA_boot", {
    expect_s3_class(alba_boot, "lfqBoot")
    expect_equal(names(alba_boot),
                 c("bootRaw", "seed", "Rn_max", "maxDen", "median", "CI",
                   "multiCI", "misc"))
    expect_equal(nrow(alba_boot$bootRaw), 200)
    expect_equal(names(alba_boot$bootRaw), c("Linf", "K", "t_anchor", "phiL", "agemax"))
})

test_that("bootstrap plots are drawn", {
    expect_equal(count_plot_pages(univariate_density(alba_boot)), 1)
    expect_equal(count_plot_pages(LinfK_scatterhist(alba_boot, phi.contour = TRUE)), 1)
    expect_equal(count_plot_pages(plotBoot(alba, alba_boot)), 1)
    expect_equal(count_plot_pages({
        plot(c(5, 20), c(0.1, 2), t = "n")
        add_phiprime()
    }), 1)
})

test_that("vbgfCI_time returns CI limits that enclose the max. density curve", {
    res <- NULL
    pages <- count_plot_pages(res <- vbgfCI_time(alba_boot, agemax = 2, CI = c(95, 50)))
    expect_equal(pages, 1)
    expect_equal(names(res), c("limCI", "inCI", "density", "max_dens"))
    expect_equal(names(res$limCI), c("CI95", "CI50"))
    expect_equal(nrow(res$inCI), nrow(alba_boot$bootRaw))
    ## fewer resamples in the 50 % than in the 95 % region
    expect_lt(sum(res$inCI$CI50), sum(res$inCI$CI95))
    ci <- res$limCI$CI95
    ok <- is.finite(ci$min) & is.finite(ci$max)
    expect_true(all(ci$min[ok] <= ci$max[ok]))
    expect_equal(names(res$max_dens), c("Linf", "K", "t_anchor"))
})
