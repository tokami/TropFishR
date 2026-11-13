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

## Millar's selectivity models: predict and select_Millar_curve (issue #43)

test_that("predict.select_Millar reproduces the fitted retention", {
    o <- select_Millar(gillnet, x0 = NULL, rtype = "lognorm", plot = FALSE)
    L <- seq(40, 90, 5)
    ref <- outer(L, o$meshSizes, rtypes_Millar("lognorm"), o$par)
    expect_equal(unname(predict(o, lengths = L)), ref)
    expect_equal(dim(predict(o)), dim(gillnet$CatchPerNet_mat))
    expect_equal(max(predict(o, lengths = L, standardise = TRUE)), 1)
    ## a new mesh size keeps the first fitted mesh size as reference
    expect_equal(unname(predict(o, lengths = L, meshSizes = 17.8)[, 1]),
                 unname(predict(o, lengths = L)[, 7]))
    expect_false(isTRUE(all.equal(predict(o, lengths = L, meshSizes = 17.8)[, 1],
                                  predict(o, lengths = L)[, 1])))
    expect_error(predict(o, meshSizes = 16, rel.power = c(1, 1)), "one value per mesh")
})

test_that("select_Millar_curve with the estimates equals predict for all curve types", {
    L <- seq(40, 90, 5)
    fits <- list(
        select_Millar(gillnet, x0 = NULL, rtype = "norm.loc", plot = FALSE),
        select_Millar(gillnet, x0 = NULL, rtype = "norm.sca", plot = FALSE),
        select_Millar(gillnet, x0 = NULL, rtype = "lognorm", plot = FALSE),
        select_Millar(gillnet, x0 = c(55, 4, 65, 4, 3), rtype = "binorm.sca", plot = FALSE),
        select_Millar(gillnet, x0 = c(4, 0.2, 4.2, 0.1, 2), rtype = "bilognorm", plot = FALSE),
        select_Millar(gillnet, x0 = c(30, 0.15), rtype = "gamma", plot = FALSE))
    for(o in fits){
        expect_equal(select_Millar_curve(L, o$meshSizes, o$rtype, o$estimates),
                     predict(o, lengths = L), info = o$rtype)
        ## one mesh size with the reference mesh size
        expect_equal(select_Millar_curve(L, 16.6, o$rtype, o$estimates[, "par"],
                                         refMesh = 13.5)[, 1],
                     predict(o, lengths = L)[, 6], info = o$rtype)
    }
})

test_that("select_Millar_curve: mode and SD of published parameters", {
    L <- seq(30, 90, 0.01)
    for(rt in c("norm.loc", "norm.sca", "lognorm", "gamma")){
        sel <- select_Millar_curve(L, c(10, 12), rt, par = c(50, 5))
        ## maximum 1 at the mode of the reference mesh size
        expect_equal(L[which.max(sel[, 1])], 50, tolerance = 1e-3, info = rt)
        expect_equal(max(sel[, 1]), 1, tolerance = 1e-6, info = rt)
        ## mode proportional to the mesh size
        expect_equal(L[which.max(sel[, 2])], 60, tolerance = 1e-3, info = rt)
    }
    ## SD of the lognormal and gamma curves (as densities)
    for(rt in c("lognorm", "gamma")){
        th <- millarPar(rt, c(50, 5), refMesh = 10)
        d <- if(rt == "lognorm") function(x) dlnorm(x, th[1], th[2]) else
                 function(x) dgamma(x, shape = th[1], scale = th[2] * 10)
        m1 <- integrate(function(x) x * d(x), 0, Inf)$value
        m2 <- integrate(function(x) x^2 * d(x), 0, Inf)$value
        expect_equal(sqrt(m2 - m1^2), 5, tolerance = 1e-6, info = rt)
    }
    ## logistic curve: 25, 50 and 75 % at L50 - SR/2, L50, L50 + SR/2
    sel <- select_Millar_curve(c(28, 30, 32), 1, "tt.logistic", par = c(30, 4))
    expect_equal(unname(sel[, 1]), c(0.25, 0.5, 0.75))
    expect_error(select_Millar_curve(L, 10, "lognorm", par = 50), "2 values")
    expect_error(select_Millar_curve(L, 10, "foo", par = c(50, 5)), "not recognised")
    expect_error(select_Millar_curve(L, 10, "binorm.sca", par = c(50, 5, 60, 5, 1)),
                 "proportion")
})

test_that("select_Millar fits the gamma curve", {
    ## catches of 8 gillnets with a gamma selection curve (shape 40, mode
    ## 55 cm for mesh size 13.5)
    k <- 40; a <- 55 / ((k - 1) * 13.5)
    L <- seq(40.5, 89.5, 1)
    sim <- list(midLengths = L, meshSizes = gillnet$meshSizes,
                CatchPerNet_mat = 1e5 * dnorm(L, 65, 10) *
                    outer(L, gillnet$meshSizes, rtypes_Millar("gamma"), c(k, a)))
    g <- select_Millar(sim, x0 = c(30, 0.15), rtype = "gamma", plot = FALSE)
    expect_equal(g$par, c(k, a), tolerance = 1e-3)
    expect_equal(unname(g$estimates[, "par"]), c(55, sqrt(k) * a * 13.5),
                 tolerance = 1e-3)
    expect_true(all(is.finite(g$estimates[, "s.e."])))
})

test_that("plot.select_Millar uses predict", {
    o <- select_Millar(gillnet, x0 = NULL, rtype = "lognorm", plot = FALSE)
    expect_equal(count_plot_pages(plot(o, plotlens = seq(40, 90, 1))), 1)
})

## comparison of Millar's selectivity models (issue #39)

test_that("select_Millar_compare computes AICc and QAICc with a common c_hat", {
    fits <- list(
        lognorm = select_Millar(gillnet, x0 = NULL, rtype = "lognorm", plot = FALSE),
        norm.loc = select_Millar(gillnet, x0 = NULL, rtype = "norm.loc", plot = FALSE),
        bilognorm = select_Millar(gillnet, x0 = c(4, 0.2, 4.2, 0.1, 2),
                                  rtype = "bilognorm", plot = FALSE))
    tab <- select_Millar_compare(fits)
    n <- sum(gillnet$CatchPerNet_mat > 0)
    ## c_hat from the most general model (bilognorm)
    c_hat <- unname(fits$bilognorm$out["Deviance", 1] / fits$bilognorm$out["d.o.f.", 1])
    expect_equal(tab$c_hat, rep(c_hat, 3))
    ## by hand: -2 logL = 2 * value (optim minimises the negative logL)
    qaicc <- function(f, K) 2 * f$value / c_hat + 2 * K + 2 * K * (K + 1) / (n - K - 1)
    aicc <- function(f, K) 2 * f$value + 2 * K + 2 * K * (K + 1) / (n - K - 1)
    row <- function(m) tab[tab$model == m, ]
    expect_equal(row("lognorm")$QAICc, qaicc(fits$lognorm, 3))
    expect_equal(row("bilognorm")$QAICc, qaicc(fits$bilognorm, 6))
    expect_equal(row("norm.loc")$AICc, aicc(fits$norm.loc, 2))
    ## sorted, best model first
    expect_equal(tab$model, c("bilognorm", "lognorm", "norm.loc"))
    expect_equal(tab$delta, tab$QAICc - min(tab$QAICc))
    ## differences of the QAIC follow the differences of the deviance
    expect_equal(row("norm.loc")$QAICc - row("lognorm")$QAICc,
                 (row("norm.loc")$deviance - row("lognorm")$deviance) / c_hat)
    ## arguments instead of a list, names from rtype
    tab2 <- select_Millar_compare(fits$lognorm, fits$norm.loc)
    expect_setequal(tab2$model, c("lognorm", "norm.loc"))
    ## given c_hat; c_hat = 1 gives the AICc
    tab1 <- select_Millar_compare(fits, c_hat = 1)
    expect_equal(tab1$QAICc, tab1$AICc)
    ## Pearson
    tabP <- select_Millar_compare(fits, c_hat_method = "pearson")
    expect_equal(tabP$c_hat[1], unname(fits$bilognorm$out["Pearson.chisq", 1] /
                                       fits$bilognorm$out["d.o.f.", 1]))
})

test_that("select_Millar_compare requires fits of the same data", {
    a <- select_Millar(gillnet, x0 = NULL, rtype = "lognorm", plot = FALSE)
    other <- gillnet
    other$CatchPerNet_mat[1, 1] <- other$CatchPerNet_mat[1, 1] + 10
    b <- select_Millar(other, x0 = NULL, rtype = "lognorm", plot = FALSE)
    expect_error(select_Millar_compare(a, b), "same data")
    expect_error(select_Millar_compare(a, list(1)), "select_Millar")
    expect_error(select_Millar_compare(a, c_hat = -1), "c_hat")
})

test_that("logLik.select_Millar allows AIC and BIC", {
    a <- select_Millar(gillnet, x0 = NULL, rtype = "lognorm", plot = FALSE)
    ll <- logLik(a)
    expect_equal(as.numeric(ll), -a$value)
    expect_equal(attr(ll, "df"), 2)
    expect_equal(attr(ll, "nobs"), sum(gillnet$CatchPerNet_mat > 0))
    expect_equal(AIC(a), 2 * a$value + 2 * 2)
    expect_equal(BIC(a), 2 * a$value + 2 * log(sum(gillnet$CatchPerNet_mat > 0)))
})
