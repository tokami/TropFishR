# Empirical M formula by Then et al. (2015) with uncertainty

Estimates the natural mortality rate (M) for each resample of a
bootstrapped growth analysis with the growth-based estimator of Then et
al. (2015) ("Then_growth" in
[`M_empirical`](https://tokami.github.io/TropFishR/reference/M_empirical.md)),
including the uncertainty of the estimator itself: for each resample,
the model of Then et al. (2015) with coefficients from a bootstrap of
their data is used for each resample.

## Usage

``` r
MstochThen(boot, resid = FALSE, seed = NULL, CI = 95)
```

## Arguments

- boot:

  an object of class \`lfqBoot\` with bootstrapped growth parameters
  (e.g. from
  [`ELEFAN_GA_boot`](https://tokami.github.io/TropFishR/reference/ELEFAN_GA_boot.md))

- resid:

  logical; should the residual error of the model of Then et al. (2015)
  be included (Default: \`FALSE\`)? See details.

- seed:

  optional; seed for the random numbers. Resample \`x\` uses
  \`set.seed(seed + x)\`. Default (\`NULL\`): the seeds of the resamples
  in \`boot\$seed\`.

- CI:

  percentage for confidence intervals (Default: 95)

## Value

The `lfqBoot` object with the additional column `M_Then` (in `bootRaw`,
`maxDen`, `median` and `CI`).

## Details

The model \\M = a K^b L\_\infty^c\\ of Then et al. (2015) was refitted
(with `nls`) to 2000 bootstrap samples of their data of 230 stocks. The
data themselves are not included in TropFishR, only the resulting 2000
sets of coefficients and the residuals of the fit to all data. For each
resample, one of the coefficient sets is selected at random and M is
predicted with the growth parameters of the resample. This accounts for
the uncertainty of the estimated coefficients (i.e. of the mean M of
species with the given growth parameters). With \`resid = TRUE\`, the
scatter of individual species around the model is included as well: the
prediction is multiplied with a randomly selected residual of the fit to
all data (on the log scale, \\M\_{obs}/M\_{fit}\\), which results in
much wider confidence intervals. The results are reproducible, because
each resample uses its own seed (see \`seed\`).

## References

Then, A. Y., J. M. Hoenig, N. G. Hall, D. A. Hewitt. 2015. Evaluating
the predictive performance of empirical estimators of natural mortality
rate using information on over 200 fish species. ICES J. Mar. Sci. 72:
82-92.

## Examples

``` r
data(alba_boot)
res <- MstochThen(alba_boot)
res$maxDen
#>       Linf          K   t_anchor       phiL     agemax     M_Then 
#> 10.0096156  2.7150279  0.4908732  2.4346091  2.0000000  3.1399201 
res$CI
#>         Linf        K  t_anchor     phiL agemax   M_Then
#> lo  9.171726 1.547109 0.2486596 2.114423      1 1.758586
#> up 11.478216 4.165708 0.5486528 2.739438      2 6.588096

## including the residual error of the model
res <- MstochThen(alba_boot, resid = TRUE)
res$CI
#>         Linf        K  t_anchor     phiL agemax     M_Then
#> lo  9.171726 1.547109 0.2486596 2.114423      1  0.4997372
#> up 11.478216 4.165708 0.5486528 2.739438      2 12.5964694
```
