# Compare Millar's selectivity models

Compares several selection curves fitted with
[`select_Millar`](http://tokami.github.io/TropFishR/reference/select_Millar.md)
to the same data by the (quasi-) Akaike information criterion corrected
for small sample sizes (AICc and QAICc).

## Usage

``` r
select_Millar_compare(..., c_hat = NULL, c_hat_method = "deviance")
```

## Arguments

- ...:

  results of
  [`select_Millar`](http://tokami.github.io/TropFishR/reference/select_Millar.md)
  (class `"select_Millar"`) fitted to the same data, or one list of such
  results. The names of the arguments (or of the list) are used as model
  names, otherwise the curve types (`rtype`).

- c_hat:

  optional; overdispersion factor. By default (`NULL`), it is estimated
  from the model with the most parameters (see Details).

- c_hat_method:

  method to estimate `c_hat`: `"deviance"` (default; deviance / degrees
  of freedom) or `"pearson"` (Pearson chi-square statistic / degrees of
  freedom).

## Value

A data frame with one row per model, sorted by the QAICc, with the
columns `model`, `rtype`, number of parameters (`npar`), `deviance`,
degrees of freedom (`dof`), `c_hat`, `AICc`, `QAICc`, and the difference
to the lowest QAICc (`delta`).

## Details

Catches of gillnet and hook experiments are often overdispersed, i.e.
more variable than expected from the multinomial model of
[`select_Millar`](http://tokami.github.io/TropFishR/reference/select_Millar.md).
The QAICc accounts for this by dividing the log-likelihood by the
overdispersion factor `c_hat` (Burnham and Anderson 2002):

QAICc = -2 logL / c_hat + 2 K + 2 K (K + 1) / (n - K - 1),

where logL is the log-likelihood of the model, K the number of
parameters (plus 1 for `c_hat` if `c_hat > 1`), and n the number of
observations, here the number of non-zero cells of the catch matrix
(length classes x mesh sizes). The models can only be compared with the
same `c_hat`, which is therefore estimated from the most general model
(with the most parameters; if several models have the most parameters,
the one with the lowest deviance) and set to 1 if it is smaller than 1.
The AICc is the same criterion without overdispersion (`c_hat = 1`, K
without the additional parameter).

The log-likelihood is the multinomial log-likelihood of the catches
conditional on the total catch per length class, without constant terms.
Therefore, only differences between models fitted to the same data are
meaningful. The model with the lowest QAICc (`delta = 0`) is the best
supported model; models with `delta` \< 2 have similar support (Burnham
and Anderson 2002).

## References

Burnham, K. P., Anderson, D. R., 2002. Model selection and multimodel
inference: a practical information-theoretic approach. 2nd edition.
Springer, New York.

Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook
selectivity using log-linear models. *ICES Journal of Marine Science:
Journal du Conseil*, 54(3), 471-477.

## See also

[`select_Millar`](http://tokami.github.io/TropFishR/reference/select_Millar.md)

## Examples

``` r
data(gillnet)
fits <- list(
  norm.loc = select_Millar(gillnet, x0 = NULL, rtype = "norm.loc", plot = FALSE),
  norm.sca = select_Millar(gillnet, x0 = NULL, rtype = "norm.sca", plot = FALSE),
  lognorm = select_Millar(gillnet, x0 = NULL, rtype = "lognorm", plot = FALSE),
  binorm.sca = select_Millar(gillnet, x0 = c(55, 4, 65, 4, 3),
                             rtype = "binorm.sca", plot = FALSE),
  bilognorm = select_Millar(gillnet, x0 = c(4, 0.2, 4.2, 0.1, 2),
                            rtype = "bilognorm", plot = FALSE))
select_Millar_compare(fits)
#>        model      rtype npar deviance dof    c_hat     AICc    QAICc      delta
#> 1  bilognorm  bilognorm    5 545.5352  72 7.576878 22124.56 2931.706  0.0000000
#> 2 binorm.sca binorm.sca    5 550.4941  72 7.576878 22129.52 2932.361  0.6544772
#> 3    lognorm    lognorm    2 704.2809  75 7.576878 22276.67 2945.845 14.1390214
#> 4   norm.sca   norm.sca    2 772.7613  75 7.576878 22345.15 2954.883 23.1770946
#> 5   norm.loc   norm.loc    2 862.9241  75 7.576878 22435.31 2966.783 35.0768303

## AIC of a single fit (without overdispersion)
AIC(fits$lognorm)
#> [1] 22276.52
```
