# Log-likelihood of Millar's selectivity model

Log-likelihood of a selection curve fitted with
[`select_Millar`](https://tokami.github.io/TropFishR/reference/select_Millar.md),
so that `AIC` and `BIC` can be used. The log-likelihood is the
multinomial log-likelihood conditional on the total catch per length
class, without constant terms; the number of observations is the number
of non-zero cells of the catch matrix. For overdispersed data use
[`select_Millar_compare`](https://tokami.github.io/TropFishR/reference/select_Millar_compare.md).

## Usage

``` r
# S3 method for class 'select_Millar'
logLik(object, ...)
```

## Arguments

- object:

  a list of the class `"select_Millar"`.

- ...:

  unused.

## Value

An object of class `"logLik"`.
