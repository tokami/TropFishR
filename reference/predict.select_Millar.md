# Predict the retention of Millar's selectivity model

Relative retention of a fitted selectivity model
([`select_Millar`](http://tokami.github.io/TropFishR/reference/select_Millar.md))
for any lengths and mesh sizes.

## Usage

``` r
# S3 method for class 'select_Millar'
predict(
  object,
  lengths = object$midLengths,
  meshSizes = object$meshSizes,
  rel.power = NULL,
  standardise = FALSE,
  ...
)
```

## Arguments

- object:

  a list of the class `"select_Millar"` with the results of
  [`select_Millar`](http://tokami.github.io/TropFishR/reference/select_Millar.md).

- lengths:

  lengths for which the retention is predicted (default: the length
  classes of the fitted data).

- meshSizes:

  mesh sizes for which the retention is predicted (default: the mesh
  sizes of the fitted data). The parameters of the model refer to the
  first (smallest) mesh size of the fitted data, which stays the
  reference mesh size for new mesh sizes.

- rel.power:

  optional; relative fishing power of the mesh sizes (one value per mesh
  size). By default (`NULL`), the relative power of the fit is used if
  `meshSizes` are the fitted mesh sizes, otherwise 1.

- standardise:

  logical; should the retention be divided by its maximum (Default:
  FALSE)?

- ...:

  unused.

## Value

A matrix with the relative retention with one row per length and one
column per mesh size.

## Details

The retention curves have a maximum of 1 (except bi-modal curves and
`"tt.logistic"`), multiplied with the relative power of the mesh size.
For `"tt.logistic"`, the retention of the control gear (mesh size equal
to the reference mesh size) is `1 - p` and of the experimental gear `p`
times the logistic curve.

## See also

[`select_Millar`](http://tokami.github.io/TropFishR/reference/select_Millar.md),
[`select_Millar_curve`](http://tokami.github.io/TropFishR/reference/select_Millar_curve.md)

## Examples

``` r
data(gillnet)
output <- select_Millar(gillnet, x0 = NULL, rtype = "lognorm", plot = FALSE)

## retention at the fitted mesh sizes
predict(output, lengths = seq(40, 90, 10))
#>            13.5           14         14.8         15.4         15.9
#> 40 4.137726e-04 6.244458e-05 2.335651e-06 1.660280e-07 1.658797e-08
#> 50 5.060169e-01 2.691024e-01 6.897279e-02 1.941783e-02 5.866897e-03
#> 60 5.435442e-01 8.089745e-01 9.991528e-01 8.660670e-01 6.463011e-01
#> 70 1.031329e-02 3.664257e-02 1.710401e-01 3.836639e-01 6.149459e-01
#> 80 1.687224e-05 1.273807e-04 1.881027e-03 9.614958e-03 2.988319e-02
#> 90 5.914385e-09 8.681135e-08 3.540442e-06 3.742166e-05 2.085826e-04
#>            16.6         17.8           19
#> 40 5.779795e-10 1.363332e-12 2.409029e-15
#> 50 9.090006e-04 2.404794e-05 4.070291e-07
#> 60 3.389027e-01 6.461861e-02 6.929337e-03
#> 70 9.039608e-01 9.155153e-01 4.676279e-01
#> 80 1.072831e-01 4.616077e-01 9.114487e-01
#> 90 1.646035e-03 2.536969e-02 1.650984e-01

## retention of a new mesh size
predict(output, lengths = seq(40, 90, 10), meshSizes = 16)
#>              16
#> 40 1.036140e-08
#> 50 4.553416e-03
#> 60 5.989843e-01
#> 70 6.621611e-01
#> 80 3.664257e-02
#> 90 2.868230e-04
```
