# Millar's selectivity types

This function returns a function corresponding to the type of curve
which was selected to represent the selectivity of nets or hooks.

## Usage

``` r
rtypes_Millar(rtype)
```

## Source

https://www.stat.auckland.ac.nz/~millar/selectware/

## Arguments

- rtype:

  a character string indicating which method for the estimation of
  selection curves should be used: `"norm.loc"` for normal with common
  spread method, `"norm.sca"` for normal with variable spread method,
  `"lognorm"` for lognormal method, `"binorm.sca"` for bi-normal method,
  `"bilognorm"` for bi-lognormal method, `"tt.logistic"` for control and
  logistic method, `"gamma"` for gamma method.

## Details

Function adapted from the selectivity functions provided by Prof. Dr.
Russell Millar (https://www.stat.auckland.ac.nz/~millar/). Until now
following curves are incorporated: `"norm.loc"` for a normal curve with
common spread, `"norm.sca"` for a normal curve with variable spread,
`"lognorm"` for a lognormal curve, `"binorm.sca"` for a bi-normal curve,
`"bilognorm"` for a bi-lognormal curve, `"tt.logistic"` for a control
and logistic curve, `"gamma"` for a gamma curve.

The returned function has the arguments `classes`, `meshSizes` and `th`
(parameters on the internal scale, see
[`select_Millar`](https://tokami.github.io/TropFishR/reference/select_Millar.md)).
The first value of `meshSizes` is the reference mesh size to which the
parameters refer. Use
[`predict.select_Millar`](https://tokami.github.io/TropFishR/reference/predict.select_Millar.md)
or
[`select_Millar_curve`](https://tokami.github.io/TropFishR/reference/select_Millar_curve.md)
instead of this internal function.

## References

Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook
selectivity using log-linear models. *ICES Journal of Marine Science:
Journal du Conseil*, 54(3), 471-477.
