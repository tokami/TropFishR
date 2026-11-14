# Predict gillnet selectivity (old Millar method)

Predict gillnet selectivity (old Millar method)

## Usage

``` r
rcurves_Millar(type, meshsizes, rel, pars, plotlens)
```

## Source

https://www.stat.auckland.ac.nz/~millar/selectware/

## Arguments

- type:

  ("norm.loc", "norm.sca", "gamma", "lognorm")

- meshsizes:

  mesh sizes

- rel:

  relative powers (one for each mesh size)

- pars:

  selection curve parameters

- plotlens:

  vector of new mesh sizes for selectivity prediction

## Value

selectivities

## References

Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook
selectivity using log-linear models. *ICES Journal of Marine Science:
Journal du Conseil*, 54(3):471-477
