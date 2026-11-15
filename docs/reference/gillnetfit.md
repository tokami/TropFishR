# Millar's original gillnet selectivity fitting function

Function to estimate selectivity parameters from experimental data. This
function is applied within
[`select_Millar`](http://tokami.github.io/TropFishR/reference/select_Millar.md)
to derive starting parameters.
[`select_Millar`](http://tokami.github.io/TropFishR/reference/select_Millar.md)
is the recommended function for selectivity estimation.

## Usage

``` r
gillnetfit(
  data,
  meshsizes,
  rtype = "norm.loc",
  rel.power = NULL,
  plotlens = NULL,
  details = FALSE
)
```

## Source

https://www.stat.auckland.ac.nz/~millar/selectware/

## Arguments

- data:

  matrix with the number of individuals caught with each sized mesh
  (`CatchPerNet_mat`).

- meshsizes:

  vector with meshSizes in increasing order (`meshSizes`),

- rtype:

  A character string indicating which method for estimating selection
  curves should be used: `"norm.loc"` for a normal curve with common
  spread, `"norm.sca"` for a normal curve with variable spread,
  `"lognorm"` for a lognormal curve, `"gamma"` for a gamma curve.

- rel.power:

  A string indicating the relative power of different meshSizes, must
  have same length as `meshSizes` (Default: `rel.power = NULL`).

- plotlens:

  lengths which should be used for graphical output, for more detailed
  curves. Default : NULL

- details:

  logical; should details be included in the output?

## Value

list of fitted parameters

## References

Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook
selectivity using log-linear models. *ICES Journal of Marine Science:
Journal du Conseil*, 54(3):471-477

## Examples

``` r
data(gillnet)

dat <- matrix(c(gillnet$midLengths, gillnet$CatchPerNet_mat),
         byrow = FALSE, ncol=(dim(gillnet$CatchPerNet_mat)[2]+1))

gillnetfit(data = dat, meshsizes = gillnet$meshSizes)
#> $fit.type
#> [1] "norm.loc "
#> 
#> $gear.pars
#>                      estimate        s.e.
#> k                    4.044451 0.005411702
#> sigma                5.147949 0.067096964
#> mode(mesh1)         54.600084 0.073057978
#> std_dev(all meshes)  5.147949 0.067096964
#> 
#> $fit.stats
#>      model_dev Pearson chi-sq            dof       null_dev 
#>       862.9241      2580.9355        75.0000     50908.7399 
#> 

```
