# Trammel net data

Data of an experiment with several trammel nets with different mesh
sizes. Can be used for function
[`select_Millar`](http://tokami.github.io/TropFishR/reference/select_Millar.md).

## Usage

``` r
data(trammelnet)
```

## Format

A list consiting of:

- `$midLengths` the midlengths of size classes,

- `$meshSizes` the meshsizes,

- `$catchPerNet_mat` a matrix with the numbers in catch of the
  corresponding mesh sizes (same order),

## Source

Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook
selectivity using log-linear models. *ICES Journal of Marine Science:
Journal du Conseil*, 54(3), 471-477.

## Examples

``` r
data(trammelnet)
str(trammelnet)
#> List of 3
#>  $ midLengths     : num [1:30] 17 18 19 20 21 22 23 24 25 26 ...
#>  $ meshSizes      : num [1:6] 6.1 7.6 7.9 9.1 10.6 13
#>  $ CatchPerNet_mat: num [1:30, 1:6] 2 1 3 3 3 1 1 0 2 0 ...
summary(trammelnet)
#>                 Length Class  Mode   
#> midLengths       30    -none- numeric
#> meshSizes         6    -none- numeric
#> CatchPerNet_mat 180    -none- numeric
```
