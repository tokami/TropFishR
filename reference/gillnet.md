# Gillnet data

Data of an experiment with several gillnets with different mesh sizes.
Data can be analysed with function
[`select_Millar`](http://tokami.github.io/TropFishR/reference/select_Millar.md).

## Usage

``` r
data(gillnet)
```

## Format

A list consiting of:

- `midLengths` the midlengths of size classes,

- `meshSizes` the meshsizes,

- `catchPerNet_mat` a matrix with the numbers in catch of the
  corresponding mesh sizes (same order),

## Source

Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook
selectivity using log-linear models. *ICES Journal of Marine Science:
Journal du Conseil*, 54(3), 471-477.

Holt, S. J. 1963. A method for determining gear selectivity and its
application. *ICNAF Special Publication*, 5: 106-115.

## Examples

``` r
data(gillnet)
str(gillnet)
#> List of 3
#>  $ midLengths     : num [1:11] 52.5 54.5 56.5 58.5 60.5 62.5 64.5 66.5 68.5 70.5 ...
#>  $ meshSizes      : num [1:8] 13.5 14 14.8 15.4 15.9 16.6 17.8 19
#>  $ CatchPerNet_mat: num [1:11, 1:8] 52 102 295 309 118 79 27 14 8 7 ...
summary(gillnet)
#>                 Length Class  Mode   
#> midLengths      11     -none- numeric
#> meshSizes        8     -none- numeric
#> CatchPerNet_mat 88     -none- numeric
```
