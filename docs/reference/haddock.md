# Haddock data

Data of a covered codend experimental catch of the haddock
(*Melanogrammus aeglefinus*). Can be used for function
[`select_Millar`](http://tokami.github.io/TropFishR/reference/select_Millar.md).

## Usage

``` r
data(haddock)
```

## Format

A list consisting of:

- `midLengths` the midlengths of size classes,

- `numCodend` the number of fish retained in codend,

- `numCover` the number of fish retained in cover,

## Source

Millar, R. B., Holst, R., 1997. Estimation of gillnet and hook
selectivity using log-linear models. *ICES Journal of Marine Science:
Journal du Conseil*, 54(3), 471-477.

## Examples

``` r
data(haddock)
str(haddock)
#> List of 4
#>  $ midLengths: num [1:24] 24 25 26 27 28 29 30 31 32 33 ...
#>  $ meshSizes : num [1:2] 0 130
#>  $ numCover  : num [1:24] 1 1 3 14 30 49 60 50 70 108 ...
#>  $ numCodend : num [1:24] 0 0 0 1 5 19 29 51 91 120 ...
summary(haddock)
#>            Length Class  Mode   
#> midLengths 24     -none- numeric
#> meshSizes   2     -none- numeric
#> numCover   24     -none- numeric
#> numCodend  24     -none- numeric
```
