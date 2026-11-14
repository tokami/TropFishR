# Hake data

This dataset contains length-frequency data and biological
characteristics about hake (*Merluccius merluccius*) and its fisheries
off Senegal. It can be used for
[`VPA`](http://tokami.github.io/TropFishR/reference/VPA.md) or
[`predict_mod`](http://tokami.github.io/TropFishR/reference/predict_mod.md).

## Usage

``` r
data(hake)
```

## Format

A list consisting of: 1. a vector with midlengths of size classes, 2. a
vector with catch in numbers, 3. K value, 4. Linf value, 5. M value, 6.
a value, 7. b value, 8. a vector with fishing mortalities, and 9. a
vector with average value of fish per kg

## Source

Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
assessment. Part 1. Manual. *FAO Fisheries Technical Paper*, (306.1,
Rev. 2). 407 p.

## Examples

``` r

data(hake)
str(hake)
#> List of 9
#>  $ midLengths: num [1:14] 9 15 21 27 33 39 45 51 57 63 ...
#>  $ catch     : num [1:14] 1823 14463 25227 8134 3889 ...
#>  $ K         : num 0.1
#>  $ Linf      : num 130
#>  $ M         : num 0.28
#>  $ a         : num 1e-05
#>  $ b         : num 3
#>  $ FM        : num [1:14] 0.04 0.39 1.07 0.65 0.49 0.59 0.65 0.39 0.29 0.31 ...
#>  $ meanValue : num [1:14] 1 1 1 1.5 1.5 2 2 2.5 2.5 2.5 ...
summary(hake)
#>            Length Class  Mode   
#> midLengths 14     -none- numeric
#> catch      14     -none- numeric
#> K           1     -none- numeric
#> Linf        1     -none- numeric
#> M           1     -none- numeric
#> a           1     -none- numeric
#> b           1     -none- numeric
#> FM         14     -none- numeric
#> meanValue  14     -none- numeric

```
