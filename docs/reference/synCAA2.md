# Synthetic Catch-at-age data II

Synthetic Catch-at-age data from Sparre & Venema (1998). Can be used for
the estimation of the instantaneous mortality rate (Z) by means of the
cumulative catch curve
([`catchCurve`](https://tokami.github.io/TropFishR/reference/catchCurve.md)).

## Usage

``` r
data(synCAA2)
```

## Format

A list consisting of:

- `midAge` a vector of the mid ages of the age groups,

- `catch` a matrix with the catches for different years.

## Source

Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
assessment. Part 1. Manual. *FAO Fisheries Technical Paper*, (306.1,
Rev. 2). 407 p.

## Examples

``` r
data(synCAA2)
str(synCAA2)
#> List of 2
#>  $ age  : num [1:8] 0 1 2 3 4 5 6 7
#>  $ catch: num [1:8] 488 612 601 237 62.3 15.7 4.7 1.4
summary(synCAA2)
#>       Length Class  Mode   
#> age   8      -none- numeric
#> catch 8      -none- numeric

```
