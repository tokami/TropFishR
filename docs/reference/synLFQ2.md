# Synthetic length frequency data II

Synthetic length frequency data from Sparre & Venema (1998). Can be used
for the estimation of the instantaneous total mortality rate (Z) by
means of
[`Z_BevertonHolt`](http://tokami.github.io/TropFishR/reference/Z_BevertonHolt.md).

## Usage

``` r
data(synLFQ2)
```

## Format

A list consisting of:

- **dates** dates of sampling times (class Date),

- `midLengths`: midpoints of the length classes,

- `Linf`: infinite length for investigated species in cm \[cm\],

- `K`: growth coefficent for investigated species per year \[1/year\],

- `catch`: matrix with catches/counts per length class (row) and
  sampling date (column).

## Source

Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
assessment. Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev.
2). 407 p.

## Examples

``` r
data(synLFQ2)
str(synLFQ2)
#> List of 5
#>  $ dates     : Date[1:3], format: "1960-01-01" "1970-01-01" ...
#>  $ midLengths: num [1:11] 47.5 52.5 57.5 62.5 67.5 72.5 77.5 82.5 87.5 92.5 ...
#>  $ Linf      : num 100
#>  $ K         : num 0.3
#>  $ catch     : num [1:11, 1:3] 256 237 211 187 161 138 113 87 62 36 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "catch1960" "catch1970" "catch1980"
summary(synLFQ2)
#>            Length Class  Mode   
#> dates       3     Date   numeric
#> midLengths 11     -none- numeric
#> Linf        1     -none- numeric
#> K           1     -none- numeric
#> catch      33     -none- numeric

```
