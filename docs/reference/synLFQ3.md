# Synthetic length frequency data III

Synthetic length frequency data from Sparre & Venema (1998). Can be used
for the estimation of the instantaneous total mortality rate (Z) by
means of
[`Z_BevertonHolt`](https://tokami.github.io/TropFishR/reference/Z_BevertonHolt.md).

## Usage

``` r
data(synLFQ3)
```

## Format

A list consisting of:

- `midLengths`: midpoints of the length classes,

- `Linf`: infinite length for investigated species in cm \[cm\],

- `K`: growth coefficent for investigated species per year \[1/year\],

- `t0`: theoretical time zero, at which individuals of this species
  hatch,

- `catch`: a vector with catches.

## Source

Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
assessment. Part 1. Manual. *FAO Fisheries Technical Paper*, (306.1,
Rev. 2). 407 p.

## Examples

``` r
data(synLFQ3)
str(synLFQ3)
#> List of 5
#>  $ midLengths: num [1:22] 4 6 8 10 12 14 16 18 20 22 ...
#>  $ Linf      : num 50
#>  $ K         : num 0.3
#>  $ t0        : num 0
#>  $ catch     : num [1:22] 37 56 86 129 188 258 319 352 351 324 ...
summary(synLFQ3)
#>            Length Class  Mode   
#> midLengths 22     -none- numeric
#> Linf        1     -none- numeric
#> K           1     -none- numeric
#> t0          1     -none- numeric
#> catch      22     -none- numeric

```
