# Synthetic Catch-at-age data I

Synthetic Catch-at-age data from Sparre & Venema (1998). Can be used for
the estimation of the instantaneous total mortality rate (Z) by means of
[`Z_BevertonHolt`](https://tokami.github.io/TropFishR/reference/Z_BevertonHolt.md).

## Usage

``` r
data(synCAA1)
```

## Format

A list consisting of:

- `midAge`: a vector of the mid ages of the age groups,

- `catch`: a matrix with the catches for different years.

## Source

Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
assessment. Part 1. Manual. *FAO Fisheries Technical Paper*, (306.1,
Rev. 2). 407 p.

## Examples

``` r
data(synCAA1)
str(synCAA1)
#> List of 2
#>  $ midAge: num [1:19] 2.25 2.75 3.25 3.75 4.25 4.75 5.25 5.75 6.25 6.75 ...
#>  $ catch : num [1:19, 1:3] 390 289 214 150 118 97 64 48 35 26 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:3] "catch1960" "catch1970" "catch1980"
summary(synCAA1)
#>        Length Class  Mode   
#> midAge 19     -none- numeric
#> catch  57     -none- numeric

```
