# Synthetical catch per unit of effort (CPUE) dataset

Synthetical CPUE dataset from Exercise 4.3 in Sparre & Venema (1999).
Can be used to estimate the instantaneous total mortality rate (Z) by
means of the function
[`Z_CPUE`](http://tokami.github.io/TropFishR/reference/Z_CPUE.md).

## Usage

``` r
data(synCPUE)
```

## Format

A dataframe consisting of following columns:

- `cohort`: name of cohort, e.g. 1982 S, meaning summer cohort of 1982,

- `age`: age of cohorts,

- `CPUE`: catch per unit of effort of cohorts.

## Source

Sparre, P., Venema, S.C., 1999. Introduction to tropical fish stock
assessment. Part 2. Excercises. *FAO Fisheries Technical Paper*, (306.2,
Rev. 2). 94 p.

## Examples

``` r
data(synCPUE)
head(synCPUE)
#>   cohort  age CPUE
#> 1  1983S 0.64  182
#> 2  1982A 1.14  111
#> 3  1982S 1.64   67
#> 4  1981A 2.14   40
#> 5  1981S 2.64   24
#> 6  1980A 3.14   15
str(synCPUE)
#> 'data.frame':    6 obs. of  3 variables:
#>  $ cohort: Factor w/ 6 levels "1980A","1981A",..: 6 4 5 2 3 1
#>  $ age   : num  0.64 1.14 1.64 2.14 2.64 3.14
#>  $ CPUE  : num  182 111 67 40 24 15

```
