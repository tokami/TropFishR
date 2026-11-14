# Beverton & Holt's Z-Equations

A method to estimate the instantaneous total mortality rate (Z) based on
a method derived by Beverton and Holt (1956).

## Usage

``` r
Z_BevertonHolt(param, catch_columns = NA, Lprime_tprime)
```

## Arguments

- param:

  a list consisting of following parameters:

  - `midLengths` or `age`: midpoints of length groups (length-frequency
    data) or ages (age composition data),

  - `Linf`: infinite length for investigated species in cm \[cm\],

  - `K`: growth coefficent for investigated species per year \[1/year\],

  - `t0`: theoretical time zero, at which individuals of this species
    hatch,

  - `catch`: catch as vector, or a matrix with catches of subsequent
    years;

- catch_columns:

  optional; in case catch is a matrix or data.frame, a number or vector
  indicating which column(s) of the matrix should be analysed (Default:
  `NA`).

- Lprime_tprime:

  length or age prime, above which all fish are under full exploitation
  as mid length or age class.

## Value

A list with the input parameters and following objects:

- **tmean** or **Lmean**: mean age or length of fish,

- **tprime** or **Lprime**: some age or length for which all fish of
  that length and longer are under full exploitation,

- **Z**: total mortality.

## Details

The first length group or age class within the list object `midLengths`
or `age` will be used as the Lprim or tprime (length of recruitment to
fishery).

## References

Beverton R.J.H and S.J. Holt, 1956. A review of methods of estimating
mortality rates in exploited fish populations, with special reference to
sources of bias in catch sampling. *Rapp.P.-v.Reun.CIEM*, 140:67-83

Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
assessment. Part 1. Manual. *FAO Fisheries Technical Paper*, (306.1,
Rev. 2). 407 p.

## Examples

``` r
# based on length-frequency data
data(synLFQ2)
Z_BevertonHolt(synLFQ2, catch_columns = 2, Lprime_tprime = 47.5)
#> $dates
#> [1] "1960-01-01" "1970-01-01" "1980-01-01"
#> 
#> $midLengths
#>  [1] 47.5 52.5 57.5 62.5 67.5 72.5 77.5 82.5 87.5 92.5 97.5
#> 
#> $Linf
#> [1] 100
#> 
#> $K
#> [1] 0.3
#> 
#> $catch
#>       catch1960 catch1970 catch1980
#>  [1,]       256       268       212
#>  [2,]       237       226       161
#>  [3,]       211       180       116
#>  [4,]       187       141        79
#>  [5,]       161       105        52
#>  [6,]       138        76        31
#>  [7,]       113        50        17
#>  [8,]        87        30         8
#>  [9,]        62        15         3
#> [10,]        36         6         1
#> [11,]        12         1         0
#> 
#> $Lmean
#> [1] 58.898
#> 
#> $Lprime
#> [1] 45
#> 
#> $Z
#> [1] 0.8872215
#> 

# based on age composition data
data(synCAA1)
Z_BevertonHolt(synCAA1, catch_columns = 3, Lprime_tprime = 2.5)
#> $midAge
#>  [1]  2.25  2.75  3.25  3.75  4.25  4.75  5.25  5.75  6.25  6.75  7.25  7.75
#> [13]  8.25  8.75  9.25  9.75 10.25 10.75 11.25
#> 
#> $catch
#>       catch1960 catch1970 catch1980
#>  [1,]       390       399       308
#>  [2,]       289       255       169
#>  [3,]       214       162        92
#>  [4,]       150       103        52
#>  [5,]       118        66        28
#>  [6,]        97        42        15
#>  [7,]        64        27         8
#>  [8,]        48        17         6
#>  [9,]        35        11         2
#> [10,]        26         7        NA
#> [11,]        19         4        NA
#> [12,]        14         3        NA
#> [13,]        11         2        NA
#> [14,]         8        NA        NA
#> [15,]         6        NA        NA
#> [16,]         4        NA        NA
#> [17,]         3        NA        NA
#> [18,]         2        NA        NA
#> [19,]         2        NA        NA
#> 
#> $tmean
#> [1] 2.839706
#> 
#> $tprime
#> [1] 2.25
#> 
#> $Z
#> [1] 1.695761
#> 
```
