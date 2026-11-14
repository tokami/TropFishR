# Estimate Z from CPUE data

Method to estimate the instantaneous total mortality rate (Z) from catch
per unit of effort (CPUE) data according to standard, Heincke's, or
Robson & Chapman's method.

## Usage

``` r
Z_CPUE(param, method = "standard", omit_age1 = FALSE)
```

## Arguments

- param:

  a list consisting of following parameters:

  - `cohort`: a vector with with a cohort label,

  - `age`: a vector with ages,

  - `CPUE`: a vector with CPUE values;

- method:

  a character string indicating which assessment method should be used:
  `"standard"`, `"Heincke"`, or `"RobsonChapman"`.

- omit_age1:

  logical; if `TRUE` the first age group is omitted (Default `FALSE`).

## Value

A list with input parameters and a Z value or matrix depending on the
method.

## Details

In Heincke's and RobsonChapman's method age groups older than 4 are
lumped, because age groups older than 3 or 4 are said to be hard to
seperate (Ricker, 1975). Sparre and Venema (1998) recommend to omit the
first age group in case it is not fully exploited by the fishery.

## References

Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
assessment. Part 1. Manual. *FAO Fisheries Technical Paper*, (306.1,
Rev. 2). 407 p.

Sparre, P., Venema, S.C., 1999. Introduction to tropical fish stock
assessment. Part 2. Excercises. *FAO Fisheries Technical Paper*, (306.2,
Rev. 2). 94 p.

Ricker, W.E., 1975. Computation and interpretation of biological
statistics of fish populations. *Bull.Fish.Res.Board Can.*, (191):382 p.

## Examples

``` r
# load data
data(synCPUE)

# run model with standard method
Z_CPUE(synCPUE, method = "standard")
#> $cohort
#> [1] 1983S 1982A 1982S 1981A 1981S 1980A
#> Levels: 1980A 1981A 1981S 1982A 1982S 1983S
#> 
#> $age
#> [1] 0.64 1.14 1.64 2.14 2.64 3.14
#> 
#> $CPUE
#> [1] 182 111  67  40  24  15
#> 
#> $Z_mat
#>   cohort 1982A 1982S 1981A 1981S 1980A
#> 1  1983S  0.99  1.00  1.01  1.01  1.00
#> 2  1982A    NA  1.01  1.02  1.02  1.00
#> 3  1982S    NA    NA  1.03  1.03  1.00
#> 4  1981A    NA    NA    NA  1.02  0.98
#> 5  1981S    NA    NA    NA    NA  0.94
#> 

# run model with Heincke's method
Z_CPUE(synCPUE, method = "Heincke")
#> $cohort
#> [1] 1983S 1982A 1982S 1981A 1981S 1980A
#> Levels: 1980A 1981A 1981S 1982A 1982S 1983S
#> 
#> $age
#> [1] 0.64 1.14 1.64 2.14 2.64 3.14
#> 
#> $CPUE
#> [1] 182 111  67  40  24  15
#> 
#> $Z
#> [1] 0.5354233
#> 

# run model with Robson and Chapman's method
Z_CPUE(synCPUE, method = "RobsonChapman", omit_age1 = TRUE)
#> $cohort
#> [1] 1983S 1982A 1982S 1981A 1981S 1980A
#> Levels: 1980A 1981A 1981S 1982A 1982S 1983S
#> 
#> $age
#> [1] 0.64 1.14 1.64 2.14 2.64 3.14
#> 
#> $CPUE
#> [1] 182 111  67  40  24  15
#> 
#> $Z
#> [1] 0.651055
#> 
```
