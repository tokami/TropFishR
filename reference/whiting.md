# Whiting data

A dataset of North Sea whiting *Merlangius merlangus* caught during the
period 1974-1980. Can be used for
[`catchCurve`](http://tokami.github.io/TropFishR/reference/catchCurve.md)
and [`VPA`](http://tokami.github.io/TropFishR/reference/VPA.md).

## Usage

``` r
data(whiting)
```

## Format

A list consisting of:

- `age`: a vector with age groups,

- `M`: natural mortality rate,

- `a`: length-weight relationship coefficent (W = a \* L^b),

- `b`: length-weight relationship coefficent (W = a \* L^b),

- `catch`:a matrix with catches from 1974 to 1980.

## Source

ICES, 1981. Report of the *Ad hoc* working group on the use of effort
data in assessment, Copenhagen, 2-6 March 1981. *ICES C.M.* 1981/G:5
(mimeo)

## Examples

``` r
data(whiting)
str(whiting)
#> List of 5
#>  $ age  : chr [1:8] "0" "1" "2" "3" ...
#>  $ M    : num 0.2
#>  $ a    : num 0.00984
#>  $ b    : num 2.93
#>  $ catch:'data.frame':   8 obs. of  7 variables:
#>   ..$ c1974: num [1:8] 599 678 1097 275 40 ...
#>   ..$ c1975: num [1:8] 239 860 390 298 54 9 8 0
#>   ..$ c1976: num [1:8] 424 431 1071 159 75 ...
#>   ..$ c1977: num [1:8] 664 1004 532 269 32 ...
#>   ..$ c1978: num [1:8] 685 418 335 203 69 8 5 1
#>   ..$ c1979: num [1:8] 478 607 464 211 86 25 3 1
#>   ..$ c1980: num [1:8] 330 288 323 243 80 31 8 1
summary(whiting)
#>       Length Class      Mode     
#> age   8      -none-     character
#> M     1      -none-     numeric  
#> a     1      -none-     numeric  
#> b     1      -none-     numeric  
#> catch 7      data.frame list     

```
