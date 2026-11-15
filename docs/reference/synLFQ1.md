# Synthetic length-frequency data I

Synthetic length-frequency data as provided in Sparre & Venema (1998).
Can be used to apply the function
[`Bhattacharya`](http://tokami.github.io/TropFishR/reference/Bhattacharya.md)
or [`ELEFAN`](http://tokami.github.io/TropFishR/reference/ELEFAN.md).

## Usage

``` r
data(synLFQ1)
```

## Format

A list consisting of:

- **dates** dates of sampling times (class Date),

- **midLengths** midpoints of the length classes,

- **catch** matrix with catches/counts per length class (row) and
  sampling date (column).

## Source

Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
assessment. Part 1. Manual. *FAO Fisheries Technical Paper*, (306.1,
Rev. 2). 407 p.

## Examples

``` r
data(synLFQ1)
str(synLFQ1)
#> List of 3
#>  $ dates     : Date[1:6], format: "1980-10-15" "1981-04-15" ...
#>  $ midLengths: num [1:40] 12.5 13.5 14.5 15.5 16.5 17.5 18.5 19.5 20.5 21.5 ...
#>  $ catch     : num [1:40, 1:6] 1 4 11 24 38 42 33 20 7 2 ...
#>   ..- attr(*, "dimnames")=List of 2
#>   .. ..$ : NULL
#>   .. ..$ : chr [1:6] "autumn1980" "spring1981" "autumn1981" "spring1982" ...
summary(synLFQ1)
#>            Length Class  Mode   
#> dates        6    Date   numeric
#> midLengths  40    -none- numeric
#> catch      240    -none- numeric

```
