# Shrimp data

Data from the Kuwait shrimp (*Penaeus semisulcatus*) fishery. They can
be used to simulate different management scenarios by applying
prediction models
([`predict_mod`](http://tokami.github.io/TropFishR/reference/predict_mod.md)).
Fishing mortalities were estimated from catch data and biomass estimates
were obtained from the swept area method. Recuritment to fishery at the
age of one month (Tr = 1).

## Usage

``` r
data(shrimps)
```

## Format

A list consisting of: 1. age groups \[months\], 2. meanWeight, 3.
meanValue, 4.FM, and 5.M

## Source

Garcia, S. and N.P. van Zalinge, 1982. Shrimp fishing in Kuwait:
methodology for a joint analysis of the artisanal and industrial
fisheries. pp. 119-142 In: Report on the Workshop on assessment of the
shrimp stocks of the west coast of the Gulf between Iran and the Arabian
Peninsula. Fisheries development in the Gulf. Rome, FAO,
FI:DP/RAB/80/015/1, 163 p.

## Examples

``` r
data(shrimps)
str(shrimps)
#> List of 5
#>  $ age       : num [1:13] 1 2 3 4 5 6 7 8 9 10 ...
#>  $ meanWeight: num [1:13] 5.7 9.3 13 17.6 22 26.1 30.3 33.8 37 40.3 ...
#>  $ meanValue : num [1:13] 0.73 0.93 1.2 1.45 1.7 1.9 2.08 2.14 2.18 2.23 ...
#>  $ FM        : num [1:13] 1.2 1.32 1.32 1.44 1.92 1.2 1.56 1.2 1.2 1.8 ...
#>  $ M         : num 3
summary(shrimps)
#>            Length Class  Mode   
#> age        13     -none- numeric
#> meanWeight 13     -none- numeric
#> meanValue  13     -none- numeric
#> FM         13     -none- numeric
#> M           1     -none- numeric

```
