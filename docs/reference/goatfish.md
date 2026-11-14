# Yellowstriped goatfish data

Data of Yellowstriped goatfish (*Upeneus vittatus*) from Manila Bay,
Philippines. Can be used for the estimation of the instantaneous
mortality rate (Z) by means of
[`catchCurve`](https://tokami.github.io/TropFishR/reference/catchCurve.md).

## Usage

``` r
data(goatfish)
```

## Format

A list consisting of:

- `midLengths`: mid points of length classes,

- `catch`: a vector with catches in numbers,

- `Linf`: infinite length in cm \[cm\],

- `K`: growth coefficent per year \[1/year\].

## Source

Ziegler, B., 1979. Growth and mortality rates of some fishes of Manila
Bay, Philippines, as estimated from analysis of length-frequencies.
Thesis. Kiel University, 115 p.

## Examples

``` r
data(goatfish)
str(goatfish)
#> List of 4
#>  $ midLengths: num [1:18] 6.5 7.5 8.5 9.5 10.5 11.5 12.5 13.5 14.5 15.5 ...
#>  $ catch     : num [1:18] 3 143 271 318 416 488 614 613 493 278 ...
#>  $ Linf      : num 23.1
#>  $ K         : num 0.59
summary(goatfish)
#>            Length Class  Mode   
#> midLengths 18     -none- numeric
#> catch      18     -none- numeric
#> Linf        1     -none- numeric
#> K           1     -none- numeric

```
