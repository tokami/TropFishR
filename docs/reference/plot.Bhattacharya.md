# Bhattacharya plot

This function plots the seperated frequency distributions and selected
regression lines of
[Bhattacharya](https://tokami.github.io/TropFishR/reference/Bhattacharya.md)
method.

## Usage

``` r
# S3 method for class 'Bhattacharya'
plot(x, analysisPlot = TRUE, ...)
```

## Arguments

- x:

  a list of the class `"Bhattacharya"` containing the results of
  Bhattacharya's method.

- analysisPlot:

  logical; indicating wheter the anaylsis graph with the regression
  lines should be created

- ...:

  additional options of the
  [`plot`](https://rdrr.io/r/graphics/plot.default.html) function

## Details

This function plots the results of the Bhattacharya method.

## References

Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
assessment. Part 1. Manual. *FAO Fisheries Technical Paper*, (306.1,
Rev. 2). 407 p.

## Examples

``` r
if (FALSE) { # \dontrun{
 data(synLFQ1)
 output <- Bhattacharya(param = synLFQ1)
 plot(output)
} # }
```
