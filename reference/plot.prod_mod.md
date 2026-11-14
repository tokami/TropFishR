# Plotting production models

This function plots CPUE and yield values against fishing effort
resulting from the production models
([prod_mod](http://tokami.github.io/TropFishR/reference/prod_mod.md)).

## Usage

``` r
# S3 method for class 'prod_mod'
plot(x, display_MSY = TRUE, ...)
```

## Arguments

- x:

  a object of the class `"prod_mod"`

- display_MSY:

  logical; should MSY be displayed in the graph?

- ...:

  optional parameters of plot function

## References

Fox, W. W. Jr., 1970. An exponential surplus-yield model for optimizing
exploited fish populations. *Trans.Am.Fish.Soc.*, 99:80-88

Graham, M., 1935. Modern theory of exploiting a fishery and application
to North Sea trawling. *J.Cons.CIEM*, 10(3):264-274

Schaefer, M., 1954. Some aspects of the dynamics of populations
important to the management of the commercial marine fisheries.
*Bull.I-ATTC/Bol. CIAT*, 1(2):27-56

Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
assessment. Part 1. Manual. FAO Fisheries Technical Paper, (306.1, Rev.
2). 407 p.

## Examples

``` r
data(trawl_fishery_Java)
output <-  prod_mod(data = trawl_fishery_Java)
plot(output, display_years = TRUE)

```
