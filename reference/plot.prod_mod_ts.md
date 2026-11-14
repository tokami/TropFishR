# Plotting time series production models

This function plots objects of the class "prod_mod_ts".

## Usage

``` r
# S3 method for class 'prod_mod_ts'
plot(x, correlation_plots = FALSE, ...)
```

## Arguments

- x:

  a object of the class "prod_mod_ts",

- correlation_plots:

  logical; indicating if correlation plots should be displayed

- ...:

  additional parameters of the
  [plot](https://rdrr.io/r/graphics/plot.default.html) function

## References

Dharmendra, D., Solmundsson, J., 2005. Stock assessment of the offshore
Mauritian banks using dynamic biomass models and analysis of length
frequency of the Sky Emperor (*Lethrinus mahsena*). Fisheries Training
Program The United Nations University, 61

## Examples

``` r
data(emperor)
output <-  prod_mod_ts(emperor, method  = "Fox")

plot(output, correlation_plots = TRUE)

```
