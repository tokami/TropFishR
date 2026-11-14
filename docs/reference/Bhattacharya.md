# Bhattacharya's method

Find relative frequencies and frequency distribution of distinct cohorts
in the observed length frequency distribution by resolving it into
Gaussian components.

## Usage

``` r
Bhattacharya(param, n_rnorm = 1000, savePlots = FALSE)
```

## Arguments

- param:

  a list consisting of following parameters:

  - `midLengths` midpoints of the length class as vector,

  - `catch` a vector with the catch per length class or a matrix with
    catches per length class of subsequent years;

- n_rnorm:

  number of fish of a cohort used to define its length range: the range
  spans from the expected smallest to the expected largest of `n_rnorm`
  lengths from the normal distribution of the cohort (Blom's
  approximation, quantiles `0.625/(n_rnorm+0.25)` and
  `(n_rnorm-0.375)/(n_rnorm+0.25)`). The numbers of the cohort are
  calculated up to the largest length. The default is 1000.

- savePlots:

  logical; indicating whether the analyis graphs should be recorded

## Value

A list with the input parameters and

- **regression_lines** dataframe with intercept, slope, start and end
  points of the regression lines,

- **Lmean_SD_list** dataframe with the mean length (Lmean), standard
  deviation (SD), and seperation index (SI) for each cohort,

- **bhat_results** dataframe with the results of the Bhattacharya
  method,

- **distributions** list with the x and y values of selected
  distributions,

- **cohort_plots** list with analysis plots (when savePlots = TRUE).

## Details

This method includes the
[`identify`](https://rdrr.io/r/graphics/identify.html) function, which
allows to choose points on a graphical device manually. The points of
each regression line have to lie on a falling line (negative slope),
otherwise the selection is rejected. Start with the distribution on the
left: the points of a distribution only lie on a falling line after the
distributions to the left of it have been removed. To stop this process
please press right mouse click on the graph, and in case you are working
on a windows maschine click on "Stop". An error will be caused if the
graphical device is closed manually. After you have selected the points
for regression analysis you will be asked if you want to redo the
selection or if you want to continue. Please enter in the Console "y"
for continuation if you are satisfied with your selection and the
corresponding Gaussian distribution or enter "redo" if you want to
repeat the selection procedure. This function allows a maximum of 12
cohorts or seperate distributions in one sample. Please find more
details in the Vignette of this package or in the FAO manual by Sparre
and Venema (1998) (**Section 3.4.1, p. 80**).

## References

Bhattacharya, C.G., 1967. A simple method of resolution of a
distribution into Gaussian components, *Biometrics*, 23:115-135

Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
assessment. Part 1. Manual. *FAO Fisheries Technical Paper*, (306.1,
Rev. 2). 407 p.

## Examples

``` r
# The following example requires to choose the points for the regression
# lines one cohort after the other in the same run, starting on the left
# (each selection requires that the cohorts to the left are already removed):
# first cohort:   point  2 to  8, then type 'y' in the console
# second cohort:  point 12 to 17, then type 'y'
# third cohort:   point 19 to 23, then type 'y'
# fourth cohort:  point 26 to 30, then type 'y'
# and stop the selection (right click or Escape).

# \donttest{
 data(synLFQ1)
 Bhattacharya(param = synLFQ1)
#> Interactive session needed for Bhattacharya.
#> NULL
# }
```
