# Add phi prime contours to a Linf-K plot

Adds isolines of the growth performance index \\\phi' = \log\_{10}(K) +
2 \log\_{10}(L\_\infty)\\ to an existing plot with Linf on the x-axis
and K on the y-axis (e.g.
[LinfK_scatterhist](http://tokami.github.io/TropFishR/reference/LinfK_scatterhist.md)).

## Usage

``` r
add_phiprime(gridsize = 20, ...)
```

## Arguments

- gridsize:

  default = 20

- ...:

  Additional arguments passed to
  [`contour`](https://rdrr.io/r/graphics/contour.html).

## Value

plot

## Examples

``` r
plot(c(5, 20), c(0.1, 2), t = "n", xlab = "Linf", ylab = "K")
add_phiprime()

```
