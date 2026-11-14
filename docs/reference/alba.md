# Length-frequency data of the clam Abra alba

Length-frequency data of the clam Abra alba as presented by Brey et al.
(1988). Includes 7 approximately bi-monthly samplings of A. alba binned
into 14 length classes. The data is of class "lfq", which can be used
used e.g. in functions estimating growth parameters:
[`ELEFAN`](https://tokami.github.io/TropFishR/reference/ELEFAN.md),
[`ELEFAN_SA`](https://tokami.github.io/TropFishR/reference/ELEFAN_SA.md)
or
[`ELEFAN_GA`](https://tokami.github.io/TropFishR/reference/ELEFAN_GA.md).

## Usage

``` r
data(alba)
```

## Format

A list of class lfq consisting of:

- **dates** dates of sampling times (class Date),

- **midLengths** midpoints of the length classes,

- **catch** matrix with catches/counts per length class (row) and
  sampling date (column).

## Source

Brey, T., Soriano, M., and Pauly, D. 1988. Electronic length frequency
analysis: a revised and expanded user's guide to ELEFAN 0, 1 and 2.
ICLARM Contrib. 261.

## Examples

``` r
data(alba)

# plot raw catch frequencies
plot(alba, Fname = "catch")


# plot restructured frequencies
alba <- lfqRestructure(alba, MA=5)
plot(alba, Fname = "rcounts")


# growth parameters: see the examples of ELEFAN_SA and ELEFAN_GA
```
