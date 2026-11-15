# Raise length-frequency samples to the total catch

The function raises the sampled numbers of fish per length class
(\`lfq\$catch\`) to the total catch, using the total landings in weight
and the weight of the samples. This is required for methods that need
absolute numbers, such as the length-based virtual population analysis
([`VPA`](http://tokami.github.io/TropFishR/reference/VPA.md)). Methods
that only use the relative length composition (e.g. ELEFAN,
[`catchCurve`](http://tokami.github.io/TropFishR/reference/catchCurve.md))
and the resampling
([`lfqResample`](http://tokami.github.io/TropFishR/reference/lfqResample.md))
work with the sampled numbers.

## Usage

``` r
lfqRaise(
  lfq,
  landings,
  sampleWeight = NULL,
  by = "sample",
  a = NULL,
  b = NULL,
  LW_unit = "kg"
)
```

## Arguments

- lfq:

  a list of the class "lfq" consisting of following parameters:

  - **midLengths** midpoints of the length classes,

  - **dates** dates of sampling times (class Date),

  - **catch** matrix with the sampled numbers of fish per length class
    (rows) and sampling date (columns),

  - **par** (optional) a list with the parameters of the length-weight
    relationship (`a` and `b`), required if `sampleWeight` is not
    provided.

- landings:

  total landings (or catch) in kg: one value per sample
  (`by = "sample"`) or per year (`by = "year"`). With `by = "year"`, the
  values are matched to the years by their names if named, otherwise
  they are used in the order of the years.

- sampleWeight:

  optional; the weight of the samples in kg (one value per sample or
  year, as `landings`). By default (`sampleWeight = NULL`), the weight
  is estimated from the sampled numbers and the length-weight
  relationship: sum of `catch * a * midLengths^b`. Use the measured
  weight of the samples if available.

- by:

  `"sample"` (default) to raise each sample (column of `lfq$catch`)
  separately, or `"year"` to raise all samples of a year with one factor
  (landings of the year / weight of all samples of the year).

- a, b:

  optional; parameters of the length-weight relationship (W = a \* L^b).
  By default taken from `lfq$par`.

- LW_unit:

  unit of the weight in the length-weight relationship (`a`): `"kg"`
  (default) or `"g"`.

## Value

The lfq list with the raised catch (numbers of fish) and the element
`raising`, a data frame with the landings, the sample weight, whether
the sample weight was estimated, and the raising factor per sample or
year.

## Details

The raising factor of each sample (or year) is the landings divided by
the weight of the sample(s). Raising each sample separately
(`by = "sample"`) gives each sample the weight of its landings, e.g.
months with large landings get more weight than months with small
landings, independent of the number of fish measured.

Keep the original lfq data with the sampled numbers for all other
analyses. In particular, the bootstrap functions
([`lfqResample`](http://tokami.github.io/TropFishR/reference/lfqResample.md),
[`ELEFAN_GA_boot`](http://tokami.github.io/TropFishR/reference/ELEFAN_GA_boot.md),
...) require the sampled numbers, because they resample the measured
fish. Instead of raising the data beforehand, the landings can be passed
directly to [`VPA`](http://tokami.github.io/TropFishR/reference/VPA.md)
(argument `landings`), which raises the data (or each resample with
`boot`) internally.

## Examples

``` r
data(synLFQ7)
## length-weight relationship (W in g)
synLFQ7$par <- list(a = 0.015, b = 3)
## landings in kg at the sampling dates (one value per sample)
landings <- seq(20000, 40000, length.out = length(synLFQ7$dates))
lfqR <- lfqRaise(synLFQ7, landings = landings, LW_unit = "g")
lfqR$raising
#>       sample landings sampleWeight estimated   factor
#> 1 2016-01-01 20000.00     14836.23      TRUE 1.348051
#> 2 2016-03-02 22857.14     15289.33      TRUE 1.494973
#> 3 2016-04-02 25714.29     15659.33      TRUE 1.642106
#> 4 2016-07-02 28571.43     15974.05      TRUE 1.788615
#> 5 2016-08-01 31428.57     16554.03      TRUE 1.898545
#> 6 2016-10-01 34285.71     14380.58      TRUE 2.384167
#> 7 2016-11-01 37142.86     14515.66      TRUE 2.558813
#> 8 2016-12-01 40000.00     14859.87      TRUE 2.691813

## one raising factor per year
lfqY <- lfqRaise(synLFQ7, landings = c("2016" = 240000), by = "year",
                 LW_unit = "g")
lfqY$raising
#>   year landings sampleWeight estimated factor
#> 1 2016   240000     122069.1      TRUE 1.9661
```
