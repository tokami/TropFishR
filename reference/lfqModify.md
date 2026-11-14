# Modify lfq data for further analysis

Modify length-frequency (LFQ) data. Allows to summarise catch matrix of
LFQ data to one column per year. This is required for e.g.
[`catchCurve`](http://tokami.github.io/TropFishR/reference/catchCurve.md).
Allows to change bin size of LFQ data. Allows to add plus group to catch
matrix.

## Usage

``` r
lfqModify(
  lfq,
  par = NULL,
  bin_size = NA,
  aggregate = NA,
  vectorise_catch = FALSE,
  plus_group = FALSE,
  minDate = NA,
  maxDate = NA,
  years = NA,
  Lmin = NA,
  Lmax = NA,
  lfq2 = NA
)
```

## Arguments

- lfq:

  lfq object with dates, midLengths, and catch

- par:

  growth parameters as resulting from e.g.
  [`ELEFAN`](http://tokami.github.io/TropFishR/reference/ELEFAN.md)

- bin_size:

  Bin size for length frequencies (in cm)

- aggregate:

  Factor to aggregate catch per year (`"year"`), per quarter
  (`"quarter"`), or per month (`"month"`). By default data is not
  aggregated (`NA`).

- vectorise_catch:

  logical; indicating if the catch matrix should be summarised to yearly
  vectors (default: FALSE).

- plus_group:

  logical or numeric; should a plus group be created? If yes you will be
  asked to insert the length for the plus group in the console (default:
  FALSE). Instead of inserting the length of the plus group via the
  console, the value can be inserted, e.g. plus_group = 85.5.

- minDate:

  minimum date to subset lfq data

- maxDate:

  maximum date to subset lfq data

- years:

  numeric with year(s) to subset lfq data

- Lmin:

  minimum length to subset lfq data

- Lmax:

  maximum length to subset lfq data

- lfq2:

  optional second lfq object which will be merged with lfq. This might
  be interesting for fleet specific lfq objects. Default: NA. Be aware
  that catches are combined without weighting!

## Value

lfq object with rearranged catch matrix (yearly sums) and growth
parameters if provided.

## Examples

``` r
data(synLFQ4)

## summarise catch matrix per year
lfq_sum <- lfqModify(synLFQ4, vectorise_catch = TRUE)

## change bin size
lfq_bin <- lfqModify(synLFQ4, bin_size = 4)

## add plus_group
lfq_plus <- lfqModify(synLFQ4, plus_group = 85.5)
```
