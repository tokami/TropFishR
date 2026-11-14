# TropFishR

[![License: GPL
v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](https://www.gnu.org/licenses/gpl-3.0)
[![CRAN
version](https://www.r-pkg.org/badges/version/TropFishR)](https://cran.r-project.org/package=TropFishR)
[![CRAN
Downloads](https://cranlogs.r-pkg.org/badges/TropFishR)](https://cran.r-project.org/package=TropFishR)
[![GitHub
release](https://img.shields.io/github/release/tokami/TropFishR.svg)](https://github.com/tokami/TropFishR/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.495176.svg)](https://doi.org/10.5281/zenodo.495176)

TropFishR is a collection of fish stock assessment methods for
length-frequency and age-based data, aimed at data-poor fisheries. It
implements the methods of the FAO manual “Introduction to tropical fish
stock assessment” by Sparre and Venema (1998, 1999), and the examples in
the help files reproduce the manual’s worked examples. In addition, it
includes newer length-based methods, such as electronic length-frequency
analysis with simulated annealing and a genetic algorithm (`ELEFAN_SA`,
`ELEFAN_GA`; Taylor and Mildenberger 2017) and bootstrapped ELEFAN for
the uncertainty of the growth parameters (`ELEFAN_GA_boot`,
`ELEFAN_SA_boot`; Schwamborn et al. 2019).

**Website:** <https://tokami.github.io/TropFishR/>

## Installation

Install the released version from CRAN:

``` r
install.packages("TropFishR")
```

Or the development version from GitHub:

``` r
# install.packages("remotes")
remotes::install_github("tokami/TropFishR")
```

## Getting started

Functions take a length-frequency data object (class `lfq`), add their
results to it and return it, so that an assessment can be built up step
by step: growth (e.g. `ELEFAN_GA`), mortality (`M_empirical`,
`catchCurve`), stock size (`VPA`) and yield per recruit (`predict_mod`).
For example, to estimate the growth parameters of a synthetic data set:

``` r
library(TropFishR)
data(synLFQ7)

lfq <- lfqModify(synLFQ7, bin_size = 4)
res <- ELEFAN_GA(lfq, MA = 7, seasonalised = TRUE, maxiter = 50)
res$par
plot(res, Fname = "rcounts")
```

The articles on the website explain each step in detail:

- [Single-species fish stock
  assessment](https://tokami.github.io/TropFishR/articles/tutorial.html):
  a complete length-based assessment with the main functions of
  TropFishR.
- [Length-frequency
  data](https://tokami.github.io/TropFishR/articles/lfqData.html): what
  LFQ data are and how to import raw length measurements (e.g. from
  spreadsheets) into TropFishR.
- [ELEFAN
  functions](https://tokami.github.io/TropFishR/articles/Using_TropFishR_ELEFAN_functions.html):
  the ELEFAN methods in detail, with recommendations for their use.
- [Bootstrapped
  ELEFAN](https://tokami.github.io/TropFishR/articles/bootstrap_ELEFAN.html):
  how to estimate the uncertainty of the growth parameters and propagate
  it through the assessment.

Changes between versions are listed in the
[changelog](https://tokami.github.io/TropFishR/news/).

## Citation

Please cite TropFishR as described by `citation("TropFishR")`.

## Questions and issues

Please report bugs and ask questions on [GitHub
issues](https://github.com/tokami/TropFishR/issues).

## References

1.  Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
    assessment. Part 1. Manual. FAO Fisheries Technical Paper, (306.1,
    Rev. 2). 407p.
    [link](https://openknowledge.fao.org/server/api/core/bitstreams/bc7c37b6-30df-49c0-b5b4-8367a872c97e/content)
2.  Sparre, P., Venema, S.C., 1999. Introduction to tropical fish stock
    assessment. Part 2. Exercises. FAO Fisheries Technical Paper,
    (306.2, Rev. 2). 94p.
    [link](https://www.fao.org/4/W5448E/W5448E00.htm)
3.  Mildenberger, T. K., Taylor, M. H., Wolff, M., 2017. TropFishR: an R
    package for fisheries analysis with length-frequency data. Methods
    in Ecology and Evolution, 8: 1520-1527.
    [doi:10.1111/2041-210X.12791](https://doi.org/10.1111/2041-210X.12791)
4.  Taylor, M. H., Mildenberger, T. K., 2017. Extending electronic
    length frequency analysis in R. Fisheries Management and Ecology,
    24: 330-338.
    [doi:10.1111/fme.12232](https://doi.org/10.1111/fme.12232)
5.  Schwamborn, R., Mildenberger, T. K., Taylor, M. H., 2019. Assessing
    sources of uncertainty in length-based estimates of body growth in
    populations of fishes and macroinvertebrates with bootstrapped
    ELEFAN. Ecological Modelling, 393: 37-51.
    [doi:10.1016/j.ecolmodel.2018.12.001](https://doi.org/10.1016/j.ecolmodel.2018.12.001)
