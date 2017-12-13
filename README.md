TropFishR :fishing_pole_and_fish:
=====

[![License: GPL v3](https://img.shields.io/badge/License-GPL%20v3-blue.svg)](http://www.gnu.org/licenses/gpl-3.0)
[![Build Status](https://api.travis-ci.org/tokami/TropFishR.svg?branch=master)](https://travis-ci.org/tokami/TropFishR)
[![CRAN version](http://www.r-pkg.org/badges/version/TropFishR)](https://cran.r-project.org/package=TropFishR) 
[![CRAN Downloads](http://cranlogs.r-pkg.org/badges/TropFishR)](https://cran.r-project.org/package=TropFishR)
[![GitHub release](https://img.shields.io/github/release/tokami/TropFishR.svg)](https://github.com/tokami/TropFishR/releases)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.495176.svg)](https://doi.org/10.5281/zenodo.495176)

## Package description
   	
TropFishR is a collection of fisheries models based on the FAO Manual 
"Introduction to tropical fish stock assessment" by Sparre and Venema 
(1998, 1999). Not only scientists working in the tropics will benefit from 
this new toolbox. The methods work with age based or length-frequency data 
and assist in the assessment of data poor fish stocks. Overall, the package 
comes with 30 functions, 19 data sets and 10 s3 methods. All objects are 
documented and provide examples that allow reproducing the examples from 
the FAO manual. 

    
## News
You can find detailed descriptions of new features, bug fixes, other changes of specific package versions [here](https://rawgit.com/tokami/TropFishR/master/inst/doc/news.html).

     
## Installation
Download the released version from CRAN:

```R
install.packages(“TropFishR”)
```

Or the development version from github:

```R
# install.packages(devtools)
devtools::install_github(“tokami/TropFishR”)
```

## Citation
Please use the R command `citation("TropFishR")` to receive information on
how to cite this package.


## Vignettes
A [tutorial](https://rawgit.com/tokami/TropFishR/master/inst/doc/tutorial.html)
demonstrates the use of some of the main functions of TropFishR for a 
single-species stock assessment with length-frequency data. The [second vignette](https://rawgit.com/tokami/TropFishR/master/inst/doc/lfqData.html) gives a brief description of LFQ data and illustrates how files with raw length measurements (e.g. excel files) can be imported into R and trimmed for the use with TropFishR.


## Questions / Issues
In case you have questions or find bugs, please report on 
[TropFishR/issues](https://github.com/tokami/TropFishR/issues). If you want to be updated with the development of the package or want to discuss with TropFishR users and developers, follow the project on [ResearchGate](https://www.researchgate.net/project/TropFishR).


## References
  1. Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock 
  assessment. Part 1. Manual. FAO Fisheries Technical Paper, 
  (306.1, Rev. 2). 407 p. [link](http://www.fao.org/docrep/w5449e/w5449e00.htm)
  2. Sparre, P., Venema, S.C., 1999. Introduction to tropical fish stock 
  assessment. Part 2. Excercises. FAO Fisheries Technical Paper, 
  (306.2, Rev. 2). 94 p. [link](http://www.fao.org/docrep/w5448e/w5448e00.htm)
  3. Mildenberger, T. K., Taylor, M. H. and Wolff, M., 2017. TropFishR: an
  R package for fisheries analysis with length-frequency data. Methods in 
  Ecology and Evolution, 8: 1520-1527. doi:10.1111/2041-210X.12791 
  [link](https://doi.org/10.1111/2041-210X.12791)
  4. Taylor, M. H., and Mildenberger, T. K., 2017. Extending electronic 
  length frequency analysis in R. Fisheries Management and Ecology, 24:330-338. 
  doi:10.1111/fme.12232 [link](https://doi.org/10.1111/fme.12232)