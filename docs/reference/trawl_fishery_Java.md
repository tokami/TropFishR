# Data from the trawl fishery off the North coast of Java

Times series of catch and effort data from the trawl fishery off the
North coast of Java. This dataset can be used for the estimation of
maximum sustainable yield by means of the production models
([`prod_mod`](https://tokami.github.io/TropFishR/reference/prod_mod.md)
and
[`prod_mod_ts`](https://tokami.github.io/TropFishR/reference/prod_mod_ts.md)).

## Usage

``` r
data(trawl_fishery_Java)
```

## Format

A dataframe consisting of: 1. **year** a vector with years, 2. **Y**
yield \[1000 tons\], and 3. **f** fishing effort \[number of standard
vessels\].

## Source

Dwiponggo, A., 1979. Review of the demersal resources and fisheries in
the Java Sea. IPFC:RRD/II/79/Inf.12. Paper presented at the SCORRAD
Meeting, 1979, Hong Kong.

## Examples

``` r
data(trawl_fishery_Java)
str(trawl_fishery_Java)
#> 'data.frame':    9 obs. of  3 variables:
#>  $ year: num  1969 1970 1971 1972 1973 ...
#>  $ Y   : num  50 49 47.5 45 51 56 66 58 52
#>  $ f   : num  623 628 520 513 661 ...
summary(trawl_fishery_Java)
#>       year            Y               f         
#>  Min.   :1969   Min.   :45.00   Min.   : 513.0  
#>  1st Qu.:1971   1st Qu.:49.00   1st Qu.: 623.0  
#>  Median :1973   Median :51.00   Median : 661.0  
#>  Mean   :1973   Mean   :52.72   Mean   : 923.2  
#>  3rd Qu.:1975   3rd Qu.:56.00   3rd Qu.:1158.0  
#>  Max.   :1977   Max.   :66.00   Max.   :1970.0  
```
