# Emperor data

This dataset contains information about sky emperor (*Lethrinus
mahsena*) and its fisheries of offshore Mauritius banks (Nazareth
banks). It can be used for production models
([`prod_mod`](http://tokami.github.io/TropFishR/reference/prod_mod.md)
and
[`prod_mod_ts`](http://tokami.github.io/TropFishR/reference/prod_mod_ts.md)).

## Usage

``` r
data(emperor)
```

## Format

A dataframe containing: 1. **year** years, 2. **Y** yield, 3. **CPUE**
CPUE or abundance index

## Source

Dharmendra, D., Solmundsson, J., 2005. Stock assessment of the offshore
Mauritian banks using dynamic biomass models and analysis of length
frequency of the Sky Emperor (*Lethrinus mahsena*). Fisheries Training
Program The United Nations University, 61

## Examples

``` r
data(emperor)
str(emperor)
#> 'data.frame':    16 obs. of  3 variables:
#>  $ year: int  1989 1990 1991 1992 1993 1994 1995 1996 1997 1998 ...
#>  $ Y   : num  837 914 793 952 1358 ...
#>  $ CPUE: num  75.7 78.5 81.3 78.2 66.1 66.5 64.3 52.5 66.1 81.4 ...
summary(emperor)
#>       year            Y               CPUE      
#>  Min.   :1989   Min.   : 468.0   Min.   :52.50  
#>  1st Qu.:1993   1st Qu.: 899.2   1st Qu.:66.40  
#>  Median :1996   Median :1083.0   Median :77.20  
#>  Mean   :1996   Mean   :1109.2   Mean   :76.67  
#>  3rd Qu.:2000   3rd Qu.:1360.0   3rd Qu.:82.10  
#>  Max.   :2004   Max.   :1720.0   Max.   :99.20  

```
