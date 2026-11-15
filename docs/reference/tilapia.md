# Tilapia data

A list with characteristics of experimental catches of tilapia *Tilapia
esculenta* with gillnets of two mesh sizes. The Results of the
experiment are the midlengths of size classes, the number of fish caught
with net 1 & 2, and the meshsizes of both nets. The data can be analysed
with the function
[`select`](http://tokami.github.io/TropFishR/reference/select.md).

## Usage

``` r
data(tilapia)
```

## Format

A list consiting of:

- `midLengths` the midlengths of size classes,

- `numNet1` the number of fish caught with net 1,

- `numNet2` the number of fish caught with net 2,

- `msNet1` the meshsize of net 1,

- `msNet2` and the meshsize of net 2.

## Source

Garrod, D.J., 1961. The selection characteristics of nylon gill nets for
*Tilapia esculenta* Graham. J.Cons.CIEM, 26(2):191-203

## Examples

``` r
data(tilapia)
str(tilapia)
#> List of 4
#>  $ midLengths     : num [1:10] 18.5 19.5 20.5 21.5 22.5 23.5 24.5 25.5 26.5 27.5
#>  $ type           : chr "gillnet"
#>  $ meshSizes      : num [1:2] 8.1 9.1
#>  $ CatchPerNet_mat: num [1:10, 1:2] 7 90 199 182 119 29 17 3 0 0 ...
summary(tilapia)
#>                 Length Class  Mode     
#> midLengths      10     -none- numeric  
#> type             1     -none- character
#> meshSizes        2     -none- numeric  
#> CatchPerNet_mat 20     -none- numeric  
```
