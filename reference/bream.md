# bream data

Data of a covered codend experimental catch of the species Threadfin
bream (*Nemipterus japonicus*) in South China Sea. This data can be
analysed with the function
[`select`](http://tokami.github.io/TropFishR/reference/select.md).

## Usage

``` r
data(bream)
```

## Format

A list consiting of:

- `midLengths` the midlengths of size classes,

- `numCodend` the number of fish retained in codend,

- `numCover` the number of fish retained in cover,

- `msCodend` the meshsize of codend.

## Source

Jones, R., 1976. Mesh regulation in the demersal fisheries of the South
China Sea area. Manila, South China Sea Fisheries Development and
Coordinating Programme, SCS/76/WP/34: 75 p.

## Examples

``` r
data(bream)
str(bream)
#> List of 4
#>  $ midLengths     : num [1:10] 9.5 10.5 11.5 12.5 13.5 14.5 15.5 16.5 17.5 18.5
#>  $ type           : chr "trawl_net"
#>  $ meshSizes      : num [1:2] 0 4
#>  $ CatchPerNet_mat: num [1:10, 1:2] 1 6 7 4 5 13 8 3 0 1 ...
summary(bream)
#>                 Length Class  Mode     
#> midLengths      10     -none- numeric  
#> type             1     -none- character
#> meshSizes        2     -none- numeric  
#> CatchPerNet_mat 20     -none- numeric  
```
