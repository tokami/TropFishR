# Stock simulation

This function estimates stock size, biomass and yield of a stock from
fishing mortality per age class or length group. This function is
embedded in the Thompson and Bell model (prediction model:
[`predict_mod`](http://tokami.github.io/TropFishR/reference/predict_mod.md)).

## Usage

``` r
stock_sim(param, age_unit = "year", stock_size_1 = NA, plus_group = NA)
```

## Arguments

- param:

  a list consisting of following parameters:

  - `age` or `midLengths`: midpoints of length classes (length-frequency
    data) or ages (age composition data),

  - `meanWeight`: mean weight in kg per age class or length group,

  - `meanValue`: mean value per kg fish per age class or length group,

  - `FM`: fishing mortality rates per age class or length group,

  - `M` or `Z`: natural or total instantaneous mortality rate.

- age_unit:

  indicates if the age groups are per month (`"month"`) or per year
  (`"year"`). Default: `"year"`

- stock_size_1:

  stock size of smallest age/length group

- plus_group:

  indicates the age group (row index), which should be turned into a
  plus group (i.e. all groups above are comprised in one group). Only
  for age-based data; for length-based data, the length class including
  Linf (or the last length class) is the plus group, and length classes
  above Linf are `NA`.

## Value

A list with the input parameters and following list objects:

- **dt**: delta t,

- **N**: population numbers,

- **dead**: number of deaths due to natural mortality,

- **C**: catch,

- **Y**: yield,

- **B**: biomass,

- **V**: value,

- **totals**: summarised output:

  - **totC** total catch,

  - **totY** total yield,

  - **totV** total value,

  - **meanB** mean biomass.

  ,

## Details

better to treat last group always as a plus group... if stock size 1 not
provided assumes 1000 as intital population size make sure that FM is
also in same unit as the classes, e.g. when classes in months than also
FM has to be provided in 1/months

## References

Garcia, S. and N.P. van Zalinge, 1982. Shrimp fishing in Kuwait:
methodology for a joint analysis of the artisanal and industrial
fisheries. pp. 119-142 In: Report on the Workshop on assessment of the
shrimp stocks of the west coast of the Gulf between Iran and the Arabian
Peninsula. Fisheries development in the Gulf. Rome, FAO,
FI:DP/RAB/80/015/1, 163 p.

Millar, R. B., & Holst, R. (1997). Estimation of gillnet and hook
selectivity using log-linear models. *ICES Journal of Marine Science:
Journal du Conseil*, 54(3):471-477

Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
assessment. Part 1. Manual. *FAO Fisheries Technical Paper*, (306.1,
Rev. 2). 407 p.

## Examples

``` r
# age-based stock simulation
data(shrimps)

# option 1: without plus group
stock_sim(shrimps, age_unit = "month")
#> $age
#>  [1]  1  2  3  4  5  6  7  8  9 10 11 12 13
#> 
#> $meanWeight
#>  [1]  5.70  9.30 13.00 17.60 22.00 26.10 30.30 33.80 37.00 40.30 43.11 44.70
#> [13]    NA
#> 
#> $meanValue
#>  [1] 0.73 0.93 1.20 1.45 1.70 1.90 2.08 2.14 2.18 2.23 2.24 2.27   NA
#> 
#> $FM
#>  [1] 1.20 1.32 1.32 1.44 1.92 1.20 1.56 1.20 1.20 1.80 2.76 2.52   NA
#> 
#> $M
#> [1] 3
#> 
#> $dt
#>  [1] 0.08333333 0.08333333 0.08333333 0.08333333 0.08333333 0.08333333
#>  [7] 0.08333333 0.08333333 0.08333333 0.08333333 0.08333333 0.08333333
#> [13]         NA
#> 
#> $N
#>  [1] 1000.000000  704.688090  491.644197  343.008517  236.927759  157.237166
#>  [7]  110.803158   75.774004   53.397038   37.628257   25.222975   15.607558
#> [13]    9.852796
#> 
#> $dead
#>  [1] 295.311910 213.043892 148.635680 106.080759  79.690592  46.434008
#>  [7]  35.029154  22.376966  15.768781  12.405282   9.615417   5.754762
#> [13]         NA
#> 
#> $C
#>  [1] 84.374832 65.096745 45.416458 34.404570 31.098768 13.266859 11.983658
#>  [8]  6.393419  4.505366  4.651981  4.607387  2.627174        NA
#> 
#> $Y
#>  [1] 480.9365 605.3997 590.4140 605.5204 684.1729 346.2650 363.1048 216.0976
#>  [9] 166.6985 187.4748 198.6245 117.4347       NA
#> 
#> $B
#>  [1] 4809.3654 5503.6339 5367.3996 5046.0037 4276.0806 3462.6503 2793.1141
#>  [8] 2160.9756 1666.9855 1249.8322  863.5846  559.2127        NA
#> 
#> $V
#>  [1]  351.0837  563.0217  708.4967  878.0046 1163.0939  657.9036  755.2581
#>  [8]  462.4488  363.4028  418.0689  444.9188  266.5767        NA
#> 
#> $totals
#>       totC     totY     totV   meanB
#> 1 308.4272 4562.143 7032.278 3146.57
#> 

# option 2: with plus group
stock_sim(param = shrimps, age_unit = "month", plus_group = 11)
#> $age
#>  [1] "1"   "2"   "3"   "4"   "5"   "6"   "7"   "8"   "9"   "10"  "11+"
#> 
#> $meanWeight
#>  [1]  5.70  9.30 13.00 17.60 22.00 26.10 30.30 33.80 37.00 40.30 43.11
#> 
#> $meanValue
#>  [1] 0.73 0.93 1.20 1.45 1.70 1.90 2.08 2.14 2.18 2.23 2.24
#> 
#> $FM
#>  [1] 1.20 1.32 1.32 1.44 1.92 1.20 1.56 1.20 1.20 1.80 2.76
#> 
#> $M
#>  [1] 3 3 3 3 3 3 3 3 3 3 3
#> 
#> $dt
#>  [1] 0.08333333 0.08333333 0.08333333 0.08333333 0.08333333 0.08333333
#>  [7] 0.08333333 0.08333333 0.08333333 0.08333333 0.08333333
#> 
#> $N
#>  [1] 1000.00000  704.68809  491.64420  343.00852  236.92776  157.23717
#>  [7]  110.80316   75.77400   53.39704   37.62826   25.22297
#> 
#> $dead
#>  [1] 295.31191 213.04389 148.63568 106.08076  79.69059  46.43401  35.02915
#>  [8]  22.37697  15.76878  12.40528  25.22297
#> 
#> $C
#>  [1] 84.374832 65.096745 45.416458 34.404570 31.098768 13.266859 11.983658
#>  [8]  6.393419  4.505366  4.651981 12.086009
#> 
#> $Y
#>  [1] 480.9365 605.3997 590.4140 605.5204 684.1729 346.2650 363.1048 216.0976
#>  [9] 166.6985 187.4748 521.0278
#> 
#> $B
#>  [1] 4809.365 5503.634 5367.400 5046.004 4276.081 3462.650 2793.114 2160.976
#>  [9] 1666.985 1249.832 2265.338
#> 
#> $V
#>  [1]  351.0837  563.0217  708.4967  878.0046 1163.0939  657.9036  755.2581
#>  [8]  462.4488  363.4028  418.0689 1167.1024
#> 
#> $totals
#>       totC     totY     totV    meanB
#> 1 313.2787 4767.112 7487.885 3509.216
#> 

# length-based stock simulation
data(hake)

stock_sim(param = hake, stock_size_1 = 98919.3)
#> $midLengths
#>  [1]  9 15 21 27 33 39 45 51 57 63 69 75 81 87
#> 
#> $catch
#>  [1]  1823 14463 25227  8134  3889  2959  1871   653   322   228   181    96
#> [13]    16    46
#> 
#> $K
#> [1] 0.1
#> 
#> $Linf
#> [1] 130
#> 
#> $M
#> [1] 0.28
#> 
#> $a
#> [1] 1e-05
#> 
#> $b
#> [1] 3
#> 
#> $FM
#>  [1] 0.04 0.39 1.07 0.65 0.49 0.59 0.65 0.39 0.29 0.31 0.40 0.39 0.11 0.28
#> 
#> $meanValue
#>  [1] 1.0 1.0 1.0 1.5 1.5 2.0 2.0 2.5 2.5 2.5 3.0 3.0 3.0 3.0
#> 
#> $N
#>  [1] 98919.30000 84400.27953 59398.43096 27519.82078 15877.10223  9812.31928
#>  [7]  5476.83521  2798.52315  1673.30744  1043.61635   611.97048   309.54971
#> [13]   146.59769    90.72902
#> 
#> $dead
#>  [1] 14519.02047 25001.84857 31878.61018 11642.71855  6064.78295  4335.48407
#>  [7]  2678.31205  1125.21572   629.69109   431.64587   302.42077   162.95202
#> [13]    55.86866          NA
#> 
#> $C
#>  [1]  1814.87756 14553.31484 25266.75029  8137.38393  3859.40733  2940.15587
#>  [7]  1871.93853   654.97631   320.36915   226.79698   177.89457    94.85267
#> [13]    15.75783    45.36451
#> 
#> $Y
#>  [1]   13.23046  491.17438 2339.95374 1601.68128 1386.95521 1744.07106
#>  [7] 1705.80399  868.83263  593.30124  567.09905  584.39967  400.15971
#> [13]   83.74356  555.73476
#> 
#> $B
#>  [1]  330.7614 1259.4215 2186.8727 2464.1250 2830.5208 2956.0526 2624.3138
#>  [8] 2227.7760 2045.8663 1829.3518 1460.9992 1026.0505  761.3051 1984.7670
#> 
#> $SSB
#>  [1] NA NA NA NA NA NA NA NA NA NA NA NA NA NA
#> 
#> $V
#>  [1]   13.23046  491.17438 2339.95374 2402.52192 2080.43282 3488.14212
#>  [7] 3411.60797 2172.08157 1483.25310 1417.74763 1753.19901 1200.47912
#> [13]  251.23068 1667.20429
#> 
#> $totals
#>       totC     totY     totV    meanB meanSSB
#> 1 59979.84 12936.14 24172.26 25988.18       0
#> 
```
