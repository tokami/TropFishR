# Empirical formulas for the estimation of natural mortality

Functions to calculate the instantaneous natural mortality rate (M)
according to 12 different empirical formulas.

## Usage

``` r
M_empirical(
  Linf = NULL,
  Winf = NULL,
  K_l = NULL,
  K_w = NULL,
  temp = NULL,
  tmax = NULL,
  tm50 = NULL,
  GSI = NULL,
  Wdry = NULL,
  Wwet = NULL,
  Bl = NULL,
  schooling = FALSE,
  method,
  boot = NULL,
  CI = 95
)
```

## Source

https://cran.r-project.org/web/packages/fishmethods/index.html

## Arguments

- Linf:

  infinite total length (TL) from a von Bertalanffy growth curve in cm.

- Winf:

  infinite weight form a von Bertalanffy growth curve in wet
  weight-grams.

- K_l:

  is the growth coefficient (per year) from a von Bertalanffy growth
  curve for length.

- K_w:

  is the growth coefficient (per year) from a von Bertalanffy growth
  curve for weight.

- temp:

  average annual temperature at the surface in degrees centigrade.

- tmax:

  the oldest age observed for the species.

- tm50:

  age when 50% of the population is mature \[year\] ("age of massive
  maturation").

- GSI:

  gonadosomatic index (wet ovary weight over wet body weight).

- Wdry:

  total dry weight in grams.

- Wwet:

  total wet weight at mean length in grams.

- Bl:

  vector with body lengths in cm for size dependent mortality estimates
  (method = "Gislason")

- schooling:

  logical; if TRUE it is accounted for the schooling behaviour of the
  species, only for Pauly's methods. Default is FALSE.

- method:

  vector of method names. Any combination of following methods can be
  employed: "AlversonCarney", "Gislason" (size dependent mortality
  estimates), "GundersonDygert", "Hoenig", "Lorenzen", "Pauly_Linf",
  "Pauly_Winf", "PetersonWroblewski", "RikhterEfanov", "Roff",
  "Then_growth", or "Then_tmax". Please refer to Details to see which
  input parameters are required by each method.

- boot:

  optional; an object of class \`lfqBoot\` with bootstrapped growth
  parameters (e.g. from
  [`ELEFAN_GA_boot`](http://tokami.github.io/TropFishR/reference/ELEFAN_GA_boot.md)).
  If provided, M is estimated for each resample with its Linf and K.
  Only the methods "Then_growth" and "Pauly_Linf" (with `temp`) are
  available. See also
  [`MstochThen`](http://tokami.github.io/TropFishR/reference/MstochThen.md).

- CI:

  percentage for confidence intervals if `boot` is provided (Default:
  95)

## Value

A matrix of M estimates. If `boot` is provided, the `lfqBoot` object
with the additional columns `M_Then` and/or `M_Pauly` (in `bootRaw`,
`maxDen`, `median` and `CI`).

## Details

Function adapted from the mortality function of the fishmethods package
by Gary A. Nelson
(https://cran.r-project.org/web/packages/fishmethods/index.html).

Depending on the method different input parameters are required:

- `"AlversonCarney"` requires `K_l` and `tmax`,

- `"Gislason"` requires `Linf`, `K_l` and `Bl`,

- `"GundersonDygert"` requires `GSI`,

- `"Hoenig"` requires `tmax`,

- `"Lorenzen"` requires `Wwet`,

- `"Pauly_Linf"` requires `Linf`, `K_l` and `temp`,

- `"Pauly_Winf"` requires `Winf`, `K_w` and `temp`,

- `"PetersonWroblewski"` requires `Wdry`,

- `"RikhterEfanov"` requires `tm50`,

- `"Roff"` requires `K_l` and `tm50`,

- `"Then_tmax"` requires `tmax`,

- `"Then_growth"` requires `Linf` and `K_l`.

If accounting for schooling behaviour M is multiplied by 0.8 according
to Pauly (1983).

## References

Alverson, D. L. and M. J. Carney. 1975. A graphic review of the growth
and decay of population cohorts. J. Cons. Int. Explor. Mer 36: 133-143.

Gislason, H., N. Daan, J. C. Rice, and J. G. Pope. 2010. Size, growth,
temperature and the natural mortality of marine fish. Fish and Fisheries
11: 149-158.

Gunderson, D. R. and P. H. Dygert. 1988. Reproductive effort as a
predictor of natural mortality rate. J. Cons. Int. Explor. Mer 44:
200-209.

Hoenig, J. M. 1983. Empirical use of longevity data to estimate
mortality rates. Fish. Bull. 82: 898-903.

Lorenzen, K. 1996. The relationship between body weight and natural
mortality in juvenile and adult fish: a comparison of natural ecosystems
and aquaculture. J. Fish. Biol. 49: 627-647.

Pauly, D. 1980. On the interrelationships between natural mortality,
growth parameters, and mean environmental temperature in 175 fish
stocks. J. Cons. Int. Explor. Mer: 175-192.

Pauly, D., 1983. Some simple methods for the assessment of tropical fish
stocks. *FAO Fish.Tech.Pap.*, (234): 52p. Issued also in French and
Spanish

Peterson, I. and J. S. Wroblewski. 1984. Mortality rate of fishes in the
pelagic ecosystem. Can. J. Fish. Aquat. Sci. 41: 1117-1120.

Rikhter, V.A., and V.N. Efanov, 1976. On one of the approaches to
estimation of natural mortality of fish populations. *ICNAF Res.Doc.*,
76/VI/8: 12p.

Roff, D. A. 1984. The evolution of life history parameters in teleosts.
Can. J. Fish. Aquat. Sci. 41: 989-1000.

Sparre, P., Venema, S.C., 1998. Introduction to tropical fish stock
assessment. Part 1. Manual. *FAO Fisheries Technical Paper*, (306.1,
Rev. 2). 407 p.

Then, A. Y., J. M. Hoenig, N. G. Hall, D. A. Hewitt. 2015. Evaluating
the predictive performance of empirical estimators of natural mortality
rate using information on over 200 fish species. ICES J. Mar. Sci. 72:
82-92.

## Examples

``` r
M_empirical(Linf = 80, K_l = 0.5, temp = 25, tmax = 30,
     method = c("Pauly_Linf","Hoenig"))
#>                                    M
#> Hoenig (1983) - Joint Equation 0.150
#> Hoenig (1983) - Fish Equation  0.139
#> Pauly (1980) - Length Equation 0.819

# with bootstrapped growth parameters
data(alba_boot)
res <- M_empirical(method = "Then_growth", boot = alba_boot)
res$maxDen
#>       Linf          K   t_anchor       phiL     agemax     M_Then 
#> 10.0096156  2.7150279  0.4908732  2.4346091  2.0000000  3.9568143 
res$CI
#>         Linf        K  t_anchor     phiL agemax   M_Then
#> lo  9.171726 1.547109 0.2486596 2.114423      1 2.513000
#> up 11.478216 4.165708 0.5486528 2.739438      2 5.512375
```
