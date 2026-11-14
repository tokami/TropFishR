# Package index

## Length-frequency data

Create, modify, restructure and resample `lfq` objects.

- [`lfqCreate()`](https://tokami.github.io/TropFishR/reference/lfqCreate.md)
  : Create lfq data from length measurements
- [`lfqModify()`](https://tokami.github.io/TropFishR/reference/lfqModify.md)
  : Modify lfq data for further analysis
- [`lfqRaise()`](https://tokami.github.io/TropFishR/reference/lfqRaise.md)
  : Raise length-frequency samples to the total catch
- [`lfqRestructure()`](https://tokami.github.io/TropFishR/reference/lfqRestructure.md)
  : Restructuring of length frequency data
- [`lfqResample()`](https://tokami.github.io/TropFishR/reference/lfqResample.md)
  : Resampling of length-frequency data
- [`plot(`*`<lfq>`*`)`](https://tokami.github.io/TropFishR/reference/plot.lfq.md)
  : Plotting of length frequency data (with VBGF curves)
- [`date2yeardec()`](https://tokami.github.io/TropFishR/reference/date2yeardec.md)
  : Date - Year conversion
- [`yeardec2date()`](https://tokami.github.io/TropFishR/reference/yeardec2date.md)
  : Year - Date conversion

## Growth

Estimate von Bertalanffy growth parameters from length-frequency,
length-at-age or tagging data.

- [`ELEFAN()`](https://tokami.github.io/TropFishR/reference/ELEFAN.md) :
  ELEFAN
- [`ELEFAN_SA()`](https://tokami.github.io/TropFishR/reference/ELEFAN_SA.md)
  : ELEFAN_SA
- [`ELEFAN_GA()`](https://tokami.github.io/TropFishR/reference/ELEFAN_GA.md)
  : ELEFAN_GA
- [`lfqFitCurves()`](https://tokami.github.io/TropFishR/reference/lfqFitCurves.md)
  : Fitting VBGF growth curves through lfq data
- [`VBGF()`](https://tokami.github.io/TropFishR/reference/VBGF.md) : Von
  Bertalanffy Growth function (VBGF)
- [`startingPoint2tanchor()`](https://tokami.github.io/TropFishR/reference/startingPoint2tanchor.md)
  : Convert FiSAT's starting point to t_anchor value
- [`powell_wetherall()`](https://tokami.github.io/TropFishR/reference/powell_wetherall.md)
  : Powell-Wetherall method
- [`Bhattacharya()`](https://tokami.github.io/TropFishR/reference/Bhattacharya.md)
  : Bhattacharya's method
- [`plot(`*`<Bhattacharya>`*`)`](https://tokami.github.io/TropFishR/reference/plot.Bhattacharya.md)
  : Bhattacharya plot
- [`growth_length_age()`](https://tokami.github.io/TropFishR/reference/growth_length_age.md)
  : Estimation of growth parameter using length-at-age data
- [`growth_tagging()`](https://tokami.github.io/TropFishR/reference/growth_tagging.md)
  : Growth from tagging data

## Bootstrapped ELEFAN

Uncertainty of growth parameters by resampling (Schwamborn et al. 2019),
and plots of the results (class `lfqBoot`).

- [`ELEFAN_GA_boot()`](https://tokami.github.io/TropFishR/reference/ELEFAN_GA_boot.md)
  : Bootstrapped ELEFAN_GA
- [`ELEFAN_SA_boot()`](https://tokami.github.io/TropFishR/reference/ELEFAN_SA_boot.md)
  : Bootstrapped ELEFAN_SA
- [`plotBoot()`](https://tokami.github.io/TropFishR/reference/plotBoot.md)
  : Plot bootstrapping results
- [`univariate_density()`](https://tokami.github.io/TropFishR/reference/univariate_density.md)
  : Univariate density plots of bootstrapping results
- [`LinfK_scatterhist()`](https://tokami.github.io/TropFishR/reference/LinfK_scatterhist.md)
  : Linf-K scatterplot of bootstrapping results
- [`add_phiprime()`](https://tokami.github.io/TropFishR/reference/add_phiprime.md)
  : Add phi prime contours to a Linf-K plot
- [`vbgfCI_time()`](https://tokami.github.io/TropFishR/reference/vbgfCI_time.md)
  : VBGF plot and CI from bootstrapping results

## Mortality

Natural, total and fishing mortality.

- [`M_empirical()`](https://tokami.github.io/TropFishR/reference/M_empirical.md)
  : Empirical formulas for the estimation of natural mortality
- [`MstochThen()`](https://tokami.github.io/TropFishR/reference/MstochThen.md)
  : Empirical M formula by Then et al. (2015) with uncertainty
- [`catchCurve()`](https://tokami.github.io/TropFishR/reference/catchCurve.md)
  : Catch curve
- [`plot(`*`<catchCurve>`*`)`](https://tokami.github.io/TropFishR/reference/plot.catchCurve.md)
  : Plotting catch curve
- [`Z_BevertonHolt()`](https://tokami.github.io/TropFishR/reference/Z_BevertonHolt.md)
  : Beverton & Holt's Z-Equations
- [`Z_CPUE()`](https://tokami.github.io/TropFishR/reference/Z_CPUE.md) :
  Estimate Z from CPUE data

## Recruitment

- [`recruitment()`](https://tokami.github.io/TropFishR/reference/recruitment.md)
  : Recruitment patterns
- [`plot(`*`<recruitment>`*`)`](https://tokami.github.io/TropFishR/reference/plot.recruitment.md)
  : Plot of recruitment patterns

## Selectivity

Trawl and gillnet selection ogives.

- [`select()`](https://tokami.github.io/TropFishR/reference/select.md) :
  Selectivity model
- [`plot(`*`<select>`*`)`](https://tokami.github.io/TropFishR/reference/plot.select.md)
  : Selectivity plot
- [`select_ogive()`](https://tokami.github.io/TropFishR/reference/select_ogive.md)
  : Selectivity patterns
- [`select_Millar()`](https://tokami.github.io/TropFishR/reference/select_Millar.md)
  : Millar's selectivity model
- [`select_Millar_curve()`](https://tokami.github.io/TropFishR/reference/select_Millar_curve.md)
  : Selection curves of Millar's selectivity models
- [`select_Millar_compare()`](https://tokami.github.io/TropFishR/reference/select_Millar_compare.md)
  : Compare Millar's selectivity models
- [`plot(`*`<select_Millar>`*`)`](https://tokami.github.io/TropFishR/reference/plot.select_Millar.md)
  : Millar's selectivity plot
- [`predict(`*`<select_Millar>`*`)`](https://tokami.github.io/TropFishR/reference/predict.select_Millar.md)
  : Predict the retention of Millar's selectivity model
- [`logLik(`*`<select_Millar>`*`)`](https://tokami.github.io/TropFishR/reference/logLik.select_Millar.md)
  : Log-likelihood of Millar's selectivity model
- [`rcurves_Millar()`](https://tokami.github.io/TropFishR/reference/rcurves_Millar.md)
  : Predict gillnet selectivity (old Millar method)
- [`gillnetfit()`](https://tokami.github.io/TropFishR/reference/gillnetfit.md)
  : Millar's original gillnet selectivity fitting function

## Stock size and prediction

Virtual population analysis, yield per recruit (Beverton-Holt,
Thompson-Bell), biological reference points and stock simulation.

- [`VPA()`](https://tokami.github.io/TropFishR/reference/VPA.md) :
  Virtual Population Analysis (VPA)
- [`plot(`*`<VPA>`*`)`](https://tokami.github.io/TropFishR/reference/plot.VPA.md)
  : VPA plot
- [`predict_mod()`](https://tokami.github.io/TropFishR/reference/predict_mod.md)
  : Prediction models
- [`plot(`*`<predict_mod>`*`)`](https://tokami.github.io/TropFishR/reference/plot.predict_mod.md)
  : Plotting prediction models yield per recruit and Thompson & Bell
- [`stock_sim()`](https://tokami.github.io/TropFishR/reference/stock_sim.md)
  : Stock simulation

## Surplus production models

- [`prod_mod()`](https://tokami.github.io/TropFishR/reference/prod_mod.md)
  : Production models
- [`plot(`*`<prod_mod>`*`)`](https://tokami.github.io/TropFishR/reference/plot.prod_mod.md)
  : Plotting production models
- [`prod_mod_ts()`](https://tokami.github.io/TropFishR/reference/prod_mod_ts.md)
  : Production models with time series fitting
- [`plot(`*`<prod_mod_ts>`*`)`](https://tokami.github.io/TropFishR/reference/plot.prod_mod_ts.md)
  : Plotting time series production models

## Data sets from Sparre & Venema

Data used in the worked examples of the FAO manual and other studies.

- [`alba`](https://tokami.github.io/TropFishR/reference/alba.md) :
  Length-frequency data of the clam Abra alba
- [`alba_boot`](https://tokami.github.io/TropFishR/reference/alba_boot.md)
  : Bootstrapped VBGF estimates for the alba lfq dataset
- [`bream`](https://tokami.github.io/TropFishR/reference/bream.md) :
  bream data
- [`emperor`](https://tokami.github.io/TropFishR/reference/emperor.md) :
  Emperor data
- [`gillnet`](https://tokami.github.io/TropFishR/reference/gillnet.md) :
  Gillnet data
- [`goatfish`](https://tokami.github.io/TropFishR/reference/goatfish.md)
  : Yellowstriped goatfish data
- [`haddock`](https://tokami.github.io/TropFishR/reference/haddock.md) :
  Haddock data
- [`hake`](https://tokami.github.io/TropFishR/reference/hake.md) : Hake
  data
- [`shrimps`](https://tokami.github.io/TropFishR/reference/shrimps.md) :
  Shrimp data
- [`tilapia`](https://tokami.github.io/TropFishR/reference/tilapia.md) :
  Tilapia data
- [`trammelnet`](https://tokami.github.io/TropFishR/reference/trammelnet.md)
  : Trammel net data
- [`trawl_fishery_Java`](https://tokami.github.io/TropFishR/reference/trawl_fishery_Java.md)
  : Data from the trawl fishery off the North coast of Java
- [`whiting`](https://tokami.github.io/TropFishR/reference/whiting.md) :
  Whiting data

## Synthetic data sets

Simulated data with known parameters, useful for testing methods.

- [`synCAA1`](https://tokami.github.io/TropFishR/reference/synCAA1.md) :
  Synthetic Catch-at-age data I
- [`synCAA2`](https://tokami.github.io/TropFishR/reference/synCAA2.md) :
  Synthetic Catch-at-age data II
- [`synCPUE`](https://tokami.github.io/TropFishR/reference/synCPUE.md) :
  Synthetical catch per unit of effort (CPUE) dataset
- [`synLFQ1`](https://tokami.github.io/TropFishR/reference/synLFQ1.md) :
  Synthetic length-frequency data I
- [`synLFQ2`](https://tokami.github.io/TropFishR/reference/synLFQ2.md) :
  Synthetic length frequency data II
- [`synLFQ3`](https://tokami.github.io/TropFishR/reference/synLFQ3.md) :
  Synthetic length frequency data III
- [`synLFQ4`](https://tokami.github.io/TropFishR/reference/synLFQ4.md) :
  Synthetic length-frequency data IV (with seasonal oscillation)
- [`synLFQ5`](https://tokami.github.io/TropFishR/reference/synLFQ5.md) :
  Synthetic length-frequency data V (without seasonal oscillation)
- [`synLFQ6`](https://tokami.github.io/TropFishR/reference/synLFQ6.md) :
  Synthetic length-frequency data VI (without seasonal oscillation)
- [`synLFQ7`](https://tokami.github.io/TropFishR/reference/synLFQ7.md) :
  Synthetic length-frequency data VII with seasonal oscillation
- [`synLFQ8`](https://tokami.github.io/TropFishR/reference/synLFQ8.md) :
  Synthetic length-frequency data VIII with variable harvest rate

## Internal helpers

Documented but not exported; used by `predict_mod`.

- [`ypr()`](https://tokami.github.io/TropFishR/reference/ypr.md) : Yield
  per recruit
- [`ypr_sel()`](https://tokami.github.io/TropFishR/reference/ypr_sel.md)
  : Yield per recruit with selection ogive
