# Package index

## Length-frequency data

Create, modify, restructure and resample `lfq` objects.

- [`lfqCreate()`](http://tokami.github.io/TropFishR/reference/lfqCreate.md)
  : Create lfq data from length measurements
- [`lfqModify()`](http://tokami.github.io/TropFishR/reference/lfqModify.md)
  : Modify lfq data for further analysis
- [`lfqRaise()`](http://tokami.github.io/TropFishR/reference/lfqRaise.md)
  : Raise length-frequency samples to the total catch
- [`lfqRestructure()`](http://tokami.github.io/TropFishR/reference/lfqRestructure.md)
  : Restructuring of length frequency data
- [`lfqResample()`](http://tokami.github.io/TropFishR/reference/lfqResample.md)
  : Resampling of length-frequency data
- [`plot(`*`<lfq>`*`)`](http://tokami.github.io/TropFishR/reference/plot.lfq.md)
  : Plotting of length frequency data (with VBGF curves)
- [`date2yeardec()`](http://tokami.github.io/TropFishR/reference/date2yeardec.md)
  : Date - Year conversion
- [`yeardec2date()`](http://tokami.github.io/TropFishR/reference/yeardec2date.md)
  : Year - Date conversion

## Growth

Estimate von Bertalanffy growth parameters from length-frequency,
length-at-age or tagging data.

- [`ELEFAN()`](http://tokami.github.io/TropFishR/reference/ELEFAN.md) :
  ELEFAN
- [`ELEFAN_SA()`](http://tokami.github.io/TropFishR/reference/ELEFAN_SA.md)
  : ELEFAN_SA
- [`ELEFAN_GA()`](http://tokami.github.io/TropFishR/reference/ELEFAN_GA.md)
  : ELEFAN_GA
- [`lfqFitCurves()`](http://tokami.github.io/TropFishR/reference/lfqFitCurves.md)
  : Fitting VBGF growth curves through lfq data
- [`VBGF()`](http://tokami.github.io/TropFishR/reference/VBGF.md) : Von
  Bertalanffy Growth function (VBGF)
- [`startingPoint2tanchor()`](http://tokami.github.io/TropFishR/reference/startingPoint2tanchor.md)
  : Convert FiSAT's starting point to t_anchor value
- [`powell_wetherall()`](http://tokami.github.io/TropFishR/reference/powell_wetherall.md)
  : Powell-Wetherall method
- [`Bhattacharya()`](http://tokami.github.io/TropFishR/reference/Bhattacharya.md)
  : Bhattacharya's method
- [`plot(`*`<Bhattacharya>`*`)`](http://tokami.github.io/TropFishR/reference/plot.Bhattacharya.md)
  : Bhattacharya plot
- [`growth_length_age()`](http://tokami.github.io/TropFishR/reference/growth_length_age.md)
  : Estimation of growth parameter using length-at-age data
- [`growth_tagging()`](http://tokami.github.io/TropFishR/reference/growth_tagging.md)
  : Growth from tagging data

## Bootstrapped ELEFAN

Uncertainty of growth parameters by resampling (Schwamborn et al. 2019),
and plots of the results (class `lfqBoot`).

- [`ELEFAN_GA_boot()`](http://tokami.github.io/TropFishR/reference/ELEFAN_GA_boot.md)
  : Bootstrapped ELEFAN_GA
- [`ELEFAN_SA_boot()`](http://tokami.github.io/TropFishR/reference/ELEFAN_SA_boot.md)
  : Bootstrapped ELEFAN_SA
- [`plotBoot()`](http://tokami.github.io/TropFishR/reference/plotBoot.md)
  : Plot bootstrapping results
- [`univariate_density()`](http://tokami.github.io/TropFishR/reference/univariate_density.md)
  : Univariate density plots of bootstrapping results
- [`LinfK_scatterhist()`](http://tokami.github.io/TropFishR/reference/LinfK_scatterhist.md)
  : Linf-K scatterplot of bootstrapping results
- [`add_phiprime()`](http://tokami.github.io/TropFishR/reference/add_phiprime.md)
  : Add phi prime contours to a Linf-K plot
- [`vbgfCI_time()`](http://tokami.github.io/TropFishR/reference/vbgfCI_time.md)
  : VBGF plot and CI from bootstrapping results

## Mortality

Natural, total and fishing mortality.

- [`M_empirical()`](http://tokami.github.io/TropFishR/reference/M_empirical.md)
  : Empirical formulas for the estimation of natural mortality
- [`MstochThen()`](http://tokami.github.io/TropFishR/reference/MstochThen.md)
  : Empirical M formula by Then et al. (2015) with uncertainty
- [`catchCurve()`](http://tokami.github.io/TropFishR/reference/catchCurve.md)
  : Catch curve
- [`plot(`*`<catchCurve>`*`)`](http://tokami.github.io/TropFishR/reference/plot.catchCurve.md)
  : Plotting catch curve
- [`Z_BevertonHolt()`](http://tokami.github.io/TropFishR/reference/Z_BevertonHolt.md)
  : Beverton & Holt's Z-Equations
- [`Z_CPUE()`](http://tokami.github.io/TropFishR/reference/Z_CPUE.md) :
  Estimate Z from CPUE data

## Recruitment

- [`recruitment()`](http://tokami.github.io/TropFishR/reference/recruitment.md)
  : Recruitment patterns
- [`plot(`*`<recruitment>`*`)`](http://tokami.github.io/TropFishR/reference/plot.recruitment.md)
  : Plot of recruitment patterns

## Selectivity

Trawl and gillnet selection ogives.

- [`select()`](http://tokami.github.io/TropFishR/reference/select.md) :
  Selectivity model
- [`plot(`*`<select>`*`)`](http://tokami.github.io/TropFishR/reference/plot.select.md)
  : Selectivity plot
- [`select_ogive()`](http://tokami.github.io/TropFishR/reference/select_ogive.md)
  : Selectivity patterns
- [`select_Millar()`](http://tokami.github.io/TropFishR/reference/select_Millar.md)
  : Millar's selectivity model
- [`select_Millar_curve()`](http://tokami.github.io/TropFishR/reference/select_Millar_curve.md)
  : Selection curves of Millar's selectivity models
- [`select_Millar_compare()`](http://tokami.github.io/TropFishR/reference/select_Millar_compare.md)
  : Compare Millar's selectivity models
- [`plot(`*`<select_Millar>`*`)`](http://tokami.github.io/TropFishR/reference/plot.select_Millar.md)
  : Millar's selectivity plot
- [`predict(`*`<select_Millar>`*`)`](http://tokami.github.io/TropFishR/reference/predict.select_Millar.md)
  : Predict the retention of Millar's selectivity model
- [`logLik(`*`<select_Millar>`*`)`](http://tokami.github.io/TropFishR/reference/logLik.select_Millar.md)
  : Log-likelihood of Millar's selectivity model
- [`rcurves_Millar()`](http://tokami.github.io/TropFishR/reference/rcurves_Millar.md)
  : Predict gillnet selectivity (old Millar method)
- [`gillnetfit()`](http://tokami.github.io/TropFishR/reference/gillnetfit.md)
  : Millar's original gillnet selectivity fitting function

## Stock size and prediction

Virtual population analysis, yield per recruit (Beverton-Holt,
Thompson-Bell), biological reference points and stock simulation.

- [`VPA()`](http://tokami.github.io/TropFishR/reference/VPA.md) :
  Virtual Population Analysis (VPA)
- [`plot(`*`<VPA>`*`)`](http://tokami.github.io/TropFishR/reference/plot.VPA.md)
  : VPA plot
- [`predict_mod()`](http://tokami.github.io/TropFishR/reference/predict_mod.md)
  : Prediction models
- [`plot(`*`<predict_mod>`*`)`](http://tokami.github.io/TropFishR/reference/plot.predict_mod.md)
  : Plotting prediction models yield per recruit and Thompson & Bell
- [`stock_sim()`](http://tokami.github.io/TropFishR/reference/stock_sim.md)
  : Stock simulation

## Surplus production models

- [`prod_mod()`](http://tokami.github.io/TropFishR/reference/prod_mod.md)
  : Production models
- [`plot(`*`<prod_mod>`*`)`](http://tokami.github.io/TropFishR/reference/plot.prod_mod.md)
  : Plotting production models
- [`prod_mod_ts()`](http://tokami.github.io/TropFishR/reference/prod_mod_ts.md)
  : Production models with time series fitting
- [`plot(`*`<prod_mod_ts>`*`)`](http://tokami.github.io/TropFishR/reference/plot.prod_mod_ts.md)
  : Plotting time series production models

## Data sets from Sparre & Venema

Data used in the worked examples of the FAO manual and other studies.

- [`alba`](http://tokami.github.io/TropFishR/reference/alba.md) :
  Length-frequency data of the clam Abra alba
- [`alba_boot`](http://tokami.github.io/TropFishR/reference/alba_boot.md)
  : Bootstrapped VBGF estimates for the alba lfq dataset
- [`bream`](http://tokami.github.io/TropFishR/reference/bream.md) :
  bream data
- [`emperor`](http://tokami.github.io/TropFishR/reference/emperor.md) :
  Emperor data
- [`gillnet`](http://tokami.github.io/TropFishR/reference/gillnet.md) :
  Gillnet data
- [`goatfish`](http://tokami.github.io/TropFishR/reference/goatfish.md)
  : Yellowstriped goatfish data
- [`haddock`](http://tokami.github.io/TropFishR/reference/haddock.md) :
  Haddock data
- [`hake`](http://tokami.github.io/TropFishR/reference/hake.md) : Hake
  data
- [`shrimps`](http://tokami.github.io/TropFishR/reference/shrimps.md) :
  Shrimp data
- [`tilapia`](http://tokami.github.io/TropFishR/reference/tilapia.md) :
  Tilapia data
- [`trammelnet`](http://tokami.github.io/TropFishR/reference/trammelnet.md)
  : Trammel net data
- [`trawl_fishery_Java`](http://tokami.github.io/TropFishR/reference/trawl_fishery_Java.md)
  : Data from the trawl fishery off the North coast of Java
- [`whiting`](http://tokami.github.io/TropFishR/reference/whiting.md) :
  Whiting data

## Synthetic data sets

Simulated data with known parameters, useful for testing methods.

- [`synCAA1`](http://tokami.github.io/TropFishR/reference/synCAA1.md) :
  Synthetic Catch-at-age data I
- [`synCAA2`](http://tokami.github.io/TropFishR/reference/synCAA2.md) :
  Synthetic Catch-at-age data II
- [`synCPUE`](http://tokami.github.io/TropFishR/reference/synCPUE.md) :
  Synthetical catch per unit of effort (CPUE) dataset
- [`synLFQ1`](http://tokami.github.io/TropFishR/reference/synLFQ1.md) :
  Synthetic length-frequency data I
- [`synLFQ2`](http://tokami.github.io/TropFishR/reference/synLFQ2.md) :
  Synthetic length frequency data II
- [`synLFQ3`](http://tokami.github.io/TropFishR/reference/synLFQ3.md) :
  Synthetic length frequency data III
- [`synLFQ4`](http://tokami.github.io/TropFishR/reference/synLFQ4.md) :
  Synthetic length-frequency data IV (with seasonal oscillation)
- [`synLFQ5`](http://tokami.github.io/TropFishR/reference/synLFQ5.md) :
  Synthetic length-frequency data V (without seasonal oscillation)
- [`synLFQ6`](http://tokami.github.io/TropFishR/reference/synLFQ6.md) :
  Synthetic length-frequency data VI (without seasonal oscillation)
- [`synLFQ7`](http://tokami.github.io/TropFishR/reference/synLFQ7.md) :
  Synthetic length-frequency data VII with seasonal oscillation
- [`synLFQ8`](http://tokami.github.io/TropFishR/reference/synLFQ8.md) :
  Synthetic length-frequency data VIII with variable harvest rate

## Internal helpers

Documented but not exported; used by `predict_mod`.

- [`ypr()`](http://tokami.github.io/TropFishR/reference/ypr.md) : Yield
  per recruit
- [`ypr_sel()`](http://tokami.github.io/TropFishR/reference/ypr_sel.md)
  : Yield per recruit with selection ogive
