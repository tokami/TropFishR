
### TropFishR (development version)

---

##### New features

- Bootstrapped ELEFAN (Schwamborn et al. 2019) to assess the uncertainty of
  the growth parameters:
  - `lfqResample` resamples length-frequency data.
  - `ELEFAN_GA_boot` and `ELEFAN_SA_boot` fit `ELEFAN_GA` or `ELEFAN_SA` to
    many resampled data sets (optionally in parallel) and return the
    estimates with maximum density, medians and confidence intervals (class
    `lfqBoot`).
  - Plots of the results: `univariate_density`, `LinfK_scatterhist` (with
    `add_phiprime` for phi' isolines), `vbgfCI_time` and `plotBoot`.
  - Example data set `alba_boot`.
  - New dependency: `ks`.
- `ELEFAN`: new argument `plot.score` (default `TRUE`), as in `ELEFAN_GA` and
  `ELEFAN_SA`, to switch off the score graph (response surface or K-scan).
  The argument `plot` only controls the plot of the restructured data with
  the growth curves.

##### Bug fixes

- `lfqCreate` with `plus_group`: the plus group no longer counts the fish in
  its own length class twice.
- `growth_tagging` with `method = "Munro"`: `K` is now the mean of the
  individual growth rates at the optimal `Linf`. Previously the coefficient
  of variation of the individual growth rates was returned as `K`.
- `stock_sim` (age-based):
  - Mean biomass is now also correct without fishing (`FM = 0`); previously
    it was 0.
  - The yield and value of the plus group now include all of its catch.
    Previously only the catch beyond the first time step of the plus group
    was used, so yield and value were too low.
  - With a plus group, `meanB` no longer depends on age classes that are
    older than the plus group.
- `predict_mod` with `type = "ThompBell"`:
  - The exploitation rates of the reference points (`E01`, `Emax`, `E05`,
    `E04`) are now calculated as `E = F / (F + M)`. Previously they could
    be greater than 1.
  - The reference points are searched over the full range of `FM_change`.
    Previously the search was limited to `F` between 0 and 5.
  - The biomass per recruit at F0.1 is now returned as `BPR_F01`. Previously
    it was returned under the duplicate name `YPR_F01` (showing up as
    `YPR_F01.1`).
  - Because of the fixes in `stock_sim`, the biomass at `F = 0` and the
    reference point F0.5 are now correct for age-based data.
- `prod_mod`: the Fox virgin biomass (`Fox_Bv`) is now `exp(c)`, the CPUE at
  zero effort of the Fox model `ln(CPUE) = c + d * f`. Previously the
  intercept `c` itself was returned.
- `yeardec2date` no longer returns the previous day for some dates because
  of floating point error (e.g. 2019-07-15 was returned as 2019-07-14).

##### Minor changes

- Added unit tests (`testthat`).


<br><br>



### TropFishR 1.6.6 (2025-04-06)

---

##### Minor changes

- Removed the dependency on the package `propagate`, which was archived on
  CRAN. The confidence band of the fitted growth curve in `growth_length_age`
  (`method = "LSM"`, `CI = TRUE`) is now calculated with the delta method.
- `growth_length_age`: the argument `do.sim` is deprecated and has no effect
  (a warning is given if `do.sim = TRUE`); the argument `nsim` was removed.


<br><br>



### TropFishR 1.6.5 (2024-10-21)

---

##### Minor changes

- Fixed CRAN notes: package authors are given as `Authors@R` (incl. ORCID),
  updated link to the FAO manual, and corrected documentation links.
- Links between the vignettes were restored.
- Shorter vignette building time: fewer iterations in the ELEFAN vignette and
  shorter `SA_time` in the `ELEFAN_SA` example.


<br><br>



### TropFishR 1.6.4 (2024-01-29)

---

##### Bug fixes

- `ELEFAN_GA` and `ELEFAN_SA` with `plot = TRUE`: the growth curves are drawn
  with the estimated parameters (`par`); previously the non-existing element
  `pars` was used.

##### Minor changes

- `Bhattacharya` requires an interactive session and returns `NULL` with a
  message otherwise.
- `catchCurve` and `powell_wetherall`: if `reg_int` is not provided (and
  `auto = FALSE` in `catchCurve`), the points for the regression line can
  only be chosen in an interactive session; otherwise `NULL` is returned with
  a message.
- Progress bars in `ELEFAN` and `predict_mod` are only shown in interactive
  sessions.
- `lfqRestructure` (and thus all ELEFAN functions) gives an informative error
  if a sample contains only zeros.
- `ELEFAN`: the score graph is only drawn if `plot = TRUE`.
- `ELEFAN_SA`: the score graph (`plot.score`) is currently not drawn, because
  the trace of the optimisation is not available from `GenSA` any more.
- `growth_length_age` (`method = "LSM"`, `CI = TRUE`) plots the confidence
  band of the growth curve again, using `propagate::predictNLS`.
- `Matrix` moved from `Depends` to `Imports`; `stats` and `utils` moved from
  `Suggests` to `Imports`.
- Updated URLs; links between vignettes were removed for CRAN.


<br><br>



### TropFishR 1.6.3 (2021-09-02)

---

#### New features

- Exclude NA in `VBGF` function with argument `na.rm = TRUE` (default is
  `FALSE`)
- SPR calculation in Thompson and Bell yield per recruit model + SPR-based
reference points


##### Bug fixes

- Correct separation index calculation in Bhattacharya's method (correct = `(L1 - L2) /
  (sd1 + sd2)/2`)
- Correct standard deviation calculation in Bhattacharya's method (correct = `sqrt(-dL/b`)
- Correct equations of Heincke and Robson & Chapman methods in the `Z_CPUE` function
- Return `score_mat` in `ELEFAN` function

##### Minor changes

- Updated vignettes
- `LazyData = false` which implies that `data("...")` has to be used to load a
  package into the R environment, e.g. `data("snyLFQ7")`
- Argument `plots_refs` in `plot.predict_mod` function to select which reference
points in the prediction models should be plotted. By default all reference
points are plotted: `plot_refs = c("F01","Fmax","F05","F04")`.
- Argument `cols_refs` in `plot.predict_mod` function to select the colour for
the reference points in the prediction models. Default colours are `cols_refs = c("goldenrod2","darkred","darkgreen","darkorange")`.
- Change naming from 'Fmsy' in prediction models to 'Fmax'.


<br><br>



### TropFishR 1.6.2 (2020-01-27)

---

##### Bug fixes
      - Use of matrices was adjusted to comply with new R matrix classes
      - ELEFAN_GA was adjusted to use new arguments of the underlying GA::ga
        function

<br><br>



### TropFishR 1.6.1

---

##### Bug fixes
      - The LBB functionality was removed from the package.


<br><br>



### TropFishR 1.6

---

##### New features
      - The Length-based Bayesian biomass estimator method (LBB) by
        Froese et al. (2018) has been added to the package.
      - Improved implementation of compilation of length-frequency
        data allowing for faster compilation of lfq data with
        lfqCreate and lfqModify (with argument bin_size) and allowing
        compilation of large data sets
      - New possibilities to easily modify length-frequency data
        objects by means of "lfqModify", including aggregating catch
        matrix per year, per quarter, or per month, subsetting lfq
        data with range of length classes or sampling dates

##### Documentation
      - A new vignette "LBBmanual" with the introduction and
        demonstration of LBB within TropFishR has been added.

##### Bug fixes
      - Seasonal growth parameters are now added correctly to the lfq
        object in lfqModify


<br><br>



### TropFishR 1.2.1

---

##### Documentation
      - A new vignette has been added to the package. The vignette "ELEFANTutorial" outlines
      all ELEFAN functions available in TropFishR in detail.


##### Bug fixes
      - The ELEFAN functions did not overwrite any element of the lfq object but instead
      concatinated parameters onto the object. This can have unintended side effects, e.g. if several
      growth parameters are saved in the lfq object and the plotting functions are called.
      Now, the application of ELEFAN functions overwrites any growth parameters in the lfq object.


<br><br>


### TropFishR 1.2

---

##### New features
      - Due to more efficient matrix computations ELEFAN is 2-4 times faster as before
      - missing seed values were added in ELEFAN_GA()
      - `plot.lfq()` allows plotting relative frequencies, this is
      in particular useful in cases where one sampling time has a large number of samples
      - `plot.catchCurve()` allows plotting results according to length rather
      than relative age. This can be done with the argument `xaxis`.
      - lfq plot can be used to plot 2 different lfq data sets; this allows to visually compare the two data sets
      - `lfqModfiy()` allows combining two different lfq data sets. This
      might be interesting when different fleets are investigated separately
      - possible to define and plot multiple regression lines
      in the catch curve analysis
      - length converted catch curve can account for the seasonalized VBGF


##### User-visible changes
      - better control over graphical devices with par(),
      when par() defined the default settings are not used
      - growth parameters can now be added to lfq lists as "par" elements, this
      is convenient as it is in line with the results of the ELEFAN methods


##### Documentation
      - new vignette with a short description of lfq data and
      how to import lfq data into R
      - updated tutorial vignette
      - more informative error and warning messages for many functions
      - document with news and changes about package version was added


##### Bug fixes
      - crash report of ELEFAN_GA with interactive sessions (Rstudio) fixed
      and more stable on windows
      - plot.catchCurve() was not displaying the regression line, but
      a straight line from the first to last point of the chosen interval
      - restructering of lfq data was not in line with FiSAT implementation
      - fixing bug in handling of leap years
      - when merging lfq data with another list using c() one has to reassign
      the class "lfq" to the merged object, this has been added in the tutorial
      - beep sound was causing R crashes on windows computers and was therefore removed
