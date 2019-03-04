### TropFishR 1.6.2
    
---

##### New features


##### User-visible changes
      - argument "x" of all ELEFAN functions (ELEFAN, ELEFAN\_GA, and
        ELEFAN\_SA) was changed to lfq
	  - Parameter "t_anchor" was renamed to "ta"
	  - All life-history parameters (e.g. Linf, K, M, a, b) are now
        stored in a sub-list of the lfq list called "par" and can be
        accessed with `lfq$par`
		

##### Documentation
      - Vignettes were updated
	  - Tutorial and thus recommended assessment steps changed


##### Bug fixes
      - "score_mat" is again reported by the ELEFAN function

      
      
<br><br>


### TropFishR 1.6.1
    
---

##### User-visible changes
      - All LBB related functions and vignettes have been removed
        again from the package due to problems stated by users

      
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

