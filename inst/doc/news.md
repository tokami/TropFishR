### TropFishR 1.1.5

##### New features

      - ELEFAN 4 times faster with more efficient calculations
      - missing seed values in ELEFAN added
      - plot.lfq() allows now to plot relative frequencies, this is
      in particular useful in cases where one sampling time has a high sample size
      - plot.catchCurve() allows now to plot results according to length rather
      than relative age. This can be done with the argument 'xaxis'.
      - lfq plot allows plotting of 2 different lfq data sets
      - lfqModfiy() allows now to combine two different lfq data sets. This
      might be interesting when different fleets are investigated separately
      - possible to define and plot multiple regression lines
      in catch curve analysis
      - length converted catch curve can now account for the seasonalized
      VBGF


##### User-visible changes

      - better control over graphical devices with par(),
      when par() defined the default settings are not used
      - growth parameters can no be added to lfq lists as "par" elements, this
      is convenient because this corresponds to the structure of the results
      from all ELEFAN methods


##### Documentation

      - new vignette with a short description of lfq data and
      how to import in R and TropFishR
      - updated tutorial vignette
      - more informative error and warning messages for many functions added

##### Bug fixes

      - crash report of ELEFAN_GA with interactive sessions (Rstudio) fixed
      and more stable on windows
      - plot.catchCurve() was not displaying the regression line, but
      a straight line from the first to last point of the chosen interval
      - restructering of lfq data was not in line with FiSAT implementation
      - fixing bug in handling of leap years
      - when merging your lfq data with another list using c() you have to reassign
      the class "lfq" to the merged object, this has been added in the tutorial

##### Deprecated & defunct



      
### TropFishR 1.1.4


    not kept track off