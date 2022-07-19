## R

The script `paralleltemp.R` contains a wrapper function to run the main Cpp code implementing the nHDP. 

The script `paralleltemp_noLR.R` contains a wrapper function to run a variation of the main Cpp code, which implements the HDP. The DP can be recovered as a special case, when there is one unique group in the HDP.
<!-- Different parameters specified in `paralleltemp` determine which method is implemented. -->

The script `nDP.R` contains the code from Zuanetti et al. (2017).

Various R functions used for the analysis of the simulations are contained in `functions_analysis.R`.

The functions to produce maps contained in the paper, highlighting the partition recovered, are contained in the script `functions_plot.R`.
