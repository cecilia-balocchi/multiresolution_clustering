# Clustering Areal Units at Multiple Levels of Resolution to Model Crime Incidence in Philadelphia

This repository contains the code for replication of the results in the paper "Clustering Areal Units at Multiple Levels of Resolution to Model Crime Incidence in Philadelphia" by Balocchi, George and Jensen ([arxiv link](https://arxiv.org/abs/2112.02059)).

## Overview

The directory `scripts` contains various R scripts to run and analyze the Synthetic Data Evaluation (Section 3 of the paper and Section S3 of the Supplementary Material) and West Philadelphia Data Evaluation (Section 4). 

The directory `src` contains the Cpp implementation of the multiresolution clustering algorithm proposed in the paper, based on the nested Hierarchical Dirichlet Process (nHDP). 
This includes the special cases of single-resolution clustering with the Dirichlet Process (DP) and the Hierarchical Dirichlet Process (HDP).

The directory `R` contains R wrapper functions that call the Cpp implementation in `src`, together with the R code for the nested Dirichlet Process (nDP) from Zuanetti et al. (2017), provided as Supplementary Code [here](https://onlinelibrary.wiley.com/doi/abs/10.1111/biom.12778). Moreover, it contains variouus functions that are used for generating, analyzing and plotting the results.

The directory `data` contains the cleaned and compiled crime density data for West Philadelphia, together with shapefiles and geographic information for the area.

The directory `results` contains some compiled summaries of the simulations performed, which can be used for easy replication of tables and figures reported in the paper. Moreover it contains the results from running the nHDP model on the West Philadelphia crime density data.

To run this code, the following packages are required:
```
list.of.packages <- c("rgeos", "spdep", "reshape2", "ggplot2", "dplyr", "gridExtra", "mapproj",
					  "mcclust", "cluster", "devtools", "compiler") 
					  # devtools only needed to install mcclust.ext
					  # compiler only needed for nDP code
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)

if("mcclust.ext" %in% installed.packages()[,"Package"]) {
	devtools::install_github("sarawade/mcclust.ext")
}
```

## Most important functions

The function to implement the nHDP sampling method is `paralleltemp`, contained in `R/paralleltemp.R`. This is a wrapper function calling the Cpp implementation in `src`.

Similarly, the function `paralleltemp_noLR` contained in `R/paralleltemp_noLR.R` is a wrapper function for the Cpp code in `src` which implements the HDP sampling method (and as a special case, the DP sampling method), by fixing the partition at the low resolution and sampling only the partition at the high resolution. 

A user could also find useful the functions contained in `R/functions_plot.R` that are used to plot the detected partition in the maps.

## Replication of paper results

### Synthetic Data Evaluation

The Synthetic Data Evaluation is composed by two parts:
- an analysis of data generated from six mixtures of normals, as described in Section 3 of the paper. The scripts related to this part are denoted with the name *simplepart*.
- an analysis of data generated from the model (equation (1) of the manuscript) as described in Section S3 of the Supplementary Material. The scripts related to this part are denoted with the name *modeldata*.

#### simplepart

The script `simplepart_simulation.R` will actually run the simulation, generating the data and comparing several methods. See `scripts/readme.md` for further details on how to replicate the simulation reported in the paper.

Using the script `simpledata_compile.R` we combined the output and results from all these simulation configurations, and saved the produced summaries in `results/summaries`.

Finally, with `simplepart_plot.R` we produce the plots reported in the paper using the saved summaries.

#### modeldata

Similarly, for the data generated from the model, the script `modeldata_simulation.R` was used to run the simulations, after generating the data with `scripts/modeldata_generate.R`. See `scripts/readme.md` for further details on how to replicate the simulation reported in the supplementary material. 

The output was combined with `modeldata_compile.R` and the produced summaries were saved in `results/summaries`.

Finally, with `modeldata_plot.R` we produce the plots reported in the paper using the saved summaries.

### Data analysis of crime density in West Philadelphia

The data (crime density in West Philadelphia) was cleaned and compiled in ... **TODO**

The analysis of the data was run using `scripts/westphilly_analysis.R` and the results were compiled and plotted using `scripts/westphilly_plot.R`.
