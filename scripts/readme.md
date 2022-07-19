# scripts

## simplepart

The script `simplepart_simulation.R` was run on a high performance cluster, from a command line script that provided several parameters determining the configuration for the simulation analysis. Specifically, it was run using the line
``` 
R --no-save --args arg1 arg2 arg3 < scripts/simplepart_simulation.R 
```
where the passed parameters determined the method and the data replicate used (`arg1`), the size of the datasets, specifically the number of high resolution units contained in each low resolution unit (`arg2`) and the cluster separation level (`arg3`). 

In our simulation we fixed the cluster separation to a moderate level (fixing `arg3` to `2`), and we considered the number of units to range between 10, 25 and 50 (respectively fixing `arg2` to `2`, `3` and `4`). Moreover, we considered three methods (the nHDP, the multi-resolution adaptation of k-means, and the nDP) and 50 data replicates. By ranging `arg1` from `1` to `150` we can produce all the required simulations results (where the first 50 implement the nHDP, the second 50 implement the adjusted k-means and the last 50 implement the nDP).
The same can be achieved by manually setting the arguments manually, or ranging them using RStudio jobs.

Using the script `simpledata_compile.R` we combined the output and results from all these simulation configurations, and saved the produced summaries in `results/summaries`.

Finally, with `simplepart_plot.R` we produce the plots reported in the paper using the saved summaries. Since the compiled summaries have been saved in `results/summaries`, this script can be run without having to replicate the full simulation analysis.

## modeldata

For the synthetic data generated from the model, we first need to save the generated data, which is done by running `scripts/modeldata_generate.R`.
Again, this can be run using a command line script passing a set of two parameters:
``` 
R --no-save --args arg1 arg2 < scripts/modeldata_generate.R 
```
where `arg1` corresponds to the data replicate (ranged from 1 to 50), and `arg2` corresponds to the configuration of hyperparameters used (1 for the configuration of 3,5,3 and 2 for 1,1,1).

The script `modeldata_simulation.R` was used to run the simulations, passing a set of four parameters: 
``` 
R --no-save --args arg1 arg2 arg3 arg4 < scripts/modeldata_simulation.R 
```
The parameters determined the method the data replicate used (`arg1`), the configuration of hyperparameters passed to the model (`arg2`), the size of the datasets, specifically the number of high resolution units contained in each low resolution unit (`arg3`) and the cluster separation level (`arg4`). 

In our simulation we fixed the cluster separation to a moderate level (fixing `arg4` to `2`), we considered the number of units to range between 10 and 50 (respectively fixing `arg3` to `2` and `4`), and we considered two configurations for the hyperparameters, fixing `arg2` to `1` and `2` to use the hyperparameters configurations 3,5,3 and 1,1,1. Moreover, we considered three methods (the nHDP, the multi-resolution adaptation of k-means, and the nDP) and 50 data replicates. By ranging `arg1` from `1` to `150` we can produce all the required simulations results (where the first 50 implement the nHDP, the second 50 implement the adjusted k-means and the last 50 implement the nDP).
The same can be achieved by manually setting the arguments manually, or ranging them using RStudio jobs.

The output was combined with `modeldata_compile.R` and the produced summaries were saved in `results/summaries`.

Finally, with `modeldata_plot.R` we produce the plots reported in the paper using the saved summaries.

## westphilly

**TODO**