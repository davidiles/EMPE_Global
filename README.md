# LaRue et al. (2024) Advances in remote sensing of emperor penguins: first multi-year time series documenting global population change

[Access this dataset on Dryad](Dataset DOI link)

# Overview

This repository contains data, code, and model output associated with the analysis of Emperor penguin population dynamics described in LaRue et al. (2024).

Raw data for the analysis are  included aerial and ground-based observations



Give a brief summary of dataset contents, contextualized in experimental procedures and results.

## File structure

- analysis
    - output
|  |  |- data_viz ** contains stuff **
|  |  |- model_checks
|  |  |- model_results
|  |  |  |- Global_Level
|  |  |  |- Regional_Level
|  |  |  |- Colony_Level
|  |  |  |- parameter_estimates.csv

|  |  |- simulation

|  |  |- model_results
|  |  |  |- 1_Global_Level
|  |  |  |  |- GLOBAL_abundance.csv
|  |  |  |  |- GLOBAL_change.csv
|  |  |  |  |- GLOBAL_trend.csv
|  |  |- tables
|  |- output_simulation/
|  |- EMPE_model_empirical.jags
|  |- EMPE_model_simulate_data.jags
|  |- script1_PrepareData.R
|  |- script2_FitModel.R
|  |- script3_Simulation.R

|- data/
|  |- colony_attributes.csv
|  |- empe_aerial_2023-05-25.csv
|  |- empe_satellite_2023-05-25.csv
|  |- fast_ice_trends.csv

## Sharing/Access information

This is a section for linking to other ways to access the data, and for linking to sources the data is derived from, if any.

Links to other publicly accessible locations of the data:
 - [http://...](http://...)

Data was derived from the following sources:
 - []()


## Code/Software

This is an optional, freeform section for describing any code in your submission and the software used to run it.

Describe any scripts, code, or notebooks (e.g., R, Python, Mathematica, MatLab) as well as the software versions (including loaded packages) that you used to run those files. If your repository contains more than one file whose relationship to other scripts is not obvious, provide information about the workflow that you used to run those scripts and notebooks.