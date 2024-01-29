# LaRue et al. (2024) Advances in remote sensing of emperor penguins: first multi-year time series documenting global population change

[Access this dataset on Dryad](Dataset DOI link)

# Overview

This repository contains data, code, and model output associated with the analysis of Emperor penguin population dynamics described in LaRue et al. (2024).

Raw data for the analysis are  included aerial and ground-based observations

Give a brief summary of dataset contents, contextualized in experimental procedures and results.

## File structure

- **analysis/** 
    - **output/** 
        - **data_viz/**    [visualizations of raw data] 
        - **model_checks/** [goodness-of-fit assessments]
        - **model_results/** 
            - **Global_Level/** [global estimates of population abundance, change, and trend]
            - **Regional_Level/** [regional estimates of population abundance, change, and trend]
            - **Colony_Level/** [colony-level estimates of population abundance, change, and trend]
            - **parameter_estimates.csv** [parameter estimates from fitted Bayesian model]
        - **simulation/** [raw results and figures related to simulations]
        - **EMPE_data_formatted.RData** [formatted data for analysis, after cleaning]
        - **fitted_model.RData** [Bayesian output in JAGS]
    - **script1_PrepareData.R** [R script to visualize, clean, and package raw data for analysis with JAGS]
    - **script2_FitModel.R** [R script to fit model, check convergence, evaluate goodness of fit, and generate results, figures, and summary statistics]
    - **script3_Simulation.R** [R script to conduct simulations; used for testing the Bayesian modeling approach]
    
- **data/** 
    - **colony_attributes.csv** [locations, names, and abbreviations for each colony]
    - **empe_aerial_2023-05-25.csv** [raw data associated with aerial counts of emperor penguin colonies]
    - **empe_satellite_2023-05-25.csv** [raw data associated with satellite surveys of emperor penguin colonies]
    - **fast_ice_trends.csv** [data related to regional trends in Antarctic fast ice, extracted from Table 1 in [Fraser et al. 2021](https://tc.copernicus.org/articles/15/5061/2021/)]

## Sharing/Access information

This is a section for linking to other ways to access the data, and for linking to sources the data is derived from, if any.

Links to other publicly accessible locations of the data:
 - [http://...](http://...)

Data was derived from the following sources:
 - []()


## Code/Software

This is an optional, freeform section for describing any code in your submission and the software used to run it.

Describe any scripts, code, or notebooks (e.g., R, Python, Mathematica, MatLab) as well as the software versions (including loaded packages) that you used to run those files. If your repository contains more than one file whose relationship to other scripts is not obvious, provide information about the workflow that you used to run those scripts and notebooks.