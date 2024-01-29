# LaRue et al. (2024): Advances in remote sensing of emperor penguins: first multi-year time series documenting global population change

[Access this dataset on Dryad](Dataset DOI link)

# Overview

This repository contains data, code, and model output associated with the global-scale analysis of Emperor penguin population dynamics described in LaRue et al. (2024), based on integrating raw data from aerial surveys with time series of circumpolar satellite surveys of known emperor penguin colonies.

The model is used to estimate an annual index of abundance at every known Emperor penguin colony in Antarctica (as of 2018), for every year between 2008 and 2018. Regional and global population indices are then calculated by summing colony-level estimates, according to regional colony membership.

Simulations are also performed to evaluate the ability of the model to accurately detect population trends, if they exist.

## File structure and code description

- **analysis/** 
    - **output/** 
        - **data_viz/** [visualizations of raw data] 
        - **model_checks/** [goodness-of-fit assessments]
        - **model_results/** 
            - **Global_Level/** [global estimates of population abundance, change, and trend]
            - **Regional_Level/** [regional estimates of population abundance, change, and trend]
            - **Colony_Level/** [colony-level estimates of population abundance, change, and trend]
            - **parameter_estimates.csv** [parameter estimates from fitted Bayesian model]
        - **simulation/** [raw results and figures related to simulations]
        - **EMPE_data_formatted.RData** [formatted data for analysis, after cleaning]
        - **fitted_model.RData** [Bayesian output in JAGS]
    - `script1_PrepareData.R` [R script to visualize, clean, and package raw data for analysis with JAGS]
    - `script2_FitModel.R` [R script to fit model, check convergence, evaluate goodness of fit, and generate results, figures, and summary statistics]
    - `script3_Simulation.R` [R script to conduct simulations; used for testing the Bayesian modeling approach]
    
- **data/** 
    - **colony_attributes.csv** [locations, names, and abbreviations for each colony]
    - **empe_aerial_2023-05-25.csv** [raw data associated with aerial counts of emperor penguin colonies]
    - **empe_satellite_2023-05-25.csv** [raw data associated with satellite surveys of emperor penguin colonies]
    - **fast_ice_trends.csv** [data related to regional trends in Antarctic fast ice, extracted from Table 1 in [Fraser et al. 2021](https://tc.copernicus.org/articles/15/5061/2021/)]

## Sharing/Access information

Data is openly accessible through Dryad (link).