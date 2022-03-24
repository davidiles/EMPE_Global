# install/load necessary packages
my.packs <- c('jagsUI',"ggplot2",'reshape2','scales','tidyverse',
              'rgeos','raster','sp','sf','ggrepel')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("~/1_Work/EMPE_Global/analysis")

rm(list=ls())

load("output_empirical/EMPE_data_prepared.RData")

# The jags script to fit the model
sink("EMPE_model_empirical.jags")
cat("
    model {
  
  # ------------------------------------
  # Population process model
  # ------------------------------------
  
  # Occupancy state in each year, at each colony
  prob_occ ~ dunif(0,1)
  for (s in 1:n_sites){
    for (t in 1:n_years){
      z_occ[s,t] ~ dbern(prob_occ)
    }
  }
  
  # SPATIAL RANDOM EFFECTS FOR r_mean[s]
  r_mean_grandmean_mu ~ dnorm(0,1)
  r_mean_grandmean_sd ~ dunif(0,2)
  r_mean_grandmean_tau <- pow(r_mean_grandmean_sd,-2)
  
  # SPATIAL RANDOM EFFECTS FOR X[s,1]
  X1_mean ~ dunif(0,50000)
  logX1_mean <- log(X1_mean)
  logX1_sd ~ dunif(0,2)
  logX1_tau <- pow(logX1_sd,-2)
  
  # Population growth process
  r_sd ~ dunif(0,2) 
  r_tau <- pow(r_sd,-2)
  
  for (s in 1:n_sites){
    
    r_mean[s] ~ dnorm(r_mean_grandmean_mu,r_mean_grandmean_tau) # Drawn from shared distribution
    log_X[s,1] ~ dnorm(logX1_mean,logX1_tau)
    N[s,1] <- exp(log_X[s,1]) * z_occ[s,1]
    
    for (t in 1:(n_years-1)){
      
      log_X[s,t+1] ~ dnorm(log_X[s,t] + r_mean[s] - 1/(2*r_tau), r_tau)
      N[s,t+1] <- exp(log_X[s,t+1]) * z_occ[s,t+1]
    }
  } 
  
  # ------------------------------------
  # Aerial observation model
  # ------------------------------------
  
  aerial_sigma ~ dunif(0,2) 
  aerial_tau <- pow(aerial_sigma,-2)
  
  for (i in 1:n_obs_aerial){
    
    # Adult counts assumed to be poisson distributed with mean centred on seasonal expected value
    log_lambda[i] ~ dnorm(log_X[aerial_site[i],aerial_year[i]] - 1/(2*aerial_tau), aerial_tau)
    adult_count[i] ~ dpois(z_occ[aerial_site[i],aerial_year[i]] * exp(log_lambda[i]))
    
    # ------------------------------------
    # FOR POSTERIOR PREDICTIVE CHECK
    # ------------------------------------
    sim_log_lambda[i] ~ dnorm(log_X[aerial_site[i],aerial_year[i]] - 1/(2*aerial_tau), aerial_tau)
    sim_adult_count[i] ~ dpois(z_occ[aerial_site[i],aerial_year[i]] * exp(sim_log_lambda[i]))
    
    # Expected values
    expected_adult_count[i] <- prob_occ * exp(log_X[aerial_site[i],aerial_year[i]] + 1/(2*aerial_tau))
    
    # Discrepancy measures
    sqE_adult_count_actual[i] <- pow(adult_count[i] - expected_adult_count[i],2)
    sqE_adult_count_sim[i] <- pow(sim_adult_count[i] - expected_adult_count[i],2)
  }
  
  # ------------------------------------
  # Satellite observation model
  # ------------------------------------
  
  # Describes proportional bias and variance in satellite counts for each level of image quality (1, 2, or 3)
  sat_CV[1] ~ dunif(0,2)
  sat_CV[2] ~ dunif(0,2)
  sat_CV[3] ~ dunif(0,2)
  
  sat_slope[1] ~ dnorm(1,25)
  sat_slope[2] ~ dnorm(1,25)
  sat_slope[3] ~ dnorm(1,25)
  
  sat_p ~ dunif(0,1)
  
  for (i in 1:n_obs_satellite){
    
    # Observation error (and bias) for satellite counts is estimated from data
    sat_mean[i] <- N[satellite_site[i],satellite_year[i]] * sat_slope[img_qual[i]]
    sat_sd[i] <- sat_mean[i] * sat_CV[img_qual[i]] + 0.001 # Tiny constant ensures tau is defined when sat_mean is zero
    sat_tau[i] <- pow(sat_sd[i],-2)
    
    sat_z[i] ~ dbern(sat_p)
    satellite[i] ~ dnorm(sat_mean[i] * sat_z[i],sat_tau[i])
    
    #------------------------------------
    # FOR POSTERIOR PREDICTIVE CHECK 
    #------------------------------------
    sim_sat_z[i] ~ dbern(sat_p)
    sim_satellite[i] ~ dnorm(sat_mean[i] * sim_sat_z[i],sat_tau[i]) # Simulate new satellite obs (based on fitted estimates of population)
    
    expected_satellite[i] <- prob_occ * sat_mean[i]
    
    # Discrepancy measures
    sqE_satellite_actual[i] <- pow(satellite[i] - expected_satellite[i],2)
    sqE_satellite_sim[i] <- pow(sim_satellite[i] - expected_satellite[i],2)
    #------------------------------------
  }
  
  # ------------------------------------
  # Derived quantities
  # ------------------------------------
  
  # Global abundance each year
  for (t in 1:n_years){
    N_global[t] <- sum(N[1:n_sites,t])
  }
  
  # Global trend (measured as least squares regression slope on log scale)
  global_trend <- inprod(log(N_global[1:n_years]),regression_weights[1,1:n_years])
  
  #------------------------------------
  # FOR POSTERIOR PREDICTIVE CHECK 
  #------------------------------------
  
  # Root mean squared error of empirical data
  RMSE_adult_count_actual <- sqrt(mean(sqE_adult_count_actual[]))
  RMSE_satellite_actual <- sqrt(mean(sqE_satellite_actual[]))
  
  # Root mean squared error of simulated data
  RMSE_adult_count_sim <- sqrt(mean(sqE_adult_count_sim[]))
  RMSE_satellite_sim <- sqrt(mean(sqE_satellite_sim[]))
  
}
    ",fill = TRUE)
sink()

out <- jags(data=jags.data,
            model.file="EMPE_model_empirical.jags",
            parameters.to.save=c(
              
              # ------------------------
              # Hyper-parameters
              # ------------------------
              "prob_occ",               # Probability colonies are "occupied"
              "r_mean_grandmean_mu",    # Hierarchical grand mean of colony-level annual growth rates
              "r_mean_grandmean_sd",    # Hierarchical sd of colony-level growth rates
              "logX1_mean",             # Hierarchical grand mean of colony-level initial abundances
              "logX1_sd",               # Hierarchical sd of colony-level initial abundances
              "r_sd",                   # Temporal variance of colony-level annual growth rates
              "aerial_sigma",           # SD of aerial observations (on log scale)
              "sat_CV",                 # Coefficient of variation in satellite observations
              "sat_slope",              # Bias in satellite observations
              "sat_p",                  # Probability satellite entirely fails to observe a colony that is present
              
              # Colony-specific mean trend
              "r_mean",
              
              # Colony-specific abundance each year
              "N",
              
              # Global abundance each year
              "N_global",
              
              # Log-linear OLS slope across the study
              "global_trend",
              
              # Discrepancy measures for posterior predictive checks
              "RMSE_adult_count_actual",
              "RMSE_satellite_actual",
              "RMSE_adult_count_sim",
              "RMSE_satellite_sim"
            ),
            
            inits = inits,
            n.chains=3,
            n.thin = 50,
            n.iter= 600000,
            n.burnin= 100000,
            parallel = TRUE)

save(out, file = "output_empirical/EMPE_out.RData")
