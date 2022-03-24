# install/load necessary packages
my.packs <- c('jagsUI')

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
  
  # Spatial random effect for X[s,1]
  X1_mean ~ dunif(0,50000)
  logX1_mean <- log(X1_mean)
  logX1_sd ~ dunif(0,2)
  logX1_tau <- pow(logX1_sd,-2)
  
  # Spatial random effect for r_mean[s]
  r_mean_grandmean_mu ~ dnorm(0,1)
  r_mean_grandmean_sd ~ dunif(0,2)
  r_mean_grandmean_tau <- pow(r_mean_grandmean_sd,-2)
  
  # Temporal variance in population growth
  r_sd ~ dunif(0,2) 
  r_tau <- pow(r_sd,-2)
  
  # Growth dynamics at each colony
  for (s in 1:n_sites){
    
    # Initial abundance at colony 
    log_X[s,1] ~ dnorm(logX1_mean,logX1_tau)
    
    # Median annual change at colony
    r_mean[s] ~ dnorm(r_mean_grandmean_mu,r_mean_grandmean_tau) 
    
    # Population growth from year to year at colony
    for (t in 1:(n_years-1)){
      log_X[s,t+1] ~ dnorm(log_X[s,t] + r_mean[s], r_tau)
    }
  } 
  
  # Probability colony is occupied each year
  prob_occ ~ dunif(0,1)
  
  # Occupancy dynamics at each colony (not a Markov process)
  for (s in 1:n_sites){
    for (t in 1:n_years){
    
      # True occupancy state at colony in each year
      z_occ[s,t] ~ dbern(prob_occ)
      
    }
  }
  
  # ------------------------------------
  # Annual population indices
  # ------------------------------------
  
  for (s in 1:n_sites){
    for (t in 1:n_years){
      N[s,t] <- exp(log_X[s,t]) * z_occ[s,t]
    }
  }
  
  # ------------------------------------
  # Aerial observation model
  # ------------------------------------
  
  aerial_sigma ~ dunif(0,2) 
  aerial_tau <- pow(aerial_sigma,-2)
  
  for (i in 1:n_obs_aerial){
    
    # Adult counts assumed to be poisson distributed with mean centred on seasonal expected value
    lambda[i] ~ dlnorm(log_X[aerial_site[i],aerial_year[i]] - 0.5*pow(aerial_sigma,2), aerial_tau)
    adult_count[i] ~ dpois(z_occ[aerial_site[i],aerial_year[i]] * lambda[i])
    
  }
  
  # ------------------------------------
  # Satellite observation model
  # ------------------------------------
  
  # Describes proportional bias (slope) and variance in satellite counts for each level of image quality (1, 2, or 3)
  for (i in 1:3){
    sat_slope[i] ~ dnorm(1,25)
    sat_CV[i] ~ dunif(0,2)
  }
  
  # Probability a satellite fails to detect a colony that is actually present
  sat_p ~ dunif(0,1)
  
  for (i in 1:n_obs_satellite){
    
    # Colony detection state (1 = colony was detected)
    sat_z[i] ~ dbern(sat_p)
    
    # Observation error (and bias) for satellite counts is estimated from data
    sat_mean[i] <- N[satellite_site[i],satellite_year[i]] * sat_slope[img_qual[i]]
    sat_sd[i] <- sat_mean[i] * sat_CV[img_qual[i]] + 0.001 # Tiny constant ensures tau is defined when N[s,y] is zero
    sat_tau[i] <- pow(sat_sd[i],-2)
    
    satellite[i] ~ dnorm(sat_mean[i]* sat_z[i],sat_tau[i])
  }
  
  # ------------------------------------
  # Global change and trend
  # ------------------------------------
  
  # Global abundance each year
  for (t in 1:n_years){
    N_global[t] <- sum(N[1:n_sites,t])
  }
  
  # Global trend (measured as least squares regression slope on log scale)
  global_trend <- inprod(log(N_global[1:n_years]),regression_weights[1,1:n_years])
  
  #------------------------------------
  # Posterior predictive check
  #------------------------------------
  
  # Posterior predictive check for aerial count data
  for (i in 1:n_obs_aerial){
  
    # Simulate aerial observations
    sim_lambda[i] ~ dlnorm(log_X[aerial_site[i],aerial_year[i]] - 0.5*pow(aerial_sigma,2), aerial_tau)
    sim_adult_count[i] ~ dpois(z_occ[aerial_site[i],aerial_year[i]] * sim_lambda[i])
    
    # Calculate discrepancy measures of actual and simulated data
    sqE_adult_count_actual[i] <- pow(adult_count[i] - N[aerial_site[i],aerial_year[i]] ,2)
    sqE_adult_count_sim[i] <- pow(sim_adult_count[i] - N[aerial_site[i],aerial_year[i]],2)
    
  }
  
  # Posterior predictive check for satellite data
  for (i in 1:n_obs_satellite){
  
    # Simulate satellite observations
    sim_sat_z[i] ~ dbern(sat_p)
    sim_satellite[i] ~ dnorm(sat_mean[i] * sim_sat_z[i],sat_tau[i]) # Simulate new satellite obs (based on fitted estimates of population)
    
    # Calculate discrepancy measures of actual and simulated data
    sqE_satellite_actual[i] <- pow(satellite[i] - N[satellite_site[i],satellite_year[i]],2)
    sqE_satellite_sim[i] <- pow(sim_satellite[i] - N[satellite_site[i],satellite_year[i]],2)
    
  }
  
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
              
              # Colony-specific mean annual differences
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
            n.thin = 5,
            n.iter= 60000,
            n.burnin= 10000,
            parallel = TRUE)

save(out, file = "output_empirical/EMPE_out.RData")
