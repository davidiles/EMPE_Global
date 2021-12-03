# install/load necessary packages
my.packs <- c('jagsUI',"ggplot2",'reshape2','scales','tidyverse',
              'rgeos','raster','sp','sf','ggrepel')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("~/1_Work/GoogleDrive_emperor_satellite/GoogleDrive_emperorsatellite/EMPE_Global/analysis")

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
  
  # Describes proportional bias and variance in satellite counts
  sat_CV ~ dunif(0,1)
  sat_slope ~ dnorm(1,25)
  sat_p ~ dunif(0,1)
  
  for (i in 1:n_obs_satellite){
    
    # Observation error (and bias) for satellite counts is estimated from data
    sat_mean[i] <- N[satellite_site[i],satellite_year[i]] * sat_slope 
    sat_sd[i] <- sat_mean[i] * sat_CV + 0.001 # Tiny constant ensures tau is defined when sat_mean is zero
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
              
              # Hyper-parameters describing population sizes and growth rates
              "r_mean_grandmean_mu",
              "r_mean_grandmean_sd",
              "logX1_mean",
              "logX1_sd",
              
              # Population dynamics parameters
              "r_mean",
              "r_sd",
              "prob_occ",
              
              # Observation parameters
              "aerial_sigma",
              "sat_CV",
              "sat_slope",
              "sat_p",
              
              # Estimates of abundance and trend
              "N_global",
              "global_trend",
              "N",
              "site_trend",
              
              # Discrepancy measures for posterior predictive checks
              "RMSE_adult_count_actual",
              "RMSE_satellite_actual",
              "RMSE_adult_count_sim",
              "RMSE_satellite_sim"
            ),
            
            inits = inits,
            n.chains=3,
            n.thin = 50,
            n.iter= 100000,
            n.burnin= 50000,
            parallel = TRUE)

save(out, file = "output_empirical/EMPE_out.RData")
