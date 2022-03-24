# install/load necessary packages
my.packs <- c('jagsUI',"ggplot2",'reshape2',
              'scales','tidyverse',
              'ggrepel','ggthemes')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("~/1_Work/EMPE_Global/analysis")

rm(list=ls())

# --------------------------------------
# Load data package and fitted model
# --------------------------------------
load("output_empirical/EMPE_data_prepared.RData") # Data
load(file = "output_empirical/EMPE_out.RData")    # Fitted model

# --------------------------------------
# - Prepare a data package to use for simulating new data
# - Use sample sizes (e.g., n_sites, n_years, n_obs, etc) from observed data to
#   ensure simulations have "realistic" data quantity and balance
# --------------------------------------

jags.data.sim <- jags.data[c("n_years","n_sites",
                             "n_obs_aerial","aerial_site","aerial_year",
                             "n_obs_satellite","satellite_site","satellite_year","img_qual",
                             "regression_weights")]

# Use mean parameter estimates from fitted model to simulate new data
# (Colony-level initial abundances and trends will be simulated based on hyper-parameters)

jags.data.sim = append(jags.data.sim,
                       out$mean[c(
                         
                         # ------------------------
                         # Hyper-parameters from fitted model
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
                         "sat_p"                   # Probability satellite entirely fails to observe a colony that is present
                         
                       )])


# --------------------------------------
# # JAGS script used to simulate data, based on jags.data.sim
# --------------------------------------

sink("EMPE_model_simulate_data.jags")
cat("
    model {
  
  # ------------------------------------
  # Population process model
  # ------------------------------------
  
  for (s in 1:n_sites){
    for (t in 1:n_years){
      z_occ[s,t] ~ dbern(prob_occ)
    }
  }
  
  r_mean_grandmean_tau <- pow(r_mean_grandmean_sd,-2)
  logX1_tau <- pow(logX1_sd,-2)
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
  
  aerial_tau <- pow(aerial_sigma,-2)
  
  for (i in 1:n_obs_aerial){
    
    # Adult counts assumed to be poisson distributed with mean centred on seasonal expected value
    log_lambda[i] ~ dnorm(log_X[aerial_site[i],aerial_year[i]] - 1/(2*aerial_tau), aerial_tau)
    adult_count[i] ~ dpois(z_occ[aerial_site[i],aerial_year[i]] * exp(log_lambda[i]))
    
  }
  
  # ------------------------------------
  # Satellite observation model
  # ------------------------------------
  
  for (i in 1:n_obs_satellite){
    
    # Observation error (and bias) for satellite counts is estimated from data
    sat_mean[i] <- N[satellite_site[i],satellite_year[i]] * sat_slope[img_qual[i]]
    sat_sd[i] <- sat_mean[i] * sat_CV[img_qual[i]] + 0.001 # Tiny constant ensures tau is defined when sat_mean is zero
    sat_tau[i] <- pow(sat_sd[i],-2)
    
    sat_z[i] ~ dbern(sat_p)
    satellite[i] ~ dnorm(sat_mean[i] * sat_z[i],sat_tau[i])

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
  
}
    ",fill = TRUE)
sink()

# --------------------------------------
# Prepare an empty file to store results of simulations
# --------------------------------------

# Create empty output file to store simulation results
if (!file.exists("./output_simulation/sim_results.RData")){
  simulation_results <- vector("list",0)
  save(simulation_results, file = "./output_simulation/sim_results.RData") 
} else{
  load(file = "./output_simulation/sim_results.RData")
}

# --------------------------------------
# Conduct repeated simulations in which new data is simulated,
# and statistical model is applied to estimate global abundance and trend
# --------------------------------------

for (sim_run in seq(1,500,1)){
  
  set.seed(sim_run)
  
  if (length(simulation_results) >= sim_run) next
  
  # Create placeholder in list for this simulation run
  load(file = "./output_simulation/sim_results.RData")
  simulation_results[[sim_run]] <- vector("list",0)
  save(simulation_results, file = "./output_simulation/sim_results.RData")
  
  # --------------------------------------
  # Simulate data
  # --------------------------------------
  
  out_sim <- jags(data=jags.data.sim,
                  model.file="EMPE_model_simulate_data.jags",
                  parameters.to.save=c(
                    
                    "N",
                    "adult_count",
                    "satellite",
                    "N_global",
                    "global_trend"
                    
                  ),
                  inits = NULL,
                  n.chains=1,
                  n.thin = 1,
                  n.iter= 2,
                  n.burnin= 1,
                  parallel = TRUE)
  
  # ------------------------------------
  # Extract simulated data and plot
  # ------------------------------------
  sim_data <- out_sim$sims.list
  
  # True (simulated) values of annual N at each colony, in this simulation run
  N_df <- sim_data$N[1,,] %>% melt() %>% rename(site = Var1, year = Var2, N = value)
  
  # Aerial observations
  aerial_df_sim <- data.frame(y = sim_data$adult_count[1,],
                              year = jags.data.sim$aerial_year,
                              site = jags.data.sim$aerial_site,
                              type = "Aerial count (adult)")
  
  # Satellite observations
  satellite_df_sim <- data.frame(y = sim_data$satellite[1,],
                                 year = jags.data.sim$satellite_year,
                                 site = jags.data.sim$satellite_site,
                                 type = "Satellite")
  # Combine into a dataframe
  obs_sim <- rbind(aerial_df_sim,satellite_df_sim)
  
  # Plot observed data at each colony (for visual confirmation that the simulation is reasonable)
  sim_data_plot <- ggplot(obs_sim) +
    geom_line(data = N_df, aes(x = year, y = N))+
    geom_point(aes(x = year, y = y, shape = type))+
    scale_shape_manual(name = 'Obs type',
                       values =c(19,4))+
    
    facet_grid(site~., scales = "free")+
    theme_bw()
  
  # Global abundance each year
  N_global_df <- data.frame(year = 1:jags.data.sim$n_years,N = out_sim$sims.list$N_global[1,] )
  
  # Plot global abundance across the simulation
  global_N_plot <- ggplot(N_global_df) +
    geom_line(data = N_global_df, aes(x = year, y = N))+
    theme_bw()
  
  # --------------------------------------
  # Fit model to simulated data
  # --------------------------------------
  
  # Data to fit
  jags.data.refit = list( n_years = jags.data.sim$n_years,
                          n_sites = jags.data.sim$n_sites,
                          
                          # aerial counts of adults
                          n_obs_aerial = length(sim_data$adult_count[1,]),
                          adult_count = sim_data$adult_count[1,],
                          aerial_site = jags.data.sim$aerial_site,
                          aerial_year = jags.data.sim$aerial_year,
                          
                          # satellite counts
                          n_obs_satellite = length(sim_data$satellite[1,]),
                          satellite = sim_data$satellite[1,],
                          satellite_site = jags.data.sim$satellite_site,
                          satellite_year = jags.data.sim$satellite_year,
                          img_qual = jags.data.sim$img_qual,
                          
                          regression_weights = jags.data.sim$regression_weights
  )
  
  z_init <- matrix(1,ncol = jags.data.refit$n_years, nrow = jags.data.refit$n_sites)
  inits <- function()list(z_occ = z_init)
  
  out_refit <- jags(data=jags.data.refit,
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
                      "global_trend"
                    ),
                    
                    inits = inits,
                    n.chains=3,
                    n.thin = 5,
                    n.iter= 60000,
                    n.burnin= 10000,
                    parallel = TRUE)
  
  # -----------------------------------------
  # Extract estimates of global trend and change estimates
  # -----------------------------------------
  
  # Percent change estimates
  percent_change_est = 100*(out_refit$sims.list$N_global[,jags.data$n_years] - out_refit$sims.list$N_global[,1])/out_refit$sims.list$N_global[,1]
  percent_change_true = 100*(out_sim$sims.list$N_global[,jags.data$n_years] - out_sim$sims.list$N_global[,1])/out_sim$sims.list$N_global[,1]
  
  # Store in dataframe
  trend_results_summary <- data.frame(seed = sim_run,
                                      max_Rhat = max(unlist(out_refit$Rhat),na.rm = TRUE),
                                      global_trend_est_mean = mean(out_refit$sims.list$global_trend),
                                      global_trend_est_q025 = quantile(out_refit$sims.list$global_trend,0.025),
                                      global_trend_est_q500 = quantile(out_refit$sims.list$global_trend,0.500),
                                      global_trend_est_q975 = quantile(out_refit$sims.list$global_trend,0.975),
                                      global_trend_true = out_sim$sims.list$global_trend,
                                      
                                      percent_change_est_mean = mean(percent_change_est),
                                      percent_change_est_q025 = quantile(percent_change_est,0.025),
                                      percent_change_est_q500 = quantile(percent_change_est,0.500),
                                      percent_change_est_q975 = quantile(percent_change_est,0.975),
                                      percent_change_true = percent_change_true)
  
  # -----------------------------------------
  # Save output
  # -----------------------------------------
  
  load(file = "./output_simulation/sim_results.RData")
  simulation_results[[sim_run]] = trend_results_summary
  save(simulation_results,file = "./output_simulation/sim_results.RData")
  
}

# -----------------------------------------------------
# Summarize simulation results
# -----------------------------------------------------

load(file = "./output_simulation/sim_results.RData")
length(simulation_results)

trend_results_summary_all = data.frame()
for (i in 1:length(simulation_results)){
  
  if (is.null(simulation_results[[i]])) next
  
  trend_results_summary_all <- rbind(trend_results_summary_all,simulation_results[[i]])
  
}
dim(trend_results_summary_all)

# ---------------------------------------
# Trend (average annual percent change from 2009 to 2018)
# ---------------------------------------

# Convert to percent change per year using 100*(exp(OLS_regression_slope)-1)
trend_results_summary_all$global_trend_est_mean = 100*(exp(trend_results_summary_all$global_trend_est_mean)-1)
trend_results_summary_all$global_trend_est_q025 = 100*(exp(trend_results_summary_all$global_trend_est_q025)-1)
trend_results_summary_all$global_trend_est_q500 = 100*(exp(trend_results_summary_all$global_trend_est_q500)-1)
trend_results_summary_all$global_trend_est_q975 = 100*(exp(trend_results_summary_all$global_trend_est_q975)-1)
trend_results_summary_all$global_trend_true = 100*(exp(trend_results_summary_all$global_trend_true)-1)

# Credible interval coverage
trend_results_summary_all$trend_cov <- trend_results_summary_all$global_trend_est_q025 <= trend_results_summary_all$global_trend_true & trend_results_summary_all$global_trend_est_q975 >= trend_results_summary_all$global_trend_true
trend_coverage <- mean(trend_results_summary_all$trend_cov) %>% round(2)
mean_trend_bias <- mean(trend_results_summary_all$global_trend_est_mean - trend_results_summary_all$global_trend_true) %>% round(1)

ylim = range(trend_results_summary_all[,c("global_trend_true","global_trend_est_q025","global_trend_est_q975")])
trend_plot <- ggplot(trend_results_summary_all,aes(x = global_trend_true, 
                                                   y = global_trend_est_mean, 
                                                   ymin = global_trend_est_q025,
                                                   ymax = global_trend_est_q975, 
                                                   col = trend_cov))+
  geom_abline(slope = 1, col = "gray90")+
  geom_errorbar(width = 0, alpha = 0.3)+
  geom_point()+
  scale_color_manual(values = c("orangered","dodgerblue"), name = "95% CRI overlaps\ntrue trend")+
  coord_cartesian(ylim = c(ylim[1],ylim[2]), xlim = c(ylim[1],ylim[2]))+
  xlab("True (simulated) global trend")+
  ylab("Estimated global trend")+
  labs(title = paste0("Simulation results"),
       subtitle = paste0("Mean bias of trend estimate = ",mean_trend_bias,"%\n95% credible interval coverage = ",trend_coverage*100,"%"))+
  theme_bw()
print(trend_plot)

pdf("./output_simulation/trend_sim_results.pdf", width = 5, height = 4)
print(trend_plot)
dev.off()

# ---------------------------------------
# Total percent change (2018 vs 2009)
# ---------------------------------------

# Credible interval coverage
trend_results_summary_all$percent_change_cov <- trend_results_summary_all$percent_change_est_q025 <= trend_results_summary_all$percent_change_true & trend_results_summary_all$percent_change_est_q975 >= trend_results_summary_all$percent_change_true
percent_change_coverage <- mean(trend_results_summary_all$percent_change_cov) %>% round(2)
mean_percent_change_bias <- mean(trend_results_summary_all$percent_change_est_mean - trend_results_summary_all$percent_change_true) %>% round(1)

# Just to help identify a useful scale
ymax = 400
ymax_scale = log(ymax/100+1)
ymin_scale = -ymax_scale
ymin = 100*(exp(ymin_scale)-1)

ylim = 400
ylim = log(ylim/100 + 1)
y_labels = data.frame(pchange = c(-80,-50,0,100,400),labels = c("-80%","-50%","0","+100%","+400%"))
y_labels$r = log(y_labels$pchange/100 + 1)

percent_change_plot <- ggplot(trend_results_summary_all,
                              aes(x = log(percent_change_true/100 + 1), 
                                  y = log(percent_change_est_mean/100 + 1), 
                                  ymin = log(percent_change_est_q025/100 + 1),
                                  ymax = log(percent_change_est_q975/100 + 1), 
                                  col = percent_change_cov))+
  geom_abline(slope = 1, col = "gray90")+
  geom_errorbar(width = 0, alpha = 0.3)+
  geom_point()+
  scale_color_manual(values = c("orangered","dodgerblue"), name = "95% CRI overlaps\ntrue trend")+
  scale_y_continuous(breaks = y_labels$r, labels = y_labels$labels)+
  scale_x_continuous(breaks = y_labels$r, labels = y_labels$labels)+
  
  coord_cartesian(ylim=c(-ylim,ylim),xlim = c(-ylim,ylim))+
  xlab("True (simulated) percent\nchange from 2009 to 2018")+
  ylab("Estimated percent\nchange from 2009 to 2018")+
  labs(title = paste0("Simulation results"),
       subtitle = paste0("Mean bias of change estimate = ",mean_percent_change_bias,"%\n95% credible interval coverage = ",percent_change_coverage*100,"%"))+
  theme_bw()
print(percent_change_plot)

pdf("./output_simulation/change_sim_results.pdf", width = 5, height = 4)
print(percent_change_plot)
dev.off()