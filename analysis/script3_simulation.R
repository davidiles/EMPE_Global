# install/load necessary packages
my.packs <- c('jagsUI',"ggplot2",'reshape2','scales','tidyverse',
              'rgeos','raster','sp','sf','ggrepel')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("~/1_Work/GoogleDrive_emperor_satellite/GoogleDrive_emperorsatellite/Global_Analysis/analysis/version_3/")

rm(list=ls())

# --------------------------------------
# Load fitted model
# --------------------------------------

load(file ="./output_empirical/empirical_output.RData")
out <- results_output$out
jags.data <- results_output$jags.data

# Create output file
if (!file.exists("./output_simulation/sim_results.RData")){
  simulation_results <- vector("list",0)
  save(simulation_results, file = "./output_simulation/sim_results.RData") 
} else{
  load(file = "./output_simulation/sim_results.RData")
}

for (sim_run in seq(1,200,1)){
  
  set.seed(sim_run)
  
  if (length(simulation_results) >= sim_run) next
  
  # Create placeholder in list for this simulation run
  load(file = "./output_simulation/sim_results.RData")
  simulation_results[[sim_run]] <- vector("list",0)
  simulation_results[[sim_run]]$seed <- sim_run
  save(simulation_results, file = "./output_simulation/sim_results.RData") 
  
  # --------------------------------------
  # Simulate data
  # --------------------------------------
  jags.data.sim <- jags.data[c("n_years","n_sites",
                               
                               "n_obs_aerial","aerial_site","aerial_year",
                               "n_obs_satellite","satellite_site","satellite_year")]
  
  # --------------------------------------
  # Simulation parameters
  # --------------------------------------
  
  # Simulate random "trends" at each colony in each simulation rep
  jags.data.sim$r_mean <- rnorm(jags.data$n_sites,runif(1,-0.2,-0.05),0.1)
  
  # Also shuffle initial abundances at colonies
  jags.data.sim$N0 <- apply(out$sims.list$N[,,1],2,mean) %>% sample(size = jags.data.sim$n_sites,replace = TRUE)
  
  jags.data.sim$r_sd <- mean(out$sims.list$r_sd)
  
  jags.data.sim$prob_occ <- mean(out$sims.list$prob_occ)
  
  jags.data.sim$aerial_sigma <- mean(out$sims.list$aerial_sigma)
  jags.data.sim$sat_CV <- mean(out$sims.list$sat_CV) 
  jags.data.sim$sat_slope <- mean(out$sims.list$sat_slope)
  jags.data.sim$sat_p <- mean(out$sims.list$sat_p)
  
  # For OLS regression on annual global index
  XX=cbind(rep(1,jags.data.sim$n_years),1:jags.data.sim$n_years)
  jags.data.sim$regression_weights <- matrix(c(0,1),1,2)%*%solve(t(XX)%*%XX)%*%t(XX)
  
  # The jags script to fit the model
  sink("EMPE_simdata.jags")
  cat("
    model {

    ###############################
    # Process model
    ###############################
    
    # Occupancy state in each year, at each colony
    for (s in 1:n_sites){
      for (t in 1:n_years){
        z_occ[s,t] ~ dbern(prob_occ)
      }
    }
    
    # Population growth process
    r_tau <- pow(r_sd,-2)
    
    for (s in 1:n_sites){
    
      log_X[s,1] <- log(N0[s])
      N[s,1] <- exp(log_X[s,1]) * z_occ[s,1]

      for (t in 1:(n_years-1)){
      
        log_X[s,t+1] ~ dnorm(log_X[s,t] + r_mean[s] + 0.5 * r_sd * r_sd, r_tau)
        N[s,t+1] <- exp(log_X[s,t+1]) * z_occ[s,t+1]
        
      }
      
    } 

    ###############################
    # Generate aerial observations
    ###############################
    
    aerial_tau <- pow(aerial_sigma,-2)
   
    for (i in 1:n_obs_aerial){
    
      # Adult counts assumed to be poisson distributed with mean centred on seasonal expected value
      log_lambda[i] ~ dnorm(log_X[aerial_site[i],aerial_year[i]] - 0.5 * aerial_sigma * aerial_sigma, aerial_tau)
      adult_count[i] ~ dpois(z_occ[aerial_site[i],aerial_year[i]] * exp(log_lambda[i]))
      
      adult_count_sim[i] ~ dpois(z_occ[aerial_site[i],aerial_year[i]] * exp(log_lambda[i]))
    }
    
    ###############################
    # Generate satellite observations
    ###############################
    
    for (i in 1:n_obs_satellite){
    
      # Observation error (and bias) for satellite counts is estimated from data
      sat_mean[i] <- N[satellite_site[i],satellite_year[i]] * sat_slope 
      sat_sd[i] <- sat_mean[i] * sat_CV + 0.001
      sat_tau[i] <- pow(sat_sd[i],-2)
      
      sat_z[i] ~ dbern(sat_p)
      satellite[i] ~ dnorm(sat_mean[i] * sat_z[i],sat_tau[i])
      
      satellite_sim[i] ~ dnorm(sat_mean[i] * sat_z[i],sat_tau[i])
    }
    
    ###############################
    # Derived quantities
    ###############################
    
    # Global abundance each year
    for (t in 1:n_years){
      N_global[t] <- sum(N[1:n_sites,t])
    }
    
    # Global trend
    global_trend <- inprod(log(N_global[1:n_years]),regression_weights[1,1:n_years])

    }
    ",fill = TRUE)
  sink()
  
  out_sim <- jags(data=jags.data.sim,
                  model.file="EMPE_simdata.jags",
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
  
  # ***************
  # True (simulated) values of annual N at each colony, in this simulation run
  # ***************
  N_df <- sim_data$N[1,,] %>% melt() %>% rename(site = Var1, year = Var2, N = value)
  
  # ***************
  # Observed values
  # ***************
  aerial_df_sim <- data.frame(y = sim_data$adult_count[1,],
                              year = jags.data.sim$aerial_year,
                              site = jags.data.sim$aerial_site,
                              type = "Aerial count (adult)")
  satellite_df_sim <- data.frame(y = sim_data$satellite[1,],
                                 year = jags.data.sim$satellite_year,
                                 site = jags.data.sim$satellite_site,
                                 type = "Satellite")
  
  obs_sim <- rbind(aerial_df_sim,satellite_df_sim)
  
  sim_data_plot <- ggplot(obs_sim) + 
    geom_line(data = N_df, aes(x = year, y = N))+
    geom_point(aes(x = year, y = y, shape = type))+
    scale_shape_manual(name = 'Obs type',
                       values =c(19,4))+
    
    facet_grid(site~., scales = "free")+
    theme_bw()
  
  N_global_df <- data.frame(year = 1:jags.data.sim$n_years,
                            N = out_sim$sims.list$N_global[1,] )
  
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
                          
                          regression_weights = jags.data.sim$regression_weights
  )
  
  #---------------------------------------------------------------------
  # Specify colony occupancy state when this is known 
  #---------------------------------------------------------------------
  
  z_init <- matrix(1,ncol = jags.data.refit$n_years, nrow = jags.data.refit$n_sites)
  
  inits <- function()list(z_occ = z_init)
  out_refit <- jags(data=jags.data.refit,
                    model.file="EMPE_model_empirical.jags",
                    parameters.to.save=c(
                      "r_mean_grandmean_mu",
                      "r_mean_grandmean_sd",
                      "logX1_mean",
                      "logX1_sd",
                      "r_sd",
                      "aerial_sigma",
                      "global_trend",
                      "N_global",
                      
                      "prob_occ",
                      "sat_CV",
                      "sat_slope",
                      "sat_p",
                      
                      "r_mean",
                      "N"
                    ),
                    inits = inits,
                    n.chains=3,
                    n.thin = 10,
                    n.iter= 100000,
                    n.burnin= 50000,
                    parallel = TRUE)
  
  # -----------------------------------------
  # Summarize results of simulation / compare fitted model estimates to 'truth'
  # -----------------------------------------
  
  # Trend estimates
  trend_results_summary <- data.frame(seed = sim_run,
                                      max_Rhat = max(unlist(out_refit$Rhat),na.rm = TRUE),
                                      global_trend_est_mean = mean(out_refit$sims.list$global_trend),
                                      global_trend_est_q025 = quantile(out_refit$sims.list$global_trend,0.025),
                                      global_trend_est_q500 = quantile(out_refit$sims.list$global_trend,0.500),
                                      global_trend_est_q975 = quantile(out_refit$sims.list$global_trend,0.975),
                                      global_trend_true = out_sim$sims.list$global_trend)
  
  N_df_est <- out_refit$sims.list$N %>% melt() %>%
    rename(samp = Var1, site = Var2, year = Var3, N_est = value) %>%
    full_join(N_df)
  N_df_est$bias <- N_df_est$N_est - N_df_est$N
  
  N_sim_global <- N_df_est %>%
    group_by(year,samp) %>%
    summarize(N_est = sum(N_est),
              N = sum(N),
              logN_bias = log(sum(N_est))-log(sum(N)))
  
  N_sim_global_summary <- N_sim_global %>%  
    group_by(year) %>%
    summarize(N = mean(N),
              N_est_mean = mean(N_est),
              N_est_q500 = quantile(N_est,0.5),
              N_est_q025 = quantile(N_est,0.025),
              N_est_q975 = quantile(N_est,0.975),
              
              mean_logN_bias = mean(logN_bias),
              median_logN_bias = median(logN_bias)) %>% 
    add_column(sim_run = sim_run)
  
  # -----------------------------------------
  # Save output
  # -----------------------------------------
  
  load(file = "./output_simulation/sim_results.RData")
  simulation_results[[sim_run]]$trend_results_summary <- trend_results_summary
  simulation_results[[sim_run]]$N_sim_global_summary <- N_sim_global_summary
  simulation_results[[sim_run]]$sat_slope_true <- jags.data.sim$sat_slope
  simulation_results[[sim_run]]$sat_slope_bias <- mean(out_refit$sims.list$sat_slope -jags.data.sim$sat_slope)
  simulation_results[[sim_run]]$sat_CV_true <- jags.data.sim$sat_CV
  simulation_results[[sim_run]]$sat_CV_bias <- mean(out_refit$sims.list$sat_CV -jags.data.sim$sat_CV)
  
  simulation_results[[sim_run]]$max_Rhat <- max(unlist(out_refit$Rhat),na.rm = TRUE)
  save(simulation_results,file = "./output_simulation/sim_results.RData")
  
}


# -----------------------------------------------------
# Visualize results of simulations
# -----------------------------------------------------

load(file = "./output_simulation/sim_results.RData")
length(simulation_results)
trend_results_summary_all <- N_sim_global_summary_all <- data.frame()

sat_slope_bias_vec <- sat_CV_bias_vec <- sat_slope_vec <- sat_CV_vec <- c()

for (i in 1:length(simulation_results)){
  
  if (is.null(simulation_results[[i]]$trend_results_summary)) next
  if(simulation_results[[i]]$max_Rhat > 1.3) next
  trend_results_summary_all <- rbind(trend_results_summary_all,simulation_results[[i]]$trend_results_summary)
  N_sim_global_summary_all <- rbind(N_sim_global_summary_all,simulation_results[[i]]$N_sim_global_summary)
  
  sat_slope_bias_vec <- c(sat_slope_bias_vec , mean(simulation_results[[i]]$sat_slope_bias))
  sat_CV_bias_vec <- c(sat_CV_bias_vec , mean(simulation_results[[i]]$sat_CV_bias))
  sat_slope_vec <- c(sat_slope_vec , mean(simulation_results[[i]]$sat_slope_true))
  sat_CV_vec <- c(sat_CV_vec , mean(simulation_results[[i]]$sat_CV_true))
  
}
nsims_converged <- length(unique(trend_results_summary_all$seed))

# Credible interval coverage
trend_results_summary_all$cov <- trend_results_summary_all$global_trend_est_q025 <= trend_results_summary_all$global_trend_true & trend_results_summary_all$global_trend_est_q975 >= trend_results_summary_all$global_trend_true

trend_coverage <- mean(trend_results_summary_all$cov) %>% round(3)
mean_trend_bias <- mean(trend_results_summary_all$global_trend_est_q500 - trend_results_summary_all$global_trend_true) %>% round(4)

trend_plot <- ggplot(trend_results_summary_all,aes(x = global_trend_true, y = global_trend_est_q500, ymin = global_trend_est_q025,ymax = global_trend_est_q975, col = cov))+
  geom_abline(slope = 1, col = "gray90")+
  geom_point(alpha = 0.3)+
  geom_errorbar(width = 0, alpha = 0.3)+
  scale_color_manual(values = c("orangered","dodgerblue"), name = "95% CRI overlaps\ntrue trend")+
  coord_cartesian(ylim = c(-0.3,0.3), xlim = c(-0.3,0.3))+
  xlab("True global trend (simulated)")+
  ylab("Estimated global trend")+
  labs(title = paste0("Simulation results (# sims = ",nsims_converged,")"),
       subtitle = paste0("Mean bias of trend estimate = ",mean_trend_bias,"\n95% credible interval coverage = ",trend_coverage))+
  theme_bw()
print(trend_plot)

pdf("./output_simulation/trend_sim_results.pdf", width = 5, height = 4)
print(trend_plot)
dev.off()

# Mean bias in logN (evaluated from median estimate in each iteration)
N_sim_global_summary_all$N_bias <- (N_sim_global_summary_all$N_est_q500 - N_sim_global_summary_all$N)/N_sim_global_summary_all$N

mean_percent_bias_N <- mean(N_sim_global_summary_all$N_bias * 100) %>% round(1)
N_sim_global_summary_all$cov <- N_sim_global_summary_all$N_est_q025 <=  N_sim_global_summary_all$N & N_sim_global_summary_all$N_est_q975 >=  N_sim_global_summary_all$N
N_coverage <- mean(N_sim_global_summary_all$cov) %>% round(2)

logN_bias_plot <- ggplot(data = N_sim_global_summary_all,aes(x = factor(year), y = mean_logN_bias))+
  geom_violin(fill = "dodgerblue",alpha = 0.3, draw_quantiles = c(0.5))+
  xlab("Year of simulation")+
  ylab("Bias in estimate of logN[t]")+
  labs(title = "Simulation results",
       subtitle = paste0("Mean bias of N estimate = ",mean_percent_bias_N,"%\n95% credible interval coverage of N[t] = ",N_coverage))+
  theme_bw()
#logN_bias_plot


mean_bias_df <- N_sim_global_summary_all %>%
  group_by(year) %>%
  summarize(N_bias = mean(N_bias))
max_limit <- max(abs(N_sim_global_summary_all$N_bias*100))
N_bias_plot <- ggplot(data = N_sim_global_summary_all,aes(x = year, y = N_bias*100, col = factor(sim_run)))+
  geom_hline(yintercept = 0, linetype = 2)+
  geom_line(data = mean_bias_df, aes(x = year, y = N_bias*100), col = "blue", size = 1.2)+
  geom_line(alpha = 0.1)+
  scale_color_manual(values = rep("blue",length(unique(N_sim_global_summary_all$sim_run))), guide = FALSE)+
  xlab("Year of simulation")+
  ylab("Bias in estimate of logN[t]")+
  labs(title = paste0("Simulation results (# sims = ",nsims_converged,")"),
       subtitle = paste0("Mean bias of N estimate = ",mean_percent_bias_N,"%\n95% credible interval coverage of N[t] = ",N_coverage))+
  theme_bw()+
  scale_y_continuous(breaks = seq(-50,50,25), limits = c(-max_limit,max_limit), labels = c("-50%","-25%","0%","+25%","+50%"))+
  scale_x_continuous(breaks = seq(1,10), labels = seq(1,10), minor_breaks = NULL)
N_bias_plot

pdf("./output_simulation/N_bias_results.pdf", width = 5, height = 4)
print(N_bias_plot)
dev.off()

hist(sat_CV_bias_vec)
mean(sat_CV_bias_vec)
mean(sat_CV_vec)
hist(sat_slope_bias_vec)
mean(sat_slope_bias_vec)
mean(sat_slope_vec)
