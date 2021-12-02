# install/load necessary packages
my.packs <- c('jagsUI',"ggplot2",'reshape2','scales','tidyverse',
              'rgeos','raster','sp','sf','ggrepel')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("~/1_Work/GoogleDrive_emperor_satellite/GoogleDrive_emperorsatellite/Global_Analysis/analysis/version_3/")

rm(list=ls())

load("EMPE_data_prepared.RData")

# The jags script to fit the model
sink("./EMPE_model_empirical.jags")
cat("
    model {

    ###############################
    # Process model
    ###############################
    
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
    X1_mean ~ dunif(1,50000)
    logX1_mean <- log(X1_mean)
    logX1_sd ~ dunif(0,2)
    logX1_tau <- pow(logX1_sd,-2)
    
    # Population growth process
    r_sd ~ dunif(0,2)      # sd of temporal growth rate
    r_tau <- pow(r_sd,-2)
    
    for (s in 1:n_sites){
    
      r_mean[s] ~ dnorm(r_mean_grandmean_mu,r_mean_grandmean_tau) # Drawn from shared distribution
      log_X[s,1] ~ dnorm(logX1_mean,logX1_tau)
      N[s,1] <- exp(log_X[s,1]) * z_occ[s,1]

      for (t in 1:(n_years-1)){
      
        log_X[s,t+1] ~ dnorm(log_X[s,t] + r_mean[s] + 0.5 * r_sd * r_sd, r_tau)
        N[s,t+1] <- exp(log_X[s,t+1]) * z_occ[s,t+1]
        
      }
      
    } 

    ###############################
    # Aerial observation model
    ###############################
    
    aerial_sigma ~ dunif(0,2)
    aerial_tau <- pow(aerial_sigma,-2)
   
    for (i in 1:n_obs_aerial){
    
      # Adult counts assumed to be poisson distributed with mean centred on seasonal expected value
      log_lambda[i] ~ dnorm(log_X[aerial_site[i],aerial_year[i]] - 0.5 * aerial_sigma * aerial_sigma, aerial_tau)
      adult_count[i] ~ dpois(z_occ[aerial_site[i],aerial_year[i]] * exp(log_lambda[i]))
      
      #------------------------------------
      # FOR POSTERIOR PREDICTIVE CHECK #1
      #------------------------------------
      sim_log_lambda_1[i] ~ dnorm(log_X[aerial_site[i],aerial_year[i]] - 0.5 * aerial_sigma * aerial_sigma, aerial_tau)
      sim_adult_count_1[i] ~ dpois(z_occ[aerial_site[i],aerial_year[i]] * exp(sim_log_lambda_1[i]))
      
      # Expected values
      expected_adult_count_1[i] <- prob_occ * exp(log_X[aerial_site[i],aerial_year[i]] - 0.5 * aerial_sigma * aerial_sigma)
      
      # Discrepancy measures
      sqE_adult_count_actual[i] <- pow(adult_count[i] - expected_adult_count_1[i],2)
      sqE_adult_count_sim_1[i] <- pow(sim_adult_count_1[i] - expected_adult_count_1[i],2)
      #------------------------------------
      
    }
    
    ###############################
    # Satellite observation model
    ###############################
    
    # Describes proportional bias and variance in satellite counts
    sat_CV ~ dunif(0,1)
    sat_slope ~ dnorm(1,25)
    sat_p ~ dunif(0,1)
    
    for (i in 1:n_obs_satellite){
    
      # Observation error (and bias) for satellite counts is estimated from data
      sat_mean[i] <- N[satellite_site[i],satellite_year[i]] * sat_slope 
      sat_sd[i] <- sat_mean[i] * sat_CV + 0.001
      sat_tau[i] <- pow(sat_sd[i],-2)
      
      sat_z[i] ~ dbern(sat_p)
      satellite[i] ~ dnorm(sat_mean[i] * sat_z[i],sat_tau[i])
      
      #------------------------------------
      # FOR POSTERIOR PREDICTIVE CHECK #1
      #------------------------------------
      sim_sat_z_1[i] ~ dbern(sat_p)
      sim_satellite_1[i] ~ dnorm(sat_mean[i] * sim_sat_z_1[i],sat_tau[i]) # Simulate new satellite obs (based on fitted estimates of population)
      
      expected_satellite_1[i] <- prob_occ * sat_mean[i]
      
      # Discrepancy measures
      sqE_satellite_actual[i] <- pow(satellite[i] - expected_satellite_1[i],2)
      sqE_satellite_sim_1[i] <- pow(sim_satellite_1[i] - expected_satellite_1[i],2)
      #------------------------------------
      
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
    
    # Site-level trends
    for (s in 1:n_sites){
      site_trend_1[s] <- inprod(log_X[s, 1:n_years],regression_weights[1,1:n_years])
      site_trend_2[s] <- inprod(log(N[s, 1:n_years] + 1),regression_weights[1,1:n_years])
    }
    
    ###############################
    # FOR POSTERIOR PREDICTIVE CHECK #2 - simulate completely new population dynamics
    ###############################
    
    # Simulate new occupancy status
    for (s in 1:n_sites){
      for (t in 1:n_years){
        sim_z_occ[s,t] ~ dbern(prob_occ)
      }
    }
    
    # Simulate new population dynamics
    for (s in 1:n_sites){
    
      sim_r_mean[s] ~ dnorm(r_mean_grandmean_mu,r_mean_grandmean_tau)
      sim_log_X[s,1] ~ dnorm(logX1_mean,logX1_tau)
      sim_N[s,1] <- exp(sim_log_X[s,1]) * sim_z_occ[s,1]

      for (t in 1:(n_years-1)){
      
        sim_log_X[s,t+1] ~ dnorm(sim_log_X[s,t] + sim_r_mean[s] + 0.5 * r_sd * r_sd, r_tau)
        sim_N[s,t+1] <- exp(sim_log_X[s,t+1]) * sim_z_occ[s,t+1]
        
      }
      
    } 

    # Simulate new aerial observations
    for (i in 1:n_obs_aerial){
  
      # Simulated data (based on new population dynamics)
      sim_log_lambda[i] ~ dnorm(sim_log_X[aerial_site[i],aerial_year[i]] - 0.5 * aerial_sigma * aerial_sigma, aerial_tau)
      sim_adult_count_2[i] ~ dpois(sim_z_occ[aerial_site[i],aerial_year[i]] * exp(sim_log_lambda[i]))
      
      # Expected values
      sim_expected_adult_count[i] <- prob_occ * exp(sim_log_X[aerial_site[i],aerial_year[i]] - 0.5 * aerial_sigma * aerial_sigma)
      
      # Discrepancy measures
      sqE_adult_count_sim_2[i] <- pow(sim_adult_count_2[i] - sim_expected_adult_count[i],2)
    }
    
    # Simulate satellite observations
    for (i in 1:n_obs_satellite){
    
      # Simulated data (based on new population dynamics)
      sim_sat_mean[i] <- sim_N[satellite_site[i],satellite_year[i]] * sat_slope 
      sim_sat_sd[i] <- sim_sat_mean[i] * sat_CV + 0.001
      sim_sat_tau[i] <- pow(sim_sat_sd[i],-2)
      sim_satellite[i] ~ dnorm(sim_sat_mean[i] * sim_sat_z_2[i],sim_sat_tau[i])
      
      sim_sat_z_2[i] ~ dbern(sat_p)
      sim_satellite_2[i] ~ dnorm(sim_sat_mean[i] * sim_sat_z_2[i],sat_tau[i]) # Simulate new satellite obs (based on fitted estimates of population)
      
      expected_satellite_2[i] <- prob_occ * sim_sat_mean[i]
      
      # Discrepancy measures
      sqE_satellite_sim_2[i] <- pow(sim_satellite_2[i] - expected_satellite_2[i],2)
    }
    
    # Sum squared errors of empirical data
    SSE_adult_count_actual <- sum(sqE_adult_count_actual[])
    SSE_satellite_actual <- sum(sqE_satellite_actual[])
    
    # Sum squared errors of simulated data (assuming fitted population dynamics)
    SSE_adult_count_sim_1 <- sum(sqE_adult_count_sim_1[])
    SSE_satellite_sim_1 <- sum(sqE_satellite_sim_1[])
    
    # Sum squared errors of simulated data (assuming new population dynamics)
    SSE_adult_count_sim_2 <- sum(sqE_adult_count_sim_2[])
    SSE_satellite_sim_2 <- sum(sqE_satellite_sim_2[])
    
    
    }
    ",fill = TRUE)
sink()

XX=cbind(rep(1,jags.data$n_years),1:jags.data$n_years)
jags.data$regression_weights <- matrix(c(0,1),1,2)%*%solve(t(XX)%*%XX)%*%t(XX)

# out <- jags(data=jags.data,
#             model.file="./EMPE_model_empirical.jags",
#             parameters.to.save=c(
#               "r_mean_grandmean_mu",
#               "r_mean_grandmean_sd",
#               "logX1_mean",
#               "logX1_sd",
#               
#               "r_sd",
#               "aerial_sigma",
#               "global_trend",
#               "N_global",
#               
#               "prob_occ",
#               "sat_CV",
#               "sat_slope",
#               "sat_p",
#               
#               "r_mean",
#               "log_X",
#               "site_trend_1",
#               "site_trend_2",
#               
#               "N",
#               
#               "SSE_adult_count_actual",
#               "SSE_satellite_actual",
#               
#               "SSE_adult_count_sim_1",
#               "SSE_satellite_sim_1",
#               
#               "SSE_adult_count_sim_2",
#               "SSE_satellite_sim_2"
#               
#               
#             ),
#             inits = inits,
#             n.chains=3,
#             n.thin = 20,
#             n.iter= 300000,
#             n.burnin= 100000,
#             parallel = TRUE)
# 
# results_output <- list(jags.data = jags.data,out = out)
# save(results_output,file = "./output_empirical/empirical_output.RData")

load(file = "./output_empirical/empirical_output.RData")
out <- results_output$out
jags.data <- results_output$jags.data

# Time required to fit
out$mcmc.info$elapsed.mins
mean(unlist(out$Rhat)>1.1, na.rm=TRUE) #proportion of Rhat values greater than 1.1
max(unlist(out$Rhat), na.rm=TRUE) # max Rhat
names(unlist(out$Rhat))[unlist(out$Rhat)>1.1] # Which parameters have not converged?

#----------------------------------------------------------
# Posterior predictive check
#----------------------------------------------------------
# Posterior predictive check #1
pval1_adult <- mean(out$sims.list$SSE_adult_count_actual > out$sims.list$SSE_adult_count_sim_1)
pval1_satellite <- mean(out$sims.list$SSE_satellite_actual > out$sims.list$SSE_satellite_sim_1)

# Posterior predictive check #2
pval2_adult <- mean(out$sims.list$SSE_adult_count_actual > out$sims.list$SSE_adult_count_sim_2)
pval2_satellite <- mean(out$sims.list$SSE_satellite_actual > out$sims.list$SSE_satellite_sim_2)

par(mfrow=c(2,2))
plot(out$sims.list$SSE_adult_count_actual ~ out$sims.list$SSE_adult_count_sim_1)
plot(out$sims.list$SSE_satellite_actual ~ out$sims.list$SSE_satellite_sim_1)
plot(out$sims.list$SSE_adult_count_actual ~ out$sims.list$SSE_adult_count_sim_1)
plot(out$sims.list$SSE_satellite_actual ~ out$sims.list$SSE_satellite_sim_1)
par(mfrow=c(1,1))

#----------------------------------------------------------
# Fig. 2: Plot of estimates at each colony in each year
#----------------------------------------------------------

N_estimates <- out$sims.list$N %>% 
  reshape2::melt() %>% 
  rename(mcmc_sample = Var1, site_number = Var2, year_number = Var3, N = value) %>%
  group_by(site_number, year_number) %>%
  summarize(N_q025 = quantile(N,0.025),
            N_q500 = quantile(N,0.500),
            N_q975 = quantile(N,0.975))
N_estimates$year <- years[N_estimates$year_number]
N_estimates$site_id <- sites[N_estimates$site_number]

p2 = ggplot(data = N_estimates) +
  geom_line(aes(x = year, y = N_q500), col = "dodgerblue")+
  geom_line(aes(x = year, y = N_q025), col = "dodgerblue", linetype = 3)+
  geom_line(aes(x = year, y = N_q975), col = "dodgerblue", linetype = 3)+
  
  geom_point(data = sat, aes(x = img_year, y = area_m2, shape = "Satellite count"))+
  geom_point(data = aer_adults, aes(x = year, y = adult_count, shape = "Aerial count (adult)"))+
  
  scale_shape_manual(name = 'Obs type', 
                     values =c('Satellite count'=4,'Aerial count (adult)'= 19))+
  
  ylab("Penguin abundance")+
  xlab("Year")+
  facet_grid(site_id~., scales = "free_y")+
  
  scale_x_continuous(breaks = years, minor_breaks = NULL)+
  #scale_y_continuous(trans = "log10")+
  
  theme_bw()+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

pdf(file = "./output_empirical/Fig2_colony_dynamics.pdf", width = 6, height = n_sites*0.8)
print(p2)
dev.off()

#----------------------------------------------------------
# Extract site-level estimates of trend, abundance, and change
#----------------------------------------------------------

# Trend in logX
site_trend_1 <- out$sims.list$site_trend_1 %>% 
  reshape2::melt() %>% 
  rename(mcmc_sample = Var1, site_number = Var2, trend = value) %>%
  group_by(site_number) %>%
  summarize(trend_1_q025 = quantile(trend,0.025),
            trend_1_q500 = quantile(trend,0.500),
            trend_1_q975 = quantile(trend,0.975))
site_trend_1$site_id <- sites[site_trend_1$site_number]

# Trend in log(N[t]+1)
site_trend_2 <- out$sims.list$site_trend_2 %>% 
  reshape2::melt() %>% 
  rename(mcmc_sample = Var1, site_number = Var2, trend = value) %>%
  group_by(site_number) %>%
  summarize(trend_2_q025 = quantile(trend,0.025),
            trend_2_q500 = quantile(trend,0.500),
            trend_2_q975 = quantile(trend,0.975))
site_trend_2$site_id <- sites[site_trend_2$site_number]
site_trend_estimates <- full_join(site_trend_1, site_trend_2)
site_trend_estimates$outlier <- abs(site_trend_estimates$trend_2_q500 - site_trend_estimates$trend_1_q500) > 0.05

# Relationship between two types of trend estimate (decided to use logX trend, if any)
ggplot(site_trend_estimates, aes(x = trend_1_q500, xmin = trend_1_q025, xmax = trend_1_q975, y = trend_2_q500, ymin = trend_2_q025, ymax = trend_2_q975))+
  geom_abline(slope = 1)+
  geom_errorbar(width = 0, col = "lightblue") +
  geom_errorbarh(height = 0, col = "lightblue") +
  
  geom_point(size = 2, col = "lightblue") +
  
  geom_text_repel(data = subset(site_trend_estimates, outlier == TRUE), aes(x = trend_1_q500, y = trend_2_q500, label = site_id),
                  col = "blue")+
  
  geom_errorbar(data = subset(site_trend_estimates, outlier == TRUE),width = 0, col = "blue") +
  geom_errorbarh(data = subset(site_trend_estimates, outlier == TRUE),height = 0, col = "blue") +
  geom_point(data = subset(site_trend_estimates, outlier == TRUE),size = 2, col = "blue") +
  
  theme_bw() +
  coord_cartesian(ylim=c(-0.5,0.5), xlim = c(-0.5,0.5))+
  xlab("Slope of log( X[t] )")+
  ylab("Slope of log( N[t]+1 )")

# N estimates
N_df <- out$sims.list$N %>%
  reshape2::melt() %>%
  rename(mcmc_sample = Var1, site_number = Var2, year_number = Var3, N = value) %>%
  group_by(site_number, year_number) %>%
  summarize(N_mean_q025 = quantile(N,0.025),
            N_mean_q500 = quantile(N,0.500),
            N_mean_q975 = quantile(N,0.975))

N_2009_df <- subset(N_df, year_number == 1) %>% rename(N_2009_q025 = N_mean_q025,
                                                   N_2009_q500 = N_mean_q500,
                                                   N_2009_q975 = N_mean_q975)

N_2018_df <- subset(N_df, year_number == 10) %>% rename(N_2018_q025 = N_mean_q025,
                                                      N_2018_q500 = N_mean_q500,
                                                      N_2018_q975 = N_mean_q975)

# Difference in abundance between 2018 and 2009
Ndiff_df <- (out$sims.list$N[,,10] - out$sims.list$N[,,1]) %>%
  reshape2::melt() %>%
  rename(mcmc_sample = Var1, site_number = Var2, Ndiff = value) %>%
  group_by(site_number) %>%
  summarize(Ndiff_q025 = quantile(Ndiff,0.025),
            Ndiff_q250 = quantile(Ndiff,0.25),
            Ndiff_q500 = quantile(Ndiff,0.500),
            Ndiff_q750 = quantile(Ndiff,0.75),
            Ndiff_q975 = quantile(Ndiff,0.975))

# Join estimates with colony coordinates
colony_estimates <- read.csv("../../../data/colony_coordinates.csv")%>%
  #st_as_sf(coords = c("long", "lat"),crs = 4269, agr = "constant", remove = FALSE) %>%
  full_join(site_trend_estimates) %>%
  full_join(N_2009_df[,-2]) %>%
  full_join(N_2018_df[,-2]) %>%
  full_join(Ndiff_df)

# Arrange by longitude
colony_estimates <- arrange(colony_estimates,long)
colony_estimates$site_id <- factor(colony_estimates$site_id,levels = colony_estimates$site_id)

#----------------------------------------------------------
# Map of colony trends
#----------------------------------------------------------

world <- map_data("world")
lim <- max(abs(colony_estimates$trend_1_q500),na.rm = TRUE)
trend_map <- ggplot(world, aes(x=long, y=lat, group=group)) +
  geom_polygon(fill = "gray95", col = "gray55",alpha=0.5) +
  scale_y_continuous(breaks=(-2:2) * 30, limits = c(-90,-60)) +
  scale_x_continuous(breaks=(-4:4) * 45) +
  coord_map("ortho", orientation=c(-90, 0, 0)) +
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank()) +

  geom_point(data = colony_estimates,aes(col = trend_1_q500,group=1, size = N_2009_q500))+
  #geom_label_repel(data = colony_estimates,aes(col = trend_1_q500,group=1, label = site_id))+
  scale_color_gradientn(colors = c("red","gray90","blue"),limits = c(-lim,lim), name = "Trend in log(X[t])")+
  scale_size_continuous(name = "N in 2009")
print(trend_map )

pdf(file = "./output_empirical/FigX_trend_map.pdf", width = 8, height = 8)
print(trend_map)
dev.off()

#----------------------------------------------------------
# Map of change (2018 relative to 2009)
#----------------------------------------------------------
lim <- max(abs(colony_estimates$Ndiff_q500),na.rm = TRUE)
change_map <- ggplot(world, aes(x=long, y=lat, group=group)) +
  geom_polygon(fill = "gray95", col = "gray55",alpha=0.5) +
  scale_y_continuous(breaks=(-2:2) * 30, limits = c(-90,-60)) +
  scale_x_continuous(breaks=(-4:4) * 45) +
  coord_map("ortho", orientation=c(-90, 0, 0)) +
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank()) +
  
  geom_point(data = colony_estimates,aes(col = Ndiff_q500 > 0,group=1, size = abs(Ndiff_q500)), alpha = 0.5)+
  #geom_label_repel(data = colony_estimates,aes(col = trend_1_q500,group=1, label = site_id))+
  scale_color_manual(values = c("red","blue"), labels = c("Decrease","Increase"), name = "Direction of Change",na.translate = FALSE)+
  scale_size_continuous(name = "Change in abundance")
print(change_map)

pdf(file = "./output_empirical/FigX_change_map.pdf", width = 8, height = 8)
print(change_map)
dev.off()

# ------------------------------------
# 2 dimensional plots of trend and change
# ------------------------------------

lim <- max(abs(colony_estimates[,c("trend_1_q025","trend_1_q975")]),na.rm = TRUE)
trend_2D <- ggplot(colony_estimates, aes(x=site_id, y = trend_1_q500, ymin = trend_1_q025, ymax = trend_1_q975)) +
  geom_hline(yintercept = 0, linetype = 2, col = "gray80")+
  geom_errorbar(width = 0)+
  geom_point()+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  coord_cartesian(ylim = c(-lim,lim))+
  ylab("Trend in logX")+
  xlab("Colony ID")
trend_2D

pdf(file = "./output_empirical/FigX_trend_2D.pdf", width = 8, height = 4)
print(trend_2D )
dev.off()
  
lim <- max(abs(colony_estimates[,c("Ndiff_q025","Ndiff_q975")]),na.rm = TRUE)
change_2D <- ggplot(colony_estimates, aes(x=site_id, y = Ndiff_q500, 
                                          ymin = Ndiff_q025, 
                                          ymax = Ndiff_q975)) +
  geom_hline(yintercept = 0, linetype = 2, col = "gray80")+
  geom_errorbar(width = 0)+
  geom_errorbar(aes(ymin = Ndiff_q250, ymax = Ndiff_q750),size = 1.2, width = 0)+
  geom_point(size = 2)+
  theme_bw()+
  scale_x_discrete(guide = guide_axis(angle = 45))+
  coord_cartesian(ylim = c(-lim,lim))+
  ylab("Change in abundance\n(2018 - 2009)")+
  xlab("Colony ID")

pdf(file = "./output_empirical/FigX_change_2D.pdf", width = 12, height = 4)
print(change_2D )
dev.off()

#----------------------------------------------------------
# Fig. XX: Global dynamics
#----------------------------------------------------------

N_global <- matrix(NA,nrow = out$mcmc.info$n.samples,ncol = jags.data$n_years)
colnames(N_global) <- year_vec
for (i in 1:out$mcmc.info$n.samples)N_global[i,] <- apply(out$sims.list$N[i,,],2,sum)

N_global_df <- N_global %>% 
  reshape2::melt() %>% 
  rename(mcmc_sample = Var1, Year = Var2, N = value)

N_global_summary <- N_global_df %>%
  group_by(Year) %>%
  summarize(N_q025 = quantile(N,0.025),
            N_q500 = quantile(N,0.500),
            N_q975 = quantile(N,0.975))

# Global trend (log-linear slope)
gt_sims <- out$sims.list$global_trend
gt <- quantile(gt_sims,c(0.025,0.5,0.975)) %>% round(3)

# Trend estimate #2 (using equation 4 in smith et al. 2014)
global_trend_2 <- 100*((N_global[,ncol(N_global)]/N_global[,1])^(1/(ncol(N_global)-1))-1)
gt2 <- round(quantile(global_trend_2,c(0.025,0.5,0.975)),1)

# Percent change from 2008 to 2018
pc_sims <- (N_global[,ncol(N_global)] - N_global[,1])/N_global[,1] * 100
pc <- quantile(pc_sims,c(0.025,0.5,0.975)) %>% round(1)

# Only plot 500 posterior samples
samples.to.plot <- round(seq(1,max(N_global_df$mcmc_sample),length.out = 3000))
tmp <- subset(N_global_df, mcmc_sample %in% samples.to.plot)

# Assuming EMPE have a generation time of 15 years, 
# the log-linear trend that would result in a 30% population decline over 3 generations is:
# log(0.7)/(GT*3) = log(0.7)/45
log(0.7)/45

P_endangered = mean(gt_sims <= log(0.7)/45) %>% round(3)

plot_global <- ggplot()+
  # Individual lines for each mcmc sample
  geom_line(data = tmp,aes(x = Year, y = N, col = factor(mcmc_sample)), alpha = 0.02)+
  scale_color_manual(values = rep("blue",length(unique(tmp$mcmc_sample))), guide = FALSE)+
  
  # Quantiles
  geom_line(data = N_global_summary, aes(x = Year, y = N_q025), linetype = 2)+
  geom_line(data = N_global_summary, aes(x = Year, y = N_q975), linetype = 2)+
  geom_line(data = N_global_summary, aes(x = Year, y = N_q500))+
  
  ylab("Global index of abundance")+
  xlab("Year")+
  ggtitle(paste0("2018 relative to 2009: ", pc[2],"% (95% CRI: ",pc[1],"% to ",pc[3],"%)\nLog-linear slope: ",gt[2]," (95% CRI: ",gt[1]," to ",gt[3],")\nProb. is slope consistent with 'Endangered' status: ",P_endangered))+
  coord_cartesian(ylim=c(min(N_global_summary$N_q025), max(N_global_summary$N_q975)))+
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())

pdf(file = "./output_empirical/Fig3_global_dynamics.pdf", width = 6, height = 4)
print(plot_global)
dev.off()


# #----------------------------------------------------------
# # Fig. 6: Satellite vs ground counts
# #----------------------------------------------------------
# 
# # Estimates of abundance at each colony, in each year
# all <- full_join(data,N_estimates)
# 
# #estimated relationship between ground count and satellite count
# n_vec = seq(min(all$N_q025,na.rm=TRUE),max(all$N_q975, na.rm=TRUE), length.out = 100)
# sat_pred = matrix(NA, nrow = out$mcmc.info$n.samples, ncol = length(n_vec))
# for (i in 1:ncol(sat_pred))sat_pred[,i] <- out$sims.list$sat_slope * n_vec[i]
# 
# #Fretwell slope: 0.933
# 
# # Plot available data
# # Comparison of estimated regression lines from current study and Fretwell et al. 2012
# 
# limits = range(all[,c("area_m2","N_q025","N_q975","adult_count")],na.rm = TRUE)
# limits = c(0,100000)
# 
# all_2 <- all %>%
#   group_by(site_number, year_number) %>%
#   summarize(adult_count = mean(adult_count),
#             area_m2 = mean(area_m2))
# p6 = ggplot() +
#   geom_abline(slope = 1, intercept = 0, col = "gray50")+
#   
#   geom_errorbarh(data = all, aes(xmin = N_q025, xmax = N_q975, y = area_m2), col = "dodgerblue", alpha = 0.2)+
#   geom_point(data = all, aes(x = N_q500, y = area_m2), col = "dodgerblue", alpha = 0.2)+
#   
#   geom_line(aes(x = n_vec, y = apply(sat_pred,2,function(x) quantile(x,0.025))), linetype = 2)+
#   geom_line(aes(x = n_vec, y = apply(sat_pred,2,function(x) quantile(x,0.975))), linetype = 2)+
#   
#   geom_line(aes(x = n_vec, y = apply(sat_pred,2,function(x) quantile(x,0.50))))+
#   
#   scale_y_continuous(limits = limits, name = expression(Satellite~area~m^2))+
#   scale_x_continuous(limits = limits, name = "Ground count")+
#   geom_point(data = all_2, aes(x = adult_count, y = area_m2))+
#   theme_bw()
# print(p6)
# 
# pdf("./output_empirical/FigX_sat_vs_ground_counts.pdf", width = 5, height = 4)
# print(p6)
# dev.off()
