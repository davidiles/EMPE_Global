# install/load necessary packages
my.packs <- c('jagsUI',"ggplot2",'reshape2','scales','tidyverse',
              'rgeos','raster','sp','sf','ggrepel','ggthemes','MCMCvis','DescTools')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)


rm(list=ls())

#----------------------------------------------------------
# Load Rose's latest results (for Ross Sea)
#----------------------------------------------------------

setwd("C:/Users/david/Documents/1_Work/EMPE_Global/from_collaborators")
colony_summary_ROSE <- read.csv("colony_estimates_summary_ROSE.csv")
Ross_colonies <- unique(colony_summary_ROSE$site_id)

region_summary_ROSE <- read.csv("RossSea_abundance_ROSE.csv")

#----------------------------------------------------------
# Load/format 'original' results
#----------------------------------------------------------

setwd("~/1_Work/EMPE_Global/analysis")
load("output_empirical/EMPE_data_prepared.RData") # Data
load(file = "output_empirical/EMPE_out.RData")    # Fitted model

colony_attributes <- subset(colony_attributes, site_id %in% Ross_colonies)

# Convert colony-level mcmc samples to dataframe
N_samples <- out$sims.list$N[,,] %>% 
  reshape2::melt() %>% 
  rename(mcmc_sample = Var1, site_number = Var2, year_number = Var3, N = value) %>%
  subset(site_number %in% colony_attributes$site_number)
N_samples$year <- years[N_samples$year_number]

# Summarize dynamics at each colony
colony_summary = N_samples %>%
  group_by(year, site_number) %>%
  summarize(N_mean = mean(N),
            N_se = sd(N),
            N_q025 = quantile(N,0.025),
            N_q05 = quantile(N,0.05),
            N_median = median(N),
            N_q95 = quantile(N,0.95),
            N_q975 = quantile(N,0.975)) %>%
  left_join(colony_attributes)

N_region <- N_samples %>% 
  group_by(mcmc_sample,year) %>%
  summarize(N = sum(N))

region_summary <- N_region %>%
  group_by(year) %>%
  summarize(N_mean = mean(N),
            N_se = sd(N),
            N_q025 = quantile(N,0.025),
            N_q05 = quantile(N,0.05),
            N_median = median(N),
            N_q95 = quantile(N,0.95),
            N_q975 = quantile(N,0.975))

# Plot dynamics within each colony
colony_plot_comparison <- ggplot(data = colony_summary_ROSE)+
  
  geom_ribbon(data = colony_summary_ROSE, aes(x = year, y = N_mean, ymin = N_q05, ymax = N_q95, fill = "New"), alpha = 0.3)+
  geom_line(data = colony_summary_ROSE, aes(x = year, y = N_mean, col = "New"))+
  
  geom_ribbon(data = colony_summary, aes(x = year, y = N_mean, ymin = N_q05, ymax = N_q95, fill = "Original"), alpha = 0.3)+
  geom_line(data = colony_summary, aes(x = year, y = N_mean, col = "Original"))+
  
  ylab("Index of abundance")+
  xlab("Year")+
  facet_wrap(site_id~., scales = "free_y")+
  scale_color_manual(values = c("dodgerblue","orangered"), name = "Estimate")+
  scale_fill_manual(values = c("dodgerblue","orangered"), name = "Estimate")+
  theme_few()

print(colony_plot_comparison)


# Plot dynamics in Ross Sea
region_plot_comparison <- ggplot()+
  
  geom_ribbon(data = region_summary_ROSE, aes(x = Year, y = N_mean, ymin = N_q05, ymax = N_q95, fill = "New"), alpha = 0.3)+
  geom_line(data = region_summary_ROSE, aes(x = Year, y = N_mean, col = "New"))+
  
  geom_ribbon(data = region_summary, aes(x = year, y = N_mean, ymin = N_q05, ymax = N_q95, fill = "Original"), alpha = 0.3)+
  geom_line(data = region_summary, aes(x = year, y = N_mean, col = "Original"))+
  
  ylab("Index of abundance")+
  xlab("Year")+
  scale_color_manual(values = c("dodgerblue","orangered"), name = "Estimate")+
  scale_fill_manual(values = c("dodgerblue","orangered"), name = "Estimate")+
  theme_few()

print(region_plot_comparison)

#----------------------------------------------------------
# Function to calculate regional trend
#----------------------------------------------------------

# Accepts an n_samp x n_year matrix of abundance estimates
change_trend_fn <- function(mat){
  n_samps = nrow(mat)
  n_years = ncol(mat)
  
  # ---------------
  # Percent change between endpoints
  # ---------------
  percent_change_samples <- 100*(mat[,n_years] - mat[,1])/mat[,1]
  percent_change_summary <- c(mean = mean(percent_change_samples),SE = sd(percent_change_samples), quantile(percent_change_samples,c(0.025,0.05,0.5,0.95,0.975)))
  
  prob_decline <- mean(percent_change_samples < 0)
  prob_30percent_decline <- mean(percent_change_samples < -30)
  prob_50percent_decline <- mean(percent_change_samples < -50)
  
  # ---------------
  # Rate of log-linear change (calculated using OLS)
  # ---------------
  XX=cbind(rep(1,n_years),1:n_years)
  regression_weights <- matrix(c(0,1),1,2)%*%solve(t(XX)%*%XX)%*%t(XX)
  OLS_regression_samples <- rep(NA,n_samps)
  for (i in 1:n_samps) OLS_regression_samples[i] <- regression_weights %*% log(mat[i,]) # OLS slope, on log scale
  
  OLS_regression_summary <- c(mean = mean(OLS_regression_samples),SE = sd(OLS_regression_samples),quantile(OLS_regression_samples,c(0.025,0.05,0.5,0.95,0.975)))
  
  # Convert to percent change per year
  OLS_regression_summary <- 100*(exp(OLS_regression_summary)-1)
  
  return(list(prob_decline = prob_decline, 
              prob_30percent_decline = prob_30percent_decline,
              prob_50percent_decline = prob_50percent_decline,
              percent_change_samples = percent_change_samples,
              OLS_regression_samples = OLS_regression_samples,
              
              percent_change_summary = percent_change_summary,
              OLS_regression_summary = OLS_regression_summary))
}

Ross_matrix <- out$sims.list$N[,colony_attributes$site_number,] %>% apply(.,c(1,3),sum)
Ross_trend <- change_trend_fn(Ross_matrix)

Ross_trend

mean(Ross_matrix[,10] < Ross_matrix[,1])
