# install/load necessary packages
my.packs <- c('jagsUI',"ggplot2",'reshape2','scales','tidyverse',
              'rgeos','raster','sp','sf','ggrepel','ggthemes','MCMCvis','DescTools')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("~/1_Work/EMPE_Global/analysis")

rm(list=ls())

#----------------------------------------------------------
# Load data and results
#----------------------------------------------------------
load("output_empirical/EMPE_data_prepared.RData") # Data
load(file = "output_empirical/EMPE_out.RData")    # Fitted model

#----------------------------------------------------------
# Output parameter estimates
#----------------------------------------------------------
parameter_estimates = out$summary[1:which(rownames(out$summary) == "sat_p"),] %>%
  as.data.frame()

write.csv(parameter_estimates, file = "output_empirical/tables/parameter_estimates.csv", row.names = TRUE)

#----------------------------------------------------------
# Functions
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


# Accepts an n_samp x n_year matrix of abundance estimates
regional_estimate_fn <- function(region_names = NA, N_samples = N_samples){
  
  tmp <- colony_attributes[,c("site_id","site_name","site_number")]
  tmp$region <- colony_attributes[,region_names]
  
  N_region <- N_samples %>% 
    left_join(tmp)%>%
    group_by(region,mcmc_sample,year) %>%
    summarize(N = sum(N))
  
  regional_abundance_summary <- N_region %>%
    group_by(region,year) %>%
    summarize(N_mean = mean(N),
              N_se = sd(N),
              N_q025 = quantile(N,0.025),
              N_q05 = quantile(N,0.05),
              N_median = median(N),
              N_q95 = quantile(N,0.95),
              N_q975 = quantile(N,0.975))
  
  n_reg <- length(unique(N_region$region))
  regional_trend_summary <- data.frame()
  regional_change_summary <- data.frame()
  regional_trend_samples <- vector(mode = "list", length = 0)
  regional_change_samples <- vector(mode = "list", length = 0)
  
  for (i in 1:n_reg){
    
    reg <- unique(N_region$region)[i]
    N_reg_matrix = N_region %>% 
      subset(region == reg) %>%
      spread(year, N) %>%
      ungroup() %>%
      dplyr::select(-region,-mcmc_sample) %>%
      as.matrix()
    
    # Estimates of change/trend
    reg_change_trend <- change_trend_fn(N_reg_matrix)
    
    regional_change_samples[[i]] <- reg_change_trend$percent_change_samples
    regional_trend_samples[[i]] <- reg_change_trend$OLS_regression_samples
    
    regional_change_summary <- rbind(regional_change_summary, 
                                     data.frame(Region = reg,
                                                Prob_Decline = reg_change_trend$prob_decline,
                                                Prob_30percent_Decline = reg_change_trend$prob_30percent_decline,
                                                Prob_50percent_Decline = reg_change_trend$prob_50percent_decline,
                                                Quantile = names(reg_change_trend$percent_change_summary),
                                                Estimate = reg_change_trend$percent_change_summary) %>%
                                       spread(Quantile, Estimate))
    
    regional_trend_summary <- rbind(regional_trend_summary, 
                                    data.frame(Region = reg,
                                               Quantile = names(reg_change_trend$OLS_regression_summary),
                                               Estimate = reg_change_trend$OLS_regression_summary) %>%
                                      spread(Quantile, Estimate))
    
    
  }
  
  # --------------------------------
  # Save summary tables
  # --------------------------------
  write.csv(regional_abundance_summary, file = paste0("output_empirical/tables/REGIONAL_abundance_",region_names,".csv"), row.names = FALSE)
  write.csv(regional_change_summary, file = paste0("output_empirical/tables/REGIONAL_change_",region_names,".csv"), row.names = FALSE)
  write.csv(regional_trend_summary, file = paste0("output_empirical/tables/REGIONAL_trend_",region_names,".csv"), row.names = FALSE)
  
  # --------------------------------
  # Generate separate plots for each region
  # --------------------------------
  region_plots <- vector(mode = "list", length = 0)
  
  for (i in 1:n_reg){
    
    reg <- unique(N_region$region)[i]
    
    reg_plot <- ggplot(subset(regional_abundance_summary, region == reg), aes(x = year, y = N_mean, ymin = N_q05, ymax = N_q95))+
      geom_ribbon(fill = "#0071fe", alpha = 0.3)+
      geom_line(col = "#0071fe")+
      ylab("Index of abundance")+
      xlab("Year")+
      ggtitle(reg)+
      theme_few()
    region_plots[[i]] <- reg_plot
    
    # --------------------------------
    # Save figures
    # --------------------------------
    tiff(filename = paste0("output_empirical/figures/REGIONAL_",region_names,"_",reg,".tif"), width = 4, height = 3, units = "in", res = 300)
    print(reg_plot)
    dev.off()
    
  }
  names(regional_change_samples) <- names(regional_trend_samples) <- names(region_plots) <- unique(N_region$region)
  
  
  
  return(list(regional_abundance_summary = regional_abundance_summary,
              regional_trend_summary = regional_trend_summary,
              regional_change_summary = regional_change_summary,
              regional_trend_samples = regional_trend_samples,
              regional_change_samples = regional_change_samples,
              region_plots = region_plots))
}


#----------------------------------------------------------
# Convert colony-level mcmc samples to dataframe
#----------------------------------------------------------

N_samples <- out$sims.list$N %>% 
  reshape2::melt() %>% 
  rename(mcmc_sample = Var1, site_number = Var2, year_number = Var3, N = value)
N_samples$year <- years[N_samples$year_number]

#----------------------------------------------------------
# Calculate change since 2009 at each colony
#----------------------------------------------------------

delta_N_samples = data.frame()
N_2009 = out$sims.list$N[,,1]

for (t in 1:jags.data$n_years){
  
  # Abundance in this year
  N_t = out$sims.list$N[,,t]
  
  # Change in abundance relative to 2009
  delta_N_t = N_t - N_2009
  
  # Long format
  delta_N_t_samples <- delta_N_t %>% 
    reshape2::melt() %>% 
    rename(mcmc_sample = Var1, site_number = Var2, delta_N = value)
  delta_N_t_samples$year <- years[t]
  delta_N_samples = rbind(delta_N_samples, delta_N_t_samples)
  
}

# Join with N_samples dataframe
N_samples = full_join(N_samples, delta_N_samples)

#----------------------------------------------------------
# Summarize dynamics at each colony
#----------------------------------------------------------

colony_summary = N_samples %>%
  group_by(year, site_number) %>%
  summarize(N_mean = mean(N),
            N_se = sd(N),
            N_q025 = quantile(N,0.025),
            N_q05 = quantile(N,0.05),
            N_median = median(N),
            N_q95 = quantile(N,0.95),
            N_q975 = quantile(N,0.975),
            
            change_since_2009_mean = mean(delta_N),
            change_since_2009_q025 = quantile(delta_N,0.025),
            change_since_2009_q05 = quantile(delta_N,0.05),
            change_since_2009_median = median(delta_N),
            change_since_2009_q95 = quantile(delta_N,0.95),
            change_since_2009_q975 = quantile(delta_N,0.975),
            
            prob_decline_since_2009 = mean(delta_N < 0)) %>%
  left_join(colony_attributes)

write.csv(colony_summary, file = "output_empirical/tables/colony_summary.csv", row.names = FALSE)

#----------------------------------------------------------
# Plot dynamics within each colony
#----------------------------------------------------------

colony_plot <- ggplot(data = colony_summary)+
  geom_ribbon(data = colony_summary, aes(x = year, y = N_mean, ymin = N_q05, ymax = N_q95),fill = "#0071fe", alpha = 0.3)+
  geom_line(data = colony_summary, aes(x = year, y = N_mean),col = "#0071fe")+
  geom_point(data = sat, aes(x = img_year, y = area_m2, shape = "Satellite count"))+
  geom_point(data = aer, aes(x = year, y = adult_count, shape = "Aerial count (adult)"))+
  
  scale_shape_manual(name = 'Obs type',
                     values =c('Satellite count'=4,'Aerial count (adult)'= 19))+
  
  ylab("Index of abundance")+
  xlab("Year")+
  facet_grid(site_id~., scales = "free_y")+
  theme_few()

pdf(file = "output_empirical/figures/colony_dynamics_fitted.pdf", width = 5, height = 50)
print(colony_plot)
dev.off()

#----------------------------------------------------------
# Plot observed aerial counts versus expected
#----------------------------------------------------------

# Match aerial counts to estimated population indices
aer_vs_expected_df = full_join(colony_summary, aer[,c("site_id","year","site_number","adult_count")]) %>%
  na.omit()

# At some sites there are many aerial observations per year.  Calculate the mean of these for plotting
aer_vs_expected_df = aer_vs_expected_df %>%
  group_by(site_id,year) # %>%
  # summarize(N_mean = mean(N_mean),
  #           adult_count = mean(adult_count))

# Percent error
sum_estimated = sum(aer_vs_expected_df$N_mean)
sum_observed = sum(aer_vs_expected_df$adult_count)
percent_error = mean(100*(sum_estimated - sum_observed)/sum_observed)

# Correlation
corr = cor.test(aer_vs_expected_df$N_mean,aer_vs_expected_df$adult_count)

lim = range(aer_vs_expected_df[,c("adult_count","N_mean")])
aerial_obs_vs_expected = ggplot(aer_vs_expected_df,aes(x = adult_count, y = N_mean)) +
  geom_abline(slope = 1, col = "gray50")+
  geom_point()+
  
  ylab("Estimated population index")+
  xlab("Observed adult count")+
  scale_y_continuous(trans = "log10", limits = lim)+
  scale_x_continuous(trans = "log10", limits = lim)+
  ggtitle(paste0("Observed adult counts vs estimated population indices\nCorrelation = ",round(corr$estimate,2)))+
  theme_few()
aerial_obs_vs_expected

#----------------------------------------------------------
# Plot predicted relationship between population index and satellite count
#----------------------------------------------------------
# Match aerial counts to estimated population indices
sat_vs_expected_df = full_join(colony_summary, sat[,c("site_id","img_year","area_m2","img_qualit")], by = c("site_id" = "site_id", "year" = "img_year")) %>%
  na.omit()

# At some sites there are many satial observations per year.  Calculate the mean of these for plotting
sat_vs_expected_df = sat_vs_expected_df %>%
  group_by(site_id,year) %>%
  summarize(N_mean = mean(N_mean),
            adult_count = mean(adult_count))

lim = c(0.1,max(sat_vs_expected_df[,c("area_m2","N_mean")]))
ggplot(sat_vs_expected_df,aes(x = area_m2, y = N_mean)) +
  geom_point()+
  geom_abline(slope = 1)+
  ylab("Estimated population index")+
  xlab("Observed satellite count")+
  scale_y_continuous(trans = "log10", limits = lim)+
  scale_x_continuous(trans = "log10", limits = lim)+
  theme_few()

#----------------------------------------------------------
# Plot magnitude of change at each colony on a map
#----------------------------------------------------------


# Change categories (>100% decrease,50-100% decrease, 0-50% decrease, 0-50% increase, 50-100% increase)
df_2009 = subset(colony_summary, year == 2009)
df_2018 = subset(colony_summary, year == 2018)
world <- map_data("world")
lim <- max(abs(df_2018$change_since_2009_mean),na.rm = TRUE)

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
  geom_point(data = df_2018,
             aes(x=lon, y=lat,group=1,
                 col = change_since_2009_mean#,size = change_since_2009_mean
                 ))+
  geom_label_repel(data = df_2018,aes(x=lon, y=lat,group=1,
                                      label = site_id,
                                      col = change_since_2009_mean
                                      ))+
  scale_color_gradientn(colors = c("red","gray90","blue"),limits = c(-lim,lim), name = "Change since 2009")+
  scale_size_continuous(name = "Change since 2009")
print(trend_map)

pdf(file = "./output_empirical/FigX_trend_map.pdf", width = 8, height = 8)
print(trend_map)
dev.off()


#----------------------------------------------------------
# Summarize regional dynamics
#----------------------------------------------------------

# Regional trends based on fast ice regions
fast_ice_reg <- regional_estimate_fn(region_names = "ice_reg", N_samples = N_samples)

# Regional trends based on pack ice regions
pack_ice_reg <- regional_estimate_fn(region_names = "p_ice_reg", N_samples = N_samples)

# Regional trends based on ccamlr regions
ccamlr_reg <- regional_estimate_fn(region_names = "ccamlr_reg", N_samples = N_samples)

#----------------------------------------------------------
# Correlation between regional sea ice trends and population trends
#----------------------------------------------------------

icetrend <- read.csv("../data/fast_ice_trends.csv")
popchange_fastice = fast_ice_reg$regional_change_summary %>% full_join(icetrend)

nameColor <- bquote(atop(Minimum~fast,
                         ice~extent~(km^2)))

sea_ice_plot <- ggplot(data = popchange_fastice,aes(x = FastIceTrend, 
                                                    y = Prob_Decline, 
                                                    col = FastIceExtent_min*1000,
                                                    label = Region)) + 
  geom_point(size = 2)+
  geom_text_repel(col = "gray70", size = 1.5, hjust = 0,direction = "x")+
  scale_color_gradientn(colors = c("gray90","blue"))+
  xlab("Fast ice trend\n(% change per year)")+
  ylab("Probability of population decline")+
  coord_cartesian(xlim = c(-3.2,3.2))+
  theme_few()+
  labs(color = nameColor)+
  theme(legend.position = "right")
print(sea_ice_plot)

# Save figure
tiff(filename = "output_empirical/figures/sea_ice_correlation.tif", width = 7, height = 3.5, units = "in", res = 300)
print(sea_ice_plot)
dev.off()

rho = DescTools::SpearmanRho(popchange_fastice$Prob_Decline,popchange_fastice$FastIceTrend,conf.level = 0.95)
rho

sea_ice_plot_black <- ggplot(data = popchange_fastice,aes(x = FastIceTrend, 
                                                          y = Prob_Decline, 
                                                          label = Region)) + 
  geom_point(size = 2)+
  geom_text_repel(col = "gray70", size = 1.5, hjust = 0,direction = "x")+
  xlab("Fast ice trend\n(% change per year)")+
  ylab("Probability of population decline")+
  coord_cartesian(xlim = c(-3.2,3.2))+
  theme_few()+
  labs(color = nameColor)+
  theme(legend.position = "right")
print(sea_ice_plot_black)

# Save figure
tiff(filename = "output_empirical/figures/sea_ice_correlation_black.tif", width = 5, height = 3.5, units = "in", res = 300)
print(sea_ice_plot_black)
dev.off()

#----------------------------------------------------------
# Summarize global dynamics
#----------------------------------------------------------

N_global <- out$sims.list$N_global 
colnames(N_global) <- year_vec
N_global <- N_global %>% 
  reshape2::melt() %>% 
  rename(mcmc_sample = Var1, Year = Var2, N = value)

N_global_matrix = N_global %>% 
  spread(Year, N) %>%
  dplyr::select(-mcmc_sample) %>%
  as.matrix()

# Estimates of abundance/change/trend
global_abundance_summary <- N_global %>%
  group_by(Year) %>%
  summarize(
    N_mean = mean(N),
    N_se = sd(N),
    N_q025 = quantile(N,0.025),
    N_q05 = quantile(N,0.05),
    N_median = quantile(N,0.5),
    N_q95 = quantile(N,0.95),
    N_q975 = quantile(N,0.975))

global_change_trend <- change_trend_fn(N_global_matrix)

global_change_summary <- data.frame(Region = "Global",
                                    Prob_Decline = global_change_trend$prob_decline,
                                    Prob_30percent_Decline = global_change_trend$prob_30percent_decline,
                                    Prob_50percent_Decline = global_change_trend$prob_50percent_decline,
                                    Quantile = names(global_change_trend$percent_change_summary),
                                    Estimate = global_change_trend$percent_change_summary) %>% 
  spread(Quantile, Estimate)

global_trend_summary <- data.frame(Region = "Global",
                                   Quantile = names(global_change_trend$OLS_regression_summary),
                                   Estimate = global_change_trend$OLS_regression_summary) %>% 
  spread(Quantile, Estimate)

write.csv(global_abundance_summary, file = "output_empirical/tables/GLOBAL_abundance.csv", row.names = FALSE)
write.csv(global_change_summary, file = "output_empirical/tables/GLOBAL_change.csv", row.names = FALSE)
write.csv(global_trend_summary, file = "output_empirical/tables/GLOBAL_trend.csv", row.names = FALSE)

global_plot <- ggplot(global_abundance_summary, aes(x = Year, y = N_mean, ymin = N_q05, ymax = N_q95))+
  geom_ribbon(fill = "#0071fe", alpha = 0.3)+
  geom_line(col = "#0071fe")+
  ylab("Index of abundance")+
  xlab("Year")+
  ggtitle("Global population")+
  scale_y_continuous(labels = comma, expand = expansion(mult = c(0.25, 0.25)))+
  theme_few()
print(global_plot)

tiff(filename = "output_empirical/figures/GLOBAL.tif", width = 4, height = 3, units = "in", res = 300)
print(global_plot)
dev.off()

#----------------------------------------------------------
# Calculate percent of population currently within an MPA
#----------------------------------------------------------

mpa_current = N_samples %>%
  left_join(colony_attributes) %>%
  group_by(year, mcmc_sample,in_MPA_current) %>%
  summarize(N = sum(N))

# Proportion of global population currently in MPA
in_current = subset(mpa_current, year == 2018 & in_MPA_current == "yes")$N
out_current = subset(mpa_current, year == 2018 & in_MPA_current == "no")$N
prop_current = in_current/(in_current+out_current)

# Proportion currently in an MPA
mean(prop_current)
quantile(prop_current, c(0.025,0.5,0.975))

# Number in an MPA
mpa_current_summary = mpa_current %>%
  group_by(year,in_MPA_current) %>%
  summarize(N_mean = mean(N),
            N_se = sd(N),
            N_q025 = quantile(N,0.025),
            N_q05 = quantile(N,0.05),
            N_median = median(N),
            N_q95 = quantile(N,0.95),
            N_q975 = quantile(N,0.975)) 

#----------------------------------------------------------
# Calculate percent of population within proposed MPAs
#----------------------------------------------------------

# Colony attributes file lists the names of the proposed MPAs.  Convert to a binary yes or no
colony_attributes$in_proposed_MPA = "yes"
colony_attributes$in_proposed_MPA[colony_attributes$MPA_proposed_name == ""] = "no"

mpa_proposed = N_samples %>%
  left_join(colony_attributes) %>%
  group_by(year, mcmc_sample,in_proposed_MPA) %>%
  summarize(N = sum(N))

# Proportion of global population in proposed MPAs (not including current MPA)
in_proposed = subset(mpa_proposed, year == 2018 & in_proposed_MPA == "yes")$N
out_proposed = subset(mpa_proposed, year == 2018 & in_proposed_MPA == "no")$N
prop_proposed = in_proposed/(in_proposed+out_proposed)

mean(prop_proposed)
quantile(prop_proposed, c(0.025,0.5,0.975))

# Number in an MPA
mpa_proposed_summary = mpa_proposed %>%
  group_by(year,in_proposed_MPA) %>%
  summarize(N_mean = mean(N),
            N_se = sd(N),
            N_q025 = quantile(N,0.025),
            N_q05 = quantile(N,0.05),
            N_median = median(N),
            N_q95 = quantile(N,0.95),
            N_q975 = quantile(N,0.975)) 

#----------------------------------------------------------
# Calculate percent of population within MPAs (current and proposed)
#----------------------------------------------------------

mpa_any = N_samples %>%
  left_join(colony_attributes) %>%
  group_by(year, mcmc_sample,in_MPA_any) %>%
  summarize(N = sum(N))

# Proportion of global population any in MPA (current and proposed)
in_any = subset(mpa_any, year == 2018 & in_MPA_any == "yes")$N
out_any = subset(mpa_any, year == 2018 & in_MPA_any == "no")$N
prop_any = in_any/(in_any+out_any)

mean(prop_any)
quantile(prop_any, c(0.025,0.5,0.975))

mpa_any_summary = mpa_any %>%
  group_by(year,in_MPA_any) %>%
  summarize(N_mean = mean(N),
            N_se = sd(N),
            N_q025 = quantile(N,0.025),
            N_q05 = quantile(N,0.05),
            N_median = median(N),
            N_q95 = quantile(N,0.95),
            N_q975 = quantile(N,0.975)) 

#----------------------------------------------------------
# Estimate change over 3 generations and IUCN thresholds
# Generation time = 16 years (Jenouvrier et al. 2014 Nature Climate Change)
# 3 generations = 48 years
#----------------------------------------------------------

annual_trend = global_change_trend$OLS_regression_samples

# Estimated percent of population remaining after 3 generations
IUCN_3generation = 100*exp(annual_trend*48)

# Probability the population will decline by more than 80% (Critically Endangered)
mean(IUCN_3generation <= 20)

# Probability the population will decline by more than 50% (Endangered)
mean(IUCN_3generation <= 50)

# Probability the population will decline by more than 30% (Vulnerable)
mean(IUCN_3generation <= 70)
