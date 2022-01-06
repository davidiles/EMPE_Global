# install/load necessary packages
my.packs <- c('jagsUI',"ggplot2",'reshape2','scales','tidyverse',
              'rgeos','raster','sp','sf','ggrepel','ggthemes','MCMCvis')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("E:/1_Work/GoogleDrive_emperor_satellite/GoogleDrive_emperorsatellite/EMPE_Global/analysis")

rm(list=ls())

#----------------------------------------------------------
# Load data and results
#----------------------------------------------------------
load("output_empirical/EMPE_data_prepared.RData") # Data
load(file = "output_empirical/EMPE_out.RData")    # Fitted model
subset(colony_attributes, !(site_id %in% sites)) # Sites not included in analysis

#----------------------------------------------------------
# Load counts of breeding adults and reproductive success at PGEO
#----------------------------------------------------------
pgeo_dat <- read.csv("../data/empe_PGEO_data.csv") # ground data at PGEO

#----------------------------------------------------------
# Convert mcmc samples to dataframe
#----------------------------------------------------------

N_samples <- out$sims.list$N %>% 
  reshape2::melt() %>% 
  rename(mcmc_sample = Var1, site_number = Var2, year_number = Var3, N = value)
N_samples$year <- years[N_samples$year_number]

#----------------------------------------------------------
# Colony-level estimates
#----------------------------------------------------------

colony_estimates <- N_samples %>%
  group_by(site_number, year) %>%
  summarize(N_q025 = quantile(N,0.025),
            N_q500 = quantile(N,0.500),
            N_q975 = quantile(N,0.975)) %>%
  left_join(., colony_attributes)
colony_estimates$site_id <- fct_reorder(colony_estimates$site_id, colony_estimates$lon)

#----------------------------------------------------------
# Restrict to PGEO
#----------------------------------------------------------

# model estimates restricted to PGEO
pgeo_est <- subset(colony_estimates, site_id == "PGEO")

# breeding pair counts at PGEO (and repro success)
pgeo_dat <- subset(pgeo_dat, Year >= 2009) %>%
  rename(Breeding_Success = 4, Breeding_Pairs = 5) %>%
  mutate(BPxBS = Breeding_Success * Breeding_Pairs)

# summarize weekly ground counts (from Celine) by mean and 1.96*SE
pgeo_weekly <- subset(aer_adults, site_id == "PGEO") %>%
  group_by(year) %>%
  summarize(mean_count = mean(adult_count),
            sd_count = sd(adult_count),
            n_count = n(),
            se_count = sd(adult_count)/sqrt(n()),
            lcl_count = mean(adult_count) - 1.96*sd(adult_count)/sqrt(n()),
            ucl_count = mean(adult_count) + 1.96*sd(adult_count)/sqrt(n()))

p1 = ggplot(data = pgeo_est) +
  
  # Estimates from Bayesian model
  geom_line(aes(x = year, y = N_q500), col = "dodgerblue")+
  geom_ribbon(aes(x = year, ymin = N_q025, ymax = N_q975), fill = "dodgerblue", alpha = 0.2)+

  # Satellite counts
  geom_point(data = subset(sat, site_id == "PGEO"), aes(x = img_year, y = area_m2, shape = "Satellite count"))+
  
  # Mean of weekly ground counts from Celine
  geom_point(data =  pgeo_weekly, aes(x = year, y = mean_count, shape = "Mean Sep-Nov ground count (Celine)"))+
  
  # Breeding Pair counts
  #geom_point(data =  pgeo_dat, aes(x = Year, y = Breeding_Pairs, shape = "Breeding Pairs"))+
  
  # Breeding Pairs x RS
  #geom_point(data =  pgeo_dat, aes(x = Year, y = BPxBS, shape = "Breeding Pairs x Breeding Success"))+
  
  scale_shape_manual(name = '', values =c('Mean Sep-Nov ground count (Celine)'= 19,'Satellite count'=4))+
  
  ylab("Count")+
  xlab("Year")+
  scale_x_continuous(breaks = years, minor_breaks = NULL)+
  theme_few()+
  theme(legend.position="top",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))
p1

p2 = ggplot(data = pgeo_est) +
  
  # Estimates from Bayesian model
  geom_line(aes(x = year, y = N_q500), col = "dodgerblue")+
  geom_ribbon(aes(x = year, ymin = N_q025, ymax = N_q975), fill = "dodgerblue", alpha = 0.2)+
  
  # Satellite counts
  #geom_point(data = subset(sat, site_id == "PGEO"), aes(x = img_year, y = area_m2, shape = "Satellite count"))+
  
  # Mean of weekly ground counts from Celine
  #geom_point(data =  pgeo_weekly, aes(x = year, y = mean_count, shape = "Mean Sep-Nov ground count (Celine)"))+
  
  # Breeding Pair counts
  geom_point(data =  pgeo_dat, aes(x = Year, y = Breeding_Pairs, shape = "Breeding Pairs"))+
  geom_line(data =  pgeo_dat, aes(x = Year, y = Breeding_Pairs, shape = "Breeding Pairs"))+
  
  # Breeding Pairs x RS
  geom_point(data =  pgeo_dat, aes(x = Year, y = BPxBS, shape = "Breeding Pairs x Breeding Success"))+
  
  scale_shape_manual(name = '', values =c('Breeding Pairs' = 11))+
  
  ylab("Count")+
  xlab("Year")+
  scale_x_continuous(breaks = years, minor_breaks = NULL)+
  theme_few()+
  theme(legend.position="top",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

p2

p3 = ggplot(data = pgeo_est) +
  
  # Estimates from Bayesian model
  geom_line(aes(x = year, y = N_q500), col = "dodgerblue")+
  geom_ribbon(aes(x = year, ymin = N_q025, ymax = N_q975), fill = "dodgerblue", alpha = 0.2)+
  
  # Satellite counts
  #geom_point(data = subset(sat, site_id == "PGEO"), aes(x = img_year, y = area_m2, shape = "Satellite count"))+
  
  # Mean of weekly ground counts from Celine
  #geom_point(data =  pgeo_weekly, aes(x = year, y = mean_count, shape = "Mean Sep-Nov ground count (Celine)"))+
  
  # Breeding Pair counts
  #geom_point(data =  pgeo_dat, aes(x = Year, y = Breeding_Pairs, shape = "Breeding Pairs"))+
  #geom_line(data =  pgeo_dat, aes(x = Year, y = Breeding_Pairs, shape = "Breeding Pairs"))+
  
  # Breeding Pairs x RS
  geom_point(data =  pgeo_dat, aes(x = Year, y = BPxBS, shape = "Breeding Pairs x Breeding Success"))+
  geom_line(data =  pgeo_dat, aes(x = Year, y = BPxBS, shape = "Breeding Pairs x Breeding Success"))+
  
  scale_shape_manual(name = '', values =c('Breeding Pairs x Breeding Success' = 14))+
  
  ylab("Count")+
  xlab("Year")+
  scale_x_continuous(breaks = years, minor_breaks = NULL)+
  theme_few()+
  theme(legend.position="top",
        legend.background = element_blank(),
        legend.box.background = element_rect(colour = "black"))

p3

#----------------------------------------------------------
# Calculate trend estimate
#----------------------------------------------------------

# Which colony is PGEO?
pgeo_site <- subset(colony_attributes, site_id == "PGEO")$site_number

pgeo_trend_est <- 100*((out$sims.list$N[,pgeo_site,jags.data$n_years-1]/out$sims.list$N[,pgeo_site,1])^(1/(jags.data$n_years-2))-1)
pgeo_trend_true <- 100*((pgeo_dat$BPxBS[nrow(pgeo_dat)]/pgeo_dat$BPxBS[1])^(1/(nrow(pgeo_dat)-1))-1)

quantile(pgeo_trend_est, c(0.025,0.5,0.975))

p4 = ggplot() +
  geom_vline(xintercept = 0, linetype = 2)+
  geom_histogram(aes(x = pgeo_trend_est), fill = "dodgerblue", alpha = 0.5)+
  
  geom_vline(xintercept = median(pgeo_trend_est), col = "dodgerblue", size = 2)+
  geom_vline(xintercept = median(pgeo_trend_est), col = "dodgerblue", size = 2)+
  
  xlab("PGEO trend estimate")+
  ylab("freq")+
  theme_few()

p4

mean(pgeo_trend_est < pgeo_trend_true)
