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

p2 = ggplot(data = colony_estimates) +
  geom_line(aes(x = year, y = N_q500), col = "dodgerblue")+
  geom_line(aes(x = year, y = N_q025), col = "dodgerblue", linetype = 3)+
  geom_line(aes(x = year, y = N_q975), col = "dodgerblue", linetype = 3)+
  
  geom_point(data = sat, aes(x = img_year, y = area_m2, shape = "Satellite count"))+
  geom_point(data = aer_adults, aes(x = year, y = adult_count, shape = "Aerial count (adult)"))+
  geom_point(aes(x = 2010,y=0), col = "transparent")+
  
  scale_shape_manual(name = 'Obs type', values =c('Satellite count'=4,'Aerial count (adult)'= 19))+
  
  ylab("Count")+
  xlab("Year")+
  facet_grid(site_id~., scales = "free_y")+
  scale_x_continuous(breaks = years, minor_breaks = NULL)+
  theme(axis.text.x = element_text(angle = 45, hjust = 1))+
  theme_few()

pdf(file = "output_empirical/Fig2_colony_dynamics.pdf", width = 6, height = n_sites*0.8)
print(p2)
dev.off()

# #----------------------------------------------------------
# # Calculate trend (log-linear slope)
# #----------------------------------------------------------
# 
# colony_trends <- out$sims.list$site_trend %>% 
#   reshape2::melt() %>% 
#   rename(mcmc_sample = Var1, site_number = Var2, trend = value) %>%
#   group_by(site_number) %>%
#   summarize(trend_q025 = quantile(trend,0.025),
#             trend_q500 = quantile(trend,0.500),
#             trend_q975 = quantile(trend,0.975)) %>%
#   left_join(colony_attributes) %>%
#   arrange(lon)
# colony_trends$site_id <- fct_reorder(colony_trends$site_id, colony_trends$lon)
# 
# lim <- max(abs(colony_trends[,c("trend_q025","trend_q975")]),na.rm = TRUE)
# trend_2D <- ggplot(colony_trends, aes(x=site_id, y = trend_q500, ymin = trend_q025, ymax = trend_q975)) +
#   geom_hline(yintercept = 0, linetype = 2, col = "gray80")+
#   geom_errorbar(width = 0)+
#   geom_point()+
#   theme_bw()+
#   scale_x_discrete(guide = guide_axis(angle = 45))+
#   coord_cartesian(ylim = c(-lim,lim))+
#   ylab("Log-linear trend")+
#   xlab("Site ID")
# trend_2D
# 
# pdf(file = "output_empirical/FigX_trend_2D.pdf", width = 8, height = 4)
# print(trend_2D )
# dev.off()

