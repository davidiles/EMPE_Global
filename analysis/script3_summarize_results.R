# install/load necessary packages
my.packs <- c('jagsUI',"ggplot2",'reshape2','scales','tidyverse',
              'rgeos','raster','sp','sf','ggrepel','ggthemes')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("~/1_Work/GoogleDrive_emperor_satellite/GoogleDrive_emperorsatellite/EMPE_Global/analysis")

rm(list=ls())

#----------------------------------------------------------
# Load data and results
#----------------------------------------------------------
load("output_empirical/EMPE_data_prepared.RData") # Data
load(file = "output_empirical/EMPE_out.RData")    # Fitted model

subset(colony_attributes, !(site_id %in% sites)) # Sites not included in analysis

#----------------------------------------------------------
# Examine convergence
#----------------------------------------------------------
mean(unlist(out$Rhat)>1.1, na.rm=TRUE) #proportion of Rhat values greater than 1.1
max(unlist(out$Rhat), na.rm=TRUE) # max Rhat
names(unlist(out$Rhat))[unlist(out$Rhat)>1.1] # Which parameters have not converged?

#----------------------------------------------------------
# Posterior predictive check
#----------------------------------------------------------

pval1_adult <- mean(out$sims.list$RMSE_adult_count_actual > out$sims.list$RMSE_adult_count_sim) %>% round(2)
pval1_satellite <- mean(out$sims.list$RMSE_satellite_actual > out$sims.list$RMSE_satellite_sim) %>% round(2)
par(mfrow=c(1,2))
lim <- range(c(out$sims.list$RMSE_adult_count_actual,out$sims.list$RMSE_adult_count_sim))
plot(out$sims.list$RMSE_adult_count_actual ~ out$sims.list$RMSE_adult_count_sim, xlim = lim, ylim = lim,
     main = paste0("Bayesian p-value = ",pval1_adult),
     xlab = "RMSE adult counts (simulated)",
     ylab = "RMSE adult counts (actual)")
abline(a = 0, b = 1)

lim <- range(c(out$sims.list$RMSE_satellite_actual,out$sims.list$RMSE_satellite_sim))
plot(out$sims.list$RMSE_satellite_actual ~ out$sims.list$RMSE_satellite_sim, xlim = lim, ylim = lim,
     main = paste0("Bayesian p-value = ",pval1_satellite),
     xlab = "RMSE satellite (simulated)",
     ylab = "RMSE satellite (actual)")
abline(a = 0, b = 1)
par(mfrow=c(1,1))

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

#----------------------------------------------------------
# Observed global dynamics
#----------------------------------------------------------

N_global <- out$sims.list$N_global
colnames(N_global) <- year_vec
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
global_trend_2 <- 100*((N_global[,jags.data$n_years]/N_global[,1])^(1/(jags.data$n_years-1))-1)
gt2 <- round(quantile(global_trend_2,c(0.025,0.5,0.975)),1)

# Percent change from 2008 to 2018
pc_sims <- (N_global[,jags.data$n_years] - N_global[,1])/N_global[,1] * 100
pc <- quantile(pc_sims,c(0.025,0.5,0.975)) %>% round(1)

# Only plot a subset of posterior samples
samples.to.plot <- round(seq(1,max(N_global_df$mcmc_sample),length.out = 3000))
tmp <- subset(N_global_df, mcmc_sample %in% samples.to.plot)

# Assuming EMPE have a generation time of 12 years, 
# the log-linear trend that would result in a 30% population decline over 3 generations is:
# log(0.7)/(GT*3 - 1)
GT <- 12
P_endangered = mean(gt_sims <= log(0.5)/(GT*3-1)) %>% round(2)
P_threatened = mean(gt_sims <= log(0.7)/(GT*3-1)) %>% round(2)

plot_global_dynamics <- ggplot()+
  
  # Individual lines for each mcmc sample
  geom_line(data = tmp,aes(x = Year, y = N, col = factor(mcmc_sample)), size = 0.2, alpha = 0.015)+
  scale_color_manual(values = rep("blue",length(unique(tmp$mcmc_sample))), guide = "none")+
  
  # Quantiles
  geom_line(data = N_global_summary, aes(x = Year, y = N_q025), linetype = 2)+
  geom_line(data = N_global_summary, aes(x = Year, y = N_q975), linetype = 2)+
  geom_line(data = N_global_summary, aes(x = Year, y = N_q500))+
  
  ylab("Index of abundance")+
  xlab("Year")+
  ggtitle(paste0("2018 relative to 2009: ", pc[2],"% (95% CRI: ",pc[1],"% to ",pc[3],"%)\nLog-linear slope: ",gt[2]," (95% CRI: ",gt[1]," to ",gt[3],")\nProb. is slope consistent with 'Threatened' status: ",P_threatened,"\nProb. is slope consistent with 'Endangered' status: ",P_endangered))+
  scale_y_continuous(limits=c(0, max(N_global_df$N)), labels = comma) +
  theme_few()

pdf(file = "output_empirical/Fig3_global_dynamics.pdf", width = 6, height = 4)
print(plot_global_dynamics)
dev.off()

#----------------------------------------------------------
# Group sites together (e.g., by sea ice region and plot dynamics)
#----------------------------------------------------------

N_ice_reg_1 <- N_samples %>% 
  left_join(colony_attributes) %>%
  group_by(ice_reg,mcmc_sample,year) %>%
  summarize(N = sum(N)) %>%
  group_by(ice_reg, year) %>%
  summarize(N_q025 = quantile(N,0.025),
            N_q500 = quantile(N,0.500),
            N_q975 = quantile(N,0.975))
  

nice_palette <- c("#016B9B","#84C5E4","#FFD146","#3D9946","#ED9484","#AF5DA4","#D6D8D8","#76C044","#07703B","#CCB776")
region_colors <- nice_palette[1:length(unique(N_ice_reg_1$ice_reg))]
plot_regional_dynamics_1 <- ggplot(N_ice_reg_1, aes(x = year, ymin = N_q025, ymax = N_q975, y = N_q500, col = ice_reg, fill = ice_reg))+
  
  # Individual lines for each mcmc sample
  geom_ribbon(alpha = 0.5, col = "transparent")+
  geom_line(size=1)+
  geom_point(aes(x = 2010,y=0), col = "transparent")+

  scale_color_manual(values = region_colors, guide = "none")+
  scale_fill_manual(values = region_colors, guide = "none")+
  
  ylab("Index of abundance")+
  xlab("Year")+
  theme_few()+
  facet_wrap(ice_reg~., scales = "free_y")

world <- map_data("world")
map_1 <- ggplot(world, aes(x=long, y=lat, group=group)) +
  geom_polygon(fill = "gray95", col = "gray55",alpha=0.5) +
  scale_y_continuous(breaks=(-2:2) * 30, limits = c(-90,-60)) +
  scale_x_continuous(breaks=(-4:4) * 45) +
  coord_map("ortho", orientation=c(-90, 0, 0)) +
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank()) +
  
  geom_point(data = colony_attributes,aes(x = lon, y = lat, col = ice_reg,group=1))+
  geom_label_repel(data = colony_attributes,aes(x = lon, y = lat, col = ice_reg,group=1, label = site_id))+
  scale_color_manual(values = region_colors, name = "Sea ice region",na.translate = FALSE)

combined_plot_1 <- cowplot::plot_grid(map_1, plot_regional_dynamics_1, nrow = 2)
pdf(file = "output_empirical/Fig5_map_region_dynamics_1.pdf", width = 12, height = 12)
print(combined_plot_1)
dev.off()

#----------------------------------------------------------
# Group sites together (e.g., by sea ice region and plot dynamics)
#----------------------------------------------------------

N_ice_reg_2 <- N_samples %>% 
  left_join(colony_attributes) %>%
  group_by(p_ice_reg,mcmc_sample,year) %>%
  summarize(N = sum(N)) %>%
  group_by(p_ice_reg, year) %>%
  summarize(N_q025 = quantile(N,0.025),
            N_q500 = quantile(N,0.500),
            N_q975 = quantile(N,0.975))


nice_palette <- c("#016B9B","#84C5E4","#FFD146","#3D9946","#ED9484","#AF5DA4","#D6D8D8","#76C044","#07703B","#CCB776")
region_colors <- nice_palette[1:length(unique(N_ice_reg_2$p_ice_reg))]
plot_regional_dynamics_2 <- ggplot(N_ice_reg_2, aes(x = year, ymin = N_q025, ymax = N_q975, y = N_q500, col = p_ice_reg, fill = p_ice_reg))+
  
  # Individual lines for each mcmc sample
  geom_ribbon(alpha = 0.5, col = "transparent")+
  geom_line(size=1)+
  geom_point(aes(x = 2010,y=0), col = "transparent")+
  
  scale_color_manual(values = region_colors, guide = "none")+
  scale_fill_manual(values = region_colors, guide = "none")+
  
  ylab("Index of abundance")+
  xlab("Year")+
  scale_y_continuous(labels = comma) +
  theme_few()+
  facet_wrap(p_ice_reg~., scales = "free_y")

world <- map_data("world")
map_2 <- ggplot(world, aes(x=long, y=lat, group=group)) +
  geom_polygon(fill = "gray95", col = "gray55",alpha=0.5) +
  scale_y_continuous(breaks=(-2:2) * 30, limits = c(-90,-60)) +
  scale_x_continuous(breaks=(-4:4) * 45) +
  coord_map("ortho", orientation=c(-90, 0, 0)) +
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank()) +
  
  geom_point(data = colony_attributes,aes(x = lon, y = lat, col = p_ice_reg,group=1))+
  geom_label_repel(data = colony_attributes,aes(x = lon, y = lat, col = p_ice_reg,group=1, label = site_id))+
  scale_color_manual(values = region_colors, name = "Sea ice region",na.translate = FALSE)

combined_plot_2 <- cowplot::plot_grid(map_2, plot_regional_dynamics_2, nrow = 2)
pdf(file = "output_empirical/Fig5_map_region_dynamics_2.pdf", width = 12, height = 12)
print(combined_plot_2)
dev.off()

#----------------------------------------------------------
# Group sites together (e.g., by sea ice region and plot dynamics)
#----------------------------------------------------------

N_ice_reg_3 <- N_samples %>% 
  left_join(colony_attributes) %>%
  group_by(Genetic.name,mcmc_sample,year) %>%
  summarize(N = sum(N)) %>%
  group_by(Genetic.name, year) %>%
  summarize(N_q025 = quantile(N,0.025),
            N_q500 = quantile(N,0.500),
            N_q975 = quantile(N,0.975))


nice_palette <- c("#016B9B","#84C5E4","#FFD146","#3D9946","#ED9484","#AF5DA4","#D6D8D8","#76C044","#07703B","#CCB776")
region_colors <- nice_palette[1:length(unique(N_ice_reg_3$Genetic.name))]
plot_regional_dynamics_3 <- ggplot(N_ice_reg_3, aes(x = year, ymin = N_q025, ymax = N_q975, y = N_q500, col = Genetic.name, fill = Genetic.name))+
  
  # Individual lines for each mcmc sample
  geom_ribbon(alpha = 0.5, col = "transparent")+
  geom_line(size=1)+
  geom_point(aes(x = 2010,y=0), col = "transparent")+
  
  scale_color_manual(values = region_colors, guide = "none")+
  scale_fill_manual(values = region_colors, guide = "none")+
  
  ylab("Index of abundance")+
  xlab("Year")+
  scale_y_continuous(labels = comma) +
  theme_few()+
  facet_wrap(Genetic.name~., scales = "free_y")

world <- map_data("world")
map_3 <- ggplot(world, aes(x=long, y=lat, group=group)) +
  geom_polygon(fill = "gray95", col = "gray55",alpha=0.5) +
  scale_y_continuous(breaks=(-2:2) * 30, limits = c(-90,-60)) +
  scale_x_continuous(breaks=(-4:4) * 45) +
  coord_map("ortho", orientation=c(-90, 0, 0)) +
  theme_bw()+
  theme(axis.text = element_blank(),
        axis.title = element_blank(),
        axis.ticks=element_blank(),
        panel.border = element_blank()) +
  
  geom_point(data = colony_attributes,aes(x = lon, y = lat, col = Genetic.name,group=1))+
  geom_label_repel(data = colony_attributes,aes(x = lon, y = lat, col = Genetic.name,group=1, label = site_id))+
  scale_color_manual(values = region_colors, name = "Genetic region",na.translate = FALSE)

combined_plot_3 <- cowplot::plot_grid(map_3, plot_regional_dynamics_3, nrow = 2)
pdf(file = "output_empirical/Fig5_map_region_dynamics_3.pdf", width = 12, height = 12)
print(combined_plot_3)
dev.off()

#----------------------------------------------------------
# Group sites together by CCAMLR region
#----------------------------------------------------------

N_ccamlr_reg <- N_samples %>% 
  left_join(colony_attributes) %>%
  group_by(ccamlr_reg,mcmc_sample,year) %>%
  summarize(N = sum(N)) %>%
  group_by(ccamlr_reg, year) %>%
  summarize(N_q025 = quantile(N,0.025),
            N_q500 = quantile(N,0.500),
            N_q975 = quantile(N,0.975))


nice_palette <- c("#016B9B","#84C5E4","#FFD146","#3D9946","#ED9484","#AF5DA4","#D6D8D8","#76C044","#07703B","#CCB776")
region_colors <- nice_palette[1:length(unique(N_ccamlr_reg$ccamlr_reg))]
plot_ccamlr_dynamics <- ggplot(N_ccamlr_reg, aes(x = year, ymin = N_q025, ymax = N_q975, y = N_q500, col = ccamlr_reg, fill = ccamlr_reg))+
  
  # Individual lines for each mcmc sample
  geom_ribbon(alpha = 0.5, col = "transparent")+
  geom_line(size=1)+
  geom_point(aes(x = 2010,y=0), col = "transparent")+
  
  scale_color_manual(values = region_colors, guide = "none")+
  scale_fill_manual(values = region_colors, guide = "none")+
  
  ylab("Index of abundance")+
  xlab("Year")+
  theme_few()+
  facet_wrap(ccamlr_reg~., scales = "free_y")

pdf(file = "output_empirical/Fig4_ccamlr_dynamics.pdf", width = 8, height = 6)
print(plot_ccamlr_dynamics)
dev.off()