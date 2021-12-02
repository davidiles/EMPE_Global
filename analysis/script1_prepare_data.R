# install/load necessary packages
my.packs <- c('jagsUI',"ggplot2",'reshape2','scales','tidyverse',
              'rgeos','raster','sp','sf','ggrepel')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("~/1_Work/GoogleDrive_emperor_satellite/GoogleDrive_emperorsatellite/Global_Analysis/analysis/version_3/")

rm(list=ls())

year_range <- 2009:2018

# Read in EMPE satellite data
sat <- read.csv("../../../data/satellite_data_2021_06_13.csv") # Previously satellite.csv
sat$site_id <- as.character(sat$site_id)
sat$area_m2 <- as.numeric(sat$area_m2)
sat <- subset(sat, !is.na(area_m2))
sat <- sat[,c("site_id","img_year","area_m2","bpresent")]
sat$bpresent[which(sat$bpresent == "NA")] <- NA
sat <- subset(sat,img_year %in% year_range)

# Read in EMPE aerial data
aer <- read.csv("../../../data/aerial_addedMAL_20210422_update.csv")
aer_adults <- subset(aer, count_type == "adults")
aer_adults$site_id <- as.character(aer_adults$site_id)
colnames(aer_adults)[which(colnames(aer_adults) == "penguin_count")] <- "adult_count"
colnames(aer_adults)[which(colnames(aer_adults) == "accuracy")] <- "adult_accuracy"
aer_adults <- aer_adults[,c("site_id","year","adult_count","adult_accuracy","count_type")]

aer_adults <- subset(aer_adults, year %in% year_range)

# ----------------------------------
# Sites to exclude (temporary)

n_counts_satellite <- sat %>%
  group_by(site_id,img_year) %>%
  summarize(n_nonzero = sum(area_m2>0)) %>%
  group_by(site_id) %>%
  summarize(enough_years= sum(n_nonzero)>=3)

sites_to_exclude <- NULL
sites_to_exclude <- subset(n_counts_satellite,!enough_years)$site_id
sites_to_exclude <- c(sites_to_exclude, "NWIS", "BARB","SANA_N","SANA_S")

aer_adults <- subset(aer_adults, !(site_id %in% sites_to_exclude))
sat <- subset(sat, !(site_id %in% sites_to_exclude))
# ----------------------------------

sites <- c(aer_adults$site_id,sat$site_id) %>% unique() %>% sort()
n_sites <- length(sites)

years <- c(aer_adults$year,sat$img_year) %>% unique() %>% sort()
year_vec <- years[1]:years[length(years)]
n_years <- length(year_vec)

aer_adults$site_number <- match(aer_adults$site_id,sites)
aer_adults$year_number <- match(aer_adults$year,years)

sat$site_number <- match(sat$site_id,sites)
sat$year_number <- match(sat$img_year,years)

#From MAPPPD - convert accuracy to lognormal precision
#accuracy_to_precision <- data.frame(accuracy = c(1,2,3,4,5), precision = c(1612.798190,407.160648,69.313726,20.419275,4.998676))

# Plot available data
p1 <- ggplot() +
  
  # Was colony present/absent?
  #geom_vline(data = sat,aes(xintercept = img_year, col = factor(bpresent)), size = 2, alpha = 0.3)+
  #scale_color_manual(values = c("orangered","dodgerblue"),na.translate = FALSE,name = "Colony Presence")+
  geom_point(data = sat, aes(x = img_year, y = area_m2, shape = "Satellite count"))+
  geom_point(data = aer_adults, aes(x = year, y = adult_count, shape = "Aerial count (adult)"))+
  
  scale_shape_manual(name = 'Obs type',
                     values =c('Satellite count'=4,'Aerial count (adult)'= 19))+
  
  #scale_x_continuous(breaks = year_vec, minor_breaks = NULL)+
  
  ylab("Penguin abundance")+
  xlab("Year")+
  facet_grid(site_id~., scales = "free_y")+
  theme_bw()

pdf("./output_empirical/Fig1_data.pdf", width = 6, height = n_sites*0.8)
print(p1)
dev.off()

#---------------------------------------------------------------------
# Specify colony occupancy state when this is known 
#---------------------------------------------------------------------

# data$z_known <- factor(data$bpresent,levels = c("no","yes")) %>% as.numeric()
# data$z_known <- data$z_known - 1
# 
# z_known <- matrix(NA,ncol = max(data$year_number), nrow = max(data$site_number))
# for (i in 1:nrow(data)){
#   
#   col <- data$year_number[i]
#   row <- data$site_number[i]
#   
#   if (!is.na(data$z_known[i])){
#     if (z_known[row,col] != data$z_known[i] & !is.na(z_known[row,col])){
#       print(paste0("Discrepancy in occupancy on row ",i))
#     }
#     z_known[row,col] <- data$z_known[i]
#   }
# }
# 
# # Determine z indices that need to be estimated
# z_ind <- z_known %>% 
#   reshape2::melt() %>%
#   rename(site_number = Var1, year_number = Var2, z = value) %>%
#   subset(is.na(z))
# 
# # Initial values for JAGS
# z_init <- z_known
# z_init[!is.na(z_init)] <- 3
# z_init[is.na(z_init)] <- 1
# z_init[z_init == 3] <- NA

#---------------------------------------------------------------------

#Package data for JAGS
jags.data <- list( n_years = n_years,
                   n_sites = n_sites,
                   
                   # aerial counts of adults
                   n_obs_aerial = length(aer_adults$adult_count),
                   adult_count = aer_adults$adult_count,
                   aerial_site = aer_adults$site_number,
                   aerial_year = aer_adults$year_number,
                   
                   # satellite counts
                   n_obs_satellite = length(sat$area_m2),
                   satellite = sat$area_m2,
                   satellite_site = sat$site_number,
                   satellite_year = sat$year_number
)




#Generate initial values


inits <- function(){list(r_mean = rep(0,jags.data$n_sites),
                         r_sd = runif(1,0,1),
                         z_occ = matrix(1,ncol = n_years, nrow = n_sites),
                         sat_z = rep(1,jags.data$n_obs_satellite),
                         sat_slope = 1,
                         sat_CV = 0.1)}


save.image("EMPE_data_prepared.RData")
