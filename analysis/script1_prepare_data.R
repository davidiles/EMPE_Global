# install/load necessary packages
my.packs <- c('jagsUI',"ggplot2",'reshape2','scales','tidyverse',
              'rgeos','raster','sp','sf','ggrepel','ggthemes')
if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("~/1_Work/EMPE_Global/analysis")

rm(list=ls())

year_range <- 2009:2018

# --------------------------------------------------
# Colony attributes
# --------------------------------------------------

colony_attributes <- read.csv("../data/colony_attributes.csv")

# --------------------------------------------------
# Read in satellite data
# --------------------------------------------------

sat <- read.csv("../data/empe_satellite_2022-01-31.csv") %>%
  mutate(area_m2 = as.numeric(area_m2)) %>%
  subset(!is.na(area_m2) & img_year %in% year_range & site_id %in% colony_attributes$site_id)

# --------------------------------------------------
# Read in aerial data
# --------------------------------------------------

aer <- read.csv("../data/empe_aerial_2022-01-31.csv") %>%
  subset(count_type == "adults") %>%
  rename(adult_count = penguin_count, adult_accuracy = accuracy) %>%
  dplyr::select(site_id,year,adult_count,adult_accuracy,count_type) %>%
  subset(!is.na(adult_count) & year %in% year_range & site_id %in% colony_attributes$site_id)

# --------------------------------------------------
# Format for plotting / summarizing
# --------------------------------------------------

sites <- c(aer$site_id,sat$site_id) %>% unique() %>% sort()
n_sites <- length(sites)

years <- c(aer$year,sat$img_year) %>% unique() %>% sort()
year_vec <- years[1]:years[length(years)]
n_years <- length(year_vec)

aer$site_number <- match(aer$site_id,sites)
aer$year_number <- match(aer$year,years)

sat$site_number <- match(sat$site_id,sites)
sat$year_number <- match(sat$img_year,years)

# Colony number identifiers
colony_attributes <- colony_attributes %>% left_join(data.frame(site_number = 1:length(sites), site_id = sites))

sat <- left_join(sat, colony_attributes) %>% mutate(site_id = fct_reorder(site_id, lon)) 
aer <- left_join(aer, colony_attributes) %>% mutate(site_id = fct_reorder(site_id, lon)) 

# Plot available data
p1 <- ggplot() +
  
  # Was colony present/absent?
  #geom_vline(data = sat,aes(xintercept = img_year, col = factor(bpresent)), size = 2, alpha = 0.3)+
  #scale_color_manual(values = c("orangered","dodgerblue"),na.translate = FALSE,name = "Colony Presence")+
  geom_point(data = sat, aes(x = img_year, y = area_m2, shape = "Satellite count"))+
  geom_point(data = aer, aes(x = year, y = adult_count, shape = "Aerial count (adult)"))+
  
  scale_shape_manual(name = 'Obs type',
                     values =c('Satellite count'=4,'Aerial count (adult)'= 19))+
  
  #scale_x_continuous(breaks = year_vec, minor_breaks = NULL)+
  
  ylab("Count")+
  xlab("Year")+
  facet_grid(site_id~., scales = "free_y")+
  theme_bw()+
  theme_few()

pdf("output_empirical/figures/raw_data.pdf", width = 6, height = n_sites*0.8)
print(p1)
dev.off()

#Package data for JAGS
jags.data <- list( n_years = n_years,
                   n_sites = n_sites,
                   
                   # aerial counts of adults
                   n_obs_aerial = length(aer$adult_count),
                   adult_count = aer$adult_count,
                   aerial_site = aer$site_number,
                   aerial_year = aer$year_number,
                   
                   # satellite counts
                   n_obs_satellite = length(sat$area_m2),
                   img_qual = sat$img_qualit,
                   satellite = sat$area_m2,
                   satellite_site = sat$site_number,
                   satellite_year = sat$year_number
)

# For calculating log-linear trend (least squares regression line)
XX=cbind(rep(1,jags.data$n_years),1:jags.data$n_years)
jags.data$regression_weights <- matrix(c(0,1),1,2)%*%solve(t(XX)%*%XX)%*%t(XX)

#Generate initial values
inits <- function(){list(r_mean = rep(0,jags.data$n_sites),
                         r_sd = runif(1,0,1),
                         z_occ = matrix(1,ncol = n_years, nrow = n_sites),
                         sat_z = rep(1,jags.data$n_obs_satellite),
                         sat_slope = c(1,1,1),
                         sat_CV = c(0.1,0.1,0.1))}

save.image("output_empirical/EMPE_data_prepared.RData")
