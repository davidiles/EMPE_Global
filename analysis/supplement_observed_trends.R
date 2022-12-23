# install/load necessary packages
my.packs <- c('tidyverse')

if (any(!my.packs %in% installed.packages()[, 'Package']))install.packages(my.packs[which(!my.packs %in% installed.packages()[, 'Package'])],dependencies = TRUE)
lapply(my.packs, require, character.only = TRUE)

setwd("~/1_Work/EMPE_Global/analysis")

rm(list=ls())

load("output_empirical/EMPE_data_prepared.RData")

min_year = min(sat$img_year)
max_year = max(sat$img_year)
n_years = length(min_year:max_year)

# Number of years of data at each colony
colony_nyear = sat %>%
  group_by(site_id) %>%
  summarize(nyear = length(unique(img_year)),
            first_year = min(img_year),
            last_year = max(img_year))

# Find colonies with data in first and final year
colonies_to_include = subset(colony_nyear, first_year == 2009 & last_year == 2018 & nyear >= 9)
sat_2 = sat %>%
  subset(site_id %in% colonies_to_include$site_id) %>%
  group_by(site_id, img_year) %>%
  summarize(index = mean(area_m2))

ggplot(data = sat_2, aes(x = img_year, y = index))+
  geom_point()+
  facet_wrap(site_id~.)+
  theme_bw()

ggplot(data = sat_2, aes(x = img_year, y = index))+
  geom_point()+
  facet_wrap(site_id~., scales = "free_y")+
  theme_bw()

# -------------------------------
# Sum up abundances each year
# -------------------------------
dat = sat_2 %>%
  group_by(img_year) %>%
  summarize(ncolony = length(unique(site_id)),
            n = n(),
            sum_n = sum(index),
            mean_per_colony = mean(index),
            se = sd(index)/sqrt(n()))
dat$ucl = dat$mean_per_colony + 1.96*dat$se
dat$lcl = dat$mean_per_colony - 1.96*dat$se

plot(mean_per_colony~img_year, data = dat, ylab = "mean N per colony", xlab = "Year", pch = 19, type = "o")

ggplot(dat)+
  geom_errorbar(aes(x = img_year, ymin = lcl, ymax = ucl),width=0)+
  geom_point(aes(x = img_year, y = mean_per_colony))+
  theme_bw()+
  xlab("Year")+
  ylab("Mean N per colony")

ggplot(dat)+
  geom_errorbar(aes(x = img_year, ymin = lcl, ymax = ucl),width=0)+
  geom_point(aes(x = img_year, y = mean_per_colony))+
  theme_bw()+
  xlab("Year")+
  ylab("Mean N per colony")


# -------------------------------
# Change relative to 2009
# -------------------------------
colonies_to_include = subset(colony_nyear, first_year == 2009)
colonies_to_include = subset(colony_nyear, first_year == 2009 & last_year == 2018 & nyear >= 9)

sat_2 = sat %>%
  subset(site_id %in% colonies_to_include$site_id) %>%
  group_by(site_id, img_year) %>%
  summarize(index = mean(area_m2))

N_2009 = sat_2 %>% 
  subset(img_year == 2009) %>%
  rename(index_2009 = index)

sat_2 <- full_join(sat_2,N_2009[,c("site_id","index_2009")])
sat_2$delta = sat_2$index - sat_2$index_2009

ggplot(sat_2)+
  geom_point(aes(x = img_year, y = delta))+
  theme_bw()+
  xlab("Year")+
  ylab("Change since 2009")+
  facet_wrap(site_id~., scales = "free_y")

change_df = sat_2 %>%
  group_by(img_year) %>%
  summarize(ncolony = length(unique(site_id)),
            n = n(),
            mean_deltaN_per_colony = mean(delta, na.rm = TRUE),
            se = sd(delta, na.rm = TRUE)/sqrt(n()))
change_df$ucl = change_df$mean_deltaN_per_colony + 1.96*change_df$se
change_df$lcl = change_df$mean_deltaN_per_colony - 1.96*change_df$se

ggplot(change_df)+
  geom_errorbar(aes(x = img_year, ymin = lcl, ymax = ucl),width=0)+
  geom_point(aes(x = img_year, y = mean_deltaN_per_colony))+
  theme_bw()+
  xlab("Year")+
  ylab("Mean change since 2009 per colony")+
  coord_cartesian(ylim=c(-4000,4000))

# -------------------------------
# GLMM
# -------------------------------

library(lme4)

sat$y <- log(sat$area_m2+1)
sat$year_adj = sat$img_year - min(sat$img_year)+1
m1 <- lmer(y ~ year_adj + site_id + (year_adj-1|site_id), data = sat)
