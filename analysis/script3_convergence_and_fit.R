# install/load necessary packages
my.packs <- c('jagsUI',"ggplot2",'reshape2','scales','tidyverse',
              'rgeos','raster','sp','sf','ggrepel','ggthemes','MCMCvis')
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
# Examine convergence
#----------------------------------------------------------
mean(unlist(out$Rhat)>1.1, na.rm=TRUE) #proportion of Rhat values greater than 1.1
max(unlist(out$Rhat), na.rm=TRUE) # max Rhat
names(unlist(out$Rhat))[unlist(out$Rhat)>1.1] # Which parameters have not converged?
MCMCtrace(out)

#----------------------------------------------------------
# Posterior predictive check
#----------------------------------------------------------

# Save figure
tiff(filename = "output_empirical/figures/Goodness-of-fit-Bayesianpval.tif", width = 7, height = 3.5, units = "in", res = 300)

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
dev.off()

#----------------------------------------------------------
# Correlation between estimates
#----------------------------------------------------------

plot(out$sims.list$prob_occ, out$sims.list$sat_p)
