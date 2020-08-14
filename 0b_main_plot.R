## new script ##################################################################
cat("\014") # clear the console
rm(list=ls()) # remove all variables
graphics.off() # close all plots
setwd("../ebola_sl/")
# install.packages("gdata", "ggplot2", "reshape2", "gridExtra", "GGally", "cowplot")
## call data and results #######################################################

source("1a_data.R") # data set
source("1b_bound.R") # boundaries of missing time points

source("2a_incubation.R") # incubation period
source("2b_distribution.R") # estimation of time lag distributions

source("13_timelines.R") # plot timelines
source("14_imputation.R") # plot imputation

# sample number of a chain
N = 1e5 # can try the minimun N = 2e3 (or less but without thinning)

# load all samples
file_samples = paste0("samples_M48", 
                      "_N", N, 
                      ".Rdata")
load(file = file_samples)

## plot (best model) ###########################################################

# the best model
flat_TRUE = 0; epsilonH_TRUE = 0; delta_TRUE = 1; 
distr_r = "geom"; distr_s = "dweibull"; 

# mcmc samples of an alternative model
samples = subset(samples_M48, 
                 flat.TRUE       == flat_TRUE     &
                   epsilonH.TRUE == epsilonH_TRUE &
                   delta.TRUE    == delta_TRUE    &
                   distr.r == distr_r &
                   distr.s == distr_s)

# mean & 95CI
t(data.frame(
  epsilon  = round(quantile(samples$epsilon *100, c(0.5, 0.025, 0.975))*10)/10, 
  epsilonH = round(quantile(samples$epsilonH*100, c(0.5, 0.025, 0.975))*10)/10, 
  R0       = round(quantile(samples$R0          , c(0.5, 0.025, 0.975))*10)/10, 
  delta    = round(quantile(samples$delta       , c(0.5, 0.025, 0.975))*1000)/1000, 
  thetaSS1 = round(quantile(samples$thetaSS1    , c(0.5, 0.025, 0.975))*10)/10, 
  thetaSS2 = round(quantile(samples$thetaSS2    , c(0.5, 0.025, 0.975))*10)/10))

source("7a_plot_mixing.R") # plot: mixing
source("7b_plot_dependence.R", echo = F) # plot: dependence

source("8a_plot_pmf_r.R") # plot: offspring distributions
source("8b_plot_pmf_s.R") # plot: serial interval distributions
source("8c_plot_pmf_s.R") # plot: serial interval distributions with data

source("9a_cal_pmf_r.R") # calculate: mean & SD of offspring distribution
source("9b_cal_pmf_s.R") # calculate: mean & SD of serial interval

source("12_infector.R") # plot: possible infectors

for (epsilonH_TRUE in c(1, 0)) { # independent / dependent epsilonH
  for (delta_TRUE in c(1, 0)) { # exponentially decreasing Rt / constant R0
    
    # source("11a_efficacy.R") # calculate: efficacy (long computation time ~ 2 hours)
    source("11b_efficacy.R") # plot: efficacy
    
  } # delta_TRUE
} # epsilonH_TRUE

## plot (alternative models) ###################################################

source("10b_sensitivity.R") # plot: comparison of parameters
source("10c_sensitivity_Rt.R") # plot: comparison of reproduction number over time
source("10d_infector_flat.R") # plot: comparison of Case 41 infector

source("10a_rank.R") # calculate: model ranking