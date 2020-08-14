## new script ##################################################################
cat("\014") # clear the console
rm(list=ls()) # remove all variables
graphics.off() # close all plots
setwd("../ebola_sl/")

## main: MCMC ##################################################################

source("1a_data.R") # data set
source("1b_bound.R") # boundaries of missing time points

source("2a_incubation.R") # incubation period
source("2b_distribution.R") # estimation of time lag distributions

# alternative models
for (flat_TRUE in c(1, 0)) { # flat proposal / pij
  for (epsilonH_TRUE in c(1, 0)) { # independent / dependent epsilonH
    for (delta_TRUE in c(1, 0)) { # exponentially decreasing Rt / constant R0
      for (distr_r in c("geom", "nbinom", "pois")) { # alternative offspring distributions
        for (distr_s in c("dweibull", "dgamma")) { # alternative serial interval distributions
          
          # sample number of a chain
          N = 1e5 # can try the minimun N = 1e3 (or less but without thinning)
          
          # counter
          print(paste0("******************************", 
                       "_N", N, 
                       "_f", flat_TRUE, 
                       "_e", epsilonH_TRUE, 
                       "_d", delta_TRUE, 
                       "_R", distr_r, 
                       "_S", distr_s, 
                       "******************************"))
          
          source("3a_pmf_r.R") # offspring distribution
          source("3b_pmf_s.R") # serial interval distribution
          
          # call the function to run one mcmc chain
          run_chain = dget("5_mcmc.R")
          run_chains = function(chain) {run_chain(N, chain, flat_TRUE, epsilonH_TRUE, delta_TRUE)}
          
          # option A: local running (recommended to work with batch submitting jobs)
          # samples = run_chains(0)
          
          # option B: parallel running (not recommended)
          source("6a_parallel.R")
          
        } # distr_s
      } # distr_r
    } # delta_TRUE
  } # epsilonH_TRUE
} # flat_TRUE

## main: summary ###############################################################

# preallocate
samples_M48 = data.frame(NULL)

for (flat_TRUE in c(1, 0)) { # flat proposal / pij
  for (epsilonH_TRUE in c(1, 0)) { # independent / dependent epsilonH
    for (delta_TRUE in c(1, 0)) { # exponentially decreasing Rt / constant R0
      for (distr_r in c("geom", "nbinom", "pois")) { # alternative offspring distributions
        for (distr_s in c("dweibull", "dgamma")) { # alternative serial interval distributions
          
          # combine chains of each alternative model
          comb_chain = dget("6b_combine_chains.R")
          samples = comb_chain(N, 1:no_cores)
          
          # combine alternative models
          samples_M48 = rbind(samples_M48, samples)
          
        } # distr_s
      } # distr_r
    } # delta_TRUE
  } # epsilonH_TRUE
} # flat_TRUE

# save all combined samples
file_samples = paste0("samples_M48", 
                      "_N", N, 
                      ".Rdata")
save(samples_M48, file = file_samples)