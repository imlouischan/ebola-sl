## MCMC ########################################################################

function (N, chain, flat_TRUE, epsilonH_TRUE, delta_TRUE) {
  
  # control randomness, choose whatever you like
  set.seed(chain + 2018)
  
  # call the function: log-posterior (or log-likelihood)
  logL = dget("4c_logL.R")
  # call the function: transmission probability
  pij = dget("4d_pij.R")
  
  # preallocate
  s.thetaSS = matrix(rep(c(NA, NA) ,N), ncol = 2)
  s.R0      = rep(NA ,N)
  s.size    = rep(NA ,N)
  s.delta   = rep(NA ,N)
  s.epsilon = matrix(rep(c(NA, NA) ,N), ncol = 2)
  s.vi      = rep(NA ,N)
  s.logL    = rep(NA ,N)
  
  # counter of acceptance
  count_acc = 0
  # timer
  timer_start = proc.time()
  
  ## start MCMC ##################################################################
  for (n in 1:N) {
    
    if (n == 1) { # initial parameters
      
      # fixed initial parameters
      thetaSS = rep(NA ,2)
      thetaSS[1] = 4
      thetaSS[2] = 15
      R0         = 1
      size       = 1
      delta      = 0.01
      epsilon = rep(NA ,2)
      epsilon[1] = 0.5
      epsilon[2] = 0.5
      
    } else { # proposal distribution
      
      # Gaussian proposal distribution
      thetaSS = rnorm(2, mean = s.thetaSS[n - 1, ], sd = c(0.5, 0.5))
      R0      = rnorm(1, mean = s.R0[n - 1],        sd = 0.5)
      size    = rnorm(1, mean = s.size[n - 1],      sd = 0.5)
      delta   = rnorm(1, mean = s.delta[n - 1],     sd = 0.01)
      epsilon = rnorm(2, mean = s.epsilon[n - 1, ], sd = c(0.1, 0.1))
      
    }
    
    # uniform prior, for ddgamma & ddweibull & dnbinom
    if (thetaSS[1] > 0 & thetaSS[1] < 20 & 
        thetaSS[2] > 0 & thetaSS[2] < 20 & 
        R0         > 0 & R0         < 48 & 
        size       > 0 & size       < 10 & 
        delta      > 0 & delta      < 1  & 
        epsilon[1] > 0 & epsilon[1] < 1  & 
        epsilon[2] > 0 & epsilon[2] < 1  ){
      
      # alternative model: proposed infector of Case 41
      if (flat_TRUE == 1) { # flat proposal (faster)
        possible <- setdiff(1:nrow(Data), c(38, 41:45, 48:49)) # 41 possible cases
        Data$vi[Data$i == 41] = sample(Data$i[possible], size = 1)
        # or, fixed infector as Case 25
        # Data$vi[Data$i == 41] = 25
      } else { # transmission prob $p_{ij}$ proposal (slower)
        Data$vi[Data$i == 41] = sample(Data$i, size = 1, prob = pij(thetaSS, R0, epsilon, Data))
      } # alternative proposal
      
      # alternative model: dependent epsilonH
      if (epsilonH_TRUE != 1) epsilon[2] = 1 - (1 - epsilon[1])^(18/22)
      
      # alternative model: constant R0
      if (!delta_TRUE) delta = 0
      
      # alternative model: not negative binomial distributed secondary case number
      if (distr_r != "nbinom") size = 1
      
      ## likelihood calculation ######################################################
      
      # log-likelihood
      logL_prop = logL(thetaSS, R0, size, delta, epsilon, Data)
      
    } else { # prior, out of range
      logL_prop = -Inf
    }
    
    ## accept or reject ############################################################
    
    # acceptance probability
    ( logA = min(0, logL_prop - s.logL[n - 1]) )
    ( logu = log(runif(1)) )
    
    # accept or reject
    if (n == 1 | logu < logA) { # initial value and accept
      s.thetaSS[n, ] = thetaSS
      s.R0[n]        = R0
      s.size[n]      = size
      s.delta[n]     = delta
      s.epsilon[n, ]   = epsilon
      s.vi[n]        = Data$vi[Data$i == 41]
      s.logL[n]      = logL_prop
      # counter for acceptance
      count_acc = count_acc + 1
    } else { # reject
      s.thetaSS[n, ] = s.thetaSS[n - 1, ]
      s.R0[n]        = s.R0[n - 1]
      s.size[n]      = s.size[n - 1]
      s.delta[n]     = s.delta[n - 1]
      s.epsilon[n, ]   = s.epsilon[n - 1, ]
      s.vi[n]        = s.vi[n - 1]
      s.logL[n]      = s.logL[n - 1]
    }
    
    # counter
    if ( (n*100)%%N == 0 ) {
      percent = paste0(n/N*100, "% (", n, "/", N, ")")
      timer = paste0(round((proc.time() - timer_start)/60), " min")
      acc_rate = paste0("accepted ", round(count_acc/n*100), "% (", count_acc, "/", n, ")")
      print( paste(percent, timer[1], acc_rate) )
    }
  }
  
  ## thinning and save ###########################################################
  
  # output sample number
  N.thin = 1000
  if (N >= N.thin) { # thinning & burn-in
    thin = seq(1, N, by = N/N.thin) + N/N.thin - 1
  } else { # no thinning & no burn-in
    thin = 1:N
  }
  
  # data frame
  samples = data.frame(thetaSS1 = s.thetaSS[thin, 1], 
                       thetaSS2 = s.thetaSS[thin, 2], 
                       R0       = s.R0[thin], 
                       size     = s.size[thin], 
                       delta    = s.delta[thin], 
                       epsilon  = s.epsilon[thin, 1], 
                       epsilonH = s.epsilon[thin, 2], 
                       vi       = s.vi[thin], 
                       logL     = s.logL[thin])
  
  # save mcmc samples
  file_samples = paste0("samples", 
                       "_N", N, 
                       "_f", flat_TRUE, 
                       "_e", epsilonH_TRUE, 
                       "_d", delta_TRUE, 
                       "_R", distr_r, 
                       "_S", distr_s,
                       "_C", chain, 
                       ".Rdata")
  save(samples, file = file_samples)
  
  # output
  return(samples)
  
}