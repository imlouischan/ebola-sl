## full log-likelihood (or posterior) ##########################################

function (thetaSS, R0, size, delta, epsilon, Data) {
  
  # call the function: likelihood of secondary case
  L_ri = dget("4a_L_ri.R")
  # call the function: likelihood of serial interval
  L_si = dget("4b_L_si.R")
  
  # preallocate
  L_ri_prop = rep(NA ,nrow(Data))
  L_si_prop = rep(NA ,nrow(Data))
  
  # each case
  for (i in 1:length(Data$i)) {
    
    # likelihood of secondary case using proposed parameters
    if ( i == 1 ) { # the index case who was not isolated
      L_ri_prop[i] = L_ri(i, thetaSS, R0, size, delta, NULL, Data)
    } else if (  Data_obs[i,3] |  Data_obs[i,5] ) { # 26 cases with known hospitalization
      L_ri_prop[i] = L_ri(i, thetaSS, R0, size, delta, epsilon[1], Data)
    } else if ( !Data_obs[i,3] & !Data_obs[i,5] ) { # 22 cases without known hospitalization
      L_ri_prop[i] = L_ri(i, thetaSS, R0, size, delta, epsilon[2], Data)
    } else { print("wrong at logL_ri"); break; }
    
    # infector of i
    vi = which(Data$i == Data$vi[i])
    # likelihood of serial interval using proposed parameters
    if ( length(vi) == 0 ) { # 9 imported cases
      L_si_prop[i] = 1
    } else if ( vi == 1 ) { # 3 cases infected by the index case who was not isolated
      L_si_prop[i] = L_si(i, thetaSS, R0, NULL, Data)
    } else if (  Data_obs[vi,3] |  Data_obs[vi,5] ) { # 21 cases infected by case with known hospitalization
      L_si_prop[i] = L_si(i, thetaSS, R0, epsilon[1], Data)
    } else if ( !Data_obs[vi,3] & !Data_obs[vi,5] ) { # 15 cases infected by case without known hospitalization
      L_si_prop[i] = L_si(i, thetaSS, R0, epsilon[2], Data)
    } else { print("wrong at logL_si"); break; }
  }
  
  # logL_r of all cases
  ( logL_r_prop = sum(log(L_ri_prop)) )
  
  # logL_s of all cases
  ( logL_s_prop = sum(log(L_si_prop)) )
  
  # log-likelihood
  ( logL_prop = logL_r_prop + logL_s_prop )
  
  # output
  if ( all(is.finite(logL_prop)) ) {
    return(logL_prop)
  } else { # uniform probability
    return(-1e10)
  }
  
}