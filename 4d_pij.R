## transmission probability $p_{ij}$ ###########################################

function (thetaSS, R0, epsilon, Data) {
  
  # call the function: likelihood of serial interval
  L_si = dget("4b_L_si.R")
  
  # preallocate
  s_hat41 = rep(NA ,nrow(Data))
  # excluding Case 41 (located at row38) and other impossible cases
  s_hat41[38] = 0    # infected himself
  s_hat41[41:45] = 0 # later illness onset
  s_hat41[48:49] = 0 # requiring long computation time (& low probabilities)
  
  for (j in setdiff(1:nrow(Data), c(38, 41:45, 48:49))){ # excluding Case 41 and 7 others
    
    # the missing infector
    Data$vi[Data$i == 41] = Data$i[j]
    
    # serial interval using proposed parameters
    if ( j == 1 ) { # the index case who was not isolated
      s_hat41[j] = L_si(i = 38, thetaSS, R0, NULL, Data)
    } else if (  Data_obs[j,3] |  Data_obs[j,5] ) { # cases with known hospitalization
      s_hat41[j] = L_si(i = 38, thetaSS, R0, epsilon[1], Data)
    } else if ( !Data_obs[j,3] & !Data_obs[j,5] ) { # cases without known hospitalization
      s_hat41[j] = L_si(i = 38, thetaSS, R0, epsilon[2], Data)
    } else { print("wrong at pij"); break; }
    
  }
  
  # transmission probability
  p_41_j = s_hat41/sum(s_hat41)
  
  # output the transmission probability
  if ( all(is.finite(p_41_j)) ) {
    return(p_41_j)
  } else { # uniform probability
    return(NULL)
  }
  
}