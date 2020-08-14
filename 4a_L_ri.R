## likelihood of offspring number ##############################################

function (i, thetaSS, R0, size, delta, epsilon, Data) {
  
  # preallocate
  L_ri_1 = 0
  L_ri_2 = 0
  
  # observed offspring number
  ri_df = as.data.frame(table(factor(Data$vi, levels = Data$i)))
  ri = ri_df$Freq
  
  # hospitalization probability
  ph = 1
  # if ( Data_obs[i,3] | Data_obs[i,5] ) { ph = 1
  # } else                               { ph = 18/22 }
  
  # categories in L_r
  # i: "xOOx", "xOMx", "OMOx", "MMOx", "OMMO(D)", "OMMM(D)", "MMMO(D)", "MMMO(R)"
  if (                       Data_obs[i,2] &  Data_obs[i,3]                 ) { # i: xOOx (1) ====
    type.i = "(1) xOOx   "
    # data of observed cases
    tauSH = Data$tauiSH[i]
    # exponentially decreasing R
    R = R0*exp(-delta*Data$tiS[i])
    # numerator and denominator
    L_ri_1 = pmf_r(ri[i], R_hat(tauSH, thetaSS, R, epsilon, ph), size)
    L_ri_2 = 1
  } else if (                  Data_obs[i,2] & !Data_obs[i,3]                 ) { # i: xOMx (2) # unknown hospitalization, except Case 1 ====
    type.i = "(2) xOMx   "
    # intervals
    a = 0
    b = Data_BH2(i) - Data$tiS[i]
    # change of variable
    tauSH = a:b
    # exponentially decreasing R
    R = R0*exp(-delta*Data$tiS[i])
    # numerator and denominator
    if (i == 1) {
      L_ri_1 = pmf_r(ri[i], R, size)
      L_ri_2 = 1
    }
    else        {
      L_ri_1 = sum( pmf_SH(tauSH)*pmf_r(ri[i], R_hat(tauSH, thetaSS, R, epsilon, ph), size) )
      L_ri_2 = sum( pmf_SH(tauSH) )
    }
  } else if ( Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]                 ) { # i: OMOx (3) ====
    type.i = "(3) OMOx   "
    # intervals
    a = 0
    b = Data_BS2(i) - Data$tiE[i]
    # change of variable
    tauES = a:b
    tauSH = Data$tiH[i] - Data$tiE[i] - tauES
    # exponentially decreasing R
    R = R0*exp(-delta*(Data$tiE[i] + tauES))
    # numerator and denominator
    L_ri_1 = sum( pmf_ES(tauES)*pmf_r(ri[i], R_hat(tauSH, thetaSS, R, epsilon, ph), size) )
    L_ri_2 = sum( pmf_ES(tauES) )
  } else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]                 ) { # i: MMOx (4) ====
    type.i = "(4) MMOx   "
    # intervals
    a = Data$tiH[i] - Data_BS2(i)
    b = Data$tiH[i] - Data_BS1(i)
    # change of variable
    tauSH = a:b
    # exponentially decreasing R
    R = R0*exp(-delta*(Data$tiH[i] - tauSH))
    # numerator and denominator
    L_ri_1 = sum( pmf_SH(tauSH)*pmf_r(ri[i], R_hat(tauSH, thetaSS, R, epsilon, ph), size) )
    L_ri_2 = sum( pmf_SH(tauSH) )
  } else if ( Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: OMMO(D) (5D) # unknown hospitalization ====
    type.i = "(5) OMMO(D)"
    # intervals
    for (tauES in 0:(Data_BS2(i) - Data$tiE[i])) {
      for (tauHD in 0:(Data$tiD[i] - Data$tiE[i] - tauES)) {
        # change of variable
        tauSH = Data$tiD[i] - tauHD - Data$tiE[i] - tauES
        # exponentially decreasing R
        R = R0*exp(-delta*(Data$tiE[i] + tauES))
        # numerator and denominator
        L_ri_1 = L_ri_1 + pmf_HD(tauHD)*pmf_ES(tauES)*pmf_r(ri[i], R_hat(tauSH, thetaSS, R, epsilon, ph), size)
        L_ri_2 = L_ri_2 + pmf_HD(tauHD)*pmf_ES(tauES)
      }
    }
  } else if ( Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] & !Data_obs[i,4]) { # i: OMMM(D) (6D) # unknown hospitalization ====
    type.i = "(6) OMMM(D)"
    # intervals
    for (tauES in 0:(Data_BS2(i) - Data$tiE[i])) {
      for (tauSH in 0:(Data_tEND - Data$tiE[i] - tauES)) {
        # exponentially decreasing R
        R = R0*exp(-delta*(Data$tiE[i] + tauES))
        # numerator and denominator
        L_ri_1 = L_ri_1 + pmf_SH(tauSH)*pmf_ES(tauES)*pmf_r(ri[i], R_hat(tauSH, thetaSS, R, epsilon, ph), size)
        L_ri_2 = L_ri_2 + pmf_SH(tauSH)*pmf_ES(tauES)
      }
    }
  } else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: MMMO(D) (7D) # unknown hospitalization ====
    type.i = "(7) MMMO(D)"
    # intervals, 1st part
    for (tauHD in 0:(Data$tiD[i] - Data_BS2(i))) {
      for (tauSH in (Data$tiD[i] - Data_BS2(i) - tauHD):(Data$tiD[i] - Data_BS1(i) - tauHD)) {
        # exponentially decreasing R
        R = R0*exp(-delta*(Data$tiD[i] - tauHD - tauSH))
        # numerator and denominator
        L_ri_1 = L_ri_1 + pmf_SH(tauSH)*pmf_HD(tauHD)*pmf_r(ri[i], R_hat(tauSH, thetaSS, R, epsilon, ph), size)
        L_ri_2 = L_ri_2 + pmf_SH(tauSH)*pmf_HD(tauHD)
      }
    }
    # intervals, 2nd part
    for (tauHD in (Data$tiD[i] - Data_BS2(i) + 1):(Data$tiD[i] - Data_BS1(i))) {
      for (tauSH in 0:(Data$tiD[i] - Data_BS1(i) - tauHD)) {
        # exponentially decreasing R
        R = R0*exp(-delta*(Data$tiD[i] - tauHD - tauSH))
        # numerator and denominator
        L_ri_1 = L_ri_1 + pmf_SH(tauSH)*pmf_HD(tauHD)*pmf_r(ri[i], R_hat(tauSH, thetaSS, R, epsilon, ph), size)
        L_ri_2 = L_ri_2 + pmf_SH(tauSH)*pmf_HD(tauHD)
      }
    }
  } else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,5]) { # i: MMMO(R) (7R) # unknown hospitalization ====
    type.i = "(7) MMMO(R)"
    # intervals, 1st part
    for (tauHR in 0:(Data$tiR[i] - Data_BS2(i))) {
      for (tauSH in (Data$tiR[i] - Data_BS2(i) - tauHR):(Data$tiR[i] - Data_BS1(i) - tauHR)) {
        # exponentially decreasing R
        R = R0*exp(-delta*(Data$tiR[i] - tauHR - tauSH))
        # numerator and denominator
        L_ri_1 = L_ri_1 + pmf_SH(tauSH)*pmf_HR(tauHR)*pmf_r(ri[i], R_hat(tauSH, thetaSS, R, epsilon, ph), size)
        L_ri_2 = L_ri_2 + pmf_SH(tauSH)*pmf_HR(tauHR)
      }
    }
    # intervals, 2nd part
    for (tauHR in (Data$tiR[i] - Data_BS2(i) + 1):(Data$tiR[i] - Data_BS1(i))) {
      for (tauSH in 0:(Data$tiR[i] - Data_BS1(i) - tauHR)) {
        # exponentially decreasing R
        R = R0*exp(-delta*(Data$tiR[i] - tauHR - tauSH))
        # numerator and denominator
        L_ri_1 = L_ri_1 + pmf_SH(tauSH)*pmf_HR(tauHR)*pmf_r(ri[i], R_hat(tauSH, thetaSS, R, epsilon, ph), size)
        L_ri_2 = L_ri_2 + pmf_SH(tauSH)*pmf_HR(tauHR)
      }
    }
  } else { print(i); break } # error(?)
  
  # resulting value
  L_ri = L_ri_1/L_ri_2
  
  # output
  return(L_ri)
  
}