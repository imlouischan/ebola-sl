## plot: efficacy ##############################################################

# mcmc samples of an alternative model
samples = subset(samples_M48, 
                 flat.TRUE       == flat_TRUE     &
                   epsilonH.TRUE == epsilonH_TRUE &
                   delta.TRUE    == delta_TRUE    &
                   distr.r == distr_r &
                   distr.s == distr_s)

source("3b_pmf_s.R") # serial interval distribution

# preallocate
type_data = rep(NA, nrow(Data))
type_col = rep(NA, nrow(Data))
ph_i = rep(NA, nrow(Data))

# preallocate
efficacy = matrix(rep(NA, nrow(samples)*nrow(Data)), ncol = nrow(Data))

## calculation #################################################################
for (n in 1:nrow(samples)) { # mcmc samples
  
  # sampled parameters
  thetaSS  = c(samples$thetaSS1[n], samples$thetaSS2[n])
  
  for (i in 1:length(Data$i)) { # 49 cases
    
    # preallocate
    efficacy_1 = 0
    efficacy_2 = 0
    
    # hospitalization probability & thickness of bar
    if ( Data_obs[i,3] | Data_obs[i,5] ) { ph = 1; ph_i[i] = 1;   epsilon = samples$epsilon[n]; 
    } else                               { ph = 1; ph_i[i] = 0.5; epsilon = samples$epsilonH[n]; }
    
    # 8 types
    if (                         Data_obs[i,2] &  Data_obs[i,3]            ) { # i: xOOx (1) ====
      type_data[i] = "(1) -OO-"
      type_col[i] = "chartreuse3"
      # data of observed cases
      tauSH = Data$tauiSH[i]
      
      # effectiveness
      epsilon_h = 1 - (1 - epsilon)^ph
      # efficacy
      epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
      # efficacy for the sample
      efficacy[n, i] = epsilon_h_i
      
    } else if (                  Data_obs[i,2] & !Data_obs[i,3]            ) { # i: xOMx (2) # unknown hospitalization, except Case 1 ====
      if (i != 1) {
        type_data[i] = "(2) -OM-"
        type_col[i] = "cornflowerblue"
        # intervals
        a = 0
        b = Data_BH2(i) - Data$tiS[i]
        # change of variable
        tauSH = a:b
        
        # effectiveness
        epsilon_h = 1 - (1 - epsilon)^ph
        # efficacy
        epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
        # efficacy for the sample
        efficacy[n, i] = sum(epsilon_h_i*pmf_SH(tauSH))/sum(pmf_SH(tauSH))
      }
    } else if ( Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: OMOx (3) ====
      type_data[i] = "(3) OMO-"
      type_col[i] = "darkgoldenrod1"
      # intervals
      a = 0
      b = Data_BS2(i) - Data$tiE[i]
      # change of variable
      tauES = a:b
      tauSH = Data$tiH[i] - Data$tiE[i] - tauES
      
      # effectiveness
      epsilon_h = 1 - (1 - epsilon)^ph
      # efficacy
      epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
      # efficacy for the sample
      efficacy[n, i] = sum(epsilon_h_i*pmf_ES(tauES))/sum(pmf_ES(tauES))
      
    } else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: MMOx (4) ====
      type_data[i] = "(4) MMO-"
      type_col[i] = "peachpuff3"
      # intervals
      a = Data$tiH[i] - Data_BS2(i)
      b = Data$tiH[i] - Data_BS1(i)
      # change of variable
      tauSH = a:b
      
      # effectiveness
      epsilon_h = 1 - (1 - epsilon)^ph
      # efficacy
      epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
      # efficacy for the sample
      efficacy[n, i] = sum(epsilon_h_i*pmf_SH(tauSH))/sum(pmf_SH(tauSH))
      
    } else if ( Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: OMMO(D) (5D) # unknown hospitalization ====
      type_data[i] = "(5) OMMO"
      type_col[i] = "mediumorchid2"
      # intervals
      for (tauES in 0:(Data_BS2(i) - Data$tiE[i])) {
        for (tauHD in 0:(Data$tiD[i] - Data$tiE[i] - tauES)) {
          # change of variable
          tauSH = Data$tiD[i] - tauHD - Data$tiE[i] - tauES
          
          # effectiveness
          epsilon_h = 1 - (1 - epsilon)^ph
          # efficacy
          epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
          # numerator and denominator
          efficacy_1 = efficacy_1 + pmf_HD(tauHD)*pmf_ES(tauES)*epsilon_h_i
          efficacy_2 = efficacy_2 + pmf_HD(tauHD)*pmf_ES(tauES)
        }
      }
      # efficacy for the sample
      efficacy[n, i] = efficacy_1/efficacy_2
      
    } else if ( Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] & !Data_obs[i,4]) { # i: OMMM(D) (6D) # unknown hospitalization ====
      type_data[i] = "(6) OMMM"
      type_col[i] = "turquoise3"
      # intervals
      for (tauES in 0:(Data_BS2(i) - Data$tiE[i])) {
        for (tauSH in 0:(Data_tEND - Data$tiE[i] - tauES)) {
          # effectiveness
          epsilon_h = 1 - (1 - epsilon)^ph
          # efficacy
          epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
          # numerator and denominator
          efficacy_1 = efficacy_1 + pmf_SH(tauSH)*pmf_ES(tauES)*epsilon_h_i
          efficacy_2 = efficacy_2 + pmf_SH(tauSH)*pmf_ES(tauES)
        }
      }
      # efficacy for the sample
      efficacy[n, i] = efficacy_1/efficacy_2
      
    } else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: MMMO(D) (7D) # unknown hospitalization ====
      type_data[i] = "(7) MMMO"
      type_col[i] = "wheat4"
      # intervals, 1st part
      for (tauHD in 0:(Data$tiD[i] - Data_BS2(i))) {
        for (tauSH in (Data$tiD[i] - Data_BS2(i) - tauHD):(Data$tiD[i] - Data_BS1(i) - tauHD)) {
          # effectiveness
          epsilon_h = 1 - (1 - epsilon)^ph
          # efficacy
          epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
          # numerator and denominator
          efficacy_1 = efficacy_1 + pmf_SH(tauSH)*pmf_HD(tauHD)*epsilon_h_i
          efficacy_2 = efficacy_2 + pmf_SH(tauSH)*pmf_HD(tauHD)
        }
      }
      # intervals, 2nd part
      for (tauHD in (Data$tiD[i] - Data_BS2(i) + 1):(Data$tiD[i] - Data_BS1(i))) {
        for (tauSH in 0:(Data$tiD[i] - Data_BS1(i) - tauHD)) {
          # effectiveness
          epsilon_h = 1 - (1 - epsilon)^ph
          # efficacy
          epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
          # numerator and denominator
          efficacy_1 = efficacy_1 + pmf_SH(tauSH)*pmf_HD(tauHD)*epsilon_h_i
          efficacy_2 = efficacy_2 + pmf_SH(tauSH)*pmf_HD(tauHD)
        }
      }
      # efficacy for the sample
      efficacy[n, i] = efficacy_1/efficacy_2
      
    } else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,5]) { # i: MMMO(R) (7R) # unknown hospitalization ====
      type_data[i] = "(8) MMMO"
      type_col[i] = "slategray2"
      # intervals, 1st part
      for (tauHR in 0:(Data$tiR[i] - Data_BS2(i))) {
        for (tauSH in (Data$tiR[i] - Data_BS2(i) - tauHR):(Data$tiR[i] - Data_BS1(i) - tauHR)) {
          # effectiveness
          epsilon_h = 1 - (1 - epsilon)^ph
          # efficacy
          epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
          # numerator and denominator
          efficacy_1 = efficacy_1 + pmf_SH(tauSH)*pmf_HR(tauHR)*epsilon_h_i
          efficacy_2 = efficacy_2 + pmf_SH(tauSH)*pmf_HR(tauHR)
        }
      }
      # intervals, 2nd part
      for (tauHR in (Data$tiR[i] - Data_BS2(i) + 1):(Data$tiR[i] - Data_BS1(i))) {
        for (tauSH in 0:(Data$tiR[i] - Data_BS1(i) - tauHR)) {
          # effectiveness
          epsilon_h = 1 - (1 - epsilon)^ph
          # efficacy
          epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
          # numerator and denominator
          efficacy_1 = efficacy_1 + pmf_SH(tauSH)*pmf_HR(tauHR)*epsilon_h_i
          efficacy_2 = efficacy_2 + pmf_SH(tauSH)*pmf_HR(tauHR)
        }
      }
      # efficacy for the sample
      efficacy[n, i] = efficacy_1/efficacy_2
      
    }
  }; print(n); # counter
}

# save efficacy
file_efficacy = paste0("efficacy/efficacy", 
                       "_N", N, 
                       "_f", flat_TRUE, 
                       "_e", epsilonH_TRUE, 
                       "_d", delta_TRUE, 
                       "_R", distr_r, 
                       "_S", distr_s,
                       ".Rdata")
save(efficacy, file = file_efficacy)