## alternative serial interval distributions ###################################

# following discrete Gamma
if (distr_s == "dgamma") {
  pmf_s = pmf_ddgamma
  cdf_s = cdf_ddgamma
}
# following discrete Weibull
if (distr_s == "dweibull") {
  pmf_s = pmf_ddweibull
  cdf_s = cdf_ddweibull
}

## serial interval distribution with isolation #################################

pmf_s_hat = function(tauSS, tauSH, thetaSS, epsilon, ph) {
  
  # effectiveness
  epsilon_h = 1 - (1 - epsilon)^ph
  
  # efficacy
  epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
  
  # serial interval distribution
  pmf_s_hat = ifelse(tauSS < tauSH, 
                                     pmf_s(tauSS, thetaSS)/(1 - epsilon_h_i),
                     (1 - epsilon_h)*pmf_s(tauSS, thetaSS)/(1 - epsilon_h_i)
  )
  
  return(pmf_s_hat)
}