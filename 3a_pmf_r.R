## alternative offspring distributions #########################################

# following geometric distribution
# pmf = p*(1-p)^x
if (distr_r == "geom") {
  pmf_r = function(ri, R, size) dgeom(ri, prob = 1/(R + 1))
}

# following Poisson distribution
if (distr_r == "pois") {
  pmf_r = function(ri, R, size) dpois(ri, lambda = R)
}

# following negative binomial (2-parameter) distribution
# if size == 1, nbinom == geom
# if size == Inf, nbinom == pois
if (distr_r == "nbinom") {
  pmf_r = function(ri, R, size) dnbinom(ri, size = size, mu = R)
}

## reproduction number with isolation ##########################################

R_hat = function(tauSH, thetaSS, R0, epsilon, ph) {
  
  # effectiveness
  epsilon_h = 1 - (1 - epsilon)^ph
  
  # efficacy
  epsilon_h_i = epsilon_h*(1 - cdf_s(tauSH - 1, thetaSS))
  
  # reproduction number
  R_hat = R0*(1 - epsilon_h_i)
  
  return(R_hat)
}