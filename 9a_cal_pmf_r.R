## calculate: mean & SD of offspring distribution ##############################

# 2 parameters
R = samples$R0
s = samples$size
delta = samples$delta

# following geometric distribution
if (distr_r == "geom") {
  var = R*(1+R)
}

# following Poisson distribution
if (distr_r == "pois") {
  var = R
}

# following negative binomial distribution
if (distr_r == "nbinom") {
  var = R + R^2/s
}

# mean & 95CI of mean & sd of R0
q_mean  = quantile(R,         c(0.025, 0.5, 0.975))
q_sd    = quantile(sqrt(var), c(0.025, 0.5, 0.975))
q_delta = quantile(delta,     c(0.025, 0.5, 0.975))

# output
print(paste("pmf_r", 
            "_N", N, 
            "_f", flat_TRUE, 
            "_e", epsilonH_TRUE, 
            "_d", delta_TRUE, 
            "_R", distr_r, 
            "_S", distr_s, 
            sep=""))
print(paste("mean with 95CI: ", 
            round(q_mean[2], 1), " ", "(", round(q_mean[1], 1), ", ", round(q_mean[3], 1), ")", 
            sep=""))
print(paste("sd   with 95CI: ", 
            round(q_sd[2], 1), " ", "(", round(q_sd[1], 1), ", ", round(q_sd[3], 1), ")", 
            sep=""))
print(paste("delta with 95CI: ", 
            round(q_delta[2], 3), " ", "(", round(q_delta[1], 3), ", ", round(q_delta[3], 3), ")", 
            sep=""))