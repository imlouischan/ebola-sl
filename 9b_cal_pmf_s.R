## calculate: mean & SD of serial interval #####################################

# 2 parameters
a = samples$thetaSS1
b = samples$thetaSS2

# following discrete Gamma
if (distr_s == "dgamma") {
  mean = a*b
  var = a*b^2
}
# following discrete Weibull
if (distr_s == "dweibull") {
  mean = b*gamma(1+1/a)
  var = b^2*(gamma(1+2/a)-gamma(1+1/a)^2)
}

# mean & 95CI of mean & sd of serial interval
q_mean = quantile(mean,      c(0.025, 0.5, 0.975))
q_sd   = quantile(sqrt(var), c(0.025, 0.5, 0.975))

# output
print(paste("pmf_s", 
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