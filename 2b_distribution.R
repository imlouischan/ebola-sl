## MLE parameter estimation of time lag distributions ##########################

## time lag
## (1) tauSH: time lags from illness onset to hospitalization
## (2) tauSD: time lags from illness onset to death
## (3) tauSR: time lags from illness onset to recovery
## (4) tauHD: time lags from hospitalization to death
## (5) tauHR: time lags from hospitalization to recovery
## (6) tauSS: serial intervals (in the presence of isolation)

## alternative distributions (pmf) and MLE
## (a) gamma distribution
## (b) weibull distribution

# alternative distributions
pmf_ddgamma   = function(tau, theta)   pgamma(tau+1, shape = theta[1], scale = theta[2]) -   pgamma(tau, shape = theta[1], scale = theta[2])
cdf_ddgamma   = function(tau, theta)   pgamma(tau+1, shape = theta[1], scale = theta[2])
pmf_ddweibull = function(tau, theta) pweibull(tau+1, shape = theta[1], scale = theta[2]) - pweibull(tau, shape = theta[1], scale = theta[2])
cdf_ddweibull = function(tau, theta) pweibull(tau+1, shape = theta[1], scale = theta[2])

pmf_alt = function(distr_index) {
  if (distr_index == 1) return(pmf_ddgamma)
  if (distr_index == 2) return(pmf_ddweibull)
}

# alternative observed time lags
taui_alt = function(taui_index) {
  if (taui_index == 1) return(tauiSH)
  if (taui_index == 2) return(tauiSD)
  if (taui_index == 3) return(tauiSR)
  if (taui_index == 4) return(tauiHD)
  if (taui_index == 5) return(tauiHR)
  if (taui_index == 6) return(tauiSS)
}
taui_name <- c("Time from illness onset to hospitalization", 
               "Time from illness onset to death", 
               "Time from illness onset to discharge", 
               "Time from hospitalization to death", 
               "Time from hospitalization to discharge", 
               "Serial interval")

# preallocate
MLE_alt   = list()
distr_MLE = list()
theta_MLE = list()

for (taui_index in 1:6) { # different observed time lags
  
  # plot name
  file_plot <- paste0("figure/2b_distribution", 
                      "_", taui_index)
  # save as eps
  setEPS()
  print( name_file <- paste0(file_plot, ".eps") )
  postscript(name_file, width = 11, height = 8.5)
  
  # preallocate
  MLE_alt[[taui_index]] = data.frame(
    dgamma   = rep(NA, 3),
    dweibull = rep(NA, 3),
    row.names = c("nlogL", "theta1", "theta2")
  )
  
  # # plot observed time lags
  # plot(prop.table(table(taui_alt(taui_index))),
  #      xlim = c(0, max(taui_alt(taui_index), 18)),
  #      cex.axis = 1.5, 
  #      xaxt = "n", 
  #      xlab = "", 
  #      ylab = "", 
  #      col = "black", lwd = 3)
  # plot frame
  par(mar=c(5, 5, 2, 5) + 0.1, mgp=c(3, 1, 0), las=0)
  plot(table(taui_alt(taui_index)),
       xlim = c(0, max(taui_alt(taui_index), 18)),
       ylim = c(0, max(prop.table(table(taui_alt(taui_index))))), 
       cex.axis = 1.5, 
       type = "n", 
       xaxt = "n", 
       xlab = "", 
       ylab = "", 
       col = "black", lwd = 3)
  axis(1, at = 0:max(taui_alt(taui_index), 18), labels = 0:max(taui_alt(taui_index), 18), cex.axis = 1.5)
  title(xlab = taui_name[taui_index], line = 2.5, cex.lab = 2)
  title(ylab = "Probability", line = 2.5, cex.lab = 2)
  
  for (distr_index in 1:2) { # alternative distributions
    
    # negative log-likelihood
    nlogL = function(theta) -sum(log(pmf_alt(distr_index)(taui_alt(taui_index), theta)))
    
    # parameter estimation using optim
    theta_optim = suppressWarnings( optim(c(1, 1), nlogL) )
    MLE_alt[[taui_index]]["nlogL", distr_index] = theta_optim$value
    MLE_alt[[taui_index]][c("theta1", "theta2"), distr_index] = theta_optim$par
    
    # plot estimated pmf
    points(0:max(taui_alt(taui_index), 18), 
           pmf_alt(distr_index)(tau = 0:max(taui_alt(taui_index), 18), theta = theta_optim$par), 
           cex = 2, 
           lwd = 3, 
           pch = distr_index, 
           col = "darkgray")
  }  
  
  # checking
  print(MLE_alt[[taui_index]])
  
  # the best distribution and estimated parameters
  distr_MLE[[taui_index]] = which.min(MLE_alt[[taui_index]]["nlogL", ])
  theta_MLE[[taui_index]] = MLE_alt[[taui_index]][c("theta1", "theta2"), distr_MLE[[taui_index]]]
  
  # best pmf
  if (taui_index == 1) pmf_SH = function(tau) pmf_alt(distr_MLE[[1]])(tau, theta_MLE[[1]])
  if (taui_index == 2) pmf_SD = function(tau) pmf_alt(distr_MLE[[2]])(tau, theta_MLE[[2]])
  if (taui_index == 3) pmf_SR = function(tau) pmf_alt(distr_MLE[[3]])(tau, theta_MLE[[3]])
  if (taui_index == 4) pmf_HD = function(tau) pmf_alt(distr_MLE[[4]])(tau, theta_MLE[[4]])
  if (taui_index == 5) pmf_HR = function(tau) pmf_alt(distr_MLE[[5]])(tau, theta_MLE[[5]])
  if (taui_index == 6) pmf_SS = function(tau) pmf_alt(distr_MLE[[6]])(tau, theta_MLE[[6]])
  
  # plot best pmf
  points(0:max(taui_alt(taui_index), 18), 
         pmf_alt(distr_MLE[[taui_index]])(tau = 0:max(taui_alt(taui_index), 18), theta = theta_MLE[[taui_index]]), 
         cex = 2, 
         lwd = 3, 
         pch = distr_MLE[[taui_index]], 
         col = "red")
  
  # observed incubation period
  par(new = T)
  hist(taui_alt(taui_index), 
       xlim = c(0, max(taui_alt(taui_index), 18)),
       ylim = c(0, max(table(taui_alt(taui_index)))), 
       breaks = seq(min(taui_alt(taui_index))-0.5, max(taui_alt(taui_index))+0.5, by = 1), 
       main = "", 
       xaxt = "n", 
       yaxt = "n", 
       xlab = "", 
       ylab = "")
  axis(4, at = 0:max(table(taui_alt(taui_index))), cex.axis = 1.5)
  mtext(side = 4, "Frequency", line = 2.5, cex = 2)
  
  # legend (Weibull is the best)
  legend(ifelse(sum( taui_index == c(1, 2, 4) ), "topright", "topleft"), 
         legend = c("data", "gamma", "Weibull (best fit)"), 
         col = c("black", "darkgray", "red"), 
         pch = c(NA, 1, 2), 
         lty = c(1, 0, 0), 
         lwd = 3, 
         cex = 2)
  
  # save as eps
  dev.off()
  
} # different observed time lags