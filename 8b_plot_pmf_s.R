## plot: serial interval distribution ##########################################
# plot name
file_plot <- paste0("figure/8b_pmf_s", 
                    "_N", N, 
                    "_f", flat_TRUE, 
                    "_e", epsilonH_TRUE, 
                    "_d", delta_TRUE, 
                    "_R", distr_r, 
                    "_S", distr_s)
# save as eps
setEPS()
print( name_file <- paste0(file_plot, ".eps") )
postscript(name_file, width = 11, height = 8.5)

source("3b_pmf_s.R") # serial interval distribution

# # observed data
# plot(prop.table(table(tauiSS)),
#      xlim = c(0, max(tauiSS)),
#      xlab = "Days", 
#      ylab = "Probability", 
#      col = "black", lwd = 10)

# preallocate
pmf_temp = matrix(rep(NA, nrow(samples)*(max(tauiSS, 22)+1)), ncol = max(tauiSS, 22)+1)
colnames(pmf_temp) = 0:max(tauiSS, 22)

# calculation
for (n in 1:nrow(samples)) { # all samples
  thetaSS = c(samples$thetaSS1[n], samples$thetaSS2[n])
  pmf_temp[n, ] = pmf_s(tau = 0:max(tauiSS, 22), thetaSS)
}

# estimated pmf
boxplot(pmf_temp, 
        # col = rgb(1, 1, 1, alpha = 0), 
        # border = "red", 
        # outcol = "darkgray", 
        outline = F, 
        lwd = 2, 
        cex = 0.5, 
        cex.axis = 1.5, 
        cex.lab = 2, 
        xlab = "", 
        ylab = "", 
        xlim = c(0, max(tauiSS, 22)), 
        at = 0:max(tauiSS, 22))
title(xlab = "Serial interval", line = 2.5, cex.lab = 2)
title(ylab = "Probability", line = 2.5, cex.lab = 2)

# save as eps
dev.off()