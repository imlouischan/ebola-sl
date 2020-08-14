## incubation period ############################################################

# results from WHO, \cite{team2014ebola}
tauES_mean = 9.1; tauES_var = 7.3^2;
# shape & scale parameters, converted using continuous Gamma
thetaES = c(tauES_mean^2 / tauES_var, tauES_var / tauES_mean)

# discretization
pmf_ES = function(tauES) pgamma(tauES+1, shape = thetaES[1], scale = thetaES[2]) - pgamma(tauES, shape = thetaES[1], scale = thetaES[2])
cdf_ES = function(tauES) pgamma(tauES+1, shape = thetaES[1], scale = thetaES[2])

## plot ########################################################################

# plot name
file_plot <- paste0("figure/2a_incubation")
# save as eps
setEPS()
print( name_file <- paste0(file_plot, ".eps") )
postscript(name_file, width = 11, height = 8.5)

# estimated pmf
par(mar=c(5, 5, 2, 5) + 0.1, mgp=c(3, 1, 0), las=0)
plot(0:max(tauiES, 18), pmf_ES(tauES = 0:max(tauiES, 18)), 
     xlim = c(0, max(tauiES, 18)), 
     ylim = c(0, max(prop.table(table(tauiES)))), 
     cex.axis = 1.5, 
     xaxt = "n", 
     xlab = "", 
     ylab = "", 
     cex = 2, 
     lwd = 3, 
     col = "blue")
axis(1, at = 0:max(tauiES, 18), labels = 0:max(tauiES, 18), cex.axis = 1.5)
title(xlab = "Incubation period", line = 2.5, cex.lab = 2)
title(ylab = "Probability", line = 2.5, cex.lab = 2)

# observed incubation period
par(new = T)
hist(tauiES, 
     xlim = c(0, max(tauiES, 18)), 
     ylim = c(0, max(table(tauiES))), 
     breaks = seq(min(tauiES)-0.5, max(tauiES)+0.5, by = 1), 
     main = "", 
     xaxt = "n", 
     yaxt = "n", 
     xlab = "", 
     ylab = "")
axis(4, at = 0:max(table(tauiES)), cex.axis = 1.5)
mtext(side = 4, "Frequency", line = 2.5, cex = 2)

# legend
legend("topleft", 
       legend = c("data", "gamma"), 
       col = c("black", "blue"), 
       pch = c(NA, 1), 
       lty = c(1, 0), 
       lwd = 3, 
       cex = 2
)

# save as eps
dev.off()