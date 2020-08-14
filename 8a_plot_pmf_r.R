## plot: offspring distribution ################################################
# plot name
file_plot <- paste0("figure/8a_pmf_r", 
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

source("3a_pmf_r.R") # offspring distribution

Data$vi[38] # not include the estimated in the data

# observed
ri_df = as.data.frame(table(factor(Data$vi, levels = Data$i)))
( ri = ri_df$Freq )

# # observed data
# plot(prop.table(table(ri)),
#      xlim = c(-0.5, max(ri)+0.5),
#      xlab = "Offspring number", 
#      ylab = "Probability", 
#      col = "black", lwd = 10)

# preallocate
pmf_temp = matrix(rep(NA, nrow(samples)*(max(ri, 9)+1)), ncol = max(ri, 9)+1)
colnames(pmf_temp) = 0:max(ri, 9)

# calculation
for (n in 1:nrow(samples)) { # all samples
  R0 = samples$R0[n]
  size = samples$size[n]
  pmf_temp[n, ] = pmf_r(ri = 0:max(ri, 9), R0, size)
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
        # xlim = c(0, max(ri)), 
        # ylim = c(0, max(pmf_temp)), 
        at = 0:max(ri, 9))
title(xlab = "Secondary case", line = 2.5, cex.lab = 2)
title(ylab = "Probability", line = 2.5, cex.lab = 2)

# save as eps
dev.off()