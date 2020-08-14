## plot: mixing ################################################################
# plot name
file_plot <- paste0("figure/7a_mixing", 
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

# alternative models
samples_index <- c("thetaSS1", "thetaSS2", "R0", "size", "delta", "epsilon", "epsilonH", "vi", "logL")
# if (distr_r != "nbinom") { samples_index = setdiff(samples_index, "size") }
# if (!delta_TRUE)         { samples_index = setdiff(samples_index, "delta") }
# if (!epsilonH_TRUE)      { samples_index = setdiff(samples_index, "epsilonH") }

par(mfrow=c(3, 3)) # subplot

for (i in 1:length(samples_index)) {
  plot(samples[, samples_index[i]], xlab = "", ylab = samples_index[i]) # plot mixing
}

par(mfrow=c(1,1)) # back to one plot

# save as eps
dev.off()