## plot: dependence ############################################################
# plot name
file_plot <- paste0("figure/7b_pairwise", 
                    "_N", N, 
                    "_f", flat_TRUE, 
                    "_e", epsilonH_TRUE, 
                    "_d", delta_TRUE, 
                    "_R", distr_r, 
                    "_S", distr_s)

# alternative models
samples_index <- c("thetaSS1", "thetaSS2", "R0", "size", "delta", "epsilon", "epsilonH")
if (distr_r != "nbinom") { samples_index = setdiff(samples_index, "size") }
if (!delta_TRUE)         { samples_index = setdiff(samples_index, "delta") }
if (!epsilonH_TRUE)      { samples_index = setdiff(samples_index, "epsilonH") }

# plot parameter dependence
library(GGally)
library(ggplot2)
plot <- GGally::ggpairs(samples[, samples_index], 
                        lower = list(continuous = GGally::wrap("points", alpha = 0.1)), 
                        upper = list(continuous = GGally::wrap("cor", col = "black", size = 8)),
                        diag  = list(continuous = GGally::wrap("densityDiag", col = "black", lwd = 1))) + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1)) + 
  theme_bw(base_size = 20)
# plot <- GGally::ggpairs(samples[, samples_index], 
#                         lower = list(continuous = GGally::wrap("points", alpha = 0.1)), 
#                         upper = list(continuous = "cor"),
#                         diag = list(continuous = "densityDiag")) + 
#   ggplot2::theme_bw()

# save as eps
print( name_file <- paste0(file_plot, ".eps") )
ggplot2::ggsave(file = name_file, plot = plot, width = 11, height = 8.5, device = cairo_ps)