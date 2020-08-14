## parameter comparison ########################################################
# plot name
file_plot <- paste0("figure/10b_sensitivity", 
                    "_N", N, 
                    "_f", flat_TRUE, 
                    # "_e", epsilonH_TRUE, 
                    # "_d", delta_TRUE, 
                    "_R", distr_r, 
                    "_S", distr_s)

# mcmc samples of two scenarios
samples_comp = subset(samples_M48, 
                      flat.TRUE       == flat_TRUE     &
                        # epsilonH.TRUE == epsilonH_TRUE &
                        # delta.TRUE    == delta_TRUE    &
                        distr.r == distr_r &
                        distr.s == distr_s)

# plots
library(ggplot2)

p_thetaSS1 = ggplot(samples_comp, aes(x = thetaSS1, linetype = delta.TRUE, col = epsilonH.TRUE)) + 
  geom_line(stat = "density", size = 1) + 
  scale_linetype_manual(breaks = c(1, 0), 
                        values = c("1" = "solid", "0" = "dashed"), 
                        labels = c("Variable Rt", "Constant R0")) + 
  scale_color_manual(breaks = c(0, 1), 
                     values = c("0" = "black", "1" = "darkgray"), 
                     labels = c("Dependent epsilonH", "Independent epsilonH")) + 
  ylab("Density") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")

p_thetaSS2 = ggplot(samples_comp, aes(x = thetaSS2, linetype = delta.TRUE, col = epsilonH.TRUE)) + 
  geom_line(stat = "density", size = 1) + 
  scale_linetype_manual(breaks = c(1, 0), 
                        values = c("1" = "solid", "0" = "dashed"), 
                        labels = c("Variable Rt", "Constant R0")) + 
  scale_color_manual(breaks = c(0, 1), 
                     values = c("0" = "black", "1" = "darkgray"), 
                     labels = c("Dependent epsilonH", "Independent epsilonH")) + 
  ylab("Density") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")

p_R0 = ggplot(samples_comp, aes(x = R0, linetype = delta.TRUE, col = epsilonH.TRUE)) + 
  geom_line(stat = "density", size = 1) + 
  scale_linetype_manual(breaks = c(1, 0), 
                        values = c("1" = "solid", "0" = "dashed"), 
                        labels = c("Variable Rt", "Constant R0")) + 
  scale_color_manual(breaks = c(0, 1), 
                     values = c("0" = "black", "1" = "darkgray"), 
                     labels = c("Dependent epsilonH", "Independent epsilonH")) + 
  ylab("Density") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")
# theme(legend.justification = c(1, 1), legend.position = c(1, 1), 
#       legend.spacing.y = unit(0, "npc"), 
#       legend.background = element_rect(fill = NA), 
#       legend.title = element_blank())

p_size = ggplot(samples_comp[samples_comp$distr.r == "nbinom", ], aes(x = size, linetype = delta.TRUE, col = epsilonH.TRUE)) + 
  geom_line(stat = "density", size = 1) + 
  scale_linetype_manual(breaks = c(1, 0), 
                        values = c("1" = "solid", "0" = "dashed"), 
                        labels = c("Variable Rt", "Constant R0")) + 
  scale_color_manual(breaks = c(0, 1), 
                     values = c("0" = "black", "1" = "darkgray"), 
                     labels = c("Dependent epsilonH", "Independent epsilonH")) + 
  ylab("Density") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")

p_delta = ggplot(samples_comp[samples_comp$delta.TRUE == 1, ], aes(x = delta, linetype = delta.TRUE, col = epsilonH.TRUE)) + 
  geom_line(stat = "density", size = 1) + 
  scale_linetype_manual(breaks = c(1, 0), 
                        values = c("1" = "solid", "0" = "dashed"), 
                        labels = c("Variable Rt", "Constant R0")) + 
  scale_color_manual(breaks = c(0, 1), 
                     values = c("0" = "black", "1" = "darkgray"), 
                     labels = c("Dependent epsilonH", "Independent epsilonH")) + 
  ylab("Density") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")

p_epsilon = ggplot(samples_comp, aes(x = epsilon, linetype = delta.TRUE, col = epsilonH.TRUE)) + 
  geom_line(stat = "density", size = 1) + 
  scale_linetype_manual(breaks = c(1, 0), 
                        values = c("1" = "solid", "0" = "dashed"), 
                        labels = c("Variable Rt", "Constant R0")) + 
  scale_color_manual(breaks = c(0, 1), 
                     values = c("0" = "black", "1" = "darkgray"), 
                     labels = c("Dependent epsilonH", "Independent epsilonH")) + 
  ylab("Density") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")

p_epsilonH = ggplot(samples_comp, aes(x = epsilonH, linetype = delta.TRUE, col = epsilonH.TRUE)) + 
  geom_line(stat = "density", size = 1) + 
  scale_linetype_manual(breaks = c(1, 0), 
                        values = c("1" = "solid", "0" = "dashed"), 
                        labels = c("Variable Rt", "Constant R0")) + 
  scale_color_manual(breaks = c(0, 1), 
                     values = c("0" = "black", "1" = "darkgray"), 
                     labels = c("Dependent epsilonH", "Independent epsilonH")) + 
  ylab("Density") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")

p_vi = ggplot(samples_comp, aes(x = reorder(vi, vi, function(x) - length(x)), linetype = delta.TRUE, fill = epsilonH.TRUE)) + 
  geom_bar(aes(y = (..count..)/(sum(..count..)/4)), position = "dodge", col = "black") + 
  scale_linetype_manual(breaks = c(1, 0), 
                        values = c("1" = "solid", "0" = "dashed"), 
                        labels = c("Variable Rt", "Constant R0"), 
                        guide = guide_legend(override.aes = list(fill = NA))) + 
  scale_fill_manual(breaks = c(0, 1), 
                    values = c("0" = "darkgray", "1" = "lightgray"), 
                    labels = c("Dependent epsilonH", "Independent epsilonH")) + 
  xlab("Infector ID") + 
  ylab("Probability") + 
  theme_bw(base_size = 20) + 
  theme(legend.justification = c(1, 1), legend.position = c(1, 1), 
        legend.direction = "vertical", legend.box = "horizontal", legend.spacing.y = unit(0, "npc"), 
        legend.background = element_rect(fill = NA), 
        legend.title = element_blank(), 
        axis.text.x = element_text(angle = 90, vjust = 0.5))

p_logL = ggplot(samples_comp, aes(x = logL, linetype = delta.TRUE, col = epsilonH.TRUE)) + 
  geom_line(stat = "density", size = 1) + 
  scale_linetype_manual(breaks = c(1, 0), 
                        values = c("1" = "solid", "0" = "dashed"), 
                        labels = c("Variable Rt", "Constant R0")) + 
  scale_color_manual(breaks = c(0, 1), 
                     values = c("0" = "black", "1" = "darkgray"), 
                     labels = c("Dependent epsilonH", "Independent epsilonH")) + 
  ylab("Density") + 
  theme_bw(base_size = 20) + 
  theme(legend.position = "none")

# multiple plots
p6 = cowplot::plot_grid(p_thetaSS1, p_thetaSS2, p_R0, p_delta, p_epsilon, p_epsilonH, 
                        nrow = 3, 
                        align = "hv")
p7 = cowplot::plot_grid(p6, p_vi, 
                        nrow = 2, 
                        rel_heights = c(3, 1))

# save as eps
print( name_file <- paste0(file_plot, ".eps") )
ggplot2::ggsave(file = name_file, plot = p7, width = 8.5, height = 11, device = cairo_ps)