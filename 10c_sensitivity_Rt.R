## plot: reproduction number (Rt) ##############################################

# preallocate
R_samples = data.frame()
for (epsilonH_TRUE in c(1, 0)) { # independent / dependent epsilonH
  
  # samples of alternative models (constant R0)
  delta_TRUE = 0
  samples = subset(samples_M48, 
                   flat.TRUE       == flat_TRUE     &
                     epsilonH.TRUE == epsilonH_TRUE &
                     delta.TRUE    == delta_TRUE    &
                     distr.r == distr_r &
                     distr.s == distr_s)
  # 95% CI
  R0_quantile = quantile(samples$R0, c(0.025, 0.5, 0.975))
  R0_mean = R0_quantile[2]
  R0_lower = R0_quantile[1]
  R0_upper = R0_quantile[3]
  
  # samples of alternative models (exponentially decreasing Rt)
  delta_TRUE = 1
  samples = subset(samples_M48, 
                   flat.TRUE       == flat_TRUE     &
                     epsilonH.TRUE == epsilonH_TRUE &
                     delta.TRUE    == delta_TRUE    &
                     distr.r == distr_r &
                     distr.s == distr_s)
  # preallocate
  t = 0:Data_tEND
  Rt_samples = matrix(NA, nrow(samples), length(t), 
                      dimnames = list(NULL, paste0("t", t)))
  # calculation
  for (n in 1:nrow(samples)) { # all samples
    R0 = samples$R0[n]
    delta = samples$delta[n]
    Rt = R0*exp(-delta*t)
    Rt_samples[n, ] = Rt
  }
  # 95% CI
  Rt_quantile = apply(Rt_samples, 2, quantile, c(0.025, 0.5, 0.975))
  Rt_mean = Rt_quantile[2, ]
  Rt_lower = Rt_quantile[1, ]
  Rt_upper = Rt_quantile[3, ]
  
  # Rt and R0
  Rt_samples = data.frame(t = t, 
                          Rt.mean = Rt_mean, Rt.lower = Rt_lower, Rt.upper = Rt_upper, 
                          R0.mean = R0_mean, R0.lower = R0_lower, R0.upper = R0_upper, 
                          epsilonH.TRUE = as.factor(epsilonH_TRUE))
  R_samples = rbind(R_samples, Rt_samples)
  
} # epsilonH_TRUE

# plot name
file_plot <- paste0("figure/10c_sensitivity", 
                    "_N", N, 
                    "_f", flat_TRUE, 
                    # "_e", epsilonH_TRUE, 
                    # "_d", delta_TRUE, 
                    "_R", distr_r, 
                    "_S", distr_s)

# plot
p_rt = ggplot(R_samples, aes(x = t, col = epsilonH.TRUE)) +
  # geom_line(aes(y = R0.lower, linetype = "R0"), size = 0.5) +
  # geom_line(aes(y = Rt.lower, linetype = "Rt"), size = 0.5) +
  # geom_line(aes(y = R0.upper, linetype = "R0"), size = 0.5) +
  # geom_line(aes(y = Rt.upper, linetype = "Rt"), size = 0.5) +
  geom_line(aes(y = R0.mean, linetype = "R0"), size = 1) + 
  geom_line(aes(y = Rt.mean, linetype = "Rt"), size = 1) + 
  scale_linetype_manual(breaks = c("Rt", "R0"), 
                        values = c("Rt" = "solid", "R0" = "dashed"), 
                        labels = c("Variable Rt", "Constant R0")) + 
  scale_color_manual(breaks = c(0, 1), 
                     values = c("0" = "black", "1" = "darkgray"), 
                     labels = c("Dependent epsilonH", "Independent epsilonH")) + 
  xlab("Day") + 
  ylab("Reproduction \n number") + 
  theme_bw(base_size = 20) + 
  theme(legend.justification = c(1, 1), legend.position = c(1, 1), 
        legend.direction = "vertical", legend.box = "horizontal", legend.spacing.y = unit(0, "npc"), 
        legend.key.width = unit(2, "lines"), 
        legend.background = element_rect(fill = NA), 
        legend.title = element_blank())

# multiple plots
p2 = cowplot::plot_grid(p_vi, p_rt, 
                        nrow = 2, 
                        align = "hv")
p8 = cowplot::plot_grid(p6, p2, 
                        nrow = 2, 
                        rel_heights = c(3, 2.5))

# save as eps
print( name_file <- paste0(file_plot, ".eps") )
ggplot2::ggsave(file = name_file, plot = p8, width = 8.5, height = 11)