## parameter comparison ########################################################
# plot name
file_plot <- paste0("figure/10d_infector_flat", 
                    "_N", N, 
                    # "_f", flat_TRUE, 
                    "_e", epsilonH_TRUE, 
                    "_d", delta_TRUE, 
                    "_R", distr_r, 
                    "_S", distr_s)

# mcmc samples of two scenarios
samples_comp = subset(samples_M48, 
                      # flat.TRUE       == flat_TRUE     &
                      epsilonH.TRUE == epsilonH_TRUE &
                        delta.TRUE    == delta_TRUE    &
                        distr.r == distr_r &
                        distr.s == distr_s)
# switch factor order
samples_comp$flat.TRUE <- factor(samples_comp$flat.TRUE, levels = c("0", "1"))

# plots
library(ggplot2)

p_vi = ggplot(samples_comp, aes(x = reorder(vi, vi, function(x) - length(x)), linetype = flat.TRUE, fill = flat.TRUE)) +
  geom_bar(aes(y = (..count..)/(sum(..count..)/2)), position = "dodge", size = 1, col = "black") +
  xlab("Infector ID") + 
  ylab("Probability") + 
  scale_linetype_manual(breaks = c(0, 1),
                        values = c("0" = "solid", "1" = "dotted"), 
                        labels = c("Pij proposal", "Flat proposal")) +
  scale_fill_manual(breaks = c(0, 1),
                    values = c("0" = "darkgray", "1" = "lightgray"),
                    labels = c("Pij proposal", "Flat proposal")) +
  theme_bw(base_size = 20) + 
  theme(legend.justification = c(1, 1), legend.position = c(1, 1), 
        legend.direction = "vertical", legend.box = "horizontal", legend.spacing.y = unit(0, "npc"), 
        legend.background = element_rect(fill = NA),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, vjust = 0.5))

# save as eps
print( name_file <- paste0(file_plot, ".eps") )
ggplot2::ggsave(file = name_file, plot = p_vi, width = 11, height = 8.5, device = cairo_ps)