## plot: unobserved infector ###################################################
# plot name
file_plot <- paste0("figure/12_infector", 
                    "_N", N, 
                    "_f", flat_TRUE, 
                    "_e", epsilonH_TRUE, 
                    "_d", delta_TRUE, 
                    "_R", distr_r, 
                    "_S", distr_s)

# plot possible infectors
library(ggplot2)
p <- 
  ggplot(samples) + 
  geom_bar(mapping = aes(x = reorder(vi, vi, function(x) - length(x)), y = ..prop.., group = 1), stat = "count") + 
  scale_y_continuous(labels = scales::percent_format()) + 
  xlab("Infector ID") + 
  ylab("Probability") + 
  theme_bw(base_size = 36) + 
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

# save as eps
print( name_file <- paste0(file_plot, ".eps") )
ggplot2::ggsave(file = name_file, plot = p, width = 11, height = 8.5)