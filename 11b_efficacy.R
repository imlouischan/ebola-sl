# load efficacy
file_efficacy = paste0("efficacy/efficacy", 
                       "_N", N, 
                       "_f", flat_TRUE, 
                       "_e", epsilonH_TRUE, 
                       "_d", delta_TRUE, 
                       "_R", distr_r, 
                       "_S", distr_s,
                       ".Rdata")
load(file = file_efficacy)

# mcmc samples of an alternative model
samples = subset(samples_M48, 
                 flat.TRUE       == flat_TRUE     &
                   epsilonH.TRUE == epsilonH_TRUE &
                   delta.TRUE    == delta_TRUE    &
                   distr.r == distr_r &
                   distr.s == distr_s)

# preallocate
type_data = rep(NA, nrow(Data))
type_col = rep(NA, nrow(Data))
ph_i = rep(NA, nrow(Data))

# 49 cases
for (i in 1:length(Data$i)) {
  
  # hospitalization probability & thickness of bar
  if ( Data_obs[i,3] | Data_obs[i,5] ) { ph = 1;     ph_i[i] = 1; 
  } else                               { ph = 18/22; ph_i[i] = 0.5; }
  
  # 8 types
  if (                         Data_obs[i,2] &  Data_obs[i,3]            ) { # i: xOOx (1) ====
    type_data[i] = "(1) -OO-"
    type_col[i] = "chartreuse3"
  } else if (                  Data_obs[i,2] & !Data_obs[i,3]            ) { # i: xOMx (2) # unknown hospitalization, except Case 1 ====
    if (i != 1) {
      type_data[i] = "(2) -OM-"
      type_col[i] = "cornflowerblue"
    }
  } else if ( Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: OMOx (3) ====
    type_data[i] = "(3) OMO-"
    type_col[i] = "darkgoldenrod1"
  } else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: MMOx (4) ====
    type_data[i] = "(4) MMO-"
    type_col[i] = "peachpuff3"
  } else if ( Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: OMMO(D) (5D) # unknown hospitalization ====
    type_data[i] = "(5) OMMO"
    type_col[i] = "mediumorchid2"
  } else if ( Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] & !Data_obs[i,4]) { # i: OMMM(D) (6D) # unknown hospitalization ====
    type_data[i] = "(6) OMMM"
    type_col[i] = "turquoise3"
  } else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: MMMO(D) (7D) # unknown hospitalization ====
    type_data[i] = "(7) MMMO"
    type_col[i] = "wheat4"
  } else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,5]) { # i: MMMO(R) (7R) # unknown hospitalization ====
    type_data[i] = "(8) MMMO"
    type_col[i] = "slategray2"
  }
}

# adding effectiveness with known and unknown hospitalization
eps <- cbind(samples$epsilon, samples$epsilonH, efficacy)

type_data <- c(NA, NA, type_data)
type_col <- c("white", "white", type_col)
ph_i <- c(1, 0.5, ph_i)
lwd_i <- c(rep(3, 2), rep(2, 49))
## plot ########################################################################
# plot name
file_plot <- paste0("figure/11_efficacy", 
                    "_N", N, 
                    "_f", flat_TRUE, 
                    "_e", epsilonH_TRUE, 
                    "_d", delta_TRUE, 
                    "_R", distr_r, 
                    "_S", distr_s)
# save as eps
setEPS()
print( name_file <- paste0(file_plot, ".eps") )
postscript(name_file, width = 8.5, height = 11)

# efficacy
boxplot(eps, 
        col = type_col, 
        lwd = lwd_i, 
        cex.axis = 1.5, 
        xlab = "", 
        ylab = "",
        ylim = c(0, 1), yaxt = "n", 
        width = ph_i, 
        outline = F, 
        horizontal = T, add = F)
# title(ylab = "Case ID (Infector ID)", line = 2.5, cex.lab = 2)
title(xlab = "Protective effect", line = 2.5, cex.lab = 2)

# Case ID & Infector ID
axis(2, at = c(1:nrow(Data))+2, las = 2, cex.axis = 1, font = 1, 
     labels = paste0(Data$i, "(", Data$vi, ")"))
# effectiveness
axis(2, at = c(1:2), las = 2, cex.axis = 0.8, font = 2, 
     labels = c("epsilon", "epsilonH"))

# effectiveness with known and unknown hospitalization
abline(v = quantile(samples$epsilon, c(0.25, 0.5, 0.75)),
       lwd = c(2, 2, 2), lty = c("dashed", "solid", "dashed"), col = "darkred")
abline(v = quantile(samples$epsilonH, c(0.25, 0.5, 0.75)),
       lwd = c(2, 2, 2), lty = c("dashed", "solid", "dashed"), col = "darkblue")

boxplot(eps, 
        col = type_col, 
        lwd = lwd_i, 
        cex.axis = 1.5, 
        xlab = "", 
        ylab = "",
        ylim = c(0, 1), yaxt = "n", 
        width = ph_i, 
        outline = F, 
        horizontal = T, add = T)

# legend
legend("topright", bg = "white", 
       legend = c("(1) -OO-", "(2) -OM-", "(3) OMO-", "(4) MMO-", 
                  "(5) OMMO", "(6) OMMM", "(7) MMMO", "(8) MMMO", 
                  "without delay"), 
       fill = c("chartreuse3", "cornflowerblue", "darkgoldenrod1", "peachpuff3", 
                "mediumorchid2", "turquoise3", "wheat4", "slategray2", 
                "white"), 
       horiz = F, cex = 1)

# save as eps
dev.off()