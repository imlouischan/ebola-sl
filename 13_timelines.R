## plot: timelines #############################################################
# http://www.statmethods.net/graphs/line.html
# http://www.statmethods.net/advgraphs/parameters.html
# http://www.statmethods.net/graphs/dot.html

# plot name
file_plot <- paste0("figure/13_timelines")
# save as eps
setEPS()
print( name_file <- paste0(file_plot, ".eps") )
postscript(name_file, width = 8.5, height = 11)

# empty plot
plot(0, 0, type = "n",
     xlim = c(0, Data_tEND), 
     ylim = c(1, nrow(Data)), 
     xlab = "", 
     ylab = "", 
     xaxt = "n", 
     yaxt = "n") 
# title(ylab = "Case ID (Infector ID)", line = 2.5, cex.lab = 2)
title(xlab = "Day", line = 2.5, cex.lab = 2)
axis(1, at = seq(0, Data_tEND, 10), 
     labels = seq(0, Data_tEND, 10), 
     cex.axis = 1.5)

# Case ID & Infector ID
axis(2, at = c(1:nrow(Data)), las = 2, cex.axis = 1,
     labels = paste(Data$i, "(", Data$vi, ")", sep=""))

# the enlarged size of estimated time points for better visualization
enlarge = 1

for (i in 1:nrow(Data)) { # 49 cases
  
  # background horizontal lines
  points(0:Data_tEND, rep(i, Data_tEND+1),
         pch = ".", col = "gray")
  
  # all observed time points
  points(Data$tiE[i], i,
         col = "darkgray", pch = 91, font = 2, cex = 1)
  points(Data$tiD[i], i,
         col = "black"   , pch = 93, font = 2, cex = 1)
  points(Data$tiR[i], i,
         col = "darkgreen"  , pch = 93, font = 2, cex = 1)
  points(Data$tiH[i], i,
         col = "blue"    , lwd = 2, pch = 1, cex = 1.5)
  points(Data$tiS[i], i,
         col = "red"     , lwd = 2, pch = 1, cex = 1.5)
  
  # estimated time points, with PMFs
  if (                         Data_obs[i,2] &  Data_obs[i,3]            ) { # i: xOOx (1) ====
    type = 1
  } else if (                  Data_obs[i,2] & !Data_obs[i,3]            ) { # i: xOMx (2) # unknown hospitalization, except Case 1 ====
    type = 2
    if (i != 1) {
      # intervals
      a = 0
      b = Data_BH2(i) - Data$tiS[i]
      # change of variable
      tauSH = a:b
      
      # possible time points
      points(Data$tiS[i] + tauSH, rep(i, length(tauSH)),
             col = "blue", font = 1, pch = 47,
             cex = enlarge)
    }
  } else if ( Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: OMOx (3) ====
    type = 3
    # intervals
    a = 0
    b = Data_BS2(i) - Data$tiE[i]
    # change of variable
    tauES = a:b
    tauSH = Data$tiH[i] - Data$tiE[i] - tauES
    
    # possible time points
    points(Data$tiE[i] + tauES, rep(i, length(tauES)),
           col = "red", font = 1, pch = 92,
           cex = enlarge)
    
  } else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: MMOx (4) ====
    type = 4
    # intervals
    a = Data$tiH[i] - Data_BS2(i)
    b = Data$tiH[i] - Data_BS1(i)
    # change of variable
    tauSH = a:b
    
    # possible time points
    points(Data$tiH[i] - tauSH, rep(i, length(tauSH)),
           col = "red", font = 1, pch = 92,
           cex = enlarge)
    
  } else if ( Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: OMMO(D) (5D) # unknown hospitalization ====
    type = 5
    # intervals
    a = Data$tiE[i]
    b = Data_BS2(i)
    c = Data$tiD[i]
    # change of variable
    tiS_hat = a:b
    tiH_hat = a:c
    
    # possible time points
    points(tiH_hat, rep(i, length(tiH_hat)),
           col = "blue", font = 1, pch = 47,
           cex = enlarge)
    points(tiS_hat, rep(i, length(tiS_hat)),
           col = "red", font = 1, pch = 92,
           cex = enlarge)
    
  } else if ( Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] & !Data_obs[i,4]) { # i: OMMM(D) (6D) # unknown hospitalization ====
    type = 6
    # intervals
    a = Data$tiE[i]
    b = Data_BS2(i)
    c = Data_tEND
    # change of variable
    tiS_hat = a:b
    tiH_hat = a:c
    
    # possible time points
    points(tiH_hat, rep(i, length(tiH_hat)),
           col = "blue", font = 1, pch = 47,
           cex = enlarge)
    points(tiS_hat, rep(i, length(tiS_hat)),
           col = "red", font = 1, pch = 92,
           cex = enlarge)
    
  } else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: MMMO(D) (7D) # unknown hospitalization ====
    type = 7
    # intervals
    a = Data_BS1(i)
    b = Data_BS2(i)
    c = Data$tiD[i]
    # change of variable
    tiS_hat = a:b
    tiH_hat = a:c
    
    # possible time points
    points(tiH_hat, rep(i, length(tiH_hat)),
           col = "blue", font = 1, pch = 47,
           cex = enlarge)
    points(tiS_hat, rep(i, length(tiS_hat)),
           col = "red", font = 1, pch = 92,
           cex = enlarge)
    
  } else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,5]) { # i: MMMO(R) (7R) # unknown hospitalization ====
    type = 8
    # intervals
    a = Data_BS1(i)
    b = Data_BS2(i)
    c = Data$tiR[i]
    # change of variable
    tiS_hat = a:b
    tiH_hat = a:c
    
    # possible time points
    points(tiH_hat, rep(i, length(tiH_hat)),
           col = "blue", font = 1, pch = 47,
           cex = enlarge)
    points(tiS_hat, rep(i, length(tiS_hat)),
           col = "red", font = 1, pch = 92,
           cex = enlarge)
    
  }
  
  # category
  cate <- c("(1) -OO-", "(2) -OM-", "(3) OMO-", "(4) MMO-", 
            "(5) OMMO", "(6) OMMM", "(7) MMMO", "(8) MMMO")
  # title
  # print( main <- paste0("Case ",i, " - ", cate[type]) )
  
}

# legend, both observed and estimated
legend("bottomright", bg = "white", 
       legend = c("Date of exposure", 
                  "Date of symptom onset (observed)", 
                  "Date of symptom onset (estimated)", 
                  "Date of hospitalization (observed)", 
                  "Date of hospitalization (estimated)", 
                  "Date of death", 
                  "Date of discharge"),
       col = c("darkgray", "red", "red", "blue", "blue", "black", "darkgreen"), 
       pch = c(91, 1, 92, 1, 47, 93, 93), 
       lwd = c(1, 2, 1, 2, 1, 1, 1), 
       lty = 0, 
       cex = 1.5
)

# save as eps
dev.off()