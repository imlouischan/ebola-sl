# http://www.sthda.com/english/wiki/ggplot2-quick-correlation-matrix-heatmap-r-software-and-data-visualization
library(ggplot2)

## plot: imputation ############################################################

plot <- p <- list()
for (i in 1:nrow(Data)) { # 49 cases
  
  # estimated time points, with PMFs
  if (                         Data_obs[i,2] &  Data_obs[i,3]            ) { # i: xOOx (1) ====
    type = 1
    p[[i]] <- NA
  } else if (                  Data_obs[i,2] & !Data_obs[i,3]            ) { # i: xOMx (2) # unknown hospitalization, except Case 1 ====
    type = 2
    if (i != 1) {
      # time lag interval
      a = 0
      b = Data_BH2(i) - Data$tiS[i]
      tauSH = a:b
      # key time events
      tiS_hat = Data$tiS[i]
      tiH_hat = Data$tiS[i] + tauSH
      # possible time points
      p[[i]] <- matrix(NA, nrow = length(tiS_hat), ncol = length(tiH_hat))
      rownames(p[[i]]) <- tiS_hat
      colnames(p[[i]]) <- tiH_hat
      for (tauSH in a:b) {
        # estimates
        tiS_hat = Data$tiS[i]
        tiH_hat = Data$tiS[i] + tauSH
        # probability
        p[[i]][as.character(tiS_hat), as.character(tiH_hat)] <- pmf_SH(tauSH)
      }
    } else { p[[i]] <- NA }
  } else if ( Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: OMOx (3) ====
    type = 3
    # time lag interval
    a = 0
    b = Data_BS2(i) - Data$tiE[i]
    tauES = a:b
    # key time events
    tiS_hat = Data$tiE[i] + tauES
    tiH_hat = Data$tiH[i]
    # possible time points
    p[[i]] <- matrix(NA, nrow = length(tiS_hat), ncol = length(tiH_hat))
    rownames(p[[i]]) <- tiS_hat
    colnames(p[[i]]) <- tiH_hat
    for (tauES in a:b) {
      # estimates
      tiS_hat = Data$tiE[i] + tauES
      tiH_hat = Data$tiH[i]
      # probability
      p[[i]][as.character(tiS_hat), as.character(tiH_hat)] <- pmf_ES(tauES)
    }
  } else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]            ) { # i: MMOx (4) ====
    type = 4
    # time lag interval
    a = Data$tiH[i] - Data_BS2(i)
    b = Data$tiH[i] - Data_BS1(i)
    tauSH = a:b
    # key time events
    tiS_hat = sort(Data$tiH[i] - tauSH) # rev
    tiH_hat = Data$tiH[i]
    # possible time points
    p[[i]] <- matrix(NA, nrow = length(tiS_hat), ncol = length(tiH_hat))
    rownames(p[[i]]) <- tiS_hat
    colnames(p[[i]]) <- tiH_hat
    for (tauSH in a:b) {
      # estimates
      tiS_hat = Data$tiH[i] - tauSH
      tiH_hat = Data$tiH[i]
      # probability
      p[[i]][as.character(tiS_hat), as.character(tiH_hat)] <- pmf_SH(tauSH)
    }
  } else if ( Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: OMMO(D) (5D) # unknown hospitalization ====
    type = 5
    # time lag intervals
    a = Data$tiE[i]
    b = Data_BS2(i)
    c = Data$tiD[i]
    # key time events
    tiS_hat = a:b
    tiH_hat = a:c
    # possible time points
    p[[i]] <- matrix(NA, nrow = length(tiS_hat), ncol = length(tiH_hat))
    rownames(p[[i]]) <- tiS_hat
    colnames(p[[i]]) <- tiH_hat
    for (tauES in 0:(Data_BS2(i) - Data$tiE[i])) {
      for (tauHD in 0:(Data$tiD[i] - Data$tiE[i] - tauES)) {
        # estimates
        tiS_hat <- Data$tiE[i] + tauES
        tiH_hat <- Data$tiD[i] - tauHD
        # probability
        p[[i]][as.character(tiS_hat), as.character(tiH_hat)] <- pmf_HD(tauHD)*pmf_ES(tauES)
      }
    }
  } else if ( Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] & !Data_obs[i,4]) { # i: OMMM(D) (6D) # unknown hospitalization ====
    type = 6
    # time lag intervals
    a = Data$tiE[i]
    b = Data_BS2(i)
    c = Data_tEND
    # key time events
    tiS_hat = a:b
    tiH_hat = a:c
    # possible time points
    p[[i]] <- matrix(NA, nrow = length(tiS_hat), ncol = length(tiH_hat))
    rownames(p[[i]]) <- tiS_hat
    colnames(p[[i]]) <- tiH_hat
    for (tauES in 0:(Data_BS2(i) - Data$tiE[i])) {
      for (tauSH in 0:(Data_tEND - Data$tiE[i] - tauES)) {
        # estimates
        tiS_hat <- Data$tiE[i] + tauES
        tiH_hat <- tiS_hat + tauSH
        # probability
        p[[i]][as.character(tiS_hat), as.character(tiH_hat)] <- pmf_SH(tauSH)*pmf_ES(tauES)
      }
    }
  } else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: MMMO(D) (7D) # unknown hospitalization ====
    type = 7
    # time lag intervals
    a = Data_BS1(i)
    b = Data_BS2(i)
    c = Data$tiD[i]
    # key time events
    tiS_hat = a:b
    tiH_hat = a:c
    # possible time points
    p[[i]] <- matrix(NA, nrow = length(tiS_hat), ncol = length(tiH_hat))
    rownames(p[[i]]) <- tiS_hat
    colnames(p[[i]]) <- tiH_hat
    for (tauHD in 0:(Data$tiD[i] - Data_BS2(i))) {
      for (tauSH in (Data$tiD[i] - Data_BS2(i) - tauHD):(Data$tiD[i] - Data_BS1(i) - tauHD)) {
        # estimates
        tiH_hat <- Data$tiD[i] - tauHD
        tiS_hat <- tiH_hat - tauSH
        # probability
        p[[i]][as.character(tiS_hat), as.character(tiH_hat)] <- pmf_SH(tauSH)*pmf_HD(tauHD)
      }
    }
    for (tauHD in (Data$tiD[i] - Data_BS2(i) + 1):(Data$tiD[i] - Data_BS1(i))) {
      for (tauSH in 0:(Data$tiD[i] - Data_BS1(i) - tauHD)) {
        # estimates
        tiH_hat <- Data$tiD[i] - tauHD
        tiS_hat <- tiH_hat - tauSH
        # probability
        p[[i]][as.character(tiS_hat), as.character(tiH_hat)] <- pmf_SH(tauSH)*pmf_HD(tauHD)
      }
    }
  } else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,5]) { # i: MMMO(R) (7R) # unknown hospitalization ====
    type = 8
    # time lag intervals
    a = Data_BS1(i)
    b = Data_BS2(i)
    c = Data$tiR[i]
    # key time events
    tiS_hat = a:b
    tiH_hat = a:c
    # possible time points
    p[[i]] <- matrix(NA, nrow = length(tiS_hat), ncol = length(tiH_hat))
    rownames(p[[i]]) <- tiS_hat
    colnames(p[[i]]) <- tiH_hat
    for (tauHR in 0:(Data$tiR[i] - Data_BS2(i))) {
      for (tauSH in (Data$tiR[i] - Data_BS2(i) - tauHR):(Data$tiR[i] - Data_BS1(i) - tauHR)) {
        # estimates
        tiH_hat <- Data$tiR[i] - tauHR
        tiS_hat <- tiH_hat - tauSH
        # probability
        p[[i]][as.character(tiS_hat), as.character(tiH_hat)] <- pmf_SH(tauSH)*pmf_HR(tauHR)
      }
    }
    for (tauHR in (Data$tiR[i] - Data_BS2(i) + 1):(Data$tiR[i] - Data_BS1(i))) {
      for (tauSH in 0:(Data$tiR[i] - Data_BS1(i) - tauHR)) {
        # estimates
        tiH_hat <- Data$tiR[i] - tauHR
        tiS_hat <- tiH_hat - tauSH
        # probability
        p[[i]][as.character(tiS_hat), as.character(tiH_hat)] <- pmf_SH(tauSH)*pmf_HR(tauHR)
      }
    }
  }
  
  # category
  cate <- c("(1) -OO-", "(2) -OM-", "(3) OMO-", "(4) MMO-", 
            "(5) OMMO", "(6) OMMM", "(7) MMMO", "(8) MMMO")
  # title
  print( main <- paste0("Case ", Data$i[i], ": ", cate[type]) )
  # normalization
  p[[i]] <- p[[i]] / sum(p[[i]], na.rm = T)
  # melt the matrix
  df <- reshape2::melt(p[[i]], 
                       na.rm = T, 
                       varnames = c("tiS", "tiH"), 
                       value.name = "prob")
  if ( nrow(df) < 10 ) {
    df$tiS <- as.character(df$tiS)
    df$tiH <- as.character(df$tiH)
  }
  # plot
  if ( type == 1 ) {
    plot[[i]] <- ggplot() + 
      labs(title = main)
  } else if ( type == 2 ) {
    if (i != 1) {
      plot[[i]] <- ggplot(data = df, aes(tiH, prob)) + 
        labs(title = main, 
             x = "Date of hospitalization", 
             y = "Probability") + 
        geom_col() + 
        # theme(aspect.ratio = 1) + 
        theme_bw()
    } else {
      plot[[i]] <- ggplot() + 
        labs(title = main)
    }
  } else if ( type == 3 | type == 4 ) {
    plot[[i]] <- ggplot(data = df, aes(tiS, prob)) + 
      labs(title = main, 
           x = "Date of symptom onset", 
           y = "Probability") + 
      geom_col() + 
      # theme(aspect.ratio = 1) + 
      theme_bw()
  } else if ( type > 4 ) {
    plot[[i]] <- ggplot(data = df, aes(tiS, tiH, fill = prob)) + 
      labs(title = main, 
           x = "Date of symptom onset", 
           y = "Date of hospitalization") + 
      scale_x_continuous(breaks = scales::pretty_breaks()) + 
      scale_y_continuous(breaks = scales::pretty_breaks()) + 
      geom_tile() + 
      coord_fixed() + 
      theme_bw()
  }
}

## plot: imputation ############################################################

# save figures
for ( type in 2:8 ) {
  
  # save eps plot
  setEPS()
  name_file <- paste0("figure/14_imputation_", type, ".eps")
  postscript(name_file, width = 11, height = 8.5)
  # plot
  if ( type == 2 ) {
    print( gridExtra::grid.arrange(
      plot[[ 3]], plot[[ 6]], plot[[ 9]], plot[[12]], plot[[16]], plot[[17]], plot[[19]], plot[[21]], 
      nrow = 4) )
  } else if ( type == 3 ) {
    print( gridExtra::grid.arrange(
      plot[[27]], plot[[31]], plot[[32]], plot[[33]], plot[[40]], plot[[45]], plot[[46]], 
      nrow = 4) )
  } else if ( type == 4 ) {
    print( gridExtra::grid.arrange(
      plot[[22]], plot[[47]], 
      nrow = 2) )
  } else if ( type == 5 ) {
    print( gridExtra::grid.arrange(
      plot[[ 2]], plot[[ 4]], plot[[23]], plot[[28]], plot[[29]], plot[[30]], plot[[37]], 
      nrow = 3) )
  } else if ( type == 6 ) {
    print( gridExtra::grid.arrange(
      plot[[48]], plot[[49]], 
      nrow = 1) )
  } else if ( type == 7 ) {
    print( gridExtra::grid.arrange(
      plot[[ 5]], plot[[11]], plot[[18]], plot[[35]], plot[[38]], plot[[39]], 
      nrow = 2) )
  } else if ( type == 8 ) {
    print(plot[[34]])
  }
  # save plot
  dev.off()
}

# # category
# c(7, 8, 10, 13, 14, 15, 20, 24, 25, 26, 36, 41, 42, 43, 44) # 15 cases
# c(3, 6, 9, 12, 16, 17, 19, 21) # 8 cases
# c(27, 31, 32, 33, 40, 45, 46) # 7 cases
# c(22, 47) # 2 cases
# c(2, 4, 23, 28, 29, 30, 37) # 7 cases
# c(48, 49) # 2 cases
# c(5, 11, 18, 35, 38, 39) # 6 cases
# c(34) # 1 case
# 
# plot[[ 3]], plot[[ 6]], plot[[ 9]], plot[[12]], plot[[16]], plot[[17]], plot[[19]], plot[[21]]
# plot[[27]], plot[[31]], plot[[32]], plot[[33]], plot[[40]], plot[[45]], plot[[46]]
# plot[[22]], plot[[47]]
# plot[[ 2]], plot[[ 4]], plot[[23]], plot[[28]], plot[[29]], plot[[30]], plot[[37]]
# plot[[48]], plot[[49]]
# plot[[ 5]], plot[[11]], plot[[18]], plot[[35]], plot[[38]], plot[[39]]
# plot[[34]]

# all cases
type <- 0
# save eps plot
setEPS()
name_file <- paste0("figure/14_imputation_", type, ".eps")
postscript(name_file, width = 8.5*2, height = 11*2)
# plot
print( gridExtra::grid.arrange(
  # ggplot()  , plot[[ 2]], plot[[ 3]], plot[[ 4]], plot[[ 5]], plot[[ 6]], ggplot()  , ggplot()  , plot[[ 9]], 
  # ggplot()  , plot[[11]], plot[[12]], ggplot()  , ggplot()  , ggplot()  , plot[[16]], plot[[17]], plot[[18]], plot[[19]], 
  # ggplot()  , plot[[21]], plot[[22]], plot[[23]], ggplot()  , ggplot()  , ggplot()  , plot[[27]], plot[[28]], plot[[29]], 
  # plot[[30]], plot[[31]], plot[[32]], plot[[33]], plot[[34]], plot[[35]], ggplot()  , plot[[37]], plot[[38]], plot[[39]], 
  # plot[[40]], ggplot()  , ggplot()  , ggplot()  , ggplot()  , plot[[45]], plot[[46]], plot[[47]], plot[[48]], plot[[49]],
  plot[[ 1]], plot[[ 2]], plot[[ 3]], plot[[ 4]], plot[[ 5]], plot[[ 6]], plot[[ 7]], plot[[ 8]], plot[[ 9]], plot[[10]], 
  plot[[11]], plot[[12]], plot[[13]], plot[[14]], plot[[15]], plot[[16]], plot[[17]], plot[[18]], plot[[19]], plot[[20]], 
  plot[[21]], plot[[22]], plot[[23]], plot[[24]], plot[[25]], plot[[26]], plot[[27]], plot[[28]], plot[[29]], plot[[30]], 
  plot[[31]], plot[[32]], plot[[33]], plot[[34]], plot[[35]], plot[[36]], plot[[37]], plot[[38]], plot[[39]], plot[[40]], 
  plot[[41]], plot[[42]], plot[[43]], plot[[44]], plot[[45]], plot[[46]], plot[[47]], plot[[48]], plot[[49]], 
  nrow = 10) )
# save plot
dev.off()