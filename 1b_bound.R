## boundaries of missing time points ###########################################

# upper bound of tiH, BH2
Data_BH2 = function(i) min(Data$tiD[i], Data$tiR[i], Data_tEND, na.rm = T)

# lower bound of tiS, BS1
Data_BS1 = function(i) {
  BS1 = max(0, Data$tiE[i], na.rm = T)
  while(length(i) > 0) { # all ancestors
    i = which(Data$i == Data$vi[i]) # the only infector of i
    BS1 = max(BS1, Data$tiE[i], Data$tiS[i], na.rm = T)
  }
  return(BS1)
}

# upper bound of tiS, BS2
Data_BS2 = function(i) {
  BS2 = min(Data$tiH[i], Data$tiD[i], Data$tiR[i], Data_tEND, na.rm = T)
  while(length(i) > 0) { # all descestors
    i = which(Data$vi %in% Data$i[i]) # all infected by i's
    BS2 = min(BS2, Data$tiE[i], Data$tiS[i], Data$tiH[i], Data$tiD[i], Data$tiR[i], na.rm = T)
  }
  return(BS2)
}