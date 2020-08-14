## likelihood of serial interval ###############################################

function (i, thetaSS, R0, epsilon, Data) {
  
  # preallocate
  L_si_1 = 0
  L_si_2 = 0
  
  # infector of i
  vi = which(Data$i == Data$vi[i])
  
  # categories in L_s
  # vi: "xOOx", "xOMx", "OMOx", "MMOx", "OMMO(D)", "OMMM(D)", "MMMO(D)", "MMMO(R)"
  #  i: "xOxx", "OMxx", "MMOx", "MMMO(D)", "MMMO(R)"
  if (length(vi) == 0) { # the index case / imported cases
    type.vi = "(0)"
    type.i = "(0)"
    # numerator and denominator
    L_si_1 = 1
    L_si_2 = 1
  } else {
    # hospitalization probability
    ph = 1
    # if ( Data_obs[vi,3] | Data_obs[vi,5] ) { ph = 1
    # } else                                 { ph = 18/22 }
    # except the index case / imported cases
    if (                        Data_obs[vi,2] &  Data_obs[vi,3]                  ) { # vi: xOOx (1) ====
      type.vi = "(1) xOOx   "
      if (                       Data_obs[i,2]                                  ) { # i: xOxx (I) ----
        type.i = "(I)    xOxx   "
        # data of observed cases
        tauSS = Data$tiS[i] - Data$tiS[vi]
        tauSH = Data$tiH[vi] - Data$tiS[vi]
        # L_s of each case
        L_si_1 = pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
        L_si_2 = 1
      }
      else if ( Data_obs[i,1] & !Data_obs[i,2]                                  ) { # i: OMxx (II) ----
        type.i = "(II)   OMxx   "
        # intervals
        a = 0
        b = Data_BS2(i) - Data$tiE[i]
        # change of variable
        tauES = a:b
        tauSS = Data$tiE[i] + tauES - Data$tiS[vi]
        tauSH = Data$tiH[vi] - Data$tiS[vi]
        # numerator and denominator
        L_si_1 = sum( pmf_ES(tauES)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph) )
        L_si_2 = sum( pmf_ES(tauES) )
      }
      # else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]                 ) { # i: MMOx (III)
      #   type.i = "(III)  MMOx   "
      # }
      else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: MMMO(D) (IV,D) ----
        type.i = "(IV,D) MMMO(D)"
        # intervals
        a = Data$tiD[i] - Data_BS2(i)
        b = Data$tiD[i] - Data$tiS[vi]
        # change of variable
        tauSD = a:b
        tauSS = Data$tiD[i] - tauSD - Data$tiS[vi]
        tauSH = Data$tiH[vi] - Data$tiS[vi]
        # numerator and denominator
        L_si_1 = sum( pmf_SD(tauSD)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph) )
        L_si_2 = sum( pmf_SD(tauSD) )
      }
      else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,5]) { # i: MMMO(R) (IV,R) ----
        type.i = "(IV,R) MMMO(R)"
        # intervals
        a = Data$tiR[i] - Data_BS2(i)
        b = Data$tiR[i] - Data$tiS[vi]
        # change of variable
        tauSR = a:b
        tauSS = Data$tiR[i] - tauSR - Data$tiS[vi]
        tauSH = Data$tiH[vi] - Data$tiS[vi]
        # numerator and denominator
        L_si_1 = sum( pmf_SR(tauSR)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph) )
        L_si_2 = sum( pmf_SR(tauSR) )
      }
      else { print(i); break }
    }
    else if (                   Data_obs[vi,2] & !Data_obs[vi,3]                  ) { # vi: xOMx (2) # unknown hospitalization, except Case 1 ====
      type.vi = "(2) xOMx   "
      if (                       Data_obs[i,2]                                  ) { # i: xOxx (I) ----
        type.i = "(I)    xOxx   "
        # intervals
        a = 0
        b = Data_BH2(vi) - Data$tiS[vi]
        # change of variable
        tauSS = Data$tiS[i] - Data$tiS[vi]
        tauSH = a:b
        # numerator and denominator
        if (vi == 1) {
          L_si_1 = pmf_s(tauSS, thetaSS)
          L_si_2 = 1
        }
        else         {
          L_si_1 = sum( pmf_SH(tauSH)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph) )
          L_si_2 = sum( pmf_SH(tauSH) )
        }
      }
      else if ( Data_obs[i,1] & !Data_obs[i,2]                                  ) { # i: OMxx (II) ----
        type.i = "(II)   OMxx   "
        if (vi == 1) {
          # intervals
          a = 0
          b = Data_BS2(i) - Data$tiE[i]
          # change of variable
          tauES = a:b
          tauSS = Data$tiE[i] + tauES - Data$tiS[vi]
          # numerator and denominator
          L_si_1 = sum( pmf_ES(tauES)*pmf_s(tauSS, thetaSS) )
          L_si_2 = sum( pmf_ES(tauES) )
        }
        else         {
          # intervals
          for (tauSH in 0:(Data_BH2(vi) - Data$tiS[vi])) {
            for (tauES in 0:(Data_BS2(i) - Data$tiE[i])) {
              # change of variable
              tauSS = Data$tiE[i] + tauES - Data$tiS[vi]
              # numerator and denominator
              L_si_1 = L_si_1 + pmf_ES(tauES)*pmf_SH(tauSH)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
              L_si_2 = L_si_2 + pmf_ES(tauES)*pmf_SH(tauSH)
            }
          }
        }
      }
      else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]                 ) { # i: MMOx (III) ----
        type.i = "(III)  MMOx   "
        if (vi == 1) {
          # intervals
          a = Data$tiH[i] - Data_BS2(i)
          b = Data$tiH[i] - Data$tiS[vi]
          # change of variable
          tauiSH = a:b
          tauSS = Data$tiH[i] - tauiSH - Data$tiS[vi]
          # numerator and denominator
          L_si_1 = sum( pmf_SH(tauiSH)*pmf_s(tauSS, thetaSS) )
          L_si_2 = sum( pmf_SH(tauiSH) )
        }
        else         {
          # intervals
          for (tauSH in 0:(Data_BH2(vi) - Data$tiS[vi])) {
            for (tauiSH in (Data$tiH[i] - Data_BS2(i)):(Data$tiH[i] - Data$tiS[vi])) {
              # change of variable
              tauSS = Data$tiH[i] - tauiSH - Data$tiS[vi]
              # numerator and denominator
              L_si_1 = L_si_1 + pmf_SH(tauiSH)*pmf_SH(tauSH)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
              L_si_2 = L_si_2 + pmf_SH(tauiSH)*pmf_SH(tauSH)
            }
          }
        }
      }
      else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: MMMO(D) (IV,D) ----
        type.i = "(IV,D) MMMO(D)"
        if (vi == 1) {
          # intervals
          a = Data$tiD[i] - Data_BS2(i)
          b = Data$tiD[i] - Data$tiS[vi]
          # change of variable
          tauSD = a:b
          tauSS = Data$tiD[i] - tauSD - Data$tiS[vi]
          # numerator and denominator
          L_si_1 = sum( pmf_SD(tauSD)*pmf_s(tauSS, thetaSS) )
          L_si_2 = sum( pmf_SD(tauSD) )
        }
        else         {
          # intervals
          for (tauSH in 0:(Data_BH2(vi) - Data$tiS[vi])) {
            for (tauSD in (Data$tiD[i] - Data_BS2(i)):(Data$tiD[i] - Data$tiS[vi])) {
              # change of variable
              tauSS = Data$tiD[i] - tauSD - Data$tiS[vi]
              # numerator and denominator
              L_si_1 = L_si_1 + pmf_SD(tauSD)*pmf_SH(tauSH)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
              L_si_2 = L_si_2 + pmf_SD(tauSD)*pmf_SH(tauSH)
            }
          }
        }
      }
      # else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,5]) { # i: MMMO(R) (IV,R)
      #   type.i = "(IV,R) MMMO(R)"
      # }
      else { print(i); break }
    }
    else if ( Data_obs[vi,1] & !Data_obs[vi,2] &  Data_obs[vi,3]                  ) { # vi: OMOx (3) ====
      type.vi = "(3) OMOx   "
      # if (                       Data_obs[i,2]                                  ) { # i: xOxx (I)
      #   type.i = "(I)    xOxx   "
      # }
      # else if ( Data_obs[i,1] & !Data_obs[i,2]                                  ) { # i: OMxx (II)
      #   type.i = "(II)   OMxx   "
      # }
      # else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]                 ) { # i: MMOx (III)
      #   type.i = "(III)  MMOx   "
      # }
      if (     !Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: MMMO(D) (IV,D) ----
        type.i = "(IV,D) MMMO(D)"
        # intervals
        for (tauES in 0:(Data_BS2(vi) - Data$tiE[vi])) {
          for (tauSD in (Data$tiD[i] - Data_BS2(i)):(Data$tiD[i] - Data$tiE[vi] - tauES)) {
            # change of variable
            tauSS = Data$tiD[i] - tauSD - Data$tiE[vi] - tauES
            tauSH = Data$tiH[vi] - Data$tiE[vi] - tauES
            # numerator and denominator
            L_si_1 = L_si_1 + pmf_SD(tauSD)*pmf_ES(tauES)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
            L_si_2 = L_si_2 + pmf_SD(tauSD)*pmf_ES(tauES)
          }
        }
      }
      # else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,5]) { # i: MMMO(R) (IV,R)
      #   type.i = "(IV,R) MMMO(R)"
      # }
      else { print(i); break }
    }
    else if (!Data_obs[vi,1] & !Data_obs[vi,2] &  Data_obs[vi,3]                  ) { # vi: MMOx (4) ====
      type.vi = "(4) MMOx   "
      if (                       Data_obs[i,2]                                  ) { # i: xOxx (I) ----
        type.i = "(I)    xOxx   "
        # intervals
        a = Data$tiH[vi] - Data_BS2(vi)
        b = Data$tiH[vi] - Data_BS1(vi)
        # change of variable
        tauSS = Data$tiS[i] - Data$tiH[vi] + (a:b)
        tauSH = a:b
        # numerator and denominator
        L_si_1 = sum( pmf_SH(tauSH)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph) )
        L_si_2 = sum( pmf_SH(tauSH) )
      }
      else if ( Data_obs[i,1] & !Data_obs[i,2]                                  ) { # i: OMxx (II) ----
        type.i = "(II)   OMxx   "
        # intervals
        for (tauSH in (Data$tiH[vi] - Data_BS1(vi)):(Data$tiH[vi] - Data_BS2(vi))) {
          for (tauES in 0:(Data_BS2(i) - Data$tiE[i])) {
            # change of variable
            tauSS = Data$tiE[i] + tauES - Data$tiH[vi] + tauSH
            # numerator and denominator
            L_si_1 = L_si_1 + pmf_ES(tauES)*pmf_SH(tauSH)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
            L_si_2 = L_si_2 + pmf_ES(tauES)*pmf_SH(tauSH)
          }
        }
      }
      # else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]                 ) { # i: MMOx (III)
      #   type.i = "(III)  MMOx   "
      # }
      else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: MMMO(D) (IV,D) ----
        type.i = "(IV,D) MMMO(D)"
        # intervals
        for (tauSH in (Data$tiH[vi] - Data_BS1(vi)):(Data$tiH[vi] - Data_BS2(vi))) {
          for (tauSD in (Data$tiD[i] - Data_BS2(i)):(Data$tiD[i] - Data$tiH[vi] - tauSH)) {
            # change of variable
            tauSS = Data$tiD[i] + tauSD - Data$tiH[vi] + tauSH
            # numerator and denominator
            L_si_1 = L_si_1 + pmf_SD(tauSD)*pmf_SH(tauSH)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
            L_si_2 = L_si_2 + pmf_SD(tauSD)*pmf_SH(tauSH)
          }
        }
      }
      # else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,5]) { # i: MMMO(R) (IV,R)
      #   type.i = "(IV,R) MMMO(R)"
      # }
      else { print(i); break }
    }
    else if ( Data_obs[vi,1] & !Data_obs[vi,2] & !Data_obs[vi,3] &  Data_obs[vi,4]) { # vi: OMMO(D) (5D) # unknown hospitalization ====
      type.vi = "(5) OMMO(D)"
      if (                       Data_obs[i,2]                                  ) { # i: xOxx (I) ----
        type.i = "(I)    xOxx   "
        # intervals
        for (tauES in 0:(Data_BS2(vi) - Data$tiE[vi])) {
          for (tauHD in 0:(Data$tiD[vi] - Data$tiE[vi] - tauES)) {
            # change of variable
            tauSS = Data$tiS[i] - Data$tiE[vi] - tauES
            tauSH = Data$tiD[vi] - tauHD - Data$tiE[vi] - tauES
            # numerator and denominator
            L_si_1 = L_si_1 + pmf_HD(tauHD)*pmf_ES(tauES)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
            L_si_2 = L_si_2 + pmf_HD(tauHD)*pmf_ES(tauES)
          }
        }
      }
      # else if ( Data_obs[i,1] & !Data_obs[i,2]                                  ) { # i: OMxx (II)
      #   type.i = "(II)   OMxx   "
      # }
      # else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]                 ) { # i: MMOx (III)
      #   type.i = "(III)  MMOx   "
      # }
      else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: MMMO(D) (IV,D) [3D] ----
        type.i = "(IV,D) MMMO(D)"
        # intervals
        for (tauES in 0:(Data_BS2(vi) - Data$tiE[vi])) {
          for (tauHD in 0:(Data$tiD[vi] - Data$tiE[vi] - tauES)) {
            for (tauSD in (Data$tiD[i] - Data_BS2(i)):(Data$tiD[i] - Data$tiE[vi] - tauES)) {
              # change of variable
              tauSS = Data$tiD[i] - tauSD - Data$tiE[vi] - tauES
              tauSH = Data$tiD[vi] - tauHD - Data$tiE[vi] - tauES
              # numerator and denominator
              L_si_1 = L_si_1 + pmf_SD(tauSD)*pmf_HD(tauHD)*pmf_ES(tauES)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
              L_si_2 = L_si_2 + pmf_SD(tauSD)*pmf_HD(tauHD)*pmf_ES(tauES)
            }
          }
        }
      }
      # else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,5]) { # i: MMMO(R) (IV,R)
      #   type.i = "(IV,R) MMMO(R)"
      # }
      else { print(i); break }
    }
    else if ( Data_obs[vi,1] & !Data_obs[vi,2] & !Data_obs[vi,3] & !Data_obs[vi,4]) { # vi: OMMM(D) (6D) # unknown hospitalization ====
      type.vi = "(6) OMMM(D)"
      # if (                       Data_obs[i,2]                                  ) { # i: xOxx (I)
      #   type.i = "(I)    xOxx   "
      # }
      # else if ( Data_obs[i,1] & !Data_obs[i,2]                                  ) { # i: OMxx (II)
      #   type.i = "(II)   OMxx   "
      # }
      # else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]                 ) { # i: MMOx (III)
      #   type.i = "(III)  MMOx   "
      # }
      if (     !Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: MMMO(D) (IV,D) [3D] ----
        type.i = "(IV,D) MMMO(D)"
        # intervals
        for (tauES in 0:(Data_BS2(vi) - Data$tiE[vi])) {
          for (tauSH in 0:(Data_tEND - Data$tiE[vi] - tauES)) {
            for (tauSD in (Data$tiD[i] - Data_BS2(i)):(Data$tiD[i] - Data$tiE[vi] - tauES)) {
              # change of variable
              tauSS = Data$tiD[i] - tauSD - Data$tiE[vi] - tauES
              # numerator and denominator
              L_si_1 = L_si_1 + pmf_SD(tauSD)*pmf_SH(tauSH)*pmf_ES(tauES)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
              L_si_2 = L_si_2 + pmf_SD(tauSD)*pmf_SH(tauSH)*pmf_ES(tauES)
            }
          }
        }
      }
      # else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,5]) { # i: MMMO(R) (IV,R)
      #   type.i = "(IV,R) MMMO(R)"
      # }
      else { print(i); break }
    }
    else if (!Data_obs[vi,1] & !Data_obs[vi,2] & !Data_obs[vi,3] &  Data_obs[vi,4]) { # vi: MMMO(D) (7D) # unknown hospitalization ====
      type.vi = "(7) MMMO(D)"
      if (                       Data_obs[i,2]                                  ) { # i: xOxx (I) ----
        type.i = "(I)    xOxx   "
        # intervals, 1st part
        for (tauHD in 0:(Data$tiD[vi] - Data_BS2(vi))) {
          for (tauSH in (Data$tiD[vi] - Data_BS2(vi) - tauHD):(Data$tiD[vi] - Data_BS1(vi) - tauHD)) {
            # change of variable
            tauSS = Data$tiS[i] - Data$tiD[vi] + tauHD + tauSH
            # numerator and denominator
            L_si_1 = L_si_1 + pmf_SH(tauSH)*pmf_HD(tauHD)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
            L_si_2 = L_si_2 + pmf_SH(tauSH)*pmf_HD(tauHD)
          }
        }
        # intervals, 2nd part
        for (tauHD in (Data$tiD[vi] - Data_BS2(vi) + 1):(Data$tiD[vi] - Data_BS1(vi))) {
          for (tauSH in 0:(Data$tiD[vi] - Data_BS1(vi) - tauHD)) {
            # change of variable
            tauSS = Data$tiS[i] - Data$tiD[vi] + tauHD + tauSH
            # numerator and denominator
            L_si_1 = L_si_1 + pmf_SH(tauSH)*pmf_HD(tauHD)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
            L_si_2 = L_si_2 + pmf_SH(tauSH)*pmf_HD(tauHD)
          }
        }
      }
      # else if ( Data_obs[i,1] & !Data_obs[i,2]                                  ) { # i: OMxx (II)
      #   type.i = "(II)   OMxx   "
      # }
      # else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]                 ) { # i: MMOx (III)
      #   type.i = "(III)  MMOx   "
      # }
      else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: MMMO(D) (IV,D) [3D] ----
        type.i = "(IV,D) MMMO(D)"
        # intervals, 1st part
        for (tauHD in 0:(Data$tiD[vi] - Data_BS2(vi))) {
          for (tauSH in (Data$tiD[vi] - Data_BS2(vi) - tauHD):(Data$tiD[vi] - Data_BS1(vi) - tauHD)) {
            for (tauSD in (Data$tiD[i] - Data_BS2(i)):(Data$tiD[i] - Data$tiD[vi] + tauHD + tauSH)) {
              # change of variable
              tauSS = Data$tiD[i] - tauSD - Data$tiD[vi] + tauHD + tauSH
              # numerator and denominator
              L_si_1 = L_si_1 + pmf_SD(tauSD)*pmf_SH(tauSH)*pmf_HD(tauHD)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
              L_si_2 = L_si_2 + pmf_SD(tauSD)*pmf_SH(tauSH)*pmf_HD(tauHD)
            }
          }
        }
        # intervals, 2nd part
        for (tauHD in (Data$tiD[vi] - Data_BS2(vi) + 1):(Data$tiD[vi] - Data_BS1(vi))) {
          for (tauSH in 0:(Data$tiD[vi] - Data_BS1(vi) - tauHD)) {
            for (tauSD in (Data$tiD[i] - Data_BS2(i)):(Data$tiD[i] - Data$tiD[vi] + tauHD + tauSH)) {
              # change of variable
              tauSS = Data$tiD[i] - tauSD - Data$tiD[vi] + tauHD + tauSH
              # numerator and denominator
              L_si_1 = L_si_1 + pmf_SD(tauSD)*pmf_SH(tauSH)*pmf_HD(tauHD)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
              L_si_2 = L_si_2 + pmf_SD(tauSD)*pmf_SH(tauSH)*pmf_HD(tauHD)
            }
          }
        }
      }
      # else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,5]) { # i: MMMO(R) (IV,R)
      #   type.i = "(IV,R) MMMO(R)"
      # }
      else { print(i); break }
    }
    else if (!Data_obs[vi,1] & !Data_obs[vi,2] & !Data_obs[vi,3] &  Data_obs[vi,5]) { # vi: MMMO(R) (7R) # unknown hospitalization ====
      type.vi = "(7) MMMO(R)"
      # if (                       Data_obs[i,2]                                  ) { # i: xOxx (I)
      #   type.i = "(I)    xOxx   "
      # }
      # else if ( Data_obs[i,1] & !Data_obs[i,2]                                  ) { # i: OMxx (II)
      #   type.i = "(II)   OMxx   "
      # }
      # else if (!Data_obs[i,1] & !Data_obs[i,2] &  Data_obs[i,3]                 ) { # i: MMOx (III)
      #   type.i = "(III)  MMOx   "
      # }
      if (     !Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,4]) { # i: MMMO(D) (IV,D) [3D] ----
        type.i = "(IV,D) MMMO(D)"
        # intervals, 1st part
        for (tauHR in 0:(Data$tiR[vi] - Data_BS2(vi))) {
          for (tauSH in (Data$tiR[vi] - Data_BS2(vi) - tauHR):(Data$tiR[vi] - Data_BS1(vi) - tauHR)) {
            for (tauSD in (Data$tiD[i] - Data_BS2(i)):(Data$tiD[i] - Data$tiR[vi] + tauHR + tauSH)) {
              # change of variable
              tauSS = Data$tiD[i] - tauSD - Data$tiR[vi] + tauHR + tauSH
              # numerator and denominator
              L_si_1 = L_si_1 + pmf_SD(tauSD)*pmf_SH(tauSH)*pmf_HR(tauHR)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
              L_si_2 = L_si_2 + pmf_SD(tauSD)*pmf_SH(tauSH)*pmf_HR(tauHR)
            }
          }
        }
        # intervals, 2nd part
        for (tauHR in (Data$tiR[vi] - Data_BS2(vi) + 1):(Data$tiR[vi] - Data_BS1(vi))) {
          for (tauSH in 0:(Data$tiR[vi] - Data_BS1(vi) - tauHR)) {
            for (tauSD in (Data$tiD[i] - Data_BS2(i)):(Data$tiD[i] - Data$tiR[vi] + tauHR + tauSH)) {
              # change of variable
              tauSS = Data$tiD[i] - tauSD - Data$tiR[vi] + tauHR + tauSH
              # numerator and denominator
              L_si_1 = L_si_1 + pmf_SD(tauSD)*pmf_SH(tauSH)*pmf_HR(tauHR)*pmf_s_hat(tauSS, tauSH, thetaSS, epsilon, ph)
              L_si_2 = L_si_2 + pmf_SD(tauSD)*pmf_SH(tauSH)*pmf_HR(tauHR)
            }
          }
        }
      }
      # else if (!Data_obs[i,1] & !Data_obs[i,2] & !Data_obs[i,3] &  Data_obs[i,5]) { # i: MMMO(R) (IV,R)
      #   type.i = "(IV,R) MMMO(R)"
      # }
      else { print(i); break }
    }
    else { print(i); break }
  }
  
  # resulting likelihood
  L_si = L_si_1/L_si_2
  
  # output each case
  return(L_si)
  
}