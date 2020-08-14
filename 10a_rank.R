## model comparison ############################################################

rank = expand.grid(
  flat.TRUE     = c(1, 0), 
  epsilonH.TRUE = c(1, 0), 
  delta.TRUE    = c(1, 0), 
  distr.r = c("geom", "pois", "nbinom"), 
  distr.s = c("dgamma", "dweibull"), 
  k        = NA,
  n        = NA,
  logMLE   = NA,
  AIC      = NA,
  AIC_diff = NA,
  BIC      = NA,
  BIC_diff = NA,
  ML       = NA,
  pM       = NA)

for (flat_TRUE in c(1, 0)) { # flat proposal / pij
  for (epsilonH_TRUE in c(1, 0)) { # independent / dependent epsilonH
    for (delta_TRUE in c(1, 0)) { # exponentially decreasing Rt / constant R0
      for (distr_r in c("geom", "pois", "nbinom")) { # alternative offspring distributions
        for (distr_s in c("dgamma", "dweibull")) { # alternative serial interval distributions
          
          # mcmc samples of an alternative model
          samples = subset(samples_M48, 
                           flat.TRUE       == flat_TRUE     &
                             epsilonH.TRUE == epsilonH_TRUE &
                             delta.TRUE    == delta_TRUE    &
                             distr.r == distr_r &
                             distr.s == distr_s)
          
          # number of parameters
          samples_index <- c("thetaSS1", "thetaSS2", "R0", "size", "delta", "epsilon", "epsilonH", "vi")
          if (distr_r != "nbinom") { samples_index = setdiff(samples_index, "size") }
          if (!delta_TRUE)         { samples_index = setdiff(samples_index, "delta") }
          if (!epsilonH_TRUE)      { samples_index = setdiff(samples_index, "epsilonH") }
          k = length(samples_index)
          # sample size
          n = nrow(Data)
          # MLE
          logMLE = max(samples$logL)
          # selection criteria, AIC
          AIC = 2*k - 2*logMLE
          # selection criteria, BIC
          BIC = log(n)*k - 2*logMLE
          
          # likelihood
          L = exp(samples$logL)
          # marginal likelihood (PHM estimator)
          ML = 1/sum(1/L) / length(L)
          
          # output
          rank[rank$flat.TRUE       == flat_TRUE     &
                 rank$epsilonH.TRUE == epsilonH_TRUE &
                 rank$delta.TRUE    == delta_TRUE    &
                 rank$distr.r == distr_r &
                 rank$distr.s == distr_s, 
               c("k", "n", "logMLE", "AIC", "BIC", "ML")] <- c(k, n, logMLE, AIC, BIC, ML)
          
        } # distr_s
      } # distr_r
    } # delta_TRUE
  } # epsilonH_TRUE
} # flat_TRUE

# subset
rank <- subset(rank, 
               flat.TRUE == 1 & 
                 epsilonH.TRUE == 0 & 
                 delta.TRUE == 1 & 
                 (distr.r == "geom" | distr.r == "pois" | distr.r == "nbinom") &
                 (distr.s == "dgamma" | distr.s == "dweibull") )

# AIC & BIC differences between models
rank$AIC_diff = rank$AIC - min(rank$AIC)
rank$BIC_diff = rank$BIC - min(rank$BIC)

# model posterior
rank$pM = rank$ML / sum(rank$ML) * 100

# round
rank$logMLE   = round(rank$logMLE  *10)/10
rank$AIC      = round(rank$AIC     *10)/10
rank$AIC_diff = round(rank$AIC_diff*10)/10
rank$BIC      = round(rank$BIC     *10)/10
rank$BIC_diff = round(rank$BIC_diff*10)/10
rank$pM       = round(rank$pM      *10)/10
# output
print(rank)

# \  flat.TRUE epsilonH.TRUE delta.TRUE distr.r  distr.s k  n logMLE   AIC AIC_diff   BIC BIC_diff           ML   pM
# 1          1             1          1    geom   dgamma 7 49 -182.3 378.6      5.0 391.9      6.9 9.380819e-94  0.0
# 2          0             1          1    geom   dgamma 7 49 -182.2 378.5      4.9 391.7      6.8 2.748573e-93  0.0
# 3          1             0          1    geom   dgamma 6 49 -183.1 378.2      4.6 389.5      4.6 2.649428e-94  0.0
# 4          0             0          1    geom   dgamma 6 49 -183.0 378.0      4.5 389.4      4.5 9.304625e-94  0.0
# 5          1             1          0    geom   dgamma 6 49 -185.2 382.5      8.9 393.8      8.9 1.392453e-94  0.0
# 6          0             1          0    geom   dgamma 6 49 -185.2 382.4      8.8 393.8      8.8 1.240621e-94  0.0
# 7          1             0          0    geom   dgamma 5 49 -185.4 380.7      7.1 390.2      5.2 3.841081e-95  0.0
# 8          0             0          0    geom   dgamma 5 49 -185.3 380.6      7.0 390.1      5.1 2.214217e-94  0.0
# 9          1             1          1    pois   dgamma 7 49 -187.3 388.6     15.0 401.8     16.9 4.345664e-95  0.0
# 10         0             1          1    pois   dgamma 7 49 -187.5 389.1     15.5 402.3     17.3 4.535290e-98  0.0
# 11         1             0          1    pois   dgamma 6 49 -188.7 389.3     15.7 400.7     15.7 1.626818e-96  0.0
# 12         0             0          1    pois   dgamma 6 49 -188.6 389.2     15.6 400.6     15.6 4.333216e-96  0.0
# 13         1             1          0    pois   dgamma 6 49 -192.3 396.5     22.9 407.9     22.9 1.237373e-97  0.0
# 14         0             1          0    pois   dgamma 6 49 -192.2 396.5     22.9 407.8     22.9 3.087676e-98  0.0
# 15         1             0          0    pois   dgamma 5 49 -192.4 394.8     21.2 404.3     19.3 1.181630e-97  0.0
# 16         0             0          0    pois   dgamma 5 49 -192.4 394.9     21.3 404.3     19.4 3.412702e-99  0.0
# 17         1             1          1  nbinom   dgamma 8 49 -182.2 380.4      6.8 395.6     10.6 8.725650e-95  0.0
# 18         0             1          1  nbinom   dgamma 8 49 -182.0 380.0      6.4 395.1     10.1 1.337303e-94  0.0
# 19         1             0          1  nbinom   dgamma 7 49 -182.4 378.8      5.2 392.0      7.1 2.179380e-93  0.0
# 20         0             0          1  nbinom   dgamma 7 49 -182.3 378.6      5.1 391.9      6.9 8.004342e-93  0.0
# 21         1             1          0  nbinom   dgamma 7 49 -184.0 381.9      8.3 395.2     10.2 1.174072e-94  0.0
# 22         0             1          0  nbinom   dgamma 7 49 -184.0 382.0      8.4 395.2     10.2 1.753628e-96  0.0
# 23         1             0          0  nbinom   dgamma 6 49 -184.0 379.9      6.3 391.3      6.3 4.344265e-94  0.0
# 24         0             0          0  nbinom   dgamma 6 49 -184.0 380.0      6.4 391.4      6.4 1.267347e-93  0.0
# 25         1             1          1    geom dweibull 7 49 -180.0 374.0      0.4 387.2      2.3 5.293435e-89 15.5
# 26         0             1          1    geom dweibull 7 49 -180.1 374.2      0.6 387.4      2.5 7.465865e-89 21.8
# 27         1             0          1    geom dweibull 6 49 -180.8 373.6      0.0 385.0      0.0 2.443391e-89  7.1
# 28         0             0          1    geom dweibull 6 49 -180.8 373.6      0.0 384.9      0.0 4.496667e-89 13.2
# 29         1             1          0    geom dweibull 6 49 -182.8 377.7      4.1 389.0      4.1 8.617506e-90  2.5
# 30         0             1          0    geom dweibull 6 49 -182.8 377.5      3.9 388.9      3.9 1.703278e-90  0.5
# 31         1             0          0    geom dweibull 5 49 -183.0 376.0      2.5 385.5      0.6 4.474077e-90  1.3
# 32         0             0          0    geom dweibull 5 49 -183.0 376.0      2.4 385.4      0.5 1.983930e-89  5.8
# 33         1             1          1    pois dweibull 7 49 -185.2 384.3     10.7 397.6     12.6 3.097612e-91  0.1
# 34         0             1          1    pois dweibull 7 49 -185.1 384.1     10.5 397.4     12.4 5.552658e-91  0.2
# 35         1             0          1    pois dweibull 6 49 -186.3 384.7     11.1 396.0     11.1 3.346501e-93  0.0
# 36         0             0          1    pois dweibull 6 49 -186.5 384.9     11.3 396.3     11.3 8.356041e-92  0.0
# 37         1             1          0    pois dweibull 6 49 -190.0 391.9     18.3 403.3     18.3 1.041746e-92  0.0
# 38         0             1          0    pois dweibull 6 49 -189.9 391.9     18.3 403.2     18.3 1.031944e-92  0.0
# 39         1             0          0    pois dweibull 5 49 -190.1 390.1     16.5 399.6     14.6 1.049775e-92  0.0
# 40         0             0          0    pois dweibull 5 49 -190.1 390.2     16.6 399.7     14.7 1.517759e-92  0.0
# 41         1             1          1  nbinom dweibull 8 49 -179.8 375.6      2.0 390.7      5.7 4.549119e-89 13.3
# 42         0             1          1  nbinom dweibull 8 49 -179.6 375.1      1.6 390.3      5.3 3.504416e-89 10.3
# 43         1             0          1  nbinom dweibull 7 49 -180.2 374.5      0.9 387.7      2.8 5.152158e-90  1.5
# 44         0             0          1  nbinom dweibull 7 49 -180.1 374.3      0.7 387.5      2.6 9.078661e-90  2.7
# 45         1             1          0  nbinom dweibull 7 49 -181.8 377.5      3.9 390.8      5.8 5.088339e-90  1.5
# 46         0             1          0  nbinom dweibull 7 49 -181.5 377.1      3.5 390.3      5.4 1.045690e-90  0.3
# 47         1             0          0  nbinom dweibull 6 49 -181.8 375.5      1.9 386.9      1.9 6.067915e-90  1.8
# 48         0             0          0  nbinom dweibull 6 49 -181.7 375.4      1.8 386.8      1.8 2.249035e-90  0.7

# # subset
# rank_sub <- subset(rank, 
#                    (distr.r == "geom" | distr.r == "nbinom") & 
#                      distr.s == "dweibull" & 
#                      delta.TRUE == 1)
# # output
# print(rank_sub)

# \  flat.TRUE epsilonH.TRUE delta.TRUE distr.r  distr.s k  n logMLE   AIC AIC_diff   BIC BIC_diff           ML   pM
# 25         1             1          1    geom dweibull 7 49 -180.0 374.0  (3) 0.4 387.2 (3)  2.3 5.293435e-89 15.5 (2)
# 26         0             1          1    geom dweibull 7 49 -180.1 374.2  (4) 0.6 387.4 (4)  2.5 7.465865e-89 21.8 (1)
# 27         1             0          1    geom dweibull 6 49 -180.8 373.6  (1) 0.0 385.0 (1)  0.0 2.443391e-89  7.1 (6)
# 28         0             0          1    geom dweibull 6 49 -180.8 373.6  (1) 0.0 384.9 (1)  0.0 4.496667e-89 13.2 (4)
# 41         1             1          1  nbinom dweibull 8 49 -179.8 375.6  (8) 2.0 390.7 (8)  5.7 4.549119e-89 13.3 (3)
# 42         0             1          1  nbinom dweibull 8 49 -179.6 375.1  (7) 1.6 390.3 (7)  5.3 3.504416e-89 10.3 (5)
# 43         1             0          1  nbinom dweibull 7 49 -180.2 374.5  (6) 0.9 387.7 (6)  2.8 5.152158e-90  1.5 (8)
# 44         0             0          1  nbinom dweibull 7 49 -180.1 374.3  (5) 0.7 387.5 (5)  2.6 9.078661e-90  2.7 (7)