## data set ####################################################################

# read in the file as a data frame
Data = gdata::read.xls("ebola_sl.xlsx")
# delete columns with NA
Data = Data[, colSums(is.na(Data)) != nrow(Data)]
# rename columns
colnames(Data) <- c("i", "vi", "tiE", "tiS", "tiH", "tiD", "tiR")

# convert into date from text
Data$tiE = as.Date(Data$tiE, "%Y-%m-%d")
Data$tiS = as.Date(Data$tiS, "%Y-%m-%d")
Data$tiH = as.Date(Data$tiH, "%Y-%m-%d")
Data$tiD = as.Date(Data$tiD, "%Y-%m-%d")
Data$tiR = as.Date(Data$tiR, "%Y-%m-%d")
# convert into numeric from date
Data$tiE = as.numeric(Data$tiE - as.Date("2014-07-07"), units="days")
Data$tiS = as.numeric(Data$tiS - as.Date("2014-07-07"), units="days")
Data$tiH = as.numeric(Data$tiH - as.Date("2014-07-07"), units="days")
Data$tiD = as.numeric(Data$tiD - as.Date("2014-07-07"), units="days")
Data$tiR = as.numeric(Data$tiR - as.Date("2014-07-07"), units="days")
# the end of the outbreak on 10 Jan 2015 (=Day 187)
Data_tEND = as.numeric(as.Date("2015-01-10") - as.Date("2014-07-07"), units="days")

# existence of time points
Data_obs = matrix(!is.na(c(Data$tiE, Data$tiS, Data$tiH, Data$tiD, Data$tiR)), ncol=5 )

# convert into numeric from factor
Data$vi = suppressWarnings( as.integer(as.character(Data$vi)) )
# Case 62 is either infected by Case 58 or 59
# assuming infected by Case 58 does not affect the results
Data$vi[Data$i == 62] = 58
# now: Data$vi[Data$i == 41] = NA; and is selected from the candidate list during MCMC

## hospitalization #############################################################
# note that Case 1 has not be isolated

## observed time lags ##########################################################

# incubation periods
Data$tauiES = Data$tiS - Data$tiE

# other time lags
Data$tauiSH = Data$tiH - Data$tiS
Data$tauiSD = Data$tiD - Data$tiS
Data$tauiSR = Data$tiR - Data$tiS
Data$tauiHD = Data$tiD - Data$tiH
Data$tauiHR = Data$tiR - Data$tiH

# serial intervals
vi = rep(NA ,nrow(Data))
for (i in 1:length(Data$i)) {
  vi_temp = which(Data$i == Data$vi[i])
  if (length(vi_temp) != 0) {
    vi[i] = vi_temp
  }
}
Data$tauiSS = Data$tiS - Data$tiS[vi]

# take non-NA values
tauiES = Data$tauiES[!is.na(Data$tauiES)]
tauiSH = Data$tauiSH[!is.na(Data$tauiSH)]
tauiSD = Data$tauiSD[!is.na(Data$tauiSD)]
tauiSR = Data$tauiSR[!is.na(Data$tauiSR)]
tauiHD = Data$tauiHD[!is.na(Data$tauiHD)]
tauiHR = Data$tauiHR[!is.na(Data$tauiHR)]
tauiSS = Data$tauiSS[!is.na(Data$tauiSS)]