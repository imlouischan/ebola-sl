## combine chains ##############################################################

function (N, chains) { 
  
  # preallocate
  samples_all = data.frame(NULL)
  
  for (chain in chains) {
    
    # load samples
    file_samples = paste0("samples", 
                          "_N", N, 
                          "_f", flat_TRUE, 
                          "_e", epsilonH_TRUE, 
                          "_d", delta_TRUE, 
                          "_R", distr_r, 
                          "_S", distr_s, 
                          "_C", chain, 
                          ".Rdata")
    load(file = file_samples)
    
    # alternative model
    samples = cbind(
      flat.TRUE     = flat_TRUE, 
      epsilonH.TRUE = epsilonH_TRUE, 
      delta.TRUE    = delta_TRUE, 
      distr.r       = distr_r,
      distr.s       = distr_s,
      samples)
    
    # factor instead of numeric
    samples$vi            = as.factor(samples$vi)
    samples$flat.TRUE     = as.factor(samples$flat.TRUE)
    samples$epsilonH.TRUE = as.factor(samples$epsilonH.TRUE)
    samples$delta.TRUE    = as.factor(samples$delta.TRUE)
    
    # combine chains
    samples_all = rbind(samples_all, samples)
    
  }
  # output samples of all chains
  return(samples_all) 
}