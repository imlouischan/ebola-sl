## parallel running on one machine #############################################
# https://www.r-bloggers.com/how-to-go-parallel-in-r-basics-tips/

# start to run
print(Sys.time())
timer_start = proc.time()

# number of local cores (as the number of chains)
no_cores = max(1, parallel::detectCores(logical = F))
# initiate cluster, type = "FORK" on Mac/Linux only
cl = parallel::makeCluster(no_cores, type = "FORK")
# check initiated cluster
print(cl)
# run
parallel_run = parallel::clusterApply(cl, 1:(1*no_cores), run_chains) 

# finish
parallel::stopCluster(cl)
print(Sys.time())
timer_stop = proc.time() - timer_start
print(timer_stop)