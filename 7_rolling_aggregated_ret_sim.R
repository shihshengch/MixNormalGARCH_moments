# 7_rolling_aggregated_ret_sim.R

# clear
rm(list = ls())

# specify
k = 2
n = 5
p = 0.025

# load MN(k)-GRACH(1,1) n-step rolling parameter estimation results
load(paste0("MNG_",k,"_roll_par_",n,"_step","_data.RData"))

# load sources
source("1_source_functions.R")

# function of simulating rolling aggregated returns
rolling.sim.R_tn = function(y) {
  r_tn = replicate(100, r_tn_sim(n, as.vector(rolling.param[y,]), 
                                 r.mat[,y], S = 1000))
  dim(r_tn) = c(n, 100000)
  R_tn = colSums(r_tn)
}

# parallel computing
no_cores = detectCores()-1
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = ls(), envir = environment())
clusterEvalQ(cl = cl, source("1_source_functions.R"))

R.mat_tn = pbsapply(1:nrow(rolling.param), 
                    function(y) rolling.sim.R_tn(y), cl = cl)

stopCluster(cl)

save(list = c("R.mat_tn", "k", "n"),
     file = paste0("MNG_",k,"_R_",n,"_step_simulation","_data.RData"))