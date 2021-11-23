# 2_MNG_k_estimation.R

# clear
rm(list = ls())

# load data
load("data.RData")

# load sources
source("1_source_functions.R")

# parallel computing
no_cores = detectCores()-1
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = ls(), envir = environment())
clusterEvalQ(cl = cl, source("1_source_functions.R"))

k = 2
result = MNG.k.est(r, k, cl)
list2env(result, envir = environment())
save(list = c("eps", "r", "c", "k", 
              "start.par", "bounds", "opt.result", "opt.par"),
     file = paste0("MNG_",k,"_data.RData"))

k = 3
result = MNG.k.est(r, k, cl)
list2env(result, envir = environment())
save(list = c("eps", "r", "c", "k", 
              "start.par", "bounds", "opt.result", "opt.par"),
     file = paste0("MNG_",k,"_data.RData"))

k = 4
result = MNG.k.est(r, k, cl)
list2env(result, envir = environment())
save(list = c("eps", "r", "c", "k", 
              "start.par", "bounds", "opt.result", "opt.par"),
     file = paste0("MNG_",k,"_data.RData"))

stopCluster(cl)