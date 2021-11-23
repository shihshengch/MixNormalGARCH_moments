# 6_rolling_parameters.R

# clear
rm(list = ls())

# specify p
k = 2
n = 10
p = 0.025

# load MN(k)-GRACH(1,1) parameter estimation results
load(paste0("MNG_",k,"_data.RData"))

# load sources
source("1_source_functions.R")

# rolling return matrix
initial = "1990/1993"
r.mat = runner(as.vector(r), k = length(r[initial]),
               f = function(x) x,
               at = seq(length(r[initial]),length(r),n))

# parallel computing
no_cores = detectCores()-1
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = ls(), envir = environment())
clusterEvalQ(cl = cl, source("1_source_functions.R"))

rolling.par = pbapply(r.mat, 2,
                      function(x) MNG.k.est(x, k, cl = NULL, TRUE),
                      cl = cl)

stopCluster(cl)

rolling.param = xts(t(rolling.par), 
                    order.by = index(r)[seq(length(r[initial]),
                                            length(r),n)])

save(list = c("rolling.param", "k", "n", "p", "r", "initial", "r.mat"),
     file = paste0("MNG_",k,"_roll_par_",n,"_step","_data.RData"))
