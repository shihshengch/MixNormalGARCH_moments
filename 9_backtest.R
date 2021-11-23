# 9_backtest.R

# clear
rm(list = ls())

# specify
k = 2
n = 10
p = 0.025

# load MN(k)-GRACH(1,1) n-step rolling parameter estimation results
load(paste0("MNG_",k,"_roll_par_",n,"_step","_data.RData"))

# load MN(k)-GRACH(1,1) n-step rolling simulated aggregated returns
load(paste0("MNG_",k,"_R_",n,"_step_simulation","_data.RData"))

# load sources
source("1_source_functions.R")

# parallel computing
no_cores = detectCores()-1
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = ls(), envir = environment())
clusterEvalQ(cl = cl, source("1_source_functions.R"))

N = 4
backtest.N.4 = backtest(p, N, n, rolling.param, 
                        r.mat, R.mat_tn, r, initial, cl)

N = 8
backtest.N.8 = backtest(p, N, n, rolling.param, 
                        r.mat, R.mat_tn, r, initial, cl)

stopCluster(cl)

ES.JSU.backtest = rbind(backtest.N.4$p_JSU, backtest.N.8$p_JSU)
ES.Sim.backtest = rbind(backtest.N.4$p_Sim, backtest.N.8$p_Sim)
rownames(ES.JSU.backtest) = paste0("N = ", c(4,8))
rownames(ES.Sim.backtest) = paste0("N = ", c(4,8))
colnames(ES.JSU.backtest) = c("LRT", "Pearson", "Nass")
colnames(ES.Sim.backtest) = c("LRT", "Pearson", "Nass")

save(list = c("ES.JSU.backtest", "ES.Sim.backtest", "k", "n", "p", "r", 
              "initial", "r.mat", "backtest.N.4", "backtest.N.8"),
     file = paste0("MNG_",k,"_backtest_",n,"_step","_data.RData"))

