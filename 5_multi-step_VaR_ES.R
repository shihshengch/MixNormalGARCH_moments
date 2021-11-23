# 5_multi-step_VaR_ES.R

# clear
rm(list = ls())

# specify k
k = 2

# load MN(k)-GRACH(1,1) parameter estimation results
load(paste0("MNG_",k,"_data.RData"))

# load sources
source("1_source_functions.R")

# specify n and p
n = c(5,10,20)
p = c(0.01, 0.025, 0.05, 0.1)

# Simulate r_tn
no_cores = detectCores()-1
cl = makeCluster(no_cores)
clusterExport(cl = cl, varlist = ls(), envir = environment())
clusterEvalQ(cl = cl, source("1_source_functions.R"))
r_tn = pbreplicate(100, r_tn_sim(max(n), param = opt.par, r, S = 1000), 
                   cl = cl)
stopCluster(cl)
dim(r_tn) = c(max(n), 100000)

# compute VaR and ES
result_VaR = result_ES = list()
for (i in 1:length(n)) {
  R_tn = colSums(r_tn[1:n[i],])
  temp_JSU = sapply(p, function(x) VaR_ES_JSU(x, n[i], opt.par, r))
  temp_Sim = sapply(p, function(x) VaR_ES_Sim(x, R_tn))
  result_VaR[[i]] = cbind(JSU = temp_JSU[1,], Sim = temp_Sim[1,])
  result_ES[[i]] = cbind(JSU = temp_JSU[2,], Sim = temp_Sim[2,])
  rownames(result_VaR[[i]]) = rownames(result_ES[[i]]) = paste("p =", p)
}

names(result_VaR) = names(result_ES) = paste("n =", n)

result_VaR
result_ES