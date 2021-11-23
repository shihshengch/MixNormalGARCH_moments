# 1_source_functions.R

# library
if (!require("MSGARCH"      )) install.packages("MSGARCH")
if (!require("optimParallel")) install.packages("optimParallel")
if (!require("xts"          )) install.packages("xts")
if (!require("tidyr"        )) install.packages("tidyr")
if (!require("ggplot2"      )) install.packages("ggplot2")
if (!require("ggthemes"     )) install.packages("ggthemes")
if (!require("expm"         )) install.packages("expm")
if (!require("parallel"     )) install.packages("parallel")
if (!require("pbapply"      )) install.packages("pbapply")
if (!require("LaplacesDemon")) install.packages("LaplacesDemon")
if (!require("SuppDists"    )) install.packages("SuppDists")
if (!require("runner"       )) install.packages("runner")
if (!require("matrixStats"  )) install.packages("matrixStats")

library(MSGARCH)   # for staring values in parameter optimization
library(optimParallel) # for parallel optimization
library(xts)       # for index()
library(tidyr)     # for gather()
library(ggplot2)   # for ggplot
library(ggthemes)  # for ggplot
library(expm)      # for matrix power %^%
library(parallel)  # for parallel computation
library(pbapply)   # for parallel computation
library(LaplacesDemon) # for rcat
library(SuppDists) # for Johnson SU distribution
library(runner)    # for rolling return matrix
library(matrixStats)   # for column quantile

########################## Codes for Chapter 3 ###########################

# MN(k)-GARCH(1,1) model parameter estimation

# MN(k)-GARCH(1,1) objects
MNG.obj = function(param) {
  
  # assign model parameters into matrix form
  # param: MNG parameters (vector): 
  #        (lambda_1,...,lambda_{k-1}, mu_1,...mu_{k-1}, 
  #           (omega,alpha,beta)_1,...(omega,alpha,beta)_k)
  
  # numbers of components
  k = (length(param)+2) / 5
  
  # mixed normal parameters
  lambda = param[1:(k-1)]
  lambda = matrix( c(lambda, 1-sum(lambda)), ncol = 1 )
  
  mu = param[k:(2*k-2)]
  mu = matrix( c(mu, -sum(lambda[-k]*mu)/lambda[k]), ncol = 1 )
  
  # GARCH parameters
  G.par = matrix( param[-(1:(2*k-2))], ncol = 3, byrow = TRUE )
  omega = G.par[,1, drop = FALSE]
  alpha = G.par[,2, drop = FALSE]
  beta  = diag(G.par[,3])
  
  # collections
  object = list(
    k = k, lambda = lambda, mu = mu, 
    omega = omega, alpha = alpha, beta = beta
  )
  
  return(object) 
}

# compute individual and overall conditional variance
MNG.eps2var = function(param, eps) {
  
  # compute from known parameters and know epsilon
  # param: MNG parameters (vector) 
  # eps: time series (vector)
  
  # time horizon
  T = length(eps)
  eps = as.vector(eps)
  
  # generate MNG object
  list2env(MNG.obj(param), envir = environment())
  
  # compute individual unconditional variance 
  # E(X_t) = (I - C_11)^(-1)(omega_tilde) 
  C_11 = alpha %*% t(lambda) + beta
  omega_tilde = omega + alpha %*% t(lambda) %*% mu^2
  E.X_t = solve(diag(k) - C_11, tol = 1e-30) %*% omega_tilde
  
  # initialize conditional variance of j-th component
  X = matrix(0, nrow = k, ncol = T)
  for (j in 1:k) {
    X[j,] = E.X_t[j]
  }
  
  # compute individual conditional variance at time t
  for (t in 2:T) {
    X[,t] = omega + alpha * eps[t-1]^2 + beta %*% X[,t-1, drop = FALSE]
  }
  
  # compute overall conditional variance at time t
  # h_t := (sigma_t)^2 = (lambda)'(X_t + mu^2)
  h = t(lambda) %*% sweep(X, 1, mu^2, "+") 
  
  return(list(X = X, h = c(h)))
}

# Extended Augmented Likelihood Estimation function
MNG.EALE = function(param, eps) {
  
  # objective function for parameter estimation
  # param: MNG parameters (vector): 
  # eps: time series (vector)
  
  # time horizon
  T = length(eps)
  eps = as.vector(eps)
  
  # generate MNG object
  list2env(MNG.obj(param), envir = environment())
  
  # check stationarity
  C_11 = alpha %*% t(lambda) + beta
  rho = max(Mod(eigen(C_11)$values))
  reject = runif(1, 1e7-50, 1e7+50)
  if (!(rho < 1)) return(reject)
  
  # compute individual conditional variance and volatility
  X = MNG.eps2var(param, eps)$X
  sigma = sqrt(X)
  
  # compute EALE
  # sum_{t=1}^{T} { log( sum_{j=1}^{k} {lambda_j * phi_{jt}}) }
  term.1 = 0
  for (t in 1:T) {
    term.1 = term.1 + log(t(lambda) %*% dnorm(eps[t], mu, sigma[,t]))
  }
  
  # sum_{j=1}^{k} { log(g[j]) }, g[j] = (prod_{t=1}^T phi_{jt})^{1/T}
  g = c()
  for (j in 1:k) {
    g[j] = exp(mean(log(dnorm(eps, mu[j], sigma[j,]))))
  }
  term.2 = sum(log(g))
  
  # sum_{j=1}^{k} { log(1+v[j]) }, 
  # v[j] = 1/T * sum_{t=1}^T {(phi_{jt} - g[j])^2}
  term.3 = 0
  for (j in 1:k) {
    term.3 = term.3 + log(1 + mean((dnorm(eps, mu[j], sigma[j,]) - 
                                      g[j])^2))
  }
  
  EALE = term.1 + term.2 - term.3
  if (is.infinite(EALE)) return(reject)
  if (is.nan(EALE)) return(reject)
  if (is.na(EALE)) return(reject)
  
  # maximize EALE = minimize negative EALE
  return(-EALE)
}

# Generate optimization arguments
MNG.bounds = function(k) {
  
  # Generate lower bound, upper bound, and initial guess
  # k: k components
  
  # mixed normal parameters
  lambda.lb  = rep(1e-6, k-1)
  lambda.ub  = 1 - lambda.lb
  
  mu.lb  = rep(-0.1, k-1)
  mu.ub  = -mu.lb
  
  # GARCH parameters
  omega.lb  = rep(1e-6, k)
  omega.ub  = 1 - omega.lb
  
  alpha.lb  = rep(1e-6, k)
  alpha.ub  = 1 - alpha.lb
  
  beta.lb  = rep(1e-6, k)
  beta.ub  = 1 - beta.lb
  
  G.lb  = rbind(omega.lb, alpha.lb, beta.lb)
  G.ub  = rbind(omega.ub, alpha.ub, beta.ub)
  
  # collections
  lb  = c(lambda.lb, mu.lb, c(G.lb))
  ub  = c(lambda.ub, mu.ub, c(G.ub))
  
  return(list(lb = lb, ub = ub))
}

# MN(k)-GARCH(1,1) parameter estimation function
MNG.k.est = function(r, k, cl, par.only = FALSE, start.par = NULL) {
  
  # estimate parameters of MN(k)-GARCH(1,1)
  # r: return series
  # k: number of components
  # cl: number of clusters for parallel computation
  # par.only: return output with only optimal parameters, default is false
  # start.par: given the starting value manually, default is null
  
  c = mean(r)
  eps = r - c
  
  if (is.null(start.par)) {
    
    # starting values using package MSGARCH
    # symmetric MN(k)-GARCH(1,1)
    spec = CreateSpec(variance.spec = list(model = rep("sGARCH", k)),
                      distribution.spec=list(distribution=rep("norm",k)),
                      switch.spec = list(do.mix = TRUE, k = k))
    # Nelder-Mead method function
    f_custom_optim = function(vPw, f_nll, spec, data, do.plm){
      out = optim(vPw, f_nll, spec = spec, data = data,
                  do.plm = do.plm, method = "Nelder-Mead")
      return(out)
    }
    sym.MNG = tryCatch(
      # L-BFGS-B method
      FitML(spec = spec, data = eps), 
      # Nelder-Mead method if L-BFGS-B error
      error = function(e) FitML(spec = spec, data = eps, 
                                ctr = list(OptimFUN = f_custom_optim))
    )
    
    # reorder starting values
    start.par = c(sym.MNG$par[-(1:(3*k))], 
                  rep(0, k-1), sym.MNG$par[1:(3*k)])
    names(start.par) = c(sprintf("lambda_%s", 1:(k-1)), 
                         sprintf("mu_%s"    , 1:(k-1)),
                         c(t(
                           sapply(c("omega_%s", "alpha_%s", "beta_%s"), 
                                  function(x) sprintf(x, 1:k))
                         )))
  } 
  
  # define lower and upper bounds
  bounds = MNG.bounds(k = k)
  
  # quasi-Newton optimization
  if (is.null(cl)) {
    opt.result = optim(par = start.par, fn = MNG.EALE, eps = eps,
                       lower = bounds$lb, upper = bounds$ub,
                       method = "L-BFGS-B", control = list(maxit = 2000))
  } else {
    opt.result = optimParallel(par = start.par, fn = MNG.EALE, eps = eps,
                               lower = bounds$lb, upper = bounds$ub,
                               method = "L-BFGS-B", 
                               control = list(maxit = 2000), 
                               parallel = list(cl = cl))
  }
  opt.par = opt.result$par
  
  if (par.only) {
    result = opt.par
  } else {
    result = list(eps = eps, r = r, c = c, k = k, 
                  start.par = start.par, bounds = bounds, 
                  opt.result = opt.result, opt.par = opt.par)
  }
  
  return(result)
}

########################## Codes for Chapter 4 ###########################

# Central Moments multi-step-ahead aggregated returns over n time periods
Moment.R = function(ith, n, param, eps) {
  
  # Compute the centered moment M_{R,n}^{(i)}
  # ith: the i-th centered moment
  # n: multi-step-ahead over n time periods
  # param: MNG parameters (vector): 
  # eps: time series (vector)
  
  # Error message if ith not in 1 to 4
  if (!(ith %in% 1:4)) {
    stop("Only the first four moments can be computed!!!")
  }
  
  # generate MNG object
  list2env(MNG.obj(param), envir = environment())
  
  # Compute known variables
  ## X_{t+1} = E_t(X_{t+1})
  eps_t = as.numeric(eps[length(eps)])
  X_t   = MNG.eps2var(param, eps)$X[,length(eps), drop = FALSE]
  X_t1  = omega + alpha * eps_t^2 + beta %*% X_t
  
  ## C_11
  C_11 = alpha %*% t(lambda) + beta
  
  ## E(X_t)
  omega_tilde = omega + alpha %*% t(lambda) %*% mu^2
  E.X_t = solve(diag(k) - C_11, tol = 1e-30) %*% omega_tilde
  
  ## E(eps2_t)
  E.eps2_t = c(t(lambda) %*% (E.X_t + mu^2))
  
  # function(i) E_t(X_{t+i})
  Et.X_ti = function(i) {
    if (i == 1) {
      Et.X_ti = X_t1
    } else {
      Et.X_ti = E.X_t + (C_11%^%(i-1)) %*% (X_t1 - E.X_t)
    }
    return(Et.X_ti)
  }
  
  # function(i) E_t(eps^2_{t+i})
  Et.eps2_ti = function(i) {
    Et.eps2_ti = t(lambda) %*% (Et.X_ti(i) + mu^2)
    return(c(Et.eps2_ti))
  }
  
  # function(i) E_t(eps^3_{t+i})
  Et.eps3_ti = function(i) {
    Et.eps3_ti = 3 * t(lambda) %*% (Et.X_ti(i) * mu) + t(lambda) %*% mu^3
    return(c(Et.eps3_ti))
  }
  
  # function(i) E_t(eps_{t+i}eps^2_{t+n})
  Et.eps_ti.eps2_tn = function(i,n) {
    Et.eps_ti.eps2_tn = t(lambda) %*% (C_11%^%{n-i-1}) %*% 
      alpha * Et.eps3_ti(i)
    return(c(Et.eps_ti.eps2_tn))
  }
  
  # vectorize matrix function
  vec = function(M) {
    dim(M) = c(dim(M)[1]*dim(M)[2], 1)
    return(M)
  }
  
  # function(i) E_t(eps^4_{t+i})
  Et.X_ti.X_ti = function(i) {
    
    C_21 = (alpha %x% omega + omega %x% alpha) %*% t(lambda) + 
      omega %x% beta +  beta %x% omega + 
      6 * (alpha %x% alpha) %*% t(lambda * mu^2) +
      (alpha %x% beta + beta %x% alpha) * c(t(lambda) %*% mu^2)
    
    C_22 = 3 * (alpha %x% alpha) %*% t(vec(diag(c(lambda)))) + 
      (alpha %*% t(lambda)) %x% beta + beta %x% (alpha %*% t(lambda)) + 
      beta %x% beta
    
    zeros = matrix(0, nrow = k, ncol = k^2)
    
    C = rbind(cbind(C_11, zeros), cbind(C_21,C_22))
    
    W = rbind(omega_tilde, omega %x% omega + 
                (alpha %x% omega + 
                   omega %x% alpha) %*% t(lambda) %*% mu^2 + 
                (alpha %x% alpha) %*% t(lambda) %*% mu^4)
    
    Et.Y_ti = matrix(0, nrow = nrow(W), ncol = i)
    Et.Y_ti[,1] = rbind(X_t1, vec(X_t1 %*% t(X_t1)))
    
    if (i > 1) {
      for (d in 2:i) {
        Et.Y_ti[,d] = W + C %*% Et.Y_ti[,d-1, drop = FALSE]
      }
    }
    
    Et.vec.X_ti.X_ti = Et.Y_ti[-(1:k), i]
    
    Et.X_ti.X_ti = matrix(Et.vec.X_ti.X_ti, nrow = k)
    
    return(Et.X_ti.X_ti)
  }
  
  # function(i) E_t(eps_{t+i}^4)
  Et.eps4_ti = function(i) {
    Et.eps4_ti = 3 * t(vec(diag(c(lambda)))) %*% vec(Et.X_ti.X_ti(i)) +
      6 * t(lambda) %*% (Et.X_ti(i) * mu^2) + t(lambda) %*% mu^4
    return(c(Et.eps4_ti))
  }
  
  # function(i) E_t(eps_{t+i} * eps_{t+n}^3)
  Et.eps_ti.eps3_tn = function(i, n) {
    Et.eps_ti.eps3_tn = 3 * t(lambda) %*% (((C_11%^%(n-i-1)) %*% alpha) 
                                           * mu) * Et.eps3_ti(i)
    return(c(Et.eps_ti.eps3_tn))
  }
  
  # function(i) E_t(eps_{t+i}^2 * eps_{t+n}^2)
  Et.eps2_ti.eps2_tn = function(i, n) {
    Et.eps2_ti.eps2_tn = Et.eps2_ti(i) * E.eps2_t + 
      t(lambda) %*% (C_11%^%(n-i-1)) %*% 
      (omega * Et.eps2_ti(i) + alpha * Et.eps4_ti(i) + 
         beta %*% (Et.X_ti.X_ti(i) + Et.X_ti(i) %*% t(mu^2)) %*% lambda -
         Et.eps2_ti(i) * E.X_t)
    return(c(Et.eps2_ti.eps2_tn))
  }
  
  # function(i) E_t(eps_{t+i} * eps_{t+j} * eps_{t+n}^2)
  Et.eps_ti.eps_tj.eps2_tn = function(i, j, n) {
    Et.eps_ti.eps_tj.eps2_tn = 3 * (t(lambda) %*% 
                                      (C_11%^%(n-j-1)) %*% alpha) %*% 
      (t(lambda) %*% (((C_11%^%{j-i-1}) %*% alpha) * mu)) * Et.eps3_ti(i)
    return(c(Et.eps_ti.eps_tj.eps2_tn))
  }
  
  # Moments
  switch(ith,
         "1" = 0,
         
         "2" = {
           V = rep(0, n)
           V[1] = Et.eps2_ti(1)
           if (n > 1) {
             for (d in 2:n) {
               V[d] = V[d-1] + Et.eps2_ti(d)
             }
           }
           return(V)
         },
         
         "3" = {
           S = rep(0,n)
           S[1] = Et.eps3_ti(1)
           if (n > 1) {
             for (d in 2:n) {
               S[d] = S[d-1] + 3 * sum(
                 sapply(1:(d-1), 
                        function(i) Et.eps_ti.eps2_tn(i, n = d)) ) + 
                 Et.eps3_ti(d)
             }
           }
           return(S)
         },
         
         "4" = {
           K = rep(0,n)
           K[1] = Et.eps4_ti(1)
           if (n > 1) {
             for (d in 2:n) {
               
               if (d == 2) {
                 sum.Et.eps_ti.eps_tj.eps2_tn = 0
               } else { 
                 sum.Et.eps_ti.eps_tj.eps2_tn = 0
                 for (r in 1:(d-2)) {
                   for (s in (r+1):(d-1)) {
                     sum.Et.eps_ti.eps_tj.eps2_tn = 
                       sum.Et.eps_ti.eps_tj.eps2_tn +
                       Et.eps_ti.eps_tj.eps2_tn(r, s , n = d)
                   }
                 }
               }
               
               K[d] = K[d-1] + 4 * sum(
                 sapply(1:(d-1), 
                        function(i) Et.eps_ti.eps3_tn(i, n = d)) ) + 
                 6 * sum(
                   sapply(1:(d-1), 
                          function(i) Et.eps2_ti.eps2_tn(i, n = d)) ) + 
                 12 * sum.Et.eps_ti.eps_tj.eps2_tn + Et.eps4_ti(d)
             }
           }
           return(K)
         }
  )
}

# multi-step-ahead conditional skewness
skew = function(n, param, eps) {
  
  # n: multi-step-ahead over n time periods
  # param: MNG parameters (vector): 
  # eps: time series (vector)
  
  Sn = Moment.R(ith = 3, n, param, eps)
  Vn = Moment.R(ith = 2, n, param, eps)
  skew = Sn / (Vn^(2/3))
  return(skew)
}

# multi-step-ahead conditional kurtosis
kurt = function(n, param, eps) {
  
  # n: multi-step-ahead over n time periods
  # param: MNG parameters (vector): 
  # eps: time series (vector)
  
  Kn = Moment.R(ith = 4, n, param, eps)
  Vn = Moment.R(ith = 2, n, param, eps)
  kurt = Kn / (Vn^2)
  return(kurt)
}

########################## Codes for Chapter 5 ###########################

# Johnson SU VaR and ES
VaR_ES_JSU = function(p, n, param, r, output = 3, ret.list = FALSE) {
  
  # p: VaR or ES level, a small number
  # n: multi-step-ahead number
  # param: model parameters of MN(k)-GARCH(1,1)
  # r: return series
  # output: 1-only VaR output, 2-only ES output, 3-VaR and ES output
  # ret.list: for output 3, if TRUE return list; otherwise vector 
  
  c = mean(r)
  eps = r - c
  
  # moments matching to Johnson SU distribution using JohnsonFit
  # first non-central moment and second,third,fourth central moment
  t = c(n*c,
        Moment.R(2, n, param, eps)[n], 
        Moment.R(3, n, param, eps)[n], 
        Moment.R(4, n, param, eps)[n])
  
  # Johnson SU parameters
  JSU = JohnsonFit(t, moment = "use")
  
  z_p = qnorm(p)
  
  VaR_JSU = JSU$xi + JSU$lambda * sinh((z_p-JSU$gamma)/JSU$delta)
  
  ES_JSU = JSU$xi + JSU$lambda/(2*p) *
    (exp((1-2*JSU$gamma*JSU$delta)/(2*JSU$delta^2)) * 
       pnorm(z_p-1/JSU$delta) -
       exp((1+2*JSU$gamma*JSU$delta)/(2*JSU$delta^2)) * 
       pnorm(z_p+1/JSU$delta))
  
  result = switch(output, 
                  "1" = VaR_JSU,
                  "2" = ES_JSU,
                  "3" = ifelse(rep(ret.list,2), 
                               list(VaR = VaR_JSU, ES = ES_JSU),
                               c(VaR = VaR_JSU, ES = ES_JSU))
  )
  
  return(result)
}

# Simulated sample paths returns
r_tn_sim = function(n, param, r, S) {
  
  # n: multi-step-ahead number
  # param: model parameters
  # r: return series
  # S: sample paths
  
  c = mean(r)
  eps = r - c
  
  # generate MNG object
  list2env(MNG.obj(param), envir = environment())
  
  # compute conditional variance
  X_t = MNG.eps2var(param, eps)$X
  
  eps_ti = matrix(0, nrow = n+1, ncol = S)
  eps_ti[1,] = eps[length(eps)]
  r_ti = matrix(0, nrow = n, ncol = S)
  j = replicate(S, rcat(n, as.vector(lambda)))
  
  X_ti = array(dim = c(n+1, k, S))
  X_ti[1,,] = X_t[,length(eps)]
  
  for (i in 2:(n+1)) {
    X_ti[i,,] = sweep(alpha %x% eps_ti[i-1,,drop=FALSE]^2,1,omega,"+") + 
      beta %*% X_ti[i-1,,]
    eps_ti[i,] = sapply(1:S, 
                        function(s) rnorm(1, mean = mu[j[i-1,s]], 
                                          sd = sqrt(X_ti[i,j[i-1,s],s])))
    r_ti[i-1,] = c + eps_ti[i,]
  }
  
  return(r_ti)
}

# Simulated VaR and ES
VaR_ES_Sim = function(p, R_tn, output = 3, ret.list = FALSE) {
  
  # p: VaR or ES level
  # R_tn: Simulated n-step-ahead aggregated returns
  
  VaR_sim = quantile(R_tn, p, name = FALSE)
  
  S = length(R_tn)
  
  ES_sim = 1/(p*S) * sum(R_tn * (R_tn < VaR_sim))
  
  result = switch(output, 
                  "1" = VaR_sim,
                  "2" = ES_sim,
                  "3" = ifelse(rep(ret.list,2), 
                               list(VaR = VaR_sim, ES = ES_sim),
                               c(VaR = VaR_sim, ES = ES_sim))
  )
  
  return(result)
}

########################## Codes for Chapter 6 ###########################

# backtesting function
backtest = function(p, N, n, rolling.param, 
                    r.mat, R.mat_tn, r, initial, cl) {
  
  # performing Kratz et al. (2018) backtesting method
  # p: VaR or ES level
  # N: number of multiple VaR levels, 4 and 8 suggested
  # n: multi-step-ahead number
  # rolling.param: parameters for each rolling period
  # r.mat: simulated returns for each rolling period
  # R.mat_tn: simulated aggregated returns for each rolling period
  # r: return series
  # initial: the initial period of rolling period
  # cl: number of clusters for parallel computation
  
  p_j = (N - (1:N) + 1) / N * p
  gamma_j = -diff(c(1, p_j, 0))
  
  # nonoverlapping aggregated returns
  ep = seq(length(r[initial]),length(r),n)
  L_t = period.apply(x = r, ep, FUN = sum)
  L_t = na.omit(lag(L_t, -1))
  
  # rolling Johnson SU VaR of multi-step-ahead aggregated returns
  VaR_JSU_t_j = t(pbsapply(1:nrow(rolling.param), function(y) 
    sapply(p_j, 
           function(x) VaR_ES_JSU(x, n, as.vector(rolling.param[y,]), 
                                  r.mat[,y], output = 1)), 
    cl = cl))
  
  VaR_JSU_t_j = xts(VaR_JSU_t_j, order.by = index(rolling.param))
  
  I_JSU_t_j = sweep(VaR_JSU_t_j[-dim(VaR_JSU_t_j)[1],], 1, L_t, ">")
  
  MX_JSU_t = rowSums(I_JSU_t_j)
  
  O_JSU_j = sapply(0:N, function(j) sum(MX_JSU_t == j))
  names(O_JSU_j) = 0:N
  
  # rolling Simulated VaR of multi-step-ahead aggregated returns
  VaR_Sim_t_j = colQuantiles(R.mat_tn, probs = p_j)
  
  VaR_Sim_t_j = xts(VaR_Sim_t_j, order.by = index(rolling.param))
  
  I_Sim_t_j = sweep(VaR_Sim_t_j[-dim(VaR_Sim_t_j)[1],], 1, L_t, ">")
  
  MX_Sim_t = rowSums(I_Sim_t_j)
  
  O_Sim_j = sapply(0:N, function(j) sum(MX_Sim_t == j))
  names(O_Sim_j) = 0:N
  
  # p.value function for test statistics
  p.value = function(O_j, gamma_j) {
    
    T = sum(O_j)
    
    # LRT
    S_LRT = 2 * sum(O_j[O_j!=0] * log(O_j[O_j!=0] / (T*gamma_j[O_j!=0])))
    p_LRT = pchisq(S_LRT, N, lower.tail = FALSE)
    
    # Pearson chi-squared test
    S_N = sum((O_j - T*gamma_j)^2/(T*gamma_j))
    p_Pearson = pchisq(S_N, N, lower.tail = FALSE)
    
    # Nass test
    E.S_N = N
    Var.S_N = 2*N - (N^2+4*N+1)/T + 1/T * sum(1/gamma_j)
    c = 2*E.S_N / Var.S_N
    nu = c * E.S_N
    p_Nass = pchisq(c*S_N, nu, lower.tail = FALSE)
    
    return(c(p_LRT, p_Pearson, p_Nass))
  }
  
  p_JSU = p.value(O_JSU_j, gamma_j)
  p_Sim = p.value(O_Sim_j, gamma_j)
  
  names(p_JSU) = names(p_Sim) = c("LRT", "Pearson", "Nass")
  
  names(VaR_JSU_t_j) = names(VaR_Sim_t_j)
  
  result = list(VaR_JSU = VaR_JSU_t_j, VaR_Sim = VaR_Sim_t_j, 
                O_JSU = O_JSU_j, O_Sim = O_Sim_j, 
                p_JSU = p_JSU  , p_Sim = p_Sim)
  
  return(result)
}