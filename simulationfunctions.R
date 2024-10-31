##################
# functions for simulations

source("esATE.R")
source("esATE_w2.R") # source the ATE estimation methods function


# given the influence function value and second phase data size
# obtain the oracle sample probability 
# phi_a is the estimated influence function value
# sn is the second phase sample size
# alpha hyperparameter that controls the probability from becoming too small
getpi <- function(phi_a, sn, alpha = 0.05){
  ## phi_a is sorted abs(phi)
  pi_a = phi_a/sum(phi_a)
  if(max(pi_a)<=(1/sn)){
    pi_a = (1-alpha)*pi_a + alpha * mean(pi_a)
    return(pi_a)
  }else{
    N = length(phi_a)
    H <- re <- c()
    k=0
    for(g in 0:(sn-1)){
      k=g+1
      H[k] = sum(phi_a[1:(N-k)])/(sn-k)
    }
    for(k in 1:sn){
      (H[k] >= phi_a[N-k])&(H[k] <= phi_a[N-k+1]) -> re[k]
    }
    g = which(re==T)
    H[g] -> H
    pi = pmin(H, phi_a)
    pi = pi/sum(pi)
    pi <- ifelse(pi>(1/sn), 1/sn, pi)
    pi = (1-alpha)*pi + alpha * mean(pi)
    return(pi)
  }
}

# given the influence function value and second phase data size
# obtain the intra-stratum sample probability for some fixed stratum strategy.
# t is the
getpi_stra <- function(t, N, n, K, alpha = 0.05){ 
  ## phi_a is sorted abs(phi)
  p = rep(1/K, K)
  phi_a = t/sum(t*p)*n/N
  if(max(phi_a)<=(1)){
    phi_a = (1-alpha)*phi_a + alpha * mean(phi_a)
    return(phi_a)
  }else{
    H <- re <- c()
    for(g in 1:(n*K/N-1)){
      H[g] = sum(t[1:(K-g)]*p[1:(K-g)])/(n/N-sum(p[(K-g+1): K]))
    }
    for(g in 1:(n*K/N-1)){
      (H[g] <= t[K-g+1]& H[g] >= t[K-g]) -> re[g]
    }
    g = which(re==T)[1]
    H[g] -> H
    pi = t/H
    pi <- ifelse(pi>1, 1, pi)
    pi = (1-alpha)*pi + alpha * mean(pi)
    return(pi)
  }
}

# Simulation data generation, generate the real first phase data
# X1 is the observed covariates (cheap)
# U is the unboserved covariates (expensive)
# Y1 and Y0 be the potential outcome
# Y is the observed outcome 
# Z is the treatment
Generatedt <- function(N){
  id <- 1:N
  X1 = runif(N, 0, 2)
  U1 = 0.5 + 0.5*X1 - 2*sin(X1) + 2*sign(sin(5*X1)) + rnorm(N, sd = 2)
  Y0 = - X1 - U1 + rnorm(N, sd = 2)
  Y1 = - X1 + 4*U1 + rnorm(N, sd = 2)
  e = plogis(1 - 0.5*X1-0.5*U1)
  Z <- rbinom(N, 1, e)
  Y = Z*Y1 + (1-Z)*Y0
  dt1 <- cbind(id, X1, U1, Y0, Y1, e, Z, Y) %>% as.data.frame()
  return(dt1)
}


# This function pick up the second phase data from first phase data
# by the design based probabilities and give an estimator of ATE.
# dt is the first phase data
# pi is the designed sampling probabilities for second phase data collection.
# method is the ATE estimation method
get_tau_hat <- function(dt, pi, method){
  dt$pi = pi
  N = nrow(dt)
  dt$V = rbinom(N, 1, dt$pi)
  dt_Z= filter(dt, dt$V==1)
  esATE_w(cbind(dt_Z$X1, dt_Z$U1), dt_Z$Y, dt_Z$Z, w = 1/dt_Z$pi, method, N) -> tau_hat
  return(tau_hat)
}

#####
# This function pick up the second phase data from other first phase data
# by the design based probabilities and give an estimator of ATE. Further 
# combine the pilot data ATE estimator with second phase data ATE estimator.
# dt_p is the pilot data
# dt1_p is the first phase data remove pilot data
# pi1_p is the designed sampling probability for second phase data collection
# method is the ATE estimation method
get_tau_hat2 <- function(dt1_p, dt_p, pi1_p, method){
  dt1_p$pi = pi1_p;
  N1_p = nrow(dt1_p);n_p = nrow(dt_p); N = n_p + N1_p
  
  # collect the second phase data and redifined as dt1_p
  dt1_p$V = rbinom(N1_p, 1, dt1_p$pi)
  dt1_p= filter(dt1_p, dt1_p$V==1) 
  
  # estimate ATE based on pilot data
  esATE(cbind(dt_p$X1, dt_p$U1), dt_p$Y, dt_p$Z, method) -> fit2
  
  # estimate ATE based on second phase data
  esATE_w(cbind(dt1_p$X1, dt1_p$U1), dt1_p$Y, dt1_p$Z, w = 1/dt1_p$pi, method, N) ->fit1_2
  
  # combined above two ATE estimators.
  (fit1_2[1]/fit1_2[2] + fit2$est/fit2$ve)/(1/fit1_2[2] + 1/fit2$ve) -> tau_hat
  (fit1_2[2]*fit2$ve)/(fit1_2[2]+fit2$ve) -> ve_hat
  
  cb_2 <- c(tau_hat, ve_hat)
  return(c(cb_2))
}

testf <- function(N, n_p, n, K, method, alpha){
  # N: the size of first phase data
  # n_p: the size of pilot data
  # n:  the size of other second phase data appart from pilot data.
  # method: is the ATE estimation method
  # alpha hyperparameter that controls the probability from becoming too small
  
  
  # get first phase data 
  Generatedt(N) -> dt1 -> dt
  
  #######
  # Outcome dependent stratification:
  # Stratify the first phase by the quantile of outcome Y, into K strata.
  # R_O = 1,...K. denote the stratum indicator. 
  te <- dt1 %>% dplyr::select(id, Y)
  te[order(te$Y), ] -> te
  R_O = rep(1:K, each = ceiling(N/K)) %>% as.factor()
  te$R_O =  R_O[1:N]
  dt1$R_O = te[order(te$id), "R_O"]
  
  ######
  # Get the true ATE value and true  influence function value (unobserved)
  # For oracle sampling design.
  x = cbind(dt1$X1, dt1$U1); y = dt1$Y; A = dt1$Z
  esATE(x, y, A, method) -> fit1
  dt1$phi = fit1$infl
  phi_bar = mean(dt1$phi)
  dt1$phi_a = abs(dt1$phi-phi_bar) # mainly for oracle sampling design
  # stratification by the true influence function (unobserved)
  # denoted as R_T = 1,...K.
  dt1 <- dt1[order(dt1$phi_a), ];
  R_T = rep(1:K, each = ceiling(N/K)) %>% as.factor()
  dt1$R_T = R_T[1:N]
  
  #####
  # uniform sampling design (Simple random sampling design)
  dt1$pi_u <- rep((n+n_p)/N, N)
  
  #####
  # Oracle sampling design (unobserved)
  getpi(dt1$phi_a, (n+n_p), alpha)*(n+n_p) -> dt1$pi_b
  
  ############
  # the pilot data collection
  pliot_ind <- sample(1:N, size = n_p, replace = F)
  #### The pilot data
  dt_p <- dt1[pliot_ind,]
  #### The other first phase data
  dt1_p <- dt1[-pliot_ind,]
  N1_p = nrow(dt1_p)
  
  
  # The pilot data collection and influence functions estimation based on pilot data
  x2 = cbind(dt_p$X1, dt_p$U1); y2 = dt_p$Y; A2 = dt_p$Z
  esATE(x2, y2, A2, method) -> fit2
  dt_p$phi_hat <- fit2$infl
  phi_hat_bar <- mean(dt_p$phi_hat)
  dt_p$phi_hat_a <- abs(dt_p$phi_hat - phi_hat_bar)
  
  ################
  #Lu's method
  #Pre fixed stratified strategy (Outcome quantile) and its intra-stratum sampling probability 
  te <- c()
  for(k in 1:K){
    te[k] <- mean((dt_p[dt_p$R_O==k,"phi_hat"]-phi_hat_bar)^2)
  }
  te2 <- te
  getpi_stra(sqrt(te2), N1_p, n, K, alpha) -> pi_O
  
  pr_O <- data.frame(R_O = 1:K, p_O = pi_O)
  merge(dt1_p, pr_O, by="R_O") -> dt1_p_O
  merge(dt_p, pr_O, by="R_O") -> dt_p_O
  
  ###################################################
  # The Adaptive stratified sampling design.
  # The clustering problem: cluster by the order of influence function value of pilot data.
  dt_p <- dt_p[order(dt_p$phi_hat_a), ]
  R_L = rep(1:K, each = ceiling(n_p/K)) %>% as.factor()
  dt_p$R_L = R_L[1:n_p]
  te <- c()
  for(k in 1:K){
    te[k] <- mean((dt_p[dt_p$R_L==k,"phi_hat"]-phi_hat_bar)^2)
  }
  te2 <- te
  getpi_stra(sqrt(te2), N1_p, n, K,  alpha) -> pi_m
  dt_p$pi_m = NA
  for(k in 1:K){
    dt_p[dt_p$R_L==k, "pi_m"] = pi_m[k]
  }
  # the classification problem: using pilot data to train a classification model.
  # and predict the stratum of each individual in first phase data.
  # Network model for classification
  digit.x = dt_p[, c("X1", "Y", "Z")]
  digit.y = dt_p[, "R_L"]
  digit.ml <- train(x=digit.x, y=digit.y,
                    method="nnet",
                    tuneGrid=expand.grid(
                      # 5个隐藏神经元
                      .size=10,
                      # 衰变率
                      .decay=0.1
                    ),
                    trControl=trainControl(method="none"),
                    # 最大权重数量
                    MaxNWts=10000,
                    # 最大迭代次数
                    maxit=500,
                    trace = FALSE)
  dt1_p$R_Z <- predict(digit.ml, newdata = dt1_p)
  pr_Z <- data.frame(R_Z = 1:K, p_Z = pi_m)
  merge(dt1_p, pr_Z, by="R_Z") -> dt1_p_Z
  
  # The random forest method for classification
  fit_f <- randomForest(R_L~X1+Y+Z , data = dt_p, proximity = T)
  dt1_p$R_F <- predict(fit_f, newdata = dt1_p)
  pr_F <- data.frame(R_F = 1:K, p_F = pi_m)
  merge(dt1_p, pr_F, by="R_F") -> dt1_p_F
  
  ######################################################
  get_tau_hat(dt1, dt1$pi_u, method) -> fit_u
  
  get_tau_hat(dt1, dt1$pi_b, method) -> fit_b
  
  get_tau_hat2(dt1_p_Z, dt_p, dt1_p_Z$p_Z, method) -> fit_Z
  
  get_tau_hat2(dt1_p_F, dt_p, dt1_p_F$p_F, method) -> fit_F
  
  get_tau_hat2(dt1_p_O, dt_p, dt1_p_O$p_O, method) -> fit_O
  
  return(c(fit_u, fit_b, fit_Z, fit_F, fit_O))
}

result_td <- function(re){
  re = cbind(apply(re, 2, mean, na.rm = T)-0.50,
             apply(re, 2, var, na.rm = T),
             (apply(re, 2, mean, na.rm = T)-0.50)^2 +apply(re, 2, var, na.rm = T) )
  colnames(re) <- c("bias", "var", "mse")
  return(re)
}
library(dplyr)
testf(N=10000, n_p=400, n=2000, K = 8, method = "reg", alpha = 0.1)

