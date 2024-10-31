#############
# Simulation codes.
rm(list = ls())
library(dplyr)
library(ggplot2)
library(randomForest)
library(MASS)
library(foreach)
require(doSNOW)
library(caret)


name = "01-01simulation"
source("simulationfunctions.R")
set.seed(1)
if(! dir.exists(paste0("result/",name))){
  dir.create(paste0("result/",name))
}

# Simulation settings 
method_set <- c("aipw", "ipw", "reg")   # ATE estimators
N_set = c(10000)                        # first phase sample size
n_p_set = c(500, 1000)                  # pilot data sample size
n_set = c(1000, 1200, 1500, 2000, 4000) # other second phase sample size
K_set = c(2, 4, 6, 8, 10)               # stratum number
alpha_set = c(0.05, 0.1, 0.2, 0.3)      # alpha control probability not too small

expand.grid(method_set, N_set, n_p_set, n_set,K_set, alpha_set)->settings
colnames(settings) <- c("method", "N", "n_p", "n", "K", "alpha")

settings[order(settings$method, settings$N, 
               settings$n_p, settings$n, 
               settings$K, settings$alpha), ] -> settings
rownames(settings) <- c(1:nrow(settings))
#########################################
for(j in 1:nrow(settings)){
  method = settings$method[j]
  N = settings$N[j]
  n_p = settings$n_p[j]
  n = settings$n[j]
  K = settings$K[j]
  alpha = settings$alpha[j]
  cat(" Method=", as.character(method), ", N=", N,", n_p=", n_p,", n=", n, ", K=", K,  ", alpha=", alpha,"\n")
  starttime <- Sys.time()
  testf(N, n_p,n,K,method,alpha) -> rr
  names(rr) <- c("SimRan", "Var", 
                 "Oracle", "Var", 
                 "AdaStrat_Network", "Var", 
                 "AdaStrat_RandomForest", "Var", 
                 "FixStrat", "Var")
  print(rr)
  endtime <- Sys.time()
  cat( "\nRunning time: ", endtime-starttime, "\n")
}







#########################################
# Repeat
cores = 40
trials =800 # repeats

for(j in 1:nrow(settings)){
  method = settings$method[j]
  N = settings$N[j]
  n_p = settings$n_p[j]
  n = settings$n[j]
  K = settings$K[j]
  alpha = settings$alpha[j]
  cat(" Method=", method, ", N=", N,", n_p=", n_p,", n=", n, ", K=", K,  ", alpha=", alpha,"\n")
  starttime <- Sys.time()
  cls <- makeSOCKcluster(cores)
  registerDoSNOW(cls)
  pb <- txtProgressBar(max=trials, style=3, char = "*",)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress=progress)
  result <-
    foreach(i=c(1: trials),
            .combine = rbind,
            .options.snow=opts,
            .errorhandling = "pass",
            .packages = c("dplyr","randomForest",
                          "nnet","caret",
                          "MASS"))%dopar%{
                            source("simulationfunctions.R")
                            return(testf(N,
                                         n_p,
                                         n,
                                         K,
                                         method,
                                         alpha))}
  stopCluster(cls)
  endtime <- Sys.time()
  cat( "\nRunning time: ", endtime-starttime, "\n")
  path = paste0("result/", name, "/", paste(N, n_p, n, K, alpha, method, sep = "_"), ".rdata")
  save(result, file = path)
}

