##########################
# classification tree plot.
# 
rm(list = ls())
set.seed(1)
source("simulationfunctions.R")
source("esATE.R")
source("esATE_w2.R")
library(dplyr)
library(randomForest)
library(nnet)
library(caret)
library(MASS)
library(foreach)
require(doSNOW)
library(rpart)
library(rpart.plot)
n_m = 1000; n = 2000; method = "aipw"; alpha = 0.3; K = 8
cat(" Method=", method, ", n_m=", n_m,", n=", n, ", K=", K, ", alpha=",alpha, "\n")
load(file = "dt.rdata")
dt1 <- dt
names(dt1) <- c("f.eid", "Ethnicity", "Sex", "Age", "Education", "PA_time", 
                "BMI", "Alcohol", "Cigarette", "Highcholesterol", "Hypertension", 
                "Diabete", "Fruit", "Vitamin", "Redmeat",
                "Calcium", "Milk", 
                "Type_PA_warking", "Type_PA_other", "Type_PA_strenuous", "Type_PA_LightDIY", 
                "hExercise", "Type_PA_None", "blood_glucose", "BMD", 
                "gpc_1", "gpc_2", "gpc_3", "gpc_4", "gpc_5", "gpc_6", "gpc_7", 
                "gpc_8", "gpc_9", "gpc_10", "gpc_11", "gpc_12", "gpc_13", "gpc_14", 
                "gpc_15", "gpc_16", "gpc_17", "gpc_18", "gpc_19", "gpc_20", "gpc_21", 
                "gpc_22", "gpc_23", "gpc_24", "gpc_25", "gpc_26", "gpc_27", "gpc_28", 
                "gpc_29", "gpc_30", "gpc_31", "gpc_32", "gpc_33", "gpc_34", "gpc_35", 
                "gpc_36", "gpc_37", "gpc_38", "gpc_39", "gpc_40")
Xvar_c <-c("Ethnicity", "Sex", "Age", "Education", "PA_time", 
           "BMI", "Alcohol", "Cigarette", 
           "Highcholesterol", "Hypertension", 
           "Diabete",
           "Vitamin", "Milk")
dt1$id = dt1$f.eid
Xvar_e <- c("gpc_1", "gpc_2", "gpc_3", "gpc_4", 
            "gpc_5", "gpc_6", "gpc_7", "gpc_8", "gpc_9", "gpc_10")
Xvar <- c(Xvar_c, Xvar_e)

Yvar <- "BMD"
Avar <- "hExercise"
X <- dt1[, Xvar] 
X_c <- dt1[, Xvar_c]
Y <- dt1[, Yvar] 
A <- dt1[, Avar]

N = nrow(dt1)

if(method == "ipw"){
  fit1 <- esATE(X, Y, A, "ipw") # unknown (True)
  fit1_c <- esATE(X_c, Y, A, "ipw") # known (bias)
}
if(method == "aipw"){
  fit1 <- esATE(X, Y, A, "aipw")# unknown (True)
  fit1_c <- esATE(X_c, Y, A, "aipw") # known (bias)
}
if(method == "reg"){
  fit1 <- esATE(X, Y, A, "reg")# unknown (True)
  fit1_c <- esATE(X_c, Y, A, "reg") # known (bias)
}


# bench mark, the ate estimator if we have the information of expensive covarites.
cat(c(fit1$est - sqrt(fit1$ve)*1.96, fit1$est + sqrt(fit1$ve)*1.96))
# biased estimator, ignoring the expensive covariates.
cat(c(fit1_c$est - sqrt(fit1_c$ve)*1.96, fit1_c$est + sqrt(fit1_c$ve)*1.96))

dt1$phi = fit1$infl
# biased influence function, though biased, 
# can be used to enhance the classification model
dt1$phi_m = fit1_c$infl 

phi_bar = mean(dt1$phi)
phi_bar_m = mean(dt1$phi_m)

dt1$phi_a = abs(dt1$phi-phi_bar)
dt1$phi_m_a = abs(dt1$phi_m-phi_bar_m)

# stratification by the true influence function (unobserved)
# denoted as R_T = 1,...K.
dt1 <- dt1[order(dt1$phi_a), ];
R_T = rep(1:K, each = ceiling(N/K)) %>% as.factor()
dt1$R_T = R_T[1:N]

#####
# uniform sampling design (Simple random sampling design)
dt1$pi_u <- rep((n+n_m)/N, N)

#####
# Oracle sampling design (unobserved)
getpi(dt1$phi_a, (n+n_m), alpha)*(n+n_m) -> dt1$pi_b

#################
# the pilot data collection
pliot_ind <- sample(1:N, size = n_m, replace = F)
#### The pilot data
dt2 <- dt1[pliot_ind,]
#### The other first phase data
dt1_2 <- dt1[-pliot_ind,]
N1_2 = nrow(dt1_2)

# The pilot data collection and influence functions estimation based on pilot data
x2 <- dt2[, Xvar]
y2 <- dt2[, Yvar]
A2 <- dt2[, Avar]
esATE(x2, y2, A2, method) -> fit2
dt2$phi_hat <- fit2$infl
phi_hat_bar <- mean(dt2$phi_hat)
dt2$phi_hat_a <- abs(dt2$phi_hat - phi_hat_bar)
dt2 <- dt2[order(dt2$phi_hat_a), ]



####################
#Lu's method
###################
#Lu's group, K==8
# (dt1$Type_PA_HeavyDIY +ifelse(dt1$edu>median(dt1$edu), 1, 0)*10+
#    ifelse(dt1$age>median(dt1$age), 1, 0)*100) %>%
#   strtoi(base = 2) ->R_O
# dt1$R_O = R_O+1
cut(dt1$BMD, breaks = quantile(dt1$BMD, probs = seq(0, 1, length.out = K + 1)),
    labels = 1:K) -> R_O
dt1$R_O = R_O
table(dt1$R_O)
merge(dt1[, c("id", "R_O")], dt2) -> dt2
merge(dt1[, c("id", "R_O")], dt1_2) -> dt1_2
dt2 <- dt2[order(dt2$phi_hat_a), ]
K_O = length(table(dt2$R_O))
te <- c()
for(k in 1:K_O){
  te[k] <- mean((dt2[dt2$R_O==k,"phi_hat"]-phi_hat_bar)^2)
}
te2 <- te
getpi_stra(sqrt(te2), N1_2, n, K_O, alpha = 0) -> pi_O
pr_O <- data.frame(R_O = 1:K_O, p_O = pi_O)
merge(dt1_2, pr_O, by="R_O", all = T) -> dt1_2_O


###################################################
# The Adaptive stratified sampling design.
# The clustering problem: cluster by the order of influence function value of pilot data.
dt2 <- dt2[order(dt2$phi_hat_a), ]
R_L = rep(1:K, each = ceiling(n_m/K)) %>% as.factor()
dt2$R_L = R_L[1:n_m]
te <- c()
for(k in 1:K){
  te[k] <- mean((dt2[dt2$R_L==k,"phi_hat_a"]-phi_hat_bar)^2)
}
te2 <- te
getpi_stra(sqrt(te2), N1_2, n, K, alpha) -> pi_m

# the classification problem: using pilot data to train a classification model.
# and predict the stratum of each individual in first phase data.
# The nnet method for classification
ind = sample(1:n_m, size = n_m, replace = F)
formula = as.formula(paste0("R_L~", paste(c(Xvar_c, Yvar, Avar), collapse = "+")))
cv_f = 5
thvar_z_hat <- thvar_f_hat <- thvar_Tr_hat <- thvar_u_hat <- c()
for(cv_i in 1:cv_f){
  cat("fold ",cv_i, ", " )
  n_m/cv_f -> cv_n
  test_ind <- ind[(cv_n*cv_i-cv_n+1):(cv_n*cv_i)]
  dt2_train <- dt2[-test_ind, ]
  dt2_test <- dt2[test_ind, ]
  fit_Tr <- rpart(formula, data = dt2_train, cp = 0.005)
  dt2_test$R_test <- predict(fit_Tr, newdata = dt2_test, type = "class")
  pr_test <- data.frame(R_test = 1:K, p_test = pi_m)
  merge(dt2_test, pr_test, by="R_test") -> dt2_test_new
  mean((dt2_test_new$phi_hat)^2/dt2_test_new$p_test)/N -> thvar_Tr_hat[cv_i]
  
  fit_f <- randomForest(formula, data = dt2_train, proximity = T)
  dt2_test$R_test <- predict(fit_f, newdata = dt2_test)
  pr_test <- data.frame(R_test = 1:K, p_test = pi_m)
  merge(dt2_test, pr_test, by="R_test") -> dt2_test_new
  mean((dt2_test_new$phi_hat)^2/dt2_test_new$p_test)/N -> thvar_f_hat[cv_i]
  
  mean((dt2_test_new$phi_hat)^2/(n/N))/N -> thvar_u_hat[cv_i]
}
((thvar_u_hat/thvar_Tr_hat) %>% mean() -> r_Tr_hat)
((thvar_u_hat/thvar_f_hat) %>% mean() -> r_f_hat)


######################
##The random foest method for classification
formula = as.formula(paste0("R_L~", paste(c(Xvar_c, Yvar, Avar), collapse = "+")))
fit_F <- randomForest(formula , data = dt2, proximity = T)
dt1_2$R_F <- predict(fit_F, newdata = dt1_2)
pr_F <- data.frame(R_F = 1:K, p_F = pi_m)
merge(dt1_2, pr_F, by="R_F") -> dt1_2_F

######################
##The random foest method for classification
# with biased influence function "phi_m" also as an predictor
formula = as.formula(paste0("R_L~", paste(c(Xvar_c, Yvar, Avar, "phi_m"), collapse = "+")))
fit_FM <- randomForest(formula , data = dt2, proximity = T)
dt1_2$R_FM <- predict(fit_FM, newdata = dt1_2)
pr_FM <- data.frame(R_FM = 1:K, p_FM = pi_m)
merge(dt1_2, pr_FM, by="R_FM") -> dt1_2_FM

#######################
##Tree method for classification
formula = as.formula(paste0("R_L~", paste(c(Xvar_c, Yvar, Avar), collapse = "+")))
fit_Tr <- rpart(formula, data = dt2, cp = 0.005)
dt1_2$R_Tr <- predict(fit_Tr, newdata = dt1_2, type = "class")
pr_Tr <- data.frame(R_Tr = 1:K, p_Tr = pi_m)
merge(dt1_2, pr_Tr, by="R_Tr") -> dt1_2_Tr

#######################
##Tree method for classification
# with biased influence function "phi_m" also as an predictor
formula = as.formula(paste0("R_L~", paste(c(Xvar_c, Yvar, Avar, "phi_m"), collapse = "+")))
fit_TrM <- rpart(formula, data = dt2, cp = 0.005)
dt1_2$R_TrM <- predict(fit_TrM, newdata = dt1_2, type = "class")
pr_TrM <- data.frame(R_TrM = 1:K, p_TrM = pi_m)
merge(dt1_2, pr_TrM, by="R_TrM") -> dt1_2_TrM

#############
sum((dt1_2$phi_a)^2/(n/N1_2))/N/N
sum((dt1_2_F$phi_a)^2/dt1_2_F$p_F)/N/N
sum((dt1_2_FM$phi_a)^2/dt1_2_FM$p_FM)/N/N
sum((dt1_2_Tr$phi_a)^2/dt1_2_Tr$p_Tr)/N/N
sum((dt1_2_TrM$phi_a)^2/dt1_2_TrM$p_TrM)/N/N
sum((dt1_2_O$phi_a)^2/dt1_2_O$p_O, na.rm = T)/N/N

((thvar_u_hat/thvar_Tr_hat) %>% mean() -> r_Tr_hat)
((thvar_u_hat/thvar_f_hat) %>% mean() -> r_f_hat)

pdf(file = paste0("classTree.pdf"))
rpart.plot(fit_Tr, type = 5,                  # 仅标记所有节点
           extra = 0,     # 节点颜色调色板
           cex = 0.7
)
dev.off()
plot(dt2$phi, dt2$phi_m)
cat(c(fit1$est - sqrt(fit1$ve)*1.96, fit1$est + sqrt(fit1$ve)*1.96))
cat(c(fit1_c$est - sqrt(fit1_c$ve)*1.96, fit1_c$est + sqrt(fit1_c$ve)*1.96))
 # 替换为您的公式和数据
 # 调整边距
