#################
# Plot the illustration density of psi under different sampling design.
rm(list = ls())
library(dplyr)
library(ggplot2)
library(randomForest)
library(MASS)
library(foreach)
require(doSNOW)
library(caret)
library(ggpubr)
source("simulationfunctions.R")
set.seed(1)

# N: the size of first phase data
# n_p: the size of pilot data
# n:  the size of other second phase data appart from pilot data.
# method: is the ATE estimation method
# alpha hyperparameter that controls the probability from becoming too small
method = "reg"
N = 30000
n_p = 1000
n = 3000
K = 5
alpha = 0.1
cat(" Method=", method, ", N=", N,", n_p=", n_p,", n=", n, ", K=", K,  ", alpha=", alpha,"\n")

# dt1 is the first phase data
# dt_p is the pilot data
# dt1_p is the first phase data remove pilot data
# method is the ATE estimation method
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

c(fit1$est - 2* sqrt(fit1$ve), fit1$est + 2* sqrt(fit1$ve))


phi_bar = mean(dt1$phi)

dt1$phi_a = abs(dt1$phi-phi_bar)

# stratification by the true influence function (unobserved)
# denoted as R_T = 1,...K.
dt1 <- dt1[order(dt1$phi_a), ];
R_T = rep(1:K, each = ceiling(N/K)) %>% as.factor()
dt1$R_T = R_T[1:N]

#####
# uniform sampling design (Simple random sampling design)
dt1$pi_u <- rep((n)/N, N)

#####
# Oracle sampling design (unobserved)
getpi(dt1$phi_a, (n), alpha)*(n) -> dt1$pi_b

############
#  the pilot data collection
pliot_ind <- sample(1:N, size = n_p, replace = F)
#### The interndata
dt_p <- dt1[pliot_ind,]
#### The otherdata
dt1_p <- dt1[-pliot_ind,]
N1_2 = nrow(dt1_p)


# The pilot data collection and influence functions estimation based on pilot data
x2 = cbind(dt_p$X1, dt_p$U1); y2 = dt_p$Y; A2 = dt_p$Z
esATE(x2, y2, A2, method) -> fit2
dt_p$phi_hat <- fit2$infl
phi_hat_bar <- mean(dt_p$phi_hat)
dt_p$phi_hat_a <- abs(dt_p$phi_hat - phi_hat_bar)

##############
#Lu's method
#Pre fixed stratified strategy (Outcome quantile) and its intra-stratum sampling probability 
te <- c()
for(k in 1:K){
  te[k] <- mean((dt_p[dt_p$R_O==k,"phi_hat"]-phi_hat_bar)^2)
}
te2 <- te
getpi_stra(sqrt(te2), N1_2, n, K, alpha) -> pi_O

pr_O <- data.frame(R_O = 1:K, p_O = pi_O)
merge(dt1_p, pr_O, by="R_O") -> dt1_p_O
merge(dt_p, pr_O, by="R_O") -> dt_p_O

###############
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
getpi_stra(sqrt(te2), N1_2, n, K, alpha) -> pi_m
dt_p$pi_m = NA
for(k in 1:K){
  dt_p[dt_p$R_L==k, "pi_m"] = pi_m[k]
}


# the classification problem: using pilot data to train a classification model.
# and predict the stratum of each individual in first phase data.
#The random foest method for classification: F method
fit_f <- randomForest(R_L~X1+Y+Z , data = dt_p, proximity = T)
dt1_p$R_F <- predict(fit_f, newdata = dt1_p)
pr_F <- data.frame(R_F = 1:K, p_F = pi_m)
merge(dt1_p, pr_F, by="R_F") -> dt1_p_F


############
# plot result
low = min(dt1_p_O$phi)-1
high = max(dt1_p_O$phi)+1
breaks = seq(low, high, by = (high-low)/50)
dt1_p_O$V_b = rbinom(length(dt1_p_O$pi_b), 1, prob = dt1_p_O$pi_b)
dt1_p_O$V_u = rbinom(length(dt1_p_O$pi_u), 1, prob = dt1_p_O$pi_u)
dt1_p_O$V_O = rbinom(length(dt1_p_O$p_O), 1, prob = dt1_p_O$p_O)
dt1_p_F$V_F = rbinom(length(dt1_p_F$p_F), 1, prob = dt1_p_F$p_F)

library(ggplot2)
theme_clean <- theme_minimal() +
  theme(
    plot.background = element_blank(),
    axis.line = element_line(colour = "gray")
  )
# 定义一组新的颜色
fill_color <- "#AED6F1"
hist_border_color <- "#3498DB"
density_line_color <- "IndianRed"

gg_1 <- ggplot(dt1_p_O[dt1_p_O$V_u==1,], aes(x=phi)) +
  geom_histogram(aes(y=..density..), color=hist_border_color, fill=fill_color, bins=length(breaks)) +
  geom_density(color=density_line_color, lwd=0.7) + 
  xlim(c(-100, 100)) +ylim(c(0, 0.027))+
  labs(title="SimRan", x = expression(psi(D, theta[0])), y = expression(paste("Density")))+theme_clean+
  annotate("text", x = -70, 
           y = 0.023, label = expression(Lambda[n]*(bold(pi)) == 0.114), 
           parse = TRUE)

gg_2 <- ggplot(dt1_p_O[dt1_p_O$V_b==1,], aes(x=phi)) +
  geom_histogram(aes(y=..density..), color=hist_border_color, fill=fill_color, bins=length(breaks)) +
  geom_density(color=density_line_color, lwd=0.7) + 
  xlim(c(-100, 100)) +ylim(c(0, 0.027))+
  labs(title="Oracle", x = expression(psi(D, theta[0])), y = expression(paste("Density")))+theme_clean+
  annotate("text", x = -70, 
           y = 0.023, label = expression(Lambda[n]*(bold(pi)) == 0.077), 
           parse = TRUE)

gg_3 <- ggplot(dt1_p_O[dt1_p_O$V_O==1,], aes(x=phi)) +
  geom_histogram(aes(y=..density..), color=hist_border_color, fill=fill_color, bins=length(breaks)) +
  geom_density(color=density_line_color, lwd=0.7) + 
  xlim(c(-100, 100)) + ylim(c(0, 0.027))+
  labs(title="FixStrat", x = expression(psi(D, theta[0])), y = expression(paste("Density")))+theme_clean+
  annotate("text", x = -70, 
           y = 0.023, label = expression(Lambda[n]*(bold(pi)) == 0.101), 
           parse = TRUE)

gg_4 <- ggplot(dt1_p_F[dt1_p_F$V_F==1,], aes(x=phi)) +
  geom_histogram(aes(y=..density..), color=hist_border_color, fill=fill_color, bins=length(breaks)) +
  geom_density(color=density_line_color, lwd=0.7) +
  xlim(c(-100, 100)) +ylim(c(0, 0.027))+
  labs(title="AdaStrat", x = expression(psi(D, theta[0])), y = expression(paste("Density")))+theme_clean+
  annotate("text", x = -70, 
           y = 0.023, label = expression(Lambda[n]*(bold(pi)) == 0.089), 
           parse = TRUE)

# 使用gridExtra包整合各个图形输出为一个2x2的图表
library(gridExtra)
pdf(file = paste0("density.pdf"),
    onefile = F, width = 8, height = 6)
grid.arrange(gg_1, gg_2, gg_3, gg_4, ncol=2)
dev.off()

