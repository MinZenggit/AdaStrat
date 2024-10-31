## Plot results 
rm(list = ls())
library(dplyr)
library(ggplot2)
name = "01-01simulation"
if(! dir.exists(paste0("result_table/",name))){
  dir.create(paste0("result_table/",name))
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



library(stringr)
isalpha <- function(str){
  str <- ifelse(grepl("[a-zA-Z]", str), NA, str) %>% unlist %>% as.numeric()
  return(str)
}
result_td <- function(re){
  re = cbind((apply(re, 2, mean, na.rm = T)-0.50)/0.50 *100,
             apply(re, 2, var, na.rm = T),
             (apply(re, 2, mean, na.rm = T)-0.50)^2 +apply(re, 2, var, na.rm = T) )
  colnames(re) <- c("bias", "var", "mse")
  return(re)
}




###################
if(! dir.exists(paste0("result_figure/",name))){
  dir.create(paste0("result_figure/",name))
}
r_td <- c()
for(j in 1:nrow(settings)){
  method = settings$method[j]
  N = settings$N[j]
  n_p = settings$n_p[j]
  n = settings$n[j]
  K = settings$K[j]
  alpha = settings$alpha[j]
  rrr <- c()
  path = paste0("result/", name, "/", paste(N, n_p, n, K, alpha, method, sep = "_"), ".rdata")
  # print(path)
  load(file = path)
  apply(result, 1, isalpha)%>%t() -> result
  result_td(result[,  c(1,3,5,7,9)]) %>% as.data.frame() -> rr
  cbind(rr,
        var_pot=colMeans(result[,c(1,3,5,7,9)+1], na.rm = T),
        method = rep(method, 5),
        N = rep(N, 5),
        n_p = rep(n_p,5),
        n = rep(n, 5),
        K = rep(K, 5),
        alpha = rep(alpha, 5),
        des = c("u", "b", "z", "f", "o"))->rrr
  colnames(rrr) <- c("Bias", "Var", "MSE",
                     "Var_est",
                     "Sampling design", "N", "n_p", "n", "K", "alpha","des")
  rbind(r_td, rrr) -> r_td
}







########################
# Plot 
#######################
method_i = "ipw"; n_p_i = 500; K_i = 2; alpha_i = 0.3
r_td %>% filter(method == method_i,n_p == n_p_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "f")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                        labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.05, 0.40))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design",
    x = "Second phase sample size", 
    y = expression("MSE of "~hat(tau)))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(a)")-> p1

method_i = "ipw"; n_p_i = 500; K_i = 4; alpha_i = 0.3
r_td %>% filter(method == method_i,n_p == n_p_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "f")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                        labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.05, 0.40))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design", 
    x = "Second phase sample size", 
    y = expression("MSE of "~hat(tau)))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(b)")-> p2


method_i = "ipw"; n_p_i = 500; K_i = 6; alpha_i = 0.3
r_td %>% filter(method == method_i,n_p == n_p_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "f")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                        labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.05, 0.40))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design", 
       x = "Second phase sample size", 
       y = expression("MSE of "~hat(tau)))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(c)")-> p3

method_i = "ipw"; n_p_i = 500; K_i = 8; alpha_i = 0.3
r_td %>% filter(method == method_i,n_p == n_p_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "f")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                        labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.05, 0.40))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design",
       x = "Second phase sample size", 
       y = expression("MSE of "~hat(tau)))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(d)")-> p4


method_i = "ipw"; n_p_i = 500; K_i = 10; alpha_i = 0.3
r_td %>% filter(method == method_i,n_p == n_p_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "f")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                        labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.05, 0.40))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design", 
       x = "Second phase sample size", 
       y = expression("MSE of "~hat(tau)))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(e)")-> p5

method_i = "ipw"; n_p_i = 1000; K_i = 2; alpha_i = 0.3
r_td %>% filter(method == method_i,n_p == n_p_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "f")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                        labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.05, 0.30))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design", 
    x = "Second phase sample size", 
    y = expression("MSE of "~hat(tau)))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(f)") -> p6

method_i = "ipw"; n_p_i = 1000; K_i = 4; alpha_i = 0.3
r_td %>% filter(method == method_i,n_p == n_p_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "f")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                        labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.05, 0.30))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design", 
    x = "Second phase sample size", 
    y = expression("MSE of "~hat(tau)))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(g)") -> p7

method_i = "ipw"; n_p_i = 1000; K_i = 6; alpha_i = 0.3
r_td %>% filter(method == method_i,n_p == n_p_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "f")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                        labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.05, 0.30))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design", 
       x = "Second phase sample size", 
       y = expression("MSE of "~hat(tau)))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(h)") -> p8

method_i = "ipw"; n_p_i = 1000; K_i = 8; alpha_i = 0.3
r_td %>% filter(method == method_i,n_p == n_p_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "f")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                        labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.05, 0.30))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design", 
       x = "Second phase sample size", 
       y = expression("MSE of "~hat(tau)))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(i)") -> p9


method_i = "ipw"; n_p_i = 1000; K_i = 10; alpha_i = 0.3
r_td %>% filter(method == method_i,n_p == n_p_i, K==K_i, alpha == alpha_i,
                des %in%c("u", "b", "o", "f")) %>%
  ggplot(aes(x = n, y = MSE, color = des, shape = des, linetype = des)) +
  geom_line()+
  geom_point(size = 3) +
  scale_color_manual(values = c("u" = "#66c2a5", "b" = "#fc8d62", "f" = "#8da0cb", "o" = "#e78ac3"),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_shape_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat")) +
  scale_linetype_manual(values = c("u" = 15, "b" = 16, "f" = 17, "o" = 18),
                     labels = c("u" ="USD"  , "b" = "Oracle" , "f" ="AdaStrat" , "o" ="FixStrat"))+
  ylim(c(0.05, 0.30))+
  labs(color = "Sampling design", shape = "Sampling design", linetype = "Sampling design", 
       x = "Second phase sample size", 
       y = expression("MSE of "~hat(tau)))+
  theme_minimal()+
  theme(axis.line = element_line(color = "gray"))+
  ggtitle("(j)")-> p10






plots <- list(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10)
library(ggpubr)
pdf(file = paste0("result_figure/", name, "/", method_i, "_",alpha_i, ".pdf"),
    onefile = F, width = 12.5, height = 5)
ggarrange(p1, p2, p3, p4, p5, p6, p7, p8, p9, p10, 
          nrow = 2, ncol = 5, common.legend = TRUE)
dev.off()
# 使用 annotate_figure 添加字母标签


#############################