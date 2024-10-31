# Estimate Average Treatment Effect (esATE)
#
# x the matrix of confounders fully observed  (n x dim(x))
# y the vector of outcome (n x 1)
# A the vector of binary treatment (n x 1)
# w is the weight for each individual
# method the estimation method for ATE
# select 1 method from c("reg","ipw","aipw")
#
# "reg": a linear regression imputation estimator of the ATE
#
# "ipw": the inverse probability of treatment weighting estimator of the ATE, where the propensity score follows a logistic regression model
#
# "aipw": the augmented inverse probability of treatment weighting estimator (Lunceford Davidian, 2004) of the ATE
#
#
# est: estimate of the ATE
#
# ve:  variance estimate for est based on the asymptotic result
#
# infl: influence function value of each individual

# This function copy some code from \url{https://github.com/shuyang1987/IntegrativeCI/tree/master}

esATE_w<-function(x, y, A, w, method, N){
  ## check argument
  method<-method[which(method%in%c("reg","ipw","aipw")==1)]
  if(length(method)==0){return(cat("Error: specify a valid argument for method"))}
  
  n<-length(y)#sample size 
  x<-as.matrix(x)
  if(method=="reg"){
    lm.y<- y
    lm.x<- x
    lm.w<- w
    dimx<-dim(lm.x)[2]
    lm.out1<-lm(lm.y[which(A==1)]~lm.x[which(A==1),],
                weights = lm.w[which(A==1)])
    reg.mu1<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(A==0)]~lm.x[which(A==0),],
                weights = lm.w[which(A==0)])
    reg.mu0<-cbind(1,lm.x)%*%lm.out0$coefficients
    reg<-mean(w*reg.mu1-w*reg.mu0)/mean(w)
    B<-apply(cbind(1,lm.x),2, weighted.mean, w = w)
    B1mat<-cbind(A,lm.x*A)
    B1<-matrix(0,1+dimx,1+dimx)
    for(jj in 1:(1+dimx)){
      for(kk in jj:(1+dimx)){
        B1[jj,kk]<-B1[kk,jj]<-weighted.mean(B1mat[,jj]*B1mat[,kk], w)
      }
    }
    B0mat<-cbind(1-A,lm.x*(1-A))
    B0<-matrix(0,1+dimx,1+dimx)
    for(jj in 1:(1+dimx)){
      for(kk in jj:(1+dimx)){
        B0[jj,kk]<-B0[kk,jj]<-weighted.mean(B0mat[,jj]*B0mat[,kk], w)
      }
    }
    b1<-B%*%ginv(B1)
    b0<-B%*%ginv(B0)
    reg.resid1<-reg.resid0<-rep(0,n)
    reg.resid1[A==1]<-(lm.out1$resid)
    reg.resid0[A==0]<-(lm.out0$resid)
    reg.adj<-apply( B1mat*matrix(b1,n,1+dimx,byrow=TRUE),1,sum )*A*reg.resid1-
      apply( B0mat*matrix(b0,n,1+dimx,byrow=TRUE),1,sum )*(1-A)*reg.resid0
    infl.reg<-(( reg.mu1-reg.mu0)-reg+reg.adj)*w
    ve_reg<-sum(infl.reg^2)/N/N
    EST<-reg
    VE<-ve_reg
    INFL<-infl.reg 
  }
  if(method=="ipw"){
    Rx<-x
    Rx<-as.matrix(Rx)
    dimx<-dim(Rx)[2]
    glm.out<-glm(A~Rx,family="binomial", weights = w)
    pshat<-glm.out$fitted.values
    psw<-mean(w*y*A/pshat-w*y*(1-A)/(1-pshat))/mean(w)
    H<-c(mean(w*(A*y*(1-pshat)/pshat+(1-A)*y*(pshat)/(1-pshat))),
         apply(w*(A*y*(1-pshat)/pshat+(1-A)*y*(pshat)/(1-pshat))*Rx , 2 , mean))/mean(w)
    
    E<-matrix(0,1+dimx,1+dimx)
    forE<-matrix(1,1+dimx,n)
    forE[2:(1+dimx),]<-t(Rx)
    for(jj in 1:(1+dimx)){
      for(kk in jj:(1+dimx)){
        E[jj,kk]<-E[kk,jj]<-mean(pshat*(1-pshat)*forE[jj,]*forE[kk,]*w)/mean(w)
      }
    }
    alpha<- -(A-pshat) * H%*%ginv(E)%*%forE
    infl.ipw<-rep(0,n);
    psi_psw<-y*A/pshat-y*(1-A)/(1-pshat)
    infl.ipw<- ( psi_psw+alpha-psw )*w
    ve_psw<-sum(infl.ipw^2)/N/N
    EST<-psw
    VE<-ve_psw
    INFL<-infl.ipw
  }
  if(method=="aipw"){
    lm.y<- y
    lm.x<- x
    lm.w<- w
    lm.out1<-lm(lm.y[which(A==1)]~lm.x[which(A==1),],
                weights = lm.w[which(A==1)])
    mu1<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(A==0)]~lm.x[which(A==0),],
                weights = lm.w[which(A==0)])
    mu0<-cbind(1,lm.x)%*%lm.out0$coefficients
    Rx<-x
    glm.out<-glm(A~Rx,family="binomial", weights = w)
    pshat<-glm.out$fitted.values
    psi.aug<-y*A/pshat-y*(1-A)/(1-pshat)+( mu1*(1-A/pshat)-mu0*(1-(1-A)/(1-pshat)))
    aipw<-mean(w*psi.aug)/mean(w)
    EST<-aipw
    infl.aug<-rep(0,n);
    infl.aug<- (psi.aug-aipw)*w
    VE <- ve_aipw <- sum(infl.aug^2)/N/N
  }
  return(c(est = EST, ve = VE))
}

