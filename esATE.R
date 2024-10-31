# Estimate Average Treatment Effect (esATE)
#
# x the matrix of confounders fully observed  (n x dim(x))
# y the vector of outcome (n x 1)
# A the vector of binary treatment (n x 1)
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

esATE<-function(x, y, A, method){
  ## check argument
  method<-method[which(method%in%c("reg","ipw","aipw")==1)]
  if(length(method)==0){return(cat("Error: specify a valid argument for method"))}
  
  n<-length(y)#sample size 
  x<-as.matrix(x)
  if(method=="reg"){
    lm.y<- y
    lm.x<- x
    lm.x<- as.matrix(lm.x)
    dimx<-dim(lm.x)[2]
    lm.out1<-lm(lm.y[which(A==1)]~lm.x[which(A==1),])
    reg.mu1<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(A==0)]~lm.x[which(A==0),])
    reg.mu0<-cbind(1,lm.x)%*%lm.out0$coefficients
    reg<-mean(reg.mu1-reg.mu0)
    B<-apply(cbind(1,lm.x),2,mean)
    B1mat<-cbind(A,lm.x*A)
    B1<-matrix(0,1+dimx,1+dimx)
    for(jj in 1:(1+dimx)){
      for(kk in jj:(1+dimx)){
        B1[jj,kk]<-B1[kk,jj]<-mean(B1mat[,jj]*B1mat[,kk])
      }
    }
    B0mat<-cbind(1-A,lm.x*(1-A))
    B0<-matrix(0,1+dimx,1+dimx)
    for(jj in 1:(1+dimx)){
      for(kk in jj:(1+dimx)){
        B0[jj,kk]<-B0[kk,jj]<-mean(B0mat[,jj]*B0mat[,kk])
      }
    }
    b1<-B%*%ginv(B1)
    b0<-B%*%ginv(B0)
    reg.resid1<-reg.resid0<-rep(0,n)
    reg.resid1[A==1]<-(lm.out1$resid)
    reg.resid0[A==0]<-(lm.out0$resid)
    reg.adj<-apply( B1mat*matrix(b1,n,1+dimx,byrow=TRUE),1,sum )*A*reg.resid1-
      apply( B0mat*matrix(b0,n,1+dimx,byrow=TRUE),1,sum )*(1-A)*reg.resid0
    infl.reg<-(( reg.mu1-reg.mu0)+reg.adj-mean(( reg.mu1-reg.mu0)))
    ve_reg<-mean(infl.reg^2)/n
    EST<-reg
    VE<-ve_reg
    INFL<-infl.reg 
  }
  if(method=="ipw"){
    Rx<-x
    Rx<-as.matrix(Rx)
    dimx<-dim(Rx)[2]
    glm.out<-glm(A~Rx,family="binomial")
    pshat<-glm.out$fitted.values
    psw<-mean(y*A/pshat)-mean(y*(1-A)/(1-pshat))
    H<-c(mean( A*y*(1-pshat)/pshat+(1-A)*y*(pshat)/(1-pshat) ),
         apply( (A*y*(1-pshat)/pshat+(1-A)*y*(pshat)/(1-pshat))*Rx , 2 , mean ) )
    
    E<-matrix(0,1+dimx,1+dimx)
    forE<-matrix(1,1+dimx,n)
    forE[2:(1+dimx),]<-t(Rx)
    for(jj in 1:(1+dimx)){
      for(kk in jj:(1+dimx)){
        E[jj,kk]<-E[kk,jj]<-mean(pshat*(1-pshat)*forE[jj,]*forE[kk,])
      }
    }
    alpha<- -(A-pshat) * H%*%ginv(E)%*%forE
    infl.ipw<-rep(0,n);
    psi_psw<-y*A/pshat-y*(1-A)/(1-pshat)
    infl.ipw<- ( psi_psw+alpha-psw )
    ve_psw<-mean(infl.ipw^2)/n
    EST<-psw
    VE<-ve_psw
    INFL<-infl.ipw
  }
  if(method=="aipw"){
    lm.y<- y
    lm.x<- x
    lm.out1<-lm(lm.y[which(A==1)]~lm.x[which(A==1),])
    mu1<-cbind(1,lm.x)%*%lm.out1$coefficients
    lm.out0<-lm(lm.y[which(A==0)]~lm.x[which(A==0),])
    mu0<-cbind(1,lm.x)%*%lm.out0$coefficients
    Rx<-x
    glm.out<-glm(A~Rx,family="binomial")
    pshat<-glm.out$fitted.values
    psi.aug<-y*A/pshat-y*(1-A)/(1-pshat)+( mu1*(1-A/pshat)-mu0*(1-(1-A)/(1-pshat)) )
    aipw<-mean(psi.aug)
    infl.aug<-rep(0,n);
    infl.aug<- (psi.aug-aipw)
    ve_aipw<- mean(infl.aug^2)/n
    EST<-aipw
    VE<-ve_aipw
    INFL<-infl.aug
  }
  return(list(est=EST, ve=VE, infl = as.numeric(INFL)))
}


