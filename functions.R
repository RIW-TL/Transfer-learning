# It includes all functions used in our paper, together with functions for other competitors such as the Trans-Lasso.

library(ncvreg)
library(mvtnorm)
library(glmnet)
library(ks)
library(bayesdistreg)
#------------------------------------------------------------
# coefficients generation
#===================================================

# s: sparsity level
# d: informative level
# p: dimension
# q: the distance between beta0 and beta_k for k notin size.A0
# K: number of sources
# size.A0: number of possible informative sources
# sig.beta: true significant coefficient
# sig.delta1: positive difference
# sig.delta2: negative difference
# exact: fixed (exact = T) / random (exact = T) setting


Coef.gen = function(s,d,q,size.A0, K, sig.beta,
                    sig.delta1, sig.delta2, p, exact){
  
  # target parameter
  beta0 = c(rep(sig.beta,s), rep(0, p - s))
  
  # source parameters
  W = rep.col(beta0, K)
  for(k in 1:K){
    
    if(k <= size.A0){
      if(exact){
        samp0<- sample((s + 1):p, d, replace=F)
        W[samp0,k] <-W[samp0,k] + rep(-sig.delta1, d)
      }else{
        samp0<- sample((s + 1):p, d, replace=F)
        W[samp0,k] <-W[samp0,k] + runif(d,0,1)
      }
    }else{
      if(exact){
        
        W[1:s,k] <- W[1:s,k] + rep(-sig.beta,s)
        samp1 <- sample((s + 1):p, q, replace = F)
        W[samp1,k] <-W[samp1,k] + rep(-sig.delta2, q)
      }else{
        W[1:s,k] <- W[1:s,k] + rep(-sig.beta,s)
        samp1 <- sample((s + 1):p, q, replace = F)
        W[samp1,k] <- W[samp1,k] + runif(q, 0, 1)
      }
    }
  }
  
  return(list(W=W, beta0=beta0))
}

# data splitting in RIW-TL & RIW-TL-U & RIW-TL-P
#==================================================

# {X01,Y01} & {X02,Y02} the splitted target data
# {X11,Y11} & {X12,Y12} & {X13,Y13} the splitted source data
# A,M: tuning parameters used in sample selection
# TU = A + M
# beta0.init & betak.init: initial estimators

fit.RIW.DS = function(X01,Y01,X02,Y02,X11,Y11,X12,Y12,X13,Y13,
                      A,M,TU){
  
  p = ncol(X01)
  K = length(X11)
  
  # Initial estimators (SCAD penalty as example)
  #---------------------------
  fit = cv.ncvreg(X01,Y01,family = "gaussian",penalty = "SCAD")
  beta0.initial.hat = coef(fit)[-1]
  
  beta1.initial.hat = matrix(0,p,K)
  for (k in 1:K) {

    # fit = cv.ncvreg(X11[[k]],Y11[[k]],family = "gaussian",penalty = "SCAD")
    fit = cv.glmnet(X11[[k]],Y11[[k]],family = "gaussian",alpha = 1)
    beta1.initial.hat[,k] = coef(fit)[-1]
  }
  
  # KDE and weighted sources data
  #-------------------------------------
  Ydata_star = Xdata_star = eta = list()
  Ydata_star_U = Xdata_star_U = list()
  error01.hat = error11.hat = list()
  
  for (k in 1:K) {
    
    error112.hat = Y12[[k]] - X12[[k]] %*% beta1.initial.hat[,k]
    kde.fit = kde(x = error112.hat[,1])
    
    error01.hat[[k]] = Y13[[k]] - X13[[k]] %*% beta0.initial.hat
    error11.hat[[k]] = Y13[[k]] - X13[[k]] %*% beta1.initial.hat[,k]
    
    # RIW-TL
    fenzi = (predict(kde.fit,x = error01.hat[[k]]) + 
               predict(kde.fit,x = -error01.hat[[k]]))/2
    fenmu = predict(kde.fit,x = error11.hat[[k]])
    weight = fenzi/fenmu;weight[which(fenzi < 0 | fenmu <0)] = 0
    
    Ydata_star[[k]] = sqrt(weight)*Y13[[k]]
    Xdata_star_matrix = matrix(0,nrow(X13[[k]]),p)
    for (i in 1:nrow(X13[[k]])) {
      
      Xdata_star_matrix[i,] = sqrt(weight[i])*X13[[k]][i,]
    }
    Xdata_star[[k]] = Xdata_star_matrix
    
    
    # RIW-TL-U
    fenzi = indicat(abs(error01.hat[[k]]),TU)/(2*TU)
    fenmu = predict(kde.fit,x = error11.hat[[k]])
    weight_U = fenzi/fenmu
    weight_U[which(fenzi < 0 | fenmu <0)] = 0
    
    Ydata_star_U[[k]] = sqrt(weight_U)*Y13[[k]]
    Xdata_star_matrix = matrix(0,nrow(X13[[k]]),p)
    for (i in 1:nrow(X13[[k]])) {
      
      Xdata_star_matrix[i,] = sqrt(weight_U[i])*X13[[k]][i,]
    }
    Xdata_star_U[[k]] = Xdata_star_matrix
    
    # contrast vectors adjusted by predictors
    eta[[k]] = X13[[k]] %*% (beta1.initial.hat[,k] - beta0.initial.hat)
  }
  
  
  # sample selection
  #-------------------------------------
  YY = Y_P = Y_U = Y02;XX = X_P = X_U = X02
  for (k in 1:K) {
    
    # RIW-TL
    id1 = which((abs(error11.hat[[k]]) <= A) & (abs(eta[[k]]) <= M))
    YY = c(YY,Ydata_star[[k]][id1])
    XX = rbind(XX,Xdata_star[[k]][id1,])
    
    # RIW-TL-U
    Y_U = c(Y_U,Ydata_star_U[[k]][id1])
    X_U = rbind(X_U,Xdata_star_U[[k]][id1,])
  }
  
  # weighted lasso optimization
  const = 0.5
  fit_F = glmnet(XX,YY,family = "gaussian",alpha = 1,
                 lambda = const*sqrt(log(p)/length(YY)))
  fit_F_U = glmnet(X_U,Y_U,family = "gaussian",alpha = 1,
                   lambda = const*sqrt(log(p)/length(Y_U)))
  
  beta_RIW_TL = coef(fit_F)[-1]
  beta_RIW_TL_U = coef(fit_F_U)[-1]
  
  list = list(beta_RIW_TL = beta_RIW_TL,
              beta_RIW_TL_U = beta_RIW_TL_U,
              SU = length(YY),
              SU_U = length(Y_U))
}


# Tuning parameter selection via cross-validation 
#================================================

# X : vector, covariate data
# y : vector, response
# n.vec : vector a sample size of target and sources

fit.RIW = function(X,y,n.vec,K,A,M,TU){
  
  # Target data
  n0 = n.vec[1]
  Xdata0 = X[1:n0,]
  Ydata0 = y[1:n0]
  
  # source data
  Ydata1 = Xdata1 = list()
  for (k in 1:K) {
    
    id.k = ind.set(n.vec,k + 1)
    Xdata1[[k]] = X[id.k,]
    Ydata1[[k]] = y[id.k]
  }
  
  # target data splitting
  X01 = Y01 = X02 = Y02 = list()
  for (j in 1:3) {
    
    id = ((j - 1)*floor(n0/3) + 1):(j*floor(n0/3))
    X01[[j]] = Xdata0[-id,]
    Y01[[j]] = Ydata0[-id]
    
    X02[[j]] = Xdata0[id,]
    Y02[[j]] = Ydata0[id]
  }
  
  # source data splitting
  X.P1 = Y.P1 = X.P2 = Y.P2 = X.P3 = Y.P3 = list()
  for (k in 1:K) {
    
    X.P1[[k]] = Xdata1[[k]][1:floor(n.vec[k+1]/3),]
    X.P2[[k]] = Xdata1[[k]][(floor(n.vec[k+1]/3) + 1):floor(2*n.vec[k+1]/3),]
    X.P3[[k]] = Xdata1[[k]][(floor(2*n.vec[k+1]/3) + 1):n.vec[k+1],]
    
    Y.P1[[k]] = Ydata1[[k]][1:floor(n.vec[k+1]/3)]
    Y.P2[[k]] = Ydata1[[k]][(floor(n.vec[k+1]/3) + 1):floor(2*n.vec[k+1]/3)]
    Y.P3[[k]] = Ydata1[[k]][(floor(2*n.vec[k+1]/3) + 1):n.vec[k+1]]
  }
  
  
  fit1 = fit.RIW.DS(X01 = X01[[1]],Y01 = Y01[[1]],
                    X02 = X02[[1]],Y02 = Y02[[1]],
                    X11 = X.P1,Y11 = Y.P1,
                    X12 = X.P2,Y12 = Y.P2,
                    X13 = X.P3,Y13 = Y.P3,
                    A,M,TU)
  
  fit2 = fit.RIW.DS(X01 = X01[[2]],Y01 = Y01[[2]],
                    X02 = X02[[2]],Y02 = Y02[[2]],
                    X11 = X.P2,Y11 = Y.P2,
                    X12 = X.P3,Y12 = Y.P3,
                    X13 = X.P1,Y13 = Y.P1,
                    A,M,TU)
  
  fit3 = fit.RIW.DS(X01 = X01[[3]],Y01 = Y01[[3]],
                    X02 = X02[[3]],Y02 = Y02[[3]],
                    X11 = X.P3,Y11 = Y.P3,
                    X12 = X.P1,Y12 = Y.P1,
                    X13 = X.P2,Y13 = Y.P2,
                    A,M,TU)
  
  
  beta_RIW_TL = (fit1$beta_RIW_TL + fit2$beta_RIW_TL + fit3$beta_RIW_TL)/3
  beta_RIW_TL_U = (fit1$beta_RIW_TL_U + fit2$beta_RIW_TL_U + fit3$beta_RIW_TL_U)/3
  
  SUR = (fit1$SU + fit2$SU + fit3$SU)/length(y)
  SUR_U = (fit1$SU_U + fit2$SU_U + fit3$SU_U)/length(y)
  
  list = list(beta_RIW_TL = beta_RIW_TL,
              beta_RIW_TL_U = beta_RIW_TL_U,
              SUR = SUR,
              SUR_U = SUR_U)
}

hyper.select.cv = function(X,y,n.vec,K,fold = 5){
  
  # Target data
  n0 = n.vec[1]
  Xdata0 = X[1:n0,]
  Ydata0 = y[1:n0]
  
  # the range of parameters (can be adjusted case-by-case)
  A.range = seq(0.5,1.5,0.5)
  M.range = seq(1,4,1)
  cv.error = cv.error.u = matrix(0,length(A.range),length(M.range))
  
  for (a in seq_along(A.range)) {
    for (m in seq_along(M.range)) {
      
      A = A.range[a]
      M = M.range[m]
      TU = A + M
      
      error.sum = error.sum.u = 0
      for (f in 1:fold) {
        
        # validation set
        val.id = ((f - 1)*floor(n0/fold) + 1):(f*floor(n0/fold))
        X.val = X[val.id,]
        y.val = y[val.id]
        
        # training set
        X.use = X[-val.id,]
        y.use = y[-val.id]
        n.vec[1] = n0 - length(val.id)
        
        fit = fit.RIW(X.use,y.use,n.vec,K,A,M,TU)
        beta_RIW_TL = fit$beta_RIW_TL
        beta_RIW_TL_U = fit$beta_RIW_TL_U
        
        # compute the cv loss
        PE.val = sum((y.val - X.val %*% beta_RIW_TL)^2)
        error.sum = error.sum + PE.val  
        
        PE.val = sum((y.val - X.val %*% beta_RIW_TL_U)^2)
        error.sum.u = error.sum.u + PE.val  
      }
      
      cv.error[a,m] = error.sum
      cv.error.u[a,m] = error.sum.u
    }
  }
  
  min.id = which(cv.error == min(cv.error),arr.ind = TRUE)
  min.id.u = which(cv.error.u == min(cv.error.u),arr.ind = TRUE)
  
  return(list(A.opt = A.range[min.id[1,1]],
              M.opt = M.range[min.id[1,2]],
              A.opt.u = A.range[min.id.u[1,1]],
              M.opt.u = M.range[min.id.u[1,2]]))
}

# Function for Trans-Lasso
#===============================================

# aggregation function
agg.fun = function(B, X.test,y.test, total.step=10, selection=F){
  if(sum(B==0)==ncol(B)*nrow(B)){
    return(rep(0,nrow(B)))
  }
  p<-nrow(B)
  K<-ncol(B)
  colnames(B)<-NULL
  if(selection){#select beta.hat with smallest prediction error
    khat<-which.min(colSums((y.test-X.test%*%B)^2))
    theta.hat<-rep(0, ncol(B))
    theta.hat[khat] <- 1
    beta=B[,khat]
    beta.ew=NULL
  }else{#Q-aggregation
    theta.hat<- exp(-colSums((y.test-X.test%*%B)^2)/2)
    theta.hat=theta.hat/sum(theta.hat)
    theta.old=theta.hat
    beta<-as.numeric(B%*%theta.hat)
    beta.ew<-beta
    # theta.old=theta.hat
    for(ss in 1:total.step){
      theta.hat<- exp(-colSums((y.test-X.test%*%B)^2)/2+colSums((as.vector(X.test%*%beta)-X.test%*%B)^2)/8)
      theta.hat<-theta.hat/sum(theta.hat)
      beta<- as.numeric(B%*%theta.hat*1/4+3/4*beta)
      if(sum(abs(theta.hat-theta.old))<10^(-3)){break}
      theta.old=theta.hat
    }
  }
  list(theta=theta.hat, beta=beta, beta.ew=beta.ew)
}


# oracle Trans-Lasso
las.kA = function(X, y, A0, n.vec, lam.const=NULL, l1=T){
  p<-ncol(X)
  size.A0<- length(A0)
  if(size.A0 > 0){
    ind.kA<- ind.set(n.vec, c(1, A0+1))
    ind.1<-1:n.vec[1]
    if(l1){
      y.A<-y[ind.kA]
    }else{ #the l0-method
      y.A<- y[ind.1]
      Sig.hat<-t(X)%*%X/nrow(X)
      for(k in 1:size.A0){
        ind.k<- ind.set(n.vec,k+1)
        lam.k <- sqrt(mean(y[ind.1]^2)/n.vec[1]+mean(y[ind.k]^2)/n.vec[k]) * sqrt(2*log(p))
        delta.hat.k<-lassoshooting(XtX=Sig.hat, 
                                   Xty=t(X[ind.k,])%*%y[ind.k]/n.vec[k+1]-t(X[1:n.vec[1],])%*%y[1:n.vec[1]]/n.vec[1],
                                   lambda=lam.k)$coef
        y.A<-c(y.A, y[ind.k]-X[ind.k,]%*%delta.hat.k)
      }
    }
    if(is.null(lam.const)){
      cv.init<-cv.glmnet(X[ind.kA,], y.A, nfolds=8, lambda=seq(1,0.1,length.out=10)*sqrt(2*log(p)/length(ind.kA)))
      lam.const <- cv.init$lambda.min/sqrt(2*log(p)/length(ind.kA))
    }
    w.kA <- as.numeric(glmnet(X[ind.kA,], y.A, lambda=lam.const*sqrt(2*log(p)/length(ind.kA)))$beta)
    w.kA<-w.kA*(abs(w.kA)>=lam.const*sqrt(2*log(p)/length(ind.kA)))
    
    # cv.delta<-cv.glmnet(x=X[ind.1,],y=y[ind.1]-X[ind.1,]%*%w.kA, lambda=seq(1,0.1,length.out=10)*sqrt(2*log(p)/length(ind.1)))
    # delta.kA<-predict(cv.delta, s='lambda.min', type='coefficients')[-1]
    delta.kA <- as.numeric(glmnet(x=X[ind.1,],y=y[ind.1]-X[ind.1,]%*%w.kA, lambda=lam.const*sqrt(2*log(p)/length(ind.1)))$beta)
    delta.kA<-delta.kA*(abs(delta.kA)>=lam.const*sqrt(2*log(p)/length(ind.1)))
    beta.kA <- w.kA + delta.kA
    lam.const=NA
  }else{
    cv.init<-cv.glmnet(X[1:n.vec[1],], y[1:n.vec[1]], nfolds=8, lambda=seq(1,0.1,length.out=20)*sqrt(2*log(p)/n.vec[1]))
    lam.const<-cv.init$lambda.min/sqrt(2*log(p)/n.vec[1])
    beta.kA <- predict(cv.init, s='lambda.min', type='coefficients')[-1]
    w.kA<-NA
  }
  list(beta.kA=as.numeric(beta.kA),w.kA=w.kA, lam.const=lam.const)
  
}

# Trans Lasso method
Trans.lasso = function(X, y, n.vec, I.til, l1=T){
  M= length(n.vec)-1
  #step 1
  X0.til<-X[I.til,] #used for aggregation
  y0.til<-y[I.til]
  X<- X[-I.til,]
  y<-y[-I.til]
  #step 2
  Rhat <- rep(0, M+1)
  p<- ncol(X)
  n.vec[1]<- n.vec[1]-length(I.til)
  ind.1<-ind.set(n.vec,1)
  for(k in 2: (M+1)){
    ind.k<-ind.set(n.vec,k)
    Xty.k <- t(X[ind.k,])%*%y[ind.k]/n.vec[k] - t(X[ind.1,])%*%y[ind.1]/ n.vec[1]
    margin.T<-sort(abs(Xty.k[,1]),decreasing=TRUE)[1:round(n.vec[1]/3)]
    Rhat[k] <-  sum(margin.T^2)
  }
  Tset<- list()
  k0=0
  kk.list<-unique(rank(Rhat[-1]))
  #cat(rank(Rhat[-1]),'\n')
  for(kk in 1:length(kk.list)){#use Rhat as the selection rule
    Tset[[k0+kk]]<- which(rank(Rhat[-1]) <= kk.list[kk])
  }
  k0=length(Tset)
  Tset<- unique(Tset)
  #cat(length(Tset),'\n')
  
  beta.T<-list()
  init.re<-las.kA(X=X, y=y, A0=NULL, n.vec=n.vec, l1=l1)
  beta.T[[1]] <- init.re$beta.kA
  beta.pool.T<-beta.T ##another method for comparison
  for(kk in 1:length(Tset)){#use pi.hat as selection rule
    T.k <- Tset[[kk]]
    re.k<- las.kA(X=X, y=y, A0=T.k, n.vec=n.vec, l1=l1, lam.const=init.re$lam.const)
    beta.T[[kk+1]] <-re.k$beta.kA
    beta.pool.T[[kk+1]]<-re.k$w.kA
  }
  beta.T<-beta.T[!duplicated((beta.T))]
  beta.T<- as.matrix(as.data.frame(beta.T))
  agg.re1 <- agg.fun(B=beta.T, X.test=X0.til, y.test=y0.til)
  beta.pool.T<-beta.pool.T[!duplicated((beta.pool.T))]
  beta.pool.T<- as.matrix(as.data.frame(beta.pool.T))
  agg.re2<-agg.fun(B=beta.pool.T, X.test=X0.til, y.test=y0.til)
  
  return(list(beta.hat=agg.re1$beta, theta.hat=agg.re1$theta, rank.pi=rank(Rhat[-1]),
              beta.pool=agg.re2$beta, theta.pool=agg.re2$theta))
}


#=====================other functions============================
# data id
ind.set = function(n.vec, k.vec){ 
  ind.re <- NULL
  for(k in k.vec){
    
    if(k==1){
      ind.re<-c(ind.re,1: n.vec[1])
    }else{
      ind.re<- c(ind.re, (sum(n.vec[1:(k-1)])+1): sum(n.vec[1:k]))
    }
    
  }
  ind.re
}

# column matrix
rep.col = function(x,n){
  matrix(rep(x,each=n), ncol=n, byrow=TRUE)
}





