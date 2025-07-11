library(ncvreg)
library(mvtnorm)
library(glmnet)
library(ks)
library(bayesdistreg)
#------------------------------------------------------------

# X : all covariates information including target and sources
# y : all responses  information including target and sources
# n.vec : vector a sample size of target and sources
# I.til : id in data splitting
# k.vec: source id
#------------------------------------------------------------

# size.A0 : number of informative sources
# sig.beta : true significant coefficient
# sig.delta1 : positive difference
# sig.delta2 : negative difference
# exact: setting fixed / random
#------------------------------------------------------------
# coefficients generation
Coef.gen = function(s,d,q, size.A0, K, sig.beta,
                    sig.delta1, sig.delta2, p, exact){
  
  beta0 = c(rep(sig.beta,s), rep(0, p - s))
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
fit.RIW.DS = function(X01,Y01,X11,Y11,X12,Y12,X13,Y13,
                      A,M,TU,K,beta0.initial){
  
  # Initial estimators
  #---------------------------
  p = ncol(X01)
  beta0.initial.hat = beta0.initial
  beta1.initial.hat = matrix(0,p,K);sigma = rep(0,K)
  for (k in 1:K) {
    
    fit = cv.ncvreg(X11[[k]],Y11[[k]],family = "gaussian",penalty = "SCAD")
    beta1.initial.hat[,k] = coef(fit)[-1]
    sigma[k] = sum((Y11[[k]] - X11[[k]] %*% matrix(beta1.initial.hat[,k]))^2)/
      (length(Y11[[k]]) - length(which(beta1.initial.hat[,k] != 0)))
  }
  
  ## sources
  #-------------------------------------
  Ydata_star = Xdata_star = eta = list()
  Ydata_star_P = Xdata_star_P = list()
  Ydata_star_U = Xdata_star_U = list()
  error01.hat = error11.hat = list()
  
  for (k in 1:K) {
    
    error01.hat[[k]] = Y13[[k]] - X13[[k]] %*% beta0.initial.hat
    error112.hat = Y12[[k]] - X12[[k]] %*% beta1.initial.hat[,k]
    
    kde.fit = kde(x = error112.hat[,1])
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
    
    # RIW-TL-P
    fenzi = dnorm(error01.hat[[k]],0,sigma[k])
    fenmu = dnorm(error11.hat[[k]],0,sigma[k])
    weight_P = fenzi/fenmu
    
    weight_P[which(fenzi < 0 | fenmu <0)] = 0
    
    Ydata_star_P[[k]] = sqrt(weight_P)*Y13[[k]]
    Xdata_star_matrix = matrix(0,nrow(X13[[k]]),p)
    for (i in 1:nrow(X13[[k]])) {
      
      Xdata_star_matrix[i,] = sqrt(weight_P[i])*X13[[k]][i,]
    }
    Xdata_star_P[[k]] = Xdata_star_matrix
    
    # RIW-U
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
    
    eta[[k]] = X13[[k]] %*% (beta1.initial.hat[,k] - beta0.initial.hat)
  }
  
  
  YY = Y_P = Y_U = Y01;XX = X_P = X_U = X01
  for (k in 1:K) {
    
    id1 = which((abs(error01.hat[[k]]) <= A) & (abs(eta[[k]]) <= M))
    YY = c(YY,Ydata_star[[k]][id1])
    XX = rbind(XX,Xdata_star[[k]][id1,])
    
    id2 = which((abs(error01.hat[[k]]) <= A) & (abs(eta[[k]]) <= M))
    Y_P = c(Y_P,Ydata_star_P[[k]][id2])
    X_P = rbind(X_P,Xdata_star_P[[k]][id2,])
    
    id3 = which((abs(error11.hat[[k]]) <= A) & (abs(eta[[k]]) <= M))
    Y_U = c(Y_U,Ydata_star_U[[k]][id3])
    X_U = rbind(X_U,Xdata_star_U[[k]][id3,])
  }
  
  fit_F = cv.glmnet(XX,YY,family = "gaussian",alpha = 1)
  fit_F_P = cv.glmnet(X_P,Y_P,family = "gaussian",alpha = 1)
  fit_F_U = cv.glmnet(X_U,Y_U,family = "gaussian",alpha = 1)
  
  beta_RIW_TL = coef(fit_F)[-1]
  beta_RIW_TL_P = coef(fit_F_P)[-1]
  beta_RIW_TL_U = coef(fit_F_U)[-1]
  
  list = list(beta_RIW_TL = beta_RIW_TL,
              beta_RIW_TL_P = beta_RIW_TL_P,
              beta_RIW_TL_U = beta_RIW_TL_U,
              SU = length(YY),
              SU_P = length(Y_P),
              SU_U = length(Y_U))
}

# main function in RIW-TL
fit.RIW = function(X,y,n.vec,K,A,M,TU,beta0.initial){
  
  # Target and source population
  n0 = n.vec[1]
  Xdata0 = X[1:n0,]
  Ydata0 = y[1:n0]
  
  Ydata1 = Xdata1 = list()
  for (k in 1:K) {
    
    id.k = ind.set(n.vec,k + 1)
    Xdata1[[k]] = X[id.k,]
    Ydata1[[k]] = y[id.k]
  }
  
  # data splitting 
  X0 = Y0 = list()
  for (j in 1:3) {
    
    X0[[j]] = Xdata0[((j - 1)*floor(n0/3) + 1):(j*floor(n0/3)),]
    Y0[[j]] = Ydata0[((j - 1)*floor(n0/3) + 1):(j*floor(n0/3))]
  }
  
  X.P1 = Y.P1 = X.P2 = Y.P2 = X.P3 = Y.P3 = list()
  for (k in 1:K) {
    
    X.P1[[k]] = Xdata1[[k]][1:floor(n.vec[k+1]/3),]
    X.P2[[k]] = Xdata1[[k]][(floor(n.vec[k+1]/3) + 1):floor(2*n.vec[k+1]/3),]
    X.P3[[k]] = Xdata1[[k]][(floor(2*n.vec[k+1]/3) + 1):n.vec[k+1],]
    
    Y.P1[[k]] = Ydata1[[k]][1:floor(n.vec[k+1]/3)]
    Y.P2[[k]] = Ydata1[[k]][(floor(n.vec[k+1]/3) + 1):floor(2*n.vec[k+1]/3)]
    Y.P3[[k]] = Ydata1[[k]][(floor(2*n.vec[k+1]/3) + 1):n.vec[k+1]]
  }
  
  
  fit1 = fit.RIW.DS(X01 = X0[[1]],Y01 = Y0[[1]],
                    X11 = X.P1,Y11 = Y.P1,
                    X12 = X.P2,Y12 = Y.P2,
                    X13 = X.P3,Y13 = Y.P3,
                    A = A,M = M,TU = TU,K = K,beta0.initial)
  
  fit2 = fit.RIW.DS(X01 = X0[[2]],Y01 = Y0[[2]],
                    X11 = X.P2,Y11 = Y.P2,
                    X12 = X.P3,Y12 = Y.P3,
                    X13 = X.P1,Y13 = Y.P1,
                    A = A,M = M,TU = TU,K = K,beta0.initial)
  
  fit3 = fit.RIW.DS(X01 = X0[[3]],Y01 = Y0[[3]],
                    X11 = X.P3,Y11 = Y.P3,
                    X12 = X.P1,Y12 = Y.P1,
                    X13 = X.P2,Y13 = Y.P2,
                    A = A,M = M,TU = TU,K = K,beta0.initial)
  
  
  beta_RIW_TL = (fit1$beta_RIW_TL + fit2$beta_RIW_TL + fit3$beta_RIW_TL)/3
  beta_RIW_TL_P = (fit1$beta_RIW_TL_P + fit2$beta_RIW_TL_P + fit3$beta_RIW_TL_P)/3
  beta_RIW_TL_U = (fit1$beta_RIW_TL_U + fit2$beta_RIW_TL_U + fit3$beta_RIW_TL_U)/3
  
  SUR = (fit1$SU + fit2$SU + fit3$SU)/length(y)
  SUR_P = (fit1$SU_P + fit2$SU_P + fit3$SU_P)/length(y)
  SUR_U = (fit1$SU_U + fit2$SU_U + fit3$SU_U)/length(y)
  
  list = list(beta_RIW_TL = beta_RIW_TL,
              beta_RIW_TL_P = beta_RIW_TL_P,
              beta_RIW_TL_U = beta_RIW_TL_U,
              SUR = SUR,
              SUR_P = SUR_P,
              SUR_U = SUR_U)
}

#==========adaptive rule in selecting tuning parameters===============
#=====================================================================

# data splitting in RIW-TL & RIW-TL-U & RIW-TL-P using adaptive rule
fit.RIW.DS.adaptive = function(X01,Y01,X11,Y11,X12,Y12,X13,Y13,K,beta0.initial){
  
  # Initial estimators
  #---------------------------
  p = ncol(X01)
  beta0.initial.hat = beta0.initial
  beta1.initial.hat = matrix(0,p,K);sigma = rep(0,K)
  for (k in 1:K) {
    
    fit = cv.ncvreg(X11[[k]],Y11[[k]],family = "gaussian",penalty = "SCAD")
    beta1.initial.hat[,k] = coef(fit)[-1]
    sigma[k] = sum((Y11[[k]] - X11[[k]] %*% matrix(beta1.initial.hat[,k]))^2)/
      (length(Y11[[k]]) - length(which(beta1.initial.hat[,k] != 0)))
  }
  
  ## sources
  #-------------------------------------
  Ydata_star = Xdata_star = eta = list()
  Ydata_star_P = Xdata_star_P = list()
  Ydata_star_U = Xdata_star_U = list()
  error01.hat = error11.hat = list()

  YY = Y_P = Y_U = Y01;XX = X_P = X_U = X01
  
  ##===============RIW-TL-adaptive=============================
  for (k in 1:K) {
    
    eta[[k]] = X13[[k]] %*% (beta1.initial.hat[,k] - beta0.initial.hat)
    
    error01.hat[[k]] = Y13[[k]] - X13[[k]] %*% beta0.initial.hat
    error112.hat = Y12[[k]] - X12[[k]] %*% beta1.initial.hat[,k]
    error11.hat[[k]] = Y13[[k]] - X13[[k]] %*% beta1.initial.hat[,k]
    
    ##------weights estimation-----------
    # optimal bandwidth 
    h_opt = sigma[k]*(3*nrow(X12[[k]])/4)^{-1/5}
  
    kde.fit = kde(x = error112.hat[,1],h = h_opt)  
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
    
    # the choice of A & M
    rho = 0.85;const = 10
    A = quantile(abs(error01.hat[[k]]),rho)
    sim.quan = sum(abs(beta0.initial.hat - beta1.initial.hat[,k]))
    M = const/(1 + sim.quan)
    
    # sample selection
    id1 = which((abs(error01.hat[[k]]) <= A) & (abs(eta[[k]]) <= M))
    
    YY = c(YY,Ydata_star[[k]][id1])
    XX = rbind(XX,Xdata_star[[k]][id1,])
  }
  
  
  ##================RIW-TL-P-adaptive===========================
  for (k in 1:K) {
    
    eta[[k]] = X13[[k]] %*% (beta1.initial.hat[,k] - beta0.initial.hat)
    
    error01.hat[[k]] = Y13[[k]] - X13[[k]] %*% beta0.initial.hat
    error112.hat = Y12[[k]] - X12[[k]] %*% beta1.initial.hat[,k]
    error11.hat[[k]] = Y13[[k]] - X13[[k]] %*% beta1.initial.hat[,k]
    
    ##-------weights estimation-------------------
    fenzi = dnorm(error01.hat[[k]],0,sigma[k])
    fenmu = dnorm(error11.hat[[k]],0,sigma[k])
    weight_P = fenzi/fenmu
    weight_P[which(fenzi < 0 | fenmu <0)] = 0
    
    Ydata_star_P[[k]] = sqrt(weight_P)*Y13[[k]]
    Xdata_star_matrix = matrix(0,nrow(X13[[k]]),p)
    for (i in 1:nrow(X13[[k]])) {
      
      Xdata_star_matrix[i,] = sqrt(weight_P[i])*X13[[k]][i,]
    }
    Xdata_star_P[[k]] = Xdata_star_matrix
    
    # the choice of A & M
    rho = 0.85;const = 10
    A = quantile(abs(error01.hat[[k]]),rho)
    sim.quan = sum(abs(beta0.initial.hat - beta1.initial.hat[,k]))
    M = const/(1 + sim.quan)
    
    # sample selection
    id2 = which((abs(error01.hat[[k]]) <= A) & (abs(eta[[k]]) <= M))
    Y_P = c(Y_P,Ydata_star_P[[k]][id2])
    X_P = rbind(X_P,Xdata_star_P[[k]][id2,])
  }
    
  
  ##================RIW-TL-U-adaptive===========================
  for (k in 1:K) {
    
    # the choice of A, M & TU
    rho = 0.85;const = 10
    A = quantile(abs(error11.hat[[k]]),rho)
    sim.quan = sum(abs(beta0.initial.hat - beta1.initial.hat[,k]))
    M = const/(1 + sim.quan)
    TU = A + M

    # weights estimation
    h_opt = sigma[k]*(3*nrow(X12[[k]])/4)^{-1/5}
    kde.fit = kde(x = error112.hat[,1],h = h_opt)
    
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
    
    id3 = which((abs(error11.hat[[k]]) <= A) & (abs(eta[[k]]) <= M))
    Y_U = c(Y_U,Ydata_star_U[[k]][id3])
    X_U = rbind(X_U,Xdata_star_U[[k]][id3,])
  }
  
  fit_F = cv.glmnet(XX,YY,family = "gaussian",alpha = 1)
  fit_F_P = cv.glmnet(X_P,Y_P,family = "gaussian",alpha = 1)
  fit_F_U = cv.glmnet(X_U,Y_U,family = "gaussian",alpha = 1)
  
  beta_RIW_TL = coef(fit_F)[-1]
  beta_RIW_TL_P = coef(fit_F_P)[-1]
  beta_RIW_TL_U = coef(fit_F_U)[-1]
  
  list = list(beta_RIW_TL = beta_RIW_TL,
              beta_RIW_TL_P = beta_RIW_TL_P,
              beta_RIW_TL_U = beta_RIW_TL_U,
              SU = length(YY),
              SU_P = length(Y_P),
              SU_U = length(Y_U))
}

# RIW-TL using adaptive rule
fit.RIW.adaptive = function(X,y,n.vec,K,beta0.initial){
  
  # Target and source population
  n0 = n.vec[1]
  Xdata0 = X[1:n0,]
  Ydata0 = y[1:n0]
  
  Ydata1 = Xdata1 = list()
  for (k in 1:K) {
    
    id.k = ind.set(n.vec,k + 1)
    Xdata1[[k]] = X[id.k,]
    Ydata1[[k]] = y[id.k]
  }
  
  # data splitting 
  X0 = Y0 = list()
  for (j in 1:3) {
    
    X0[[j]] = Xdata0[((j - 1)*floor(n0/3) + 1):(j*floor(n0/3)),]
    Y0[[j]] = Ydata0[((j - 1)*floor(n0/3) + 1):(j*floor(n0/3))]
  }
  
  X.P1 = Y.P1 = X.P2 = Y.P2 = X.P3 = Y.P3 = list()
  for (k in 1:K) {
    
    X.P1[[k]] = Xdata1[[k]][1:floor(n.vec[k+1]/3),]
    X.P2[[k]] = Xdata1[[k]][(floor(n.vec[k+1]/3) + 1):floor(2*n.vec[k+1]/3),]
    X.P3[[k]] = Xdata1[[k]][(floor(2*n.vec[k+1]/3) + 1):n.vec[k+1],]
    
    Y.P1[[k]] = Ydata1[[k]][1:floor(n.vec[k+1]/3)]
    Y.P2[[k]] = Ydata1[[k]][(floor(n.vec[k+1]/3) + 1):floor(2*n.vec[k+1]/3)]
    Y.P3[[k]] = Ydata1[[k]][(floor(2*n.vec[k+1]/3) + 1):n.vec[k+1]]
  }
  
  
  fit1 = fit.RIW.DS.adaptive(X01 = X0[[1]],Y01 = Y0[[1]],
                    X11 = X.P1,Y11 = Y.P1,
                    X12 = X.P2,Y12 = Y.P2,
                    X13 = X.P3,Y13 = Y.P3,
                    K = K,beta0.initial)
  
  fit2 = fit.RIW.DS.adaptive(X01 = X0[[2]],Y01 = Y0[[2]],
                    X11 = X.P2,Y11 = Y.P2,
                    X12 = X.P3,Y12 = Y.P3,
                    X13 = X.P1,Y13 = Y.P1,
                    K = K,beta0.initial)
  
  fit3 = fit.RIW.DS.adaptive(X01 = X0[[3]],Y01 = Y0[[3]],
                    X11 = X.P3,Y11 = Y.P3,
                    X12 = X.P1,Y12 = Y.P1,
                    X13 = X.P2,Y13 = Y.P2,
                    K = K,beta0.initial)
  
  
  beta_RIW_TL = (fit1$beta_RIW_TL + fit2$beta_RIW_TL + fit3$beta_RIW_TL)/3
  beta_RIW_TL_P = (fit1$beta_RIW_TL_P + fit2$beta_RIW_TL_P + fit3$beta_RIW_TL_P)/3
  beta_RIW_TL_U = (fit1$beta_RIW_TL_U + fit2$beta_RIW_TL_U + fit3$beta_RIW_TL_U)/3
  
  SUR = (fit1$SU + fit2$SU + fit3$SU)/length(y)
  SUR_P = (fit1$SU_P + fit2$SU_P + fit3$SU_P)/length(y)
  SUR_U = (fit1$SU_U + fit2$SU_U + fit3$SU_U)/length(y)
  
  list = list(beta_RIW_TL = beta_RIW_TL,
              beta_RIW_TL_P = beta_RIW_TL_P,
              beta_RIW_TL_U = beta_RIW_TL_U,
              SUR = SUR,
              SUR_P = SUR_P,
              SUR_U = SUR_U)
}


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

