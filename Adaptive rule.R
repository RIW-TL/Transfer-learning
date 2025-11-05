# The RIW-TL using adaptive selection rule for the tuning parameters.
#-------------------------------
library(ncvreg)
library(mvtnorm)
library(glmnet)
library(ks)
library(bayesdistreg)
#------------------------------------------------------------

# rho: quantile used in determining A
# c  : tuning parameter used in determining M_k
# the values of rho and c can be refined case-by-case

fit.RIW.DS.adaptive = function(X01,Y01,X02,Y02,
                               X11,Y11,X12,Y12,X13,Y13,
                               rho = 0.3, c = 10){
  
  p = ncol(X01)
  K = length(X11)
  
  # Initial estimators (SCAD penalty as example)
  #---------------------------
  fit = cv.ncvreg(X01,Y01,family = "gaussian",penalty = "SCAD")
  beta0.initial.hat = coef(fit)[-1]
  
  
  # beta1.initial.hat = matrix(0,p,K)
  # for (k in 1:K) {
  #   
  #   fit = cv.ncvreg(X11[[k]],Y11[[k]],family = "gaussian",penalty = "SCAD")
  #   beta1.initial.hat[,k] = coef(fit)[-1]
  # }
  beta1.initial.hat = B[,-1]
  

  #-------------------------------------
  Ydata_star = Xdata_star = eta = list()
  Ydata_star_U = Xdata_star_U = list()
  error01.hat = error11.hat = list()
  
  YY = Y_P = Y_U = Y02;XX = X_P = X_U = X02
  
  for (k in 1:K) {
    
    # KDE and weighted sources data
    error112.hat = Y12[[k]] - X12[[k]] %*% beta1.initial.hat[,k]
    kde.fit = kde(x = error112.hat[,1])
    
    error01.hat[[k]] = Y13[[k]] - X13[[k]] %*% beta0.initial.hat
    error11.hat[[k]] = Y13[[k]] - X13[[k]] %*% beta1.initial.hat[,k]
    
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
    
    # contrast vectors adjusted by predictors
    eta[[k]] = X13[[k]] %*% (beta1.initial.hat[,k] - beta0.initial.hat)
    
    # sample selection
    A = quantile(abs(error11.hat[[k]]),rho)
    M_k = c/(1 + norm(beta1.initial.hat[,k] - beta0.initial.hat,"2"))
      
    id1 = which((abs(error11.hat[[k]]) <= A) & (abs(eta[[k]]) <= M_k))
    YY = c(YY,Ydata_star[[k]][id1])
    XX = rbind(XX,Xdata_star[[k]][id1,])
  }
  
  
  # weighted lasso optimization
  const = 0.5
  fit_F = glmnet(XX,YY,family = "gaussian",alpha = 1,
                 lambda = const*sqrt(log(p)/length(YY)))
  
  beta_RIW_TL = coef(fit_F)[-1]
  
  list = list(beta_RIW_TL = beta_RIW_TL,SU = length(YY))
}

fit.RIW.adaptive = function(X,y,n.vec,K,rho = 0.3, c = 10){
  
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
  
  
  fit1 = fit.RIW.DS.adaptive(X01 = X01[[1]],Y01 = Y01[[1]],
                             X02 = X02[[1]],Y02 = Y02[[1]],
                             X11 = X.P1,Y11 = Y.P1,
                             X12 = X.P2,Y12 = Y.P2,
                             X13 = X.P3,Y13 = Y.P3,
                             B, rho, c)

  
  fit2 = fit.RIW.DS.adaptive(X01 = X01[[2]],Y01 = Y01[[2]],
                             X02 = X02[[2]],Y02 = Y02[[2]],
                             X11 = X.P2,Y11 = Y.P2,
                             X12 = X.P3,Y12 = Y.P3,
                             X13 = X.P1,Y13 = Y.P1,
                             B, rho, c)

  
  fit3 = fit.RIW.DS.adaptive(X01 = X01[[3]],Y01 = Y01[[3]],
                             X02 = X02[[3]],Y02 = Y02[[3]],
                             X11 = X.P3,Y11 = Y.P3,
                             X12 = X.P1,Y12 = Y.P1,
                             X13 = X.P2,Y13 = Y.P2,
                             B, rho, c)

  
  
  beta_RIW_TL = (fit1$beta_RIW_TL + fit2$beta_RIW_TL + fit3$beta_RIW_TL)/3
  SUR = (fit1$SU + fit2$SU + fit3$SU)/length(y)
  
  list = list(beta_RIW_TL = beta_RIW_TL,SUR = SUR)
}



