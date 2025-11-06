# This file records the codes used in analyzing the processed GTEx data.

library(glmnet)
library(ks)
library(ncvreg)

#----------------------------
load("target-list.RData")
load("source-list.RData")

#---------------------------------------------------
# p: dimension; K: number of sources; TT: number of target tissues
p = ncol(target.list[[1]]) - 1
K = length(source.list)
TT = length(target.list)
p;K;TT

# remove missing values
for (k in 1:K) {source.list[[k]] = na.omit(source.list[[k]])}
for (k in 1:TT) {target.list[[k]] = na.omit(target.list[[k]])}

#-----------------------------------------------------

RPE.result = S.result = PR.result = NR.result = matrix(0,TT,4)
rownames(RPE.result) = target.tissue.name
colnames(RPE.result) = colnames(S.result) = colnames(PR.result) = 
  colnames(NR.result) = c("RIW-TL","RIW-TL-U","Trans-Lasso","LASSO")

for (t in 1:TT) {
  
  iter = 100
  RPE = S = PR = NR = matrix(0,iter,4)
  
  for (ii in 1:iter) {
    
    # training and testing data
    #------------------------------------------
    n0 = nrow(target.list[[t]])
    train = 1:n0
    train = sample(1:n0,0.7*n0)
    test = setdiff(1:n0,train)
    
    n.vec = c();n.vec = length(train)
    X = target.list[[t]][train,-Y.id]
    y = target.list[[t]][train,Y.id]
    
    for (k in 1:K) {
      
      n.vec = c(n.vec,nrow(source.list[[k]]))
      X = rbind(X,source.list[[k]][,-Y.id])
      y = c(y,source.list[[k]][,Y.id])
    }
    
    #--------------------------------------------------------------------------
    ###########RIW-TL###################
    fit = hyper.select.cv(X,y,n.vec,K,fold = 5)
    
    A = fit$A.opt;
    M = fit$M.opt;
    TU = A + M;
    fit1 = fit.RIW(X,y,n.vec,K,A,M,TU)
    beta_RIW_TL = fit1$beta_RIW_TL
    
    ###########RIW-TL-U###################
    A = fit$A.opt.u;
    M = fit$M.opt.u;
    TU = A + M;
    fit1 = fit.RIW(X,y,n.vec,K,A,M,TU)
    beta_RIW_TL_U = fit1$beta_RIW_TL_U
    
    PE.RIW.TL = norm(X[test,] %*% beta_RIW_TL - y[test],"2")^2/length(test)
    PE.RIW.TL.U = norm(X[test,] %*% beta_RIW_TL_U - y[test],"2")^2/length(test)
    
    ##================Trans-lasso===========================
    id = 1:floor(length(train)/2)
    prop.re1 = Trans.lasso(X, y, n.vec, I.til = id,l1=T)
    id = (floor(length(train)/2) + 1):length(train)
    prop.re2 = Trans.lasso(X, y, n.vec, I.til = id,l1=T)
    beta.trans.lasso = (prop.re1$beta.hat + prop.re2$beta.hat)/2
    
    PE.trans.lasso = norm(X[test,] %*% beta.trans.lasso - y[test],"2")^2/length(test)
    
    ##======================LASSO=========================
    fit = cv.glmnet(X[1:n.vec[1],],y[1:n.vec[1]],family = "gaussian",alpha = 1)
    beta.lasso = coef(fit)[-1]
    PE.lasso = norm(X[test,] %*% beta.lasso - y[test],"2")^2/length(test)
    
    # output the relative prediction error (RPE)
    RPE[ii,] = c(c(PE.RIW.TL,PE.RIW.TL.U,PE.trans.lasso)/PE.lasso,1)
    
    # output-signal
    S[ii,] = c(length(which(beta_RIW_TL != 0))/p,
               length(which(beta_RIW_TL_U != 0))/p,
               length(which(beta.trans.lasso != 0))/p,
               length(which(beta.lasso != 0))/p)
    
    PR[ii,] = c(length(intersect(which(beta_RIW_TL != 0),which(beta.lasso != 0)))/length(which(beta.lasso != 0)),
                length(intersect(which(beta_RIW_TL_U != 0),which(beta.lasso != 0)))/length(which(beta.lasso != 0)),
                length(intersect(which(beta.trans.lasso != 0),which(beta.lasso != 0)))/length(which(beta.lasso != 0)),
                1)
    
    NR[ii,] = c(length(intersect(which(beta_RIW_TL != 0),which(beta.lasso == 0)))/length(which(beta.lasso == 0)),
                length(intersect(which(beta_RIW_TL_U != 0),which(beta.lasso == 0)))/length(which(beta.lasso == 0)),
                length(intersect(which(beta.trans.lasso != 0),which(beta.lasso == 0)))/length(which(beta.lasso == 0)),
                0)
  }
  
  RPE.result[t,] = colMeans(RPE)
  
  S.result[t,] = colMeans(S);
  PR.result[t,] = colMeans(PR);
  NR.result[t,] = colMeans(NR)
}

RPE.result
S.result;PR.result;NR.result


