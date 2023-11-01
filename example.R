library(ncvreg)
library(mvtnorm)
library(glmnet)
library(ks)
library(bayesdistreg)
#-----------------------------------------------
p = 200;s = 10;K = 10
sig.beta = 1;sig.delta1 = sig.delta2 = 0.5

n0 = 600;n1 = 600;n.vec = c(n0, rep(n1, K))
Sig.X = matrix(0,p,p)
for (i in 1:p) {
  for (j in 1:p) {
    Sig.X[i,j] = 0.5^(abs(i-j))
  }
}

# example
A = 1.5;M = 3;TU = 4.5
size.A0 = 4;d = 8
#------------------------------------------------  
coef.all = Coef.gen(s, d = d, q = 2*s, size.A0 = size.A0,  
                    K = K,sig.beta = sig.beta,
                    sig.delta1 = sig.delta1, 
                    sig.delta2 = sig.delta1, 
                    p = p, exact = F)
B = cbind(coef.all$beta0, coef.all$W)
beta0 = coef.all$beta0    

###generate the data
X <- NULL
y <- NULL
for (k in 1:(K + 1)) {
  
  # posterior drift
  # X <- rbind(X, rmvnorm(n.vec[k], rep(0, p), Sig.X))
  
  # full-distribution-shift
  if(k == 1){X = rmvt(n.vec[k],Sig.X,df = 5)}else{
    X <- rbind(X, rmvnorm(n.vec[k], rep(0, p), Sig.X))
  }
  
  ind.k <- ind.set(n.vec, k)
  epsilon.k = rnorm (n.vec[k], 0, 1)
  y <- c(y, X[ind.k,] %*% B[, k] + epsilon.k)
}    
    
#--------------------------------------------------------------------------
###########RIW-TL###################
fit = fit.RIW(X,y,n.vec,K)
beta_RIW_TL = fit$beta_RIW_TL
beta_RIW_TL_P = fit$beta_RIW_TL_P
beta_RIW_TL_U = fit$beta_RIW_TL_U

SUR = fit$SUR      
SUR_P = fit$SUR_P
SUR_U = fit$SUR_U    

