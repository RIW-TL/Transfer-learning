# It provides an example of running the above functions,
# and specifically how to use our method.
#--------------------------

library(ncvreg)
library(mvtnorm)
library(glmnet)
library(ks)
library(bayesdistreg)
#-----------------------------------------------

# parameter setting
p = 200;s = 10;K = 5;q = 2*s
sig.beta = 1;sig.delta1 = sig.delta2 = 0.5

n0 = 150;n1 = 1200;n.vec = c(n0, rep(n1, K))
Sig.X = matrix(0,p,p)
for (i in 1:p) {
  for (j in 1:p) {
    Sig.X[i,j] = 0.5^(abs(i-j))
  }
}

size.A0 = 4;d = 4

# example
#------------------------------------------------  
coef.all = Coef.gen(s, d = d, q = q, size.A0 = size.A0,  
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

###########RIW-TL based on adaptive rule####
# the parameters rho and c can also be selected by cv method
fit1 = fit.RIW.adaptive(X,y,n.vec,K,rho = 0.3, c = 10)
beta_RIW_TL_ada = fit1$beta_RIW_TL

###########Trans-Lasso#############
prop.re1 <- Trans.lasso(X, y, n.vec, I.til = 1:(n.vec[1]/2), l1 = T)
prop.re2 <- Trans.lasso(X, y, n.vec, I.til = (n.vec[1]/2 + 1):n.vec[1], l1=T)
beta.prop <- (prop.re1$beta.hat + prop.re2$beta.hat)/2

###########LASSO####################
X0 = X[1:n0,];y0 = y[1:n0]
fit = cv.glmnet(X0,y0,family = "gaussian",alpha = 1)
beta.lasso = coef(fit)[-1]


#-----------------------------------------------
sum((beta_RIW_TL - beta0)^2)
sum((beta_RIW_TL_U - beta0)^2)
sum((beta_RIW_TL_ada - beta0)^2)
sum((beta.prop - beta0)^2)
sum((beta.lasso - beta0)^2)





