This document provides the code, explanations, and relevant output results for the real data analysis.

# Part I. Data processing
```
# Some packages are needed
library(mvtnorm)
library(ks)
library(Matrix)
library(stringi)
````

```
load("data_all.RData") 
```

```
# Unify the gene_id and retain 15 digits for all genes
data.all$gene_id = stri_sub(data.all$gene_id,1,15)

# the names of tissues
ALL.tissue.name = data.all$tissue
ALL.tissue.name = ALL.tissue.name[!duplicated(ALL.tissue.name)]

# Extract the genes we focus on
MD_137 = read.csv("gene_id.csv",header = F)

# Retain only the variable data for the MD_137 gene
F.tissue = list()
gene.name = list()
for (l in 1:length(ALL.tissue.name)) {
  
  F.tissue[[l]] = subset(data.all,tissue == ALL.tissue.name[l])
  F.tissue[[l]] = subset(F.tissue[[l]],gene_id %in% MD_137$V1)
  
  # remove the missing values
  F.tissue[[l]] = data.frame(t(na.omit(t(F.tissue[[l]]))))
  
  # extract the gene names for each tissues
  gene.name[[l]] = F.tissue[[l]]$gene_id
}

# Obtain the intersect genes
intersect = Reduce(intersect,gene.name)
for (l in 1:length(ALL.tissue.name)) {
  
  F.tissue[[l]] = subset(F.tissue[[l]],gene_id %in% intersect)
}
```


```
# the location of response gene
chr_21 = subset(F.tissue[[1]],X.chr == "chr21")
response.candidate = chr_21$gene_id
r = 2
response.gene = response.candidate[r]
response.gene
```
![d6a48bd461e4caafe125ccd5c7ca1b51.png](../_resources/d6a48bd461e4caafe125ccd5c7ca1b51.png)

```
Y.id = which(intersect == response.gene)
Num.response = c()
for (l in 1:length(ALL.tissue.name)) {
  
  row = subset(F.tissue[[l]],gene_id == response.gene)
  Num.response[l] = length(which(row[1,] != "NA")) - 5
}

# get the remained tissue data after processing
delete.tissue = which(Num.response < 150)
remain.tissue = F.tissue[-delete.tissue]
remain.tissue.name = ALL.tissue.name[-delete.tissue]

# remove the irrelavant variables
for (l in 1:length(remain.tissue.name)) {
  
  remain.tissue[[l]] = remain.tissue[[l]][,-c(1:5)]
  remain.tissue[[l]] = t(remain.tissue[[l]])
  remain.tissue[[l]] = apply(remain.tissue[[l]],2,as.numeric)
  colnames(remain.tissue[[l]]) = intersect
}
```

```
# Extract the target data and souce data, respectively.
p = length(intersect) - 1
target.list = remain.tissue[7:15]
target.tissue.name = remain.tissue.name[7:15]     
source.list = remain.tissue[-c(7:15)]
source.tissue.name = remain.tissue.name[-c(7:15)]
```

# Part II. Data analysis

```
load("target-list.RData")
load("source-list.RData")
```

```
# p: dimension; K: number of sources; TT: number of target tissues
p = ncol(target.list[[1]]) - 1
K = length(source.list)
TT = length(target.list)
p;K;TT

# remove missing values
for (k in 1:K) {source.list[[k]] = na.omit(source.list[[k]])}
for (k in 1:TT) {target.list[[k]] = na.omit(target.list[[k]])}
```
![cd9cfbd5f10e31ec9ba0dcd9da66505e.png](../_resources/cd9cfbd5f10e31ec9ba0dcd9da66505e.png)

```
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
    prop.re1 = Trans.lasso(X, y, n.vec, I.til = id,l1 = T)
    id = (floor(length(train)/2) + 1):length(train)
    prop.re2 = Trans.lasso(X, y, n.vec, I.til = id,l1 = T)
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
  S.result[t,] = colMeans(S);PR.result[t,] = colMeans(PR);NR.result[t,] = colMeans(NR)
}
```

```
RPE.result
```

![51dc908359e0f29fc8d47ec060cf45ab.png](../_resources/51dc908359e0f29fc8d47ec060cf45ab.png)
