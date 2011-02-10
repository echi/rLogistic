source("call_ENLS.R")
library(glmnet)


set.seed(1)
n = 5000
p.true = 2
p.dummy = 100
p = p.true + p.dummy
X = matrix(rnorm(100),n,p.true)
beta = matrix(c(2,3),2,1)
Y = X %*% beta + matrix(rnorm(n),n,1)
X = cbind(X, matrix(rnorm(n*(p.dummy)),n,p.dummy))
alpha = 1

Xc = apply(X,2,FUN=function(x){x-mean(x)})
Yc = Y - mean(Y)

sorted = sort(abs(t(Xc) %*% Yc), index.return=T)
scores = sorted$x
ix = sorted$ix
lambdaMax = max(scores)
CS = sqrt(sum(Yc^2))*apply(Xc[,ix],2,FUN=function(x){sqrt(sum(x^2))})
lambdas = (scores + CS)/(1+CS/lambdaMax)/n

lambda = 1
beta.ENLS = ENLS(Xc, Yc, alpha, lambda, tol=1e-9)
glmnet.fit = glmnet(Xc, Yc, family="gaussian",alpha=alpha, lambda=lambda,standardize=FALSE,thresh=1e-9)
cbind(as.matrix(glmnet.fit$beta), beta.ENLS)

alpha = 0
beta.ENLS = ENLS(Xc, Yc, alpha, lambda, tol=1e-9)
beta.ridge = solve(t(Xc) %*% Xc + nrow(Xc)*lambda*(1-alpha)*diag(1,ncol(Xc),ncol(Xc)))%*%t(Xc) %*% Yc
cbind(beta.ENLS, beta.ridge)