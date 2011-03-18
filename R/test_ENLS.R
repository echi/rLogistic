source("call_ENLS.R")
library(glmnet)
library(lars)

set.seed(1)
n = 500
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

# Test to see if the R glmnet Elastic net Least squares is _not_ working by using the LASSO formulation of the Elastic-net problem

Xtilde = rbind(Xc, sqrt(lambda*(1-alpha)*n)*diag(1,p,p))
Ytilde = rbind(Yc, matrix(0,p,1))

alpha = 0.25
beta.ENLS = ENLS(Xc, Yc, alpha, lambda, tol=1e-9)
glmnet.fit = glmnet(Xc, Yc, family="gaussian",alpha=alpha, lambda=lambda,standardize=FALSE,thresh=1e-9)
glmnet.fit2 = glmnet(Xtilde, Ytilde, family="gaussian",alpha=1, lambda=lambda*alpha,standardize=FALSE,thresh=1e-9)
#lars.fit = lars(Xtilde,Ytilde,type="lasso",normalize=FALSE)
beta.ENLS2 = ENLS(Xtilde, Ytilde, alpha=1, lambda=lambda*alpha, tol=1e-9)
cbind(as.matrix(glmnet.fit$beta), as.matrix(glmnet.fit2$beta), beta.ENLS, beta.ENLS2)


lambda*(alpha*norm(glmnet.fit2$beta,"1") + 0.5*(1-alpha)*norm(glmnet.fit2$beta,"f")^2)
lambda*(alpha*norm(as.matrix(beta.ENLS2),"1") + 0.5*(1-alpha)*norm(as.matrix(beta.ENLS2),"f")^2)

lambda*(alpha*norm(glmnet.fit$beta,"1") + 0.5*(1-alpha)*norm(glmnet.fit$beta,"f")^2)
lambda*(alpha*norm(as.matrix(beta.ENLS),"1") + 0.5*(1-alpha)*norm(as.matrix(beta.ENLS),"f")^2)
