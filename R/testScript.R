source("call_logisticL2E.R")
source("logisticL2E.R")

alpha = 0.75
lambda = 0.125

n = 100
p = 200

mu = matrix(c(rep(1.1,n/2),rep(-1.1,n/2)),ncol=1)
set.seed(10)
X = matrix(rnorm(n*p),n,p)
X[,1] = X[,1] + mu
X[,17] = X[,17] - 0.75*mu
Y = matrix(1, n, 1)
Y[1:(n/2),] = 0

beta = matrix(rnorm(p),p,1)
bHx = cdL2E(X,Y,alpha,lambda,1,beta,tol=1e-9)
out = call_logisticL2E(X,Y,alpha,lambda,1,beta, tol=1e-9)
rbHx = read.table("betaHx.csv", sep=",")

k = 6
plot(bHx[,k],rbHx[k,])
abline(0,1)