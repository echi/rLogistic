source("call_logisticL2Eas.R")
source("logisticL2E.R")

alpha = 0.75
lambda = 0.075

n = 200
p = 50

mu = matrix(c(rep(1.1,n/2),rep(-1.1,n/2)),ncol=1)
set.seed(10)
X = matrix(rnorm(n*p),n,p)
X[,1] = X[,1] + mu
X[,17] = X[,17] - 0.75*mu
X = apply(X,2,FUN=function(x){x-mean(x)})

Y = matrix(1, n, 1)
Y[1:(n/2),] = 0

beta = matrix(rnorm(p),p,1)
beta = matrix(0,p,1)
system.time({bHx = cdL2E(X,Y,alpha,lambda,0,beta,tol=1e-9)})
system.time(call_logisticL2Eas(X,Y,alpha,lambda,0,beta, tol=1e-9))
rbHx = read.table("betaHx.csv", sep=",")


plot(bHx[,ncol(bHx)],rbHx[nrow(rbHx),])
abline(0,1)