source("call_ENLS.R")

maxrerr = function(b, bold) {
	mrerr = 0
	ix = which(abs(bold) > 0)
	for (j in ix) {
		mrerr = max(mrerr, abs(b[j] - bold[j])/abs(bold[j]))
	}
	return(mrerr)
}

# Make some test data
set.seed(1)
n = 500
p.true = 2
p.dummy = 98
p = p.true + p.dummy
X = matrix(rnorm(2*n*p.true),2*n,p.true)
mu = 2
X[1:n,] = X[1:n,] + mu
X[(n+1):(2*n),] = X[(n+1):(2*n),] - mu
beta = matrix(c(1,2),2,1)
X = cbind(X, matrix(rnorm(2*n*(p.dummy)),2*n,p.dummy))
X = apply(X,2,FUN=function(x){x-mean(x)})
Y = ifelse(runif(2*n) <= plogis(X[,1:p.true,drop=F] %*% beta), 1, 0)
Y = 2*Y-1

## MM algorithm
loss = function(u, t=0.5) {	
	loss = ((1-t)^2 + 2*(1-t)*(t-u))*ifelse(u <= t, 1, 0)
	loss = loss + (1-u)^2 * ifelse(t < u & u <= 1,1,0)
	return(loss)
}

lp = function(u, t=0.5) {
	lp = -2*((1-t)*ifelse(u <= t, 1, 0) + (1-u)*ifelse(u > t & u <= 1, 1, 0))
	return(lp)
}

updateResponse = function(y, X, b0, b) {
	u = X %*% b
	z = y*lp(y * (b0 + u))
	z = z - mean(z)
	z = u - 0.5*z
	return(z)
}

lambda2 = .001
svdX = svd(X)
P = svdX$v %*% diag(svdX$d/(svdX$d^2 + 0.5*lambda2)) %*% t(svdX$u)

b0 = 1
b = matrix(0, ncol(X))
L = c()
mrerr = c()
tol = 1e-7

for (i in 1:10000) {
	bold = c(b0, b)
	phi = Y*lp(Y *(b0 + X %*% b))
	b0 = b0 - 0.5*mean(phi)
	b = b - 0.5*P %*% phi
	L[i] = sum(loss(Y*(b0 + X%*%b)))
	mrerr[i] = maxrerr(c(b0, b), bold)
	if (mrerr[i] < tol)
		break
}



## Fit intercept
fitIntercept = function(b0, Y) {
	return(sum(loss(Y*b0)))	
}

b0 = nlminb(1, fitIntercept, Y=Y)$par
b = matrix(0, ncol(X))
L = c()
mrerr = c()

## with coordinate descent
alpha = 1
lambdaMax = max(abs(t(X) %*% (Y*lp(Y*b0))))/(2*nrow(X)*alpha)
eps = 0.0001
lambdaMin = eps*lambdaMax
nLambda = 20
lambda = 10^(seq(log10(lambdaMax), log10(lambdaMin), length.out=nLambda))
tol = 1e-9
maxiter = 10000

b0Hx = matrix(NA, 1, nLambda)
bHx = matrix(NA, p, nLambda)

for (iLambda in 1:nLambda) {
	for (i in 1:maxiter) {
		bold = c(b0, b)
		phi = Y*lp(Y *(b0 + X %*% b))
		b0 = b0 - 0.5*mean(phi)
		z = X %*% b - 0.5*(phi-mean(phi))
		b = ENLS(X, z, alpha, lambda[iLambda], b, tol=1e-3)
#		L[i] = sum(loss(Y*(b0 + X%*%b))) + lambda*(alpha*sum(abs(b)) + 0.5*(1-alpha)*sum(b^2))
		mrerr[i] = maxrerr(c(b0, b), bold)
		if (mrerr[i] < tol)
			break
	}
	b0Hx[iLambda] = b0
	bHx[,iLambda]	= b
}

f = function(beta, y, X, lambda) {
	b0 = beta[1]
	b = matrix(beta[-1],ncol=1)
	return(sum(loss(y*(b0 + X %*% b))) + 0.5*lambda*sum(beta^2))
}



par = c(rep(0, ncol(X)+1))
hhsvm.op = nlminb(par,f, y=Y, X=X, lambda=lambda2)

loss = function(u, t=0.5) {	
	loss = ((1-t)^2 + 2*(1-t)*(t-u))*ifelse(u <= t, 1, 0)
	loss = loss + (1-u)^2 * ifelse(t < u & u <= 1,1,0)
	return(loss)
}


alpha = 1

u = seq(-2, 2, 0.1)
u2 = -1.11
maj = loss(u2, t) + lp(u2, t)*(u - u2) + (u - u2)^2
