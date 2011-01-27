# Description: Suite of coordinate descent, variable selection functions.
#
# Author: Eric Chi
# Date: 08 Jun 2010
#
# Change log
# -----------
# 2010.06.09:	added nlambda parameter for MLEboth.xic to pass to glmnet
# 2010.06.09:	removed return.coef option for MLEboth.xic
# 2010.06.09:	added MLEboth.xic
# 2010.06.09:	added L2Eboth.xic
# 2010.06.09: 	put if..stop syntax in call_logistic_l2e
# 2010.06.09:	added lines to more automatically set paths depending on machine in use (i.e. MBP or SUG@R)
# 2010.06.08: 	added MLExic
# 2010.06.08: 	changed L2Ebic to L2Exic to handle 'aic' as well as 'bic'

#machine = "personal"
machine = "sugar"

pwd = getwd()
home = system("echo $HOME", intern=T)
if (machine == "personal") home = paste(home,"Dropbox",sep="/")
setwd(paste(home,"Research/Rfiles/L2eIterativeSolver/Rfunctions",sep="/"))
if (!is.loaded("softThreshold")) {
	dyn.load("softThreshold.so")	
}
if (!is.loaded("logistic_l2e")) {
	dyn.load("logistic_L2e.so")	
}
setwd(pwd)

soft.th = function(x,lambda) {
	n = length(x)
	y = matrix(0,nrow=n)
	posIX = which(x > lambda)
	y[posIX] = x[posIX] - lambda
	negIX = which(x < -lambda)
	y[negIX] = x[negIX] + lambda
	return(y)
}

call_soft.th = function(x,lambda) {
	n = length(x)
	shrunk = double(n)
	.C("softThreshold", res=as.double(shrunk), as.double(x), as.double(lambda), as.integer(n))$res
}

#void meanL2e(double *loss, double *X, int *Y, double *beta0, double *beta,
#             int *n, int *p, double *w, double *lambda, double *alpha)

call_meanL2e = function(par, X, Y, lambda=0, alpha=1, w = 1) {
		n = dim(X)[1]
		p = dim(X)[2]
		beta0 = par[1]
		beta = par[-1]
		if (length(beta) == p & length(Y) == n) {
			loss = double(1)
			.C("meanL2e", res=as.double(loss), as.double(X), as.integer(Y),
				as.double(beta0), as.double(beta), as.integer(n), as.integer(p),
				as.double(w), as.double(lambda),as.double(alpha))$res
		}
}

# void meanL2eBatch(double *loss, double *X, int *Y, double *beta0, double *beta,
#                   int *n, int *p, double *w, double *lambda, double *alpha, int *m)

# beta is a matrix, each column is a p length vector.
call_meanL2eBatch = function(X, Y, lambda, alpha, beta0, beta, w = 1) {
	if (is.matrix(X) & is.matrix(beta)) {
		n = dim(X)[1]
		p = dim(X)[2]
		m = ncol(beta)
		loss = double(m)
		.C("meanL2eBatch", res=as.double(loss), as.double(X), as.integer(Y),
			as.double(beta0), as.double(beta), as.integer(n), as.integer(p),
			as.double(w), as.double(lambda), as.double(alpha), as.integer(m))$res
	}	
}

mean.l2e = function(par, X, Y, lambda, alpha) {
  w = par[1]
  theta = matrix(par[-1],ncol=1)
  p.case = plogis(X%*%theta)
  p.cntl = 1 - p.case
  l2e = 0.5*mean( w^2 * (p.case^2 + p.cntl^2) -2 * w * (p.case^Y)*(p.cntl^(1-Y)))
  ridg = 0.5*(1-alpha)*sum(theta[-1]^2)
  lass = alpha*sum(abs(theta[-1]))
  return(l2e + lambda*(ridg+lass))
}

median.l2e = function(par,X,Y,lambda,alpha) {
  w = par[1]
  theta = matrix(par[-1],ncol=1)
  p.case = plogis(X%*%theta)
  p.cntl = 1 - p.case
  l2e = 0.5*median( w^2 * (p.case^2 + p.cntl^2) -2 * w * (p.case^Y)*(p.cntl^(1-Y)))
  ridg = 0.5*(1-alpha)*sum(theta[-1]^2)
  lass = alpha*sum(abs(theta[-1]))
  return(l2e + lambda*(ridg+lass))
}

# C implementation
# Estimated parameter beta is saved in 'mybeta.csv'

#void logistic_l2e(double *loss, double *X, int *Y, double *beta0, double *beta, int *n, int *p,
#                  double * K, double *w, double *lambda, double *alpha, int *maxiter, int *check) {

# save_flag: 0 if don't save intermediate values of beta
#            1 if save intermediate values of beta
# 

call_logistic_l2e = function(X, Y, alpha, lambda, beta0, beta, w = 1, K = 0.16, niter=500, check=F, tol=1e-9, save_flag=1, filename="bHx.csv") {
	if (!is.matrix(X)) stop("X must be a matrix.")

	n = dim(X)[1]
	p = dim(X)[2]
		
	if (length(Y) != n | length(beta) !=p) stop("length(Y) == nrow(X) & length(beta) == ncol(X)")
		
	print(paste("alpha = ",alpha, ", lambda = ",lambda,", w = ", w, ",K = ",K))
	loss = double(niter)
	.C("logistic_l2e", as.double(loss), as.double(X), as.integer(Y), as.double(beta0), as.double(beta),
		as.integer(n), as.integer(p), as.double(K), as.double(w), as.double(lambda), as.double(alpha),
		as.integer(niter),as.integer(check),as.double(tol),as.integer(save_flag), as.character(filename))
}


# Iterative L2e
logisticL2eIRLS = function(X,Y,b0,beta,K=0.16,w=1,miter=10000,tol=1e-8) {
	n = length(Y)
	p = ncol(X)
	if (nrow(X) != n | length(beta) != p | K < 0.16) stop("Error")

	X = cbind(1,X)
	M = t(X)%*%X
# An alternative is to compute an inverse once and use it several times.
# Perhaps we could put an option here when the condition number is nice.
#	M = solve(t(X) %*% X) %*% t(X)
	U = 2*Y-1
	nbeta = matrix(c(b0,beta),ncol=1)
	betaHx = nbeta
	last_loss = call_meanL2e(nbeta,X[,-1,drop=F], Y, lambda=0, alpha=0, w = w)
	for (i in 1:miter) {
		P = plogis(X%*%nbeta)
		Q = 2*P - 1
		WZ = matrix(P*(1-P)*(w*Q - U),ncol=1)
		nbeta = nbeta - (1/K)*solve(M, t(X) %*% WZ)
#		nbeta = nbeta - (1/K) * M %*% WZ
		current_loss = call_meanL2e(nbeta,X[,-1,drop=F], Y, lambda=0, alpha=0, w = w)
		if ( abs(current_loss-last_loss)/(abs(last_loss)+1) > tol) {
			betaHx = cbind(betaHx, nbeta)	
			last_loss = current_loss
		} else {
			break
		}
	}
	print(paste("IRLS number of iterations: ",i-1))
	return(betaHx)
}

read.beta = function(filename="mybeta.csv") {
	mybeta = read.table(file=filename,sep=",")
	mybeta = t(mybeta[,-dim(mybeta)[2]])
	mybeta = as.matrix(mybeta)
	mybeta = mybeta[nrow(mybeta):1,,drop=F]	
	return(mybeta)
}

# Calculate the largest value of lambda for which all regression coefficients are zero.
lambdaMax = function(XX, Z, w, alpha) {
	if (alpha <= 0 | alpha > 1) stop("alpha must be in (0,1]")
	U = 2*Z-1
	zbar = mean(Z)
	#beta0 = log(zbar/(1-zbar))
	beta0 = -log(w/(zbar-0.5*(1-w)) - 1)
	beta = matrix(0,p,1)
	P = plogis(beta0)
	Q = 2*P - 1
	G = P*(1-P)
	lmax = (max(abs((w*G/length(Z))*t(XX)%*%(w*Q-U))))/alpha
}

# This function returns the variables selected by both aic and bic. It's almost identicaly
# to L2Exic, but is twice as fast because coordinate descent is called only once for both
# criterion calculations.
L2Eboth.xic = function(lseq, XX, Z, beta0, beta, w, alpha, K, niter, tol, sigset,
	plotFlag=F, return.coef = F, warm.start=T) {
		
	p = ncol(XX)
	T = length(lseq)
	varSel = matrix(0,p+1,T)
	bHx0 = matrix(0,p+1,T)
	df = matrix(NA,1,T)
	for (lix in 1:T) {
		lambda = lseq[lix]
		system.time(call_logistic_l2e(XX, Z, alpha, lambda, beta0, beta, w=w,check=T, K = K, niter=niter, tol=tol))
		mybeta = read.beta()
		varIX = c(1,which(mybeta[-1,ncol(mybeta)] != 0)+1)
		df[lix] = length(varIX) - 1
		varSel[varIX,lix] = 1
		bHx0[,lix] = mybeta[,ncol(mybeta)]
		# Warm start
		if (warm.start) {
			beta0 = bHx0[1,lix]
			beta = bHx0[-1,lix,drop=F]
		}
	}

	P = plogis(cbind(1,XX) %*% bHx0)
	dev = -2 * (t(Z) %*% log(P) + t(1-Z)  %*% log(1-P))
	AIC = dev + 2*df
	BIC = dev + log(length(Z))*(df)

	ul2e.OsdIX = which(AIC==min(AIC))
	aicSel = setdiff(which(abs(bHx0[,ul2e.OsdIX ]) > 0),1)
	ul2e.OsdIX = which(BIC==min(BIC))
	bicSel = setdiff(which(abs(bHx0[,ul2e.OsdIX ]) > 0),1)

	if (plotFlag) {
		quartz(); par(mfrow=c(2,2))
		regpathPlot(log10(lseq), bHx0, subset = sigset) 
		regpathPlot(log10(lseq), bHx0, subset = sigset) 
		xlab = expression(log[10](lambda))
		ylab = expression(paste(L[2],"E", "BIC"))
		plot(log10(lseq), BIC, pch=16, xlab=xlab, ylab=ylab)
		axis(side=3,at=log10(lseq),labels=df, tick=F)
		ylab = expression(paste(L[2],"E", "AIC"))
		plot(log10(lseq), BIC, pch=16, xlab=xlab, ylab=ylab)
		axis(side=3,at=log10(lseq),labels=df, tick=F)
	}

	if (!return.coef)
		return(list(aicSel=aicSel, bicSel=bicSel))

	bTHx = logisticL2eIRLS(XX[,setdiff(varSel,1)-1,drop=F],Z, bHx0[1,ul2e.OsdIX], bHx0[varSel,ul2e.OsdIX],w=w)
		return(list(aicSel=aicSel, bicSel=bicSel, bTHx=bTHx))
}

L2Exic = function(lseq, XX, Z, beta0, beta, w, alpha, K, niter, tol, sigset,
	criterion = 'bic', plotFlag = F, return.coef = F) {
	
	if (is.na(match(criterion,c('bic','aic')))) stop("criterion must be either 'aic' or 'bic'")
	
	p = ncol(XX)
	T = length(lseq)
	varSel = matrix(0,p+1,T)
	bHx0 = matrix(0,p+1,T)
	df = matrix(NA,1,T)
	for (lix in 1:T) {
		lambda = lseq[lix]
		system.time(call_logistic_l2e(XX, Z, alpha, lambda, beta0, beta, w=w,check=T, K = K, niter=niter, tol=tol))
		mybeta = read.beta()
		varIX = c(1,which(mybeta[-1,ncol(mybeta)] != 0)+1)
		df[lix] = length(varIX) - 1
		varSel[varIX,lix] = 1
		bHx0[,lix] = mybeta[,ncol(mybeta)]
		# Warm start
		beta0 = bHx0[1,lix]
		beta = bHx0[-1,lix,drop=F]
	}

	P = plogis(cbind(1,XX) %*% bHx0)
	dev = -2 * (t(Z) %*% log(P) + t(1-Z)  %*% log(1-P))
	if (criterion == 'bic') {
		xIC = dev + log(length(Z))*(df)
	} else {
		xIC = dev + 2*df
	}

	ul2e.OsdIX = which(xIC==min(xIC))
	varSel = setdiff(which(abs(bHx0[,ul2e.OsdIX ]) > 0),1)

	if (plotFlag) {
		quartz(); par(mfrow=c(2,1))
		regpathPlot(log10(lseq), bHx0, subset = sigset) 
		if (criterion == "bic") {
			ylab = expression(paste(L[2],"E", "BIC"))
		} else {
			ylab = expression(paste(L[2],"E", "AIC"))		}	
		xlab = expression(log[10](lambda))
		plot(log10(lseq), xIC, pch=16, xlab=xlab, ylab=ylab)
		axis(side=3,at=log10(lseq),labels=df, tick=F)
	}

	if (!return.coef)
		return(list(varSel=varSel))	

	bTHx = logisticL2eIRLS(XX[,setdiff(varSel,1)-1,drop=F],Z, bHx0[1,ul2e.OsdIX], bHx0[varSel,ul2e.OsdIX],w=w)
	return(list(varSel=varSel,bTHx=bTHx))	

}

L2Ecv = function(lseq, XX, Z, beta0, beta, w, alpha, K, niter, tol, foldMatrix, plotFlag = F) {
	
	p = ncol(XX)
	T = length(lseq)
	varSel = matrix(0,p+1,T)
	bHx0 = matrix(0,p+1,T)
	df = matrix(NA,1,T)
	for (lix in 1:T) {
		lambda = lseq[lix]
		system.time(call_logistic_l2e(XX, Z, alpha, lambda, beta0, beta, w=w,check=T, K = K, niter=niter, tol=tol))
		mybeta = read.beta()
		varIX = c(1,which(mybeta[-1,ncol(mybeta)] != 0)+1)
		df[lix] = length(varIX) - 1
		varSel[varIX,lix] = 1
		bHx0[,lix] = mybeta[,ncol(mybeta)]
	}
	

	if (plotFlag) quartz(); regpathPlot(log10(lseq), bHx0, subset = c(2)) 

	# Cross validation: 10-fold
	uL2eCosts = matrix(NA,nrow(foldMatrix),length(lseq))
	dof = matrix(NA,nrow(foldMatrix),length(lseq))
	for (i in 1:nrow(foldMatrix)) {
		print(paste("@@@@@@ CV Fold ",i,"@@@@@@"))
		varSel = matrix(0,p+1,T)
		bHx = matrix(0,p+1,T)
		cvix = foldMatrix[i,]
		for (lix in 1:T) {
			lambda = lseq[lix]
			system.time(call_logistic_l2e(XX[-cvix,,drop=F], Z[-cvix], alpha, lambda, beta0, beta, w=w, check=T, K=K, niter=niter, tol=tol))
			mybeta = read.beta()
			varIX = c(1,which(mybeta[-1,ncol(mybeta)] != 0)+1)
			dof[i,lix] = length(varIX) - 1
			bHx[,lix] = mybeta[,ncol(mybeta)]
			uL2eCosts[i,lix] = mean.l2e(par=c(w,bHx[,lix]),cbind(1,XX[cvix,,drop=F]),Z[cvix],0,alpha=alpha)
		}
	}

	ul2e.mean = apply(uL2eCosts,2,FUN=function(x) {return(mean(x,na.rm=T))})
	ul2e.sd = apply(uL2eCosts,2,FUN=function(x) {return(sd(x,na.rm=T))})
	ul2e.minIX = which(ul2e.mean == min(ul2e.mean))
	ul2e.OsdIX = min(which(ul2e.mean < ul2e.mean[ul2e.minIX] + ul2e.sd[ul2e.minIX]))
	varSel = setdiff(which(abs(bHx0[,ul2e.OsdIX ]) > 0),1)
	bTHx = logisticL2eIRLS(XX[,setdiff(varSel,1)-1,drop=F],Z, bHx0[1,ul2e.OsdIX], bHx0[varSel,ul2e.OsdIX],w=w)

	if (plotFlag) {
		ylim = c(min(uL2eCosts, na.rm=T),max(uL2eCosts, na.rm=T))
		ylab = expression(paste(L[2],"E Loss"))
		quartz(); cvPlot(lambda=lseq, df=df, mx = ul2e.mean, sdev=ul2e.sd, ylim=ylim,ylab=ylab)
	}

	return(list(varSel=varSel,bTHx=bTHx))	
	
}

# This function returns the variables selected by both aic and bic. It's almost identicaly
# to MLExic, but is twice as fast because coordinate descent is called only once for both
# criterion calculations.
library(glmnet)
MLEboth.xic = function(XX, Z, alpha, sigset, plotFlag = F, nlambda = 20) {
	
	en.model = glmnet(x=XX, y=factor(Z), alpha = alpha, family="binomial",
		nlambda=nlambda, lambda.min=0.0001)
	en.bHx = as.matrix(en.model$beta)
	en.lseq = en.model$lambda
	
	p = ncol(XX)
	T = length(en.lseq)
	varSel = matrix(0,p+1,T)
#	beta = as.matrix(en.model$beta)
	P = plogis(en.model$a0 + XX%*%en.bHx)

	dev = -2 * (t(Z) %*% log(P) + t(1-Z)  %*% log(1-P))
	BIC = dev + log(length(Z))*(en.model$df)
	AIC = dev + 2*en.model$df

	ul2e.OsdIX = which(AIC==min(AIC))
	aicSel = which(abs(en.bHx[,ul2e.OsdIX ]) > 0)
	ul2e.OsdIX = which(BIC==min(BIC))
	bicSel = which(abs(en.bHx[,ul2e.OsdIX ]) > 0)

	if (plotFlag) {
		quartz(); par(mfrow=c(2,2))
		regpathPlot(log10(en.lseq), rbind(en.model$a0,en.bHx), subset = sigset)
		regpathPlot(log10(en.lseq), rbind(en.model$a0,en.bHx), subset = sigset) 
		xlab = expression(log[10](lambda))
		ylab = "MLE BIC"
		plot(log10(en.lseq), BIC, pch=16, xlab=xlab, ylab=ylab)
		axis(side=3,at=log10(en.lseq),labels=en.model$df, tick=F)
		ylab = "MLE AIC"		
		plot(log10(en.lseq), AIC, pch=16, xlab=xlab, ylab=ylab)
		axis(side=3,at=log10(en.lseq),labels=en.model$df, tick=F)
	}

	return(list(aicSel=aicSel, bicSel=bicSel))

}

MLExic = function(XX, Z, alpha, sigset, criterion = 'bic', plotFlag = F, return.coef = F) {
	
	if (is.na(match(criterion,c('bic','aic')))) stop("criterion must be either 'aic' or 'bic'")

	en.model = glmnet(x=XX, y=factor(Z), alpha = alpha, family="binomial", lambda.min=0.0001)
	en.bHx = as.matrix(en.model$beta)
	en.lseq = en.model$lambda
	
	p = ncol(XX)
	T = length(en.lseq)
	varSel = matrix(0,p+1,T)
	beta = as.matrix(en.model$beta)
	P = plogis(en.model$a0 + XX%*%as.matrix(en.model$beta))

	dev = -2 * (t(Z) %*% log(P) + t(1-Z)  %*% log(1-P))
	if (criterion == 'bic') {
		xIC = dev + log(length(Z))*(en.model$df)
	} else {
		xIC = dev + 2*en.model$df
	}

	ul2e.OsdIX = which(xIC==min(xIC))
	varSel = which(abs(beta[,ul2e.OsdIX ]) > 0)

	if (plotFlag) {
		quartz(); par(mfrow=c(2,1))
		regpathPlot(log10(en.lseq), rbind(en.model$a0,beta), subset = sigset) 
		if (criterion == "bic") {
			ylab = "MLE BIC"
		} else {
			ylab = "MLE AIC"		
		}	
		xlab = expression(log[10](lambda))
		plot(log10(en.lseq), xIC, pch=16, xlab=xlab, ylab=ylab)
		axis(side=3,at=log10(en.lseq),labels=en.model$df, tick=F)
	}

	if (!return.coef)
		return(list(varSel=varSel))	

	bTHx = glm(Y ~ ., family=binomial(logit), data=data.frame(Y=Z, X = XX[,varSel,drop=F]))$coef
	return(list(varSel=varSel,bTHx=bTHx))	

}


MLEcv = function(XX, Z, alpha = 1, foldMatrix = matrix(c(),0,0), plotFlag = F) {
	en.model = glmnet(x=XX, y=factor(Z), alpha = alpha, family="binomial", lambda.min=0.0001)
	en.bHx = as.matrix(en.model$beta)
	en.lseq = en.model$lambda

#	foldMatrix
#	if (!is.matrix(foldMatrix)) stop("foldMatrix is not a matrix.")
	if (plotFlag)	
		quartz(); regpathPlot(log10(en.model$lambda), rbind(en.model$a0,en.bHx), subset = c(2))

	# Cross validation: 10 fold
	print(nrow(foldMatrix))
	bino.dev = matrix(NA,nrow(foldMatrix),length(en.lseq))
	for (i in 1:nrow(foldMatrix)) {
		cvix = foldMatrix[i,]
		en.boot = glmnet(x=XX[-cvix,,drop=F],y=factor(Z[-cvix]),alpha=alpha,family="binomial",lambda=en.lseq)
		P = plogis(en.boot$a0 + XX[cvix,,drop=F]%*%as.matrix(en.boot$beta))
		bino.dev[i,1:ncol(P)] = apply(P,2,FUN=function(x){2*sum(-Z[cvix]*log(x) - (1-Z[cvix])*(log(1-x)))})
		print(paste("Fold ",i))
	}

	naFreeIX = which(apply(bino.dev,2,FUN=function(x) {!any(is.na(x))}))

	en.mean = apply(bino.dev[,naFreeIX],2,mean)
	en.sd = apply(bino.dev[,naFreeIX],2,sd)
	minIX = which(en.mean == min(en.mean))
	OsdIX = min(which(en.mean < en.mean[minIX] + en.sd[minIX]))

	en.varSel = which(abs(en.bHx[, naFreeIX[OsdIX] ]) > 0)
	en.bTHx = glm(Y ~ .,family=binomial(logit),data=data.frame(Y = Z, X = XX[,en.varSel]))

	if (plotFlag) {
		ylim = c(min(bino.dev[,naFreeIX],na.rm=T),max(bino.dev[,naFreeIX],na.rm=T))
		ylab = "Deviance Loss"
		quartz()
		cvPlot(lambda=en.lseq[naFreeIX], df=en.model$df[naFreeIX], mx = en.mean, sdev=en.sd, ylim=ylim,ylab=ylab)
	}

	return(list(varSel=en.varSel,bTHx=en.bTHx))	

}

bin.dev = function(beta0,beta,x,y,lambda,alpha) {
	p = plogis(beta0 + x%*%beta)
	bd = -2*sum(y*log(p) + (1-y)*(log(1-p)))
	en = lambda*(alpha*sum(abs(beta)) + 0.5*(1-alpha)*sum(beta^2))
	return(bd + en)
}

#@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
# The following three functions are R implementations of coordinate descent.

coordinate_descent.l2e = function(X, Y, alpha, lambda, beta, w = 1, K = 0.16, niter = 500, check=F, tol=1e-3) {

# Need to add checks on all the dimensions.
	n = length(Y)
	p = length(beta)

	Ki = 1/K;
	U = 2*Y-1
	S = apply(X,2,FUN = function(x) {return(t(x)%*%x)} )
	m = w*K/n
	mS = m * S

	betaHx = beta

	R = X %*% beta
	P = plogis(R)
	D = 2*P - 1
	WZ = P * (1-P) * (w*D - U)
	v = beta
	
	if (!check) {
		for (t in 1:niter) {
			v[1] = beta[1] - Ki * mean(WZ)
			R = R + (v[1] - beta[1])
			beta[1] = v[1]
			P = plogis(R)
			D = 2*P - 1
			WZ = P * (1-P) * (w * D - U)
			for (j in 2:p) {
				v[j] = soft.th(mS[j] * v[j] - (w/n) * t(X[,j,drop=F]) %*% WZ, lambda*alpha) / (mS[j] + lambda*(1-alpha))
				R = R + (v[j] - beta[j])*X[,j,drop=F]
				beta[j] = v[j]
				P = plogis(R)
				D = 2*P - 1
				WZ = P * (1-P) * (w * D - U)
			}				
			betaHx = cbind(betaHx, beta)
		}
	} else {
		lastCost = mean.l2e(c(w,v),X,Y,lambda,alpha)
		for (t in 1:niter) {
			v[1] = beta[1] - Ki * mean(WZ)
			R = R + (v[1] - beta[1])
			beta[1] = v[1]
			P = plogis(R)
			D = 2*P - 1
			WZ = P * (1-P) * (w * D - U)
			for (j in 2:p) {
				v[j] = soft.th(mS[j] * v[j] - (w/n) * t(X[,j,drop=F]) %*% WZ, lambda*alpha) / (mS[j] + lambda*(1-alpha))
				R = R + (v[j] - beta[j])*X[,j,drop=F]
				beta[j] = v[j]
				P = plogis(R)
				D = 2*P - 1
				WZ = P * (1-P) * (w * D - U)
			}				
			curL2eCost = mean.l2e(c(w,v),X,Y,lambda,alpha)
			if (abs(curL2eCost - lastCost)/(abs(lastCost) + 1) > tol) {
				beta = v
				betaHx = cbind(betaHx, beta)
			} else {
				print(paste("Number of iterations = ", t))
				break;	
			}
		}	
	}
	return(betaHx)
}

coordinate_descent.l2eB = function(X, Y, alpha, lambda, beta, w = 1, K = 0.16, niter = 500, check=F) {

# Need to add checks on all the dimensions.
	n = length(Y)
	p = length(beta)

	Ki = 1/K;
	U = 2*Y-1

	betaHx = beta
	minL2eCost = Inf

	R = X %*% beta
	P = plogis(R)
	D = 2*P - 1
	WZ = P * (1-P) * (w*D - U)
	v = beta

	if (!check) {
		for (t in 1:niter) {
			v[1] = beta[1] - Ki * mean(WZ)
#			plot(beta, ylim=c(-20,20),pch=16,col=rgb(0,0,1,0.5))
#			points(beta[1:2],col=rgb(1,0,0,.5),pch=16)
			R = R + (v[1] - beta[1])
			beta[1] = v[1]
			P = plogis(R)
			D = 2*P - 1
			WZ = P * (1-P) * (w * D - U)
			for (j in 2:length(v)) {
				v[j] = soft.th(beta[j] - Ki * t(X[,j,drop=F]) %*% WZ, lambda*alpha) / (1 + lambda*(1-alpha))
				R = R + (v[j] - beta[j])*X[,j,drop=F]
				beta[j] = v[j]
				P = plogis(R)
				D = 2*P - 1
				WZ = P * (1-P) * (w * D - U)
			}				
			betaHx = cbind(betaHx, beta)
		}
	} else {
		for (t in 1:niter) {
			v[1] = v[1] - Ki * mean(WZ)
			R = R + (v[1] - beta[1])
			beta[1] = v[1]
			P = plogis(R)
			D = 2*P - 1
			WZ = P * (1-P) * (w * D - U)
			for (j in 2:length(v)) {
				v[j] = soft.th(v[j] - Ki * t(X[,j,drop=F]) %*% WZ, lambda*alpha) / (1 + lambda*(1-alpha))
				R = R + (v[j] - beta[j])*X[,j,drop=F]
				beta[j] = v[j]
				P = plogis(R)
				D = 2*P - 1
				WZ = P * (1-P) * (w * D - U)
			}		
			curL2eCost = mean.l2e(c(w,v),X,Y,lambda,alpha)
			if ( curL2eCost < minL2eCost) {
				minL2eCost = curL2eCost
				beta = v
				betaHx = cbind(betaHx, beta)
			} else {
				break;	
			}
		}	
	}
	
	return(betaHx)
}
