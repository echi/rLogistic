source("call_logisticL2Eas.R")

rLogistic = function(X, Y, alpha, lambda, beta0, beta, w = 1, K = 0.16, niter=500, tol=1e-6, maxVar = 40) {

	modelID = c()
	for (i in 1:length(lambda)) {
		modelID[i] = paste("s", i-1,sep="")
	}
	model = list()
	model$a0 = numeric(length(lambda))
	names(model$a0) = modelID
	model$beta = matrix(NA, p, length(lambda))
	rownames(model$beta) = colnames(X)
	colnames(model$beta) = modelID
	model$df = rep(NA,length(lambda))
	model$lambda = lambda

	nIter = 0
	for (i in 1:length(lambda)) {
		nIter = nIter + 1
		call_logisticL2Eas(X, Y, alpha, lambda[i], beta0, beta, w=w, K=K, niter=niter, tol=tol)
		rbHx = read.table("betaHx.csv", sep=",")
		beta = t(rbHx[nrow(rbHx),-1,drop=F])
		beta0 = rbHx[nrow(rbHx),1]
		model$a0[i] = beta0
		model$beta[,i] = beta
		varSel = which(abs(beta) > 0)
		model$df[i] = length(varSel)
		if (length(varSel) > maxVar)
			break
	}

	model$a0 = model$a0[1:nIter]
	model$beta = model$beta[,1:nIter]
	model$df = model$df[1:nIter]
	model$lambda = model$lambda[1:nIter]
	model$dim = dim(model$beta)
	return(model)
}