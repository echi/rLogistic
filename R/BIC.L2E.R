source("mean.l2e.R")

# Calculate BIC for L2E logistic models
#

BIC.L2E = function(Z, Y, L2E.model, alpha, w = 1) {

# First calculate the degrees of freedom.
	activeSet = list()
	lrm.fit.L2E = list()
	dof = c()
	BIC = c()
	bic.L2E = list()

	for (i in 1:L2E.model$dim[2]) {

		activeSet[[i]] = which(abs(L2E.model$beta[,i]) > 0)
		sas = length(activeSet[[i]])
		if (sas > 0) {
			Za = Z[,activeSet[[i]],drop=F]
			beta0 = L2E.model$a0[i]
			beta = L2E.model$beta[activeSet[[i]],i]
			lambda = (1-alpha)*L2E.model$lambda[i]
			theta0 = c(w,matrix(0, nrow=length(activeSet[[i]])+1, ncol=1))
			upper = c(w, rep(Inf, p+1)); lower = c(w, rep(-Inf, p+1))
			X = cbind(1, Za)
			opt.mean.l2e = nlminb(start=theta0, mean.l2e, X=X, Y=Y, lambda=lambda, alpha=0,
			lower=lower, upper=upper)
			fit = list()
			fit$a0 = opt.mean.l2e$par[2]
			fit$beta = opt.mean.l2e$par[-c(1,2)]
			names(fit$beta) = names(activeSet[[i]])
			fit$df = length(activeSet[[i]])	
			fit$lambda = lambda
#			fit = rLogistic(Za, Y, 0, lambda, beta0, beta)
			lrm.fit.L2E[[i]] = fit
			beta0 = fit$a0
			beta = fit$beta
			linearPredictor = beta0 + Za %*% beta
			d = svd(Za)$d
			dof[i] = sum(d^2 / (d^2 + lambda/2))
		} else {
			dof[i] = 0
			linearPredictor = L2E.model$a0[i]
			nullModel = list()
			nullModel$a0 = linearPredictor
			nullModel$beta = NULL
			nullModel$df = dof[i]
			nullModel$lambda = lambda
			nullModel$dim = NULL
			lrm.fit.L2E[[i]] = nullModel
		}
		L = median(Y*linearPredictor - log(1 + exp(linearPredictor)))
		BIC[i] = -2*L + log(n)*dof[i]/n
	}
	bic.L2E$BIC = BIC
	bic.L2E$dof = dof
	bic.L2E$activeSets = activeSet
	bic.model.index = min(which(BIC == min(BIC)))
	bic.L2E$bestSet = activeSet[[bic.model.index]]
	bic.L2E$model = lrm.fit.L2E[[bic.model.index]]
	
	return(bic.L2E)
}