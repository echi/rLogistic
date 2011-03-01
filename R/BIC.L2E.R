# Calculate BIC for L2E logistic models
#

BIC.L2E = function(Z, Y, L2E.model, alpha) {

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
			fit = rLogistic(Za, Y, 0, lambda, beta0, beta)
			lrm.fit.L2E[[i]] = fit
			beta0 = fit$a0
			beta = fit$beta
			linearPredictor = beta0 + Za %*% beta
			d = svd(Za)$d
			dof[i] = sum(d^2 / (d^2 + lambda/2))
		} else {
			dof[i] = 0
			linearPredictor = rLogistic.fit$a0[i]
		}
		L = sum(Y * linearPredictor - log(1 + exp(linearPredictor)))
		BIC[i] = -2*L + log(n)*dof[i]
	}
	bic.L2E$BIC = BIC
	bic.L2E$dof = dof
	bic.L2E$activeSets = activeSet
	bic.model.index = min(which(BIC == min(BIC)))
	bic.L2E$bestSet = activeSet[[bic.model.index]]
	bic.L2E$model = lrm.fit.L2E[[bic.model.index]]
	
	return(bic.L2E)
}