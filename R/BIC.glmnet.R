# Calculate BIC for glmnet logistic models
#

BIC.glmnetB = function(Z, Y, glmnet.model, alpha, modelSet, reducer='median') {

# First calculate the degrees of freedom.
	activeSet.gn = list()
	lrm.fit.gn = list()
	dof = c()
	BIC = c()
	L = c()
	bic.gn = list()

	if (reducer == 'median') {

		for (i in modelSet) {	

			activeSet.gn[[i]] = which(abs(glmnet.model$beta[,i]) > 0)
			sas = length(activeSet.gn[[i]])
			lambda = glmnet.model$lambda[i]*(1-alpha)
			if (sas > 0) {
				beta0 = glmnet.model$a0[i]
				beta = matrix(glmnet.model$beta[activeSet.gn[[i]],i], ncol=1)
				Za = Z[,activeSet.gn[[i]],drop=F]
				ridge = ridge_logistic(Za, Y, n*lambda, beta0, beta)
				beta0 = ridge[[1]]
				beta = ridge[[2]]
				Zb = Za %*% beta
				L[i] = median(Y * Zb - log(1 + exp(Zb)))
				d = svd(Za)$d
				dof[i] = sum(d^2 / (d^2 + n*lambda))
			
				mod = list()
				mod$a0 = beta0
				mod$beta = beta
				mod$df = dof[i]
				mod$lambda = lambda
				lrm.fit.gn[[i]] = mod
			} else {
				dof[i] = 0
				linearPredictor = glmnet.model$a0[i]
				L[i] = median(Y * linearPredictor - log(1 + exp(linearPredictor)))
				nullModel = list()
				nullModel$a0 = linearPredictor
				nullModel$beta = NULL
				nullModel$df = 0
				nullModel$lambda = lambda
				lrm.fit.gn[[i]] = nullModel
			}
			BIC[i] = -2*L[i] + log(n)*dof[i]/n
		}
	} else if (reducer == 'mean') {
		for (i in modelSet) {	

			activeSet.gn[[i]] = which(abs(glmnet.model$beta[,i]) > 0)
			sas = length(activeSet.gn[[i]])
			lambda = glmnet.model$lambda[i]*(1-alpha)
			if (sas > 0) {
				beta0 = glmnet.model$a0[i]
				beta = matrix(glmnet.model$beta[activeSet.gn[[i]],i], ncol=1)
				Za = Z[,activeSet.gn[[i]],drop=F]
				ridge = ridge_logistic(Za, Y, lambda, beta0, beta)
				beta0 = ridge[[1]]
				beta = ridge[[2]]
				Zb = Za %*% beta
				L[i] = sum(Y * Zb - log(1 + exp(Zb)))
				d = svd(Za)$d
				dof[i] = sum(d^2 / (d^2 + n*lambda))
			
				mod = list()
				mod$a0 = beta0
				mod$beta = beta
				mod$df = dof[i]
				mod$lambda = lambda
				lrm.fit.gn[[i]] = mod
			} else {
				dof[i] = 0
				linearPredictor = glmnet.model$a0[i]
				L[i] = sum(Y * linearPredictor - log(1 + exp(linearPredictor)))
				nullModel = list()
				nullModel$a0 = linearPredictor
				nullModel$beta = NULL
				nullModel$df = 0
				nullModel$lambda = lambda
				lrm.fit.gn[[i]] = nullModel
			}
			BIC[i] = (-2*L[i] + log(n)*dof[i])/n
		}		
	} else {
		# Add error to generate a message that this reducer is not supported.
		
	}

	bic.gn$BIC = BIC
	bic.gn$dof = dof
	bic.gn$activeSets = activeSet.gn
	bic.model.index = min(which(BIC == min(BIC)))
	bic.gn$bestSet = activeSet.gn[[bic.model.index]]
	bic.gn$model = lrm.fit.gn[[bic.model.index]]
	
	return(bic.gn)
}