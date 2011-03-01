# Calculate BIC for glmnet logistic models
#

library(Design)

BIC.glmnet = function(Z, Y, glmnet.model, alpha, modelSet) {

# First calculate the degrees of freedom.
	activeSet.gn = list()
	lrm.fit.gn = list()
	dof = c()
	BIC = c()
	bic.gn = list()

	for (i in modelSet) {

		activeSet.gn[[i]] = which(abs(glmnet.model$beta[,i]) > 0)
		sas = length(activeSet.gn[[i]])
		if (sas > 0) {
			pmatrix = glmnet.model$lambda[i]*(1-alpha)*diag(1,sas,sas)
			Za = Z[,activeSet.gn[[i]],drop=F]
			lrm.fit.gn[[i]] = lrm(Y ~ Za, penalty.matrix= pmatrix)
			d = svd(Za)$d
			dof[i] = sum(d^2 / (d^2 + (1-alpha)*glmnet.model$lambda[i]/2))
			L = lrm.fit.gn[[i]]$dev[2]
		} else {
			dof[i] = 0
			linearPredictor = glmnet.model$a0[i]
			L = -2*sum(Y * linearPredictor - log(1 + exp(linearPredictor)))
		}
		BIC[i] = L + log(n)*dof[i]
	}
	bic.gn$BIC = BIC
	bic.gn$dof = dof
	bic.gn$activeSets = activeSet.gn
	bic.model.index = min(which(BIC == min(BIC)))
	bic.gn$bestSet = activeSet.gn[[bic.model.index]]
	bic.gn$model = lrm.fit.gn[[bic.model.index]]
	
	return(bic.gn)
}