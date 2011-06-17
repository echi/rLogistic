library(Design)

# Calculate BIC for HHSVM models
#

BIC.hhsvm = function(Z, Y, hhsvm.model, alpha) {

# First calculate the degrees of freedom.
	activeSet = list()
	lrm.fit.hhsvm = list()
	dof = c()
	BIC = c()
	L = c()
	success = c()
	bic.HHSVM = list()
	n = nrow(Z)

	for (i in 1:hhsvm.model$dim[2]) {

		activeSet[[i]] = which(abs(hhsvm.model$beta[,i]) > 0)
		sas = length(activeSet[[i]])
		if (sas > 0) {
#			lambda = n*(1-alpha)*hhsvm.model$lambda[i]
			lambda = hhsvm.model$lambda[i]
#			pmatrix = lambda*diag(1,sas,sas)
			pmatrix = (1-alpha)*lambda*diag(1,sas,sas)
			Za = Z[,activeSet[[i]],drop=F]
			lrm.fit.hhsvm[[i]] = try(lrm(Y ~ Za, penalty.matrix= pmatrix), TRUE)
			if (inherits(lrm.fit.hhsvm[[i]], "try-error")) {
				success[i] = 0
				break
			}
			success[i] = 1			
			d = svd(Za)$d
			dof[i] = sum(d^2 / (d^2 + (1-alpha)*lambda*sas/2))
#			L[i] = lrm.fit.hhsvm[[i]]$dev[2	
			linearPredictor = predict(lrm.fit.hhsvm[[i]])
			L[i] = median(Y * linearPredictor - log(1 + exp(linearPredictor)))
			L[i] = mean(Y * linearPredictor - log(1 + exp(linearPredictor)))
		} else {
			dof[i] = 0
			linearPredictor = hhsvm.model$a0[i]
			L[i] = median(Y * linearPredictor - log(1 + exp(linearPredictor)))
			L[i] = mean(Y * linearPredictor - log(1 + exp(linearPredictor)))
			success[i] = 1
		}
		BIC[i] = -2*L[i] + log(n)*dof[i]/n
	}
	bic.HHSVM$BIC = BIC
	bic.HHSVM$dof = dof
	bic.HHSVM$activeSets = activeSet
	bic.model.index = min(which(BIC == min(BIC)))
	bic.HHSVM$bestSet = activeSet[[bic.model.index]]
	bic.HHSVM$model = lrm.fit.hhsvm[[bic.model.index]]
	
	print(sucess)
	
	return(bic.HHSVM)
}