if (!is.loaded("selectVariablesSAFE")) {
	dyn.load("discardRules.so")	
} else {
	dyn.unload("discardRules.so")
	dyn.load("discardRules.so")	
}

selectVariablesSAFE = function(X, Y, alpha, lambda) {
	if (!is.matrix(X)) stop("X must be a matrix.")

	n = dim(X)[1]
	p = dim(X)[2]
		
	if (length(Y) != n) stop("length(Y) == nrow(X)")
		
	print(paste("alpha = ",alpha, ", lambda = ",lambda))
	keepSet = integer(p)
	sizeKeepSet = integer(1)
	.C("selectVariablesSAFE", as.double(X), as.double(Y), keepSet = as.integer(keepSet), sizeKeepSet = as.integer(sizeKeepSet), as.integer(n), as.integer(p), as.double(lambda), as.double(alpha))
}

selectVariablesSTRONG = function(X, Y, alpha, lambda) {
	if (!is.matrix(X)) stop("X must be a matrix.")

	n = dim(X)[1]
	p = dim(X)[2]
		
	if (length(Y) != n) stop("length(Y) == nrow(X)")
		
	print(paste("alpha = ",alpha, ", lambda = ",lambda))
	keepSet = integer(p)
	sizeKeepSet = integer(1)
	.C("selectVariablesSTRONG", as.double(X), as.double(Y), keepSet = as.integer(keepSet), sizeKeepSet = as.integer(sizeKeepSet),as.integer(n), as.integer(p), as.double(lambda), as.double(alpha))
}