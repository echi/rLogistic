if (!is.loaded("ENLS")) {
	dyn.load("ENLS.so")	
} else {
	dyn.unload("ENLS.so")
	dyn.load("ENLS.so")	
}

ENLS = function(X, Y, alpha, lambda, maxiter = 1000, tol=1e-6) {
	if (!is.matrix(X)) stop("X must be a matrix.")

	n = dim(X)[1]
	p = dim(X)[2]
		
	if (length(Y) != n) stop("length(Y) == nrow(X)")
	
	print(paste("alpha = ",alpha, ", lambda = ",lambda))
	beta = matrix(rnorm(p),p,1)
	Y = Y - mean(Y)
	X = apply(X,2,FUN=function(x){x - mean(x)})
	.C("ENLS", as.double(X), as.double(Y), beta = as.double(beta), as.integer(n), as.integer(p), as.double(lambda), as.double(alpha),
		as.integer(maxiter), as.double(tol))$beta
}