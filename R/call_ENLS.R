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

source("stepSizeCalculator.R")

call_logisticL2E = function(X, Y, alpha, lambda, beta0, beta, w = 1, niter=500, tol=1e-9, save_flag=1, filename="betaHx.csv") {
	if (!is.matrix(X)) stop("X must be a matrix.")

	n = dim(X)[1]
	p = dim(X)[2]

	K = stepSizeCalculator(w)
		
	if (length(Y) != n | length(beta) !=p) stop("length(Y) == nrow(X) & length(beta) == ncol(X)")
		
	print(paste("alpha = ",alpha, ", lambda = ",lambda,", w = ", w, ",K = ",K))
	.C("logisticL2E", as.double(X), as.integer(Y), as.double(beta0), as.double(beta),
		as.integer(n), as.integer(p), as.double(K), as.double(w), as.double(lambda), as.double(alpha),
		as.integer(niter),as.double(tol),as.integer(save_flag), as.character(filename))
}
