if (!is.loaded("logisticL2E")) {
	dyn.load("logisticL2E.so")	
} else {
	dyn.unload("logisticL2E.so")
	dyn.load("logisticL2E.so")	
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
