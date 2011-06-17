function(v1, v2) {	
	return (max(abs(v1 - v2) / (1 + abs(v1))))
}

function(X, Y, lambda, beta0, beta, maxiter=1000, tol=1e-10) {
	
	n = dim(X)[1]; p = dim(X)[2]
	X = scale(X, center=T, scale=F)
	svdX = svd(X)
	U = svdX$u
	V = svdX$v
	d = svdX$d
	d = d / (d^2 + 4*lambda)
	
	for (iter in 1:maxiter) {
		P = plogis(beta0 + X %*% beta)
		beta0_last = beta0; beta_last = beta
		# Update beta0
		beta0 = beta0 + 4*mean(Y-P)
	
		# Update beta
		z = X %*% beta + 4 * (Y - P)
		beta = V %*% (d * (t(U) %*% z))

		if (rerr(c(beta0_last, beta_last), c(beta0, beta)) < tol)
			break
	}
	return(list(b0 = beta0, beta=beta, iters=iter))
	
}