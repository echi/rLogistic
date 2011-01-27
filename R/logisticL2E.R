soft.th = function(x,lambda) {
	n = length(x)
	y = matrix(0,nrow=n)
	posIX = which(x > lambda)
	y[posIX] = x[posIX] - lambda
	negIX = which(x < -lambda)
	y[negIX] = x[negIX] + lambda
	return(y)
}

# X's first column does not have ones.
cdL2E = function(X, Y, alpha, lambda, beta0, beta, w = 1, K = 0.16, niterInner = 500, niterOuter = 500, tol=1e-3) {

	if (!is.matrix(X)) stop("X must be a matrix.")

	n = dim(X)[1]
	p = dim(X)[2]
		
	if (length(Y) != n | length(beta) !=p) stop("length(Y) == nrow(X) & length(beta) == ncol(X)")

	Ki = 1/K;
	U = 2*Y-1
	S = apply(X,2,FUN = function(x) {return(t(x)%*%x)} )
	m = w*K/n
	mS = m * S

	betaHx = c(beta0, beta)

	for (t in 1:niterOuter) {
		lastOuter = c(beta0, beta)
		Xbeta = X %*% beta
		P = plogis(beta0 + Xbeta)
		D = 2*P - 1
		WZ = P * (1-P) * (w * D - U)
		zeta = matrix(beta0 + Xbeta - Ki*WZ, ncol=1)
		for (iter in 1:niterInner) {
			lastInner = c(beta0, beta)
			# Update beta0 
			beta0 = mean(zeta - X %*% beta)
			# Update beta
			for (j in seq(p,1,-1)) {
				mask = matrix(1, p, 1)
				mask[j] = 0
				pResponse = beta0 + X %*% (mask * beta)
				pResidual = zeta - pResponse
				v = soft.th(m * t(X[,j,drop=F]) %*% pResidual, lambda*alpha)
				v = v / (mS[j] + lambda*(1-alpha))
				beta[j] = v
			}
			# Check inner loop convergence
			if (sqrt(sum(c(beta0,beta) - lastInner)^2) < tol) {
				print(paste("Number of inner iterations = ", iter))
				break
			}
		}
		betaHx = cbind(betaHx, c(beta0, beta))
		# Check outer loop convergence
		if (sqrt(sum(c(beta0,beta) - lastOuter)^2) < tol) {
			print(paste("Number of outer iterations = ", t))
			break
		}
	}
	return(betaHx)
}

# Add non-coordinate descent version in this file!