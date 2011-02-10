soft.th = function(x,lambda) {
	n = length(x)
	return(sign(x) * max(abs(x) - lambda, 0))
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
#-	m = w*K/n #-
	m = 1/n   #+
	mS = m * S
	lambda = lambda/(w*K) #+

	betaHx = c(beta0, beta)

	activeSet = 1:p
	inActiveSet = setdiff(1:p, activeSet)

	for (t in 1:niterOuter) {
		lastOuter = c(beta0, beta)
		Xbeta = X %*% beta
		P = plogis(beta0 + Xbeta)
		D = 2*P - 1
		WZ = P * (1-P) * (w * D - U)
		zeta = matrix(beta0 + Xbeta - Ki*WZ, ncol=1)

#+		beta0 = mean(zeta)  				#+
#+		zeta = zeta - beta0 				#+
		repeat {
			print("Active Set:")
			print(activeSet)
			print("Inactive Set:")
			print(inActiveSet)
			if (length(activeSet) > 0) {
				Xa = X[, activeSet, drop=F]
				betaActive = beta[activeSet]
				for (iter in 1:niterInner) {
					lastInner = c(beta0, betaActive)
					# Update beta0 
					beta0 = mean(zeta - Xa %*% betaActive)        #-
					# Update beta
					for (j in 1:length(activeSet)) {
						pResponse = beta0 + Xa[,-j,drop=F] %*% betaActive[-j] #-
#+						pResponse = Xa[,-j,drop=F] %*% betaActive[-j]         #+
						pResidual = zeta - pResponse
						v = soft.th(m * t(Xa[,j,drop=F]) %*% pResidual, lambda*alpha)
						v = v / (mS[activeSet[j]] + lambda*(1-alpha))
						betaActive[j] = v
					}
					# Check inner loop convergence
					if (sqrt(sum(c(beta0, betaActive) - lastInner)^2) < tol) {
						print(paste("Number of inner iterations = ", iter))
						break
					}
				}
				beta[activeSet] = betaActive
			}

			pResponse = beta0 + X %*% beta                #-
#+			pResponse = Xa[,-j,drop=F] %*% betaActive[-j] #+
			pResidual = zeta - pResponse
			for (j in 1:length(inActiveSet)) {
				v = soft.th(m * t(X[,inActiveSet[j],drop=F]) %*% pResidual, lambda*alpha)
				v = v / (mS[inActiveSet[j]] + lambda*(1-alpha))
				beta[inActiveSet[j]] = v
			}
			# Remove variables from the active set
			deActiveSet = activeSet[which(abs(betaActive) < 1e-16)]
			activeSet = setdiff(activeSet, deActiveSet)
			# Add variables to inactive set
			inActiveSet = union(inActiveSet, deActiveSet)
			# Check inactive variables if they need to be added.
			addActiveSet = inActiveSet[which(abs(beta[inActiveSet]) > 0)]
			print(addActiveSet)
			if (length(addActiveSet) == 0)
				break
			# Remove variables from inActiveSet
			inActiveSet = setdiff(inActiveSet, addActiveSet)
			# Add variables to activeSet
			activeSet = union(activeSet, addActiveSet)
			if (length(activeSet) == 0)
				break
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