# Description:  This function calculates the L2E logistic regression
#				coefficients (no regularization). It uses coordinate descent
#				with a fixed number of iterations. See coordinateDescentFunctions.R
#				for the coordinate descent solver that tracks the objective function
#				as a stopping rule. This function can be thought as coordinate descent
#				light.
#
#				It uses the majorization described in Eric's thesis. The initial
#				starting position is always the zero vector.
#
# Author: Eric Chi
# Date: 6/1/2010
#
# Usage:
#
#     l2e.cd(X, Y, w=1, noi=1000, K=0.16, path=T)
#
# Input:
#  X		An n-by-p matrix. Rows are observations; columns are covariates.
#			X'X should be non-singular.
#  Y 		A vector of binary response variables i.e. should be a vector of
#			zeros and ones.
#  w 		A value between zero and 1; it represents the fraction of data
#			that is not contaminant.
#  noi 		The number of iterations.
#  K		A parameter that controls the step size. See Eric's thesis for
#			discussion on choosing its value.
#  path 	A boolean variable that specifies whether the returns the regression
#			coefficients at each iteration or the final iteration.
# 
# Output:
#     betaHx is a p+1-by-noi matrix containing the regression coefficients computed at each iteration.
#     beta is a p+1 length vector of the regression coefficients at the last iteration.
#
# Modification History:


l2e.cd = function(X, Y, w=1, noi=1000, K = 0.16, path=T) {

	# Checks
	# X is a matrix
	if (!is.matrix(X)) stop("X must be a matrix.")
	n = nrow(X); p = ncol(X)
	
	# Y is zero/one vector
	if(any(is.na(match(unique(Y),c(0,1))))) stop("Y must be a numeric of zeros and ones.")

	# w is between 0 and 1
	if (w < 0 | w > 1) stop("w must be between 0 and 1")

	# X and Y have the right dimensions.
	if (length(Y) != n) stop("nrow(X) must equal length(Y)")

	beta = matrix(0,p+1,1)
	X = cbind(1,X)
	U = 2*Y - 1
	M = solve(t(X) %*% X) %*% t(X)

	if (path) { 	# Save the regression coefficients at each iteration or not.
		betaHx = beta
		for (t in 1:noi) {
			P = plogis(X %*% beta)
			W = P * (1 - P)
			D = 2*P - 1
			v = beta
			v[1] = v[1] - (1/K) * sum(W * (w*D - U)) /(2*n)
			for (j in 2:(p+1)) {
				v[j] = v[j] -(1/K) * X[,j] %*% (W * (w * D - U)) / (X[,j] %*% X[,j])
			}
			beta = v
			betaHx = cbind(betaHx, beta)	
		}
		return(betaHx)
	} else {
		for (t in 1:noi) {
			P = plogis(X %*% beta)
			W = P * (1 - P)
			D = 2*P - 1
			v = beta
			v[1] = v[1] - (1/K) * sum(W * (w*D - U)) /(2*n)
			for (j in 2:(p+1)) {
				v[j] = v[j] -(1/K) * X[,j] %*% (W * (w * D - U)) / (X[,j] %*% X[,j])
			}
			beta = v
		}	
		return(beta)
	}
}