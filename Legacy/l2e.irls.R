# Description: This function calculates the L2E logistic regression coefficients (no regularization).
#              It uses the majorization described in Eric's thesis. The initial starting position is always the zero vector.
#
# Author: Eric Chi
# Date: 5/31/2010
#
# Usage:
#
#     l2e.irls(X,Y,w=1,noi=20000, K=0.16)
#
# Input:
#  X		An n-by-p matrix. Rows are observations; columns are covariates. X'X should be non-singular.
#  Y 		A vector of binary response variables i.e. should be a vector of zeros and ones.
#  w 		A value between zero and 1; it represents the fraction of data that is not contaminant.
#  noi 		The number of iterations.
#  K		A parameter that controls the step size. See Eric's thesis for discussion on choosing its value.
#  path 	A boolean variable that specifies whether the returns the regression coefficients at each iteration
#         	or the final iteration.
# 
# Output:
#     betaHx is a p+1-by-noi matrix containing the regression coefficients computed at each iteration.
#     beta is a p+1 length vector of the regression coefficients at the last iteration.
#
# Modification History:


l2e.irls = function(X, Y, w=1, noi=20000, K = 0.16, path=T) {

	# Checks
	# X'X is non-singular
	
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
			D = 2*P - 1
			WZ = P * (1 - P) * (w * D - U)
			beta = beta - (1/K) * M %*% WZ
			betaHx = cbind(betaHx, beta)	
		}
		return(betaHx)
	} else {
		for (t in 1:noi) {
			P = plogis(X %*% beta)
			D = 2*P - 1
			WZ = P * (1 - P) * (w * D - U)
			beta = beta - (1/K) * M %*% WZ
		}	
		return(beta)
	}
}