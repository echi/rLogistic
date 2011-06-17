source("stepSizeCalculator.R")

calculateLambdaMax = function(Y,Z,alpha,w) {

	U = 2*Y - 1
	P = 0.5 + (mean(Y) - 0.5)/w
	D = 2*P-1
	WZ = P * (1-P) * (w*D- U)
	beta0 = log(P/(1-P))
	K = stepSizeCalculator(w)
	zeta = beta0 - (1/K)*WZ

	return(max(abs(t(Z) %*%(zeta - mean(zeta))))*K/(nrow(Z)*alpha))
}

calculateLambdaSequence = function(Y,Z,alpha,w,epsilon=0.05,nLambda=50) {

	lambdaMax = calculateLambdaMax(Y,Z,alpha,w)
	lambdaMin = epsilon*lambdaMax
	return(exp(seq(log(lambdaMax),log(lambdaMin),length.out=nLambda)))
}
