source("call_discardRules.R")

set.seed(1)
X = matrix(rnorm(100),50,2)
beta = matrix(4,2,1)
Y = X %*% beta + matrix(rnorm(50),50,1)
X = cbind(X, matrix(rnorm(200),50,4))
alpha = 0.95
#lambda = 180

scores = abs(t(X)%*%Y)
lambdaMax = max(scores)
normY = sqrt(sum(Y^2))

keepSet = c()
threshold = c()
for (j in 1:length(scores)) {
	varXj = sum(X[,j]^2)
	threshold[j] = alpha*lambda - normY*sqrt(varXj + (1.-alpha)*lambda)*(lambdaMax-lambda)/lambdaMax;
	if (scores[j] > threshold[j])
		keepSet = union(keepSet, j)
}

outputSAFE = selectVariablesSAFE(X,Y,alpha,lambda)
outputSAFE$keepSet[1:outputSAFE$sizeKeepSet]
outputSTRONG = selectVariablesSTRONG(X,Y,alpha,lambda)
outputSTRONG$keepSet[1:outputSTRONG$sizeKeepSet]

