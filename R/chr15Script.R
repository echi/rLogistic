source("rLogistic.R")
source("logisticL2E.R")
X = as.matrix(read.table("../Lung/X.csv.gz", sep=","))
Y = as.matrix(read.table("../Lung/Y.csv", sep=","))

alpha = 0.75

start = 5700
end = 5800

Z = X[,start:end,drop=F]
# Center SNP data
Z = X
Z = apply(Z, 2, FUN=function(x){return(x-mean(x))})

n = nrow(Z)
p = ncol(Z)

beta = matrix(0,p,1)

# Calculate lambda_max
lastHx = c()

beta0 = 0
beta = matrix(0,p,1)
U = 2*Y - 1
P = plogis(rep(1,n))
D = 2*P-1
WZ = P * (1-P) * (D- U)
zeta = beta0 - (1/0.16)*WZ

lambdaMax = max(abs(t(Z) %*%(zeta - mean(zeta))))*0.16/(n*alpha) #+ 0.015
lambdaMax = 1.5 * lambdaMax
epsilon = 0.05
lambdaMin = epsilon*lambdaMax
nSteps = 100
lambda = exp(seq(log(lambdaMax),log(lambdaMin),length.out=nSteps))

# Loop the loop
system.time({rLogistic.fit = rLogistic(Z,Y,alpha,lambda,beta0,beta)})

varSel = c()
for (i in 1:rLogistic.fit$dim[2]) {
	varSel = union(varSel, which(abs(rLogistic.fit$beta[,i]) > 0))
}

#filename = paste("alpha=",alpha,"start=",start,"end=",end)
#save(file=filename)

fMax = c()
fMin = c()
for (i in 1:rLogistic.fit$dim[2]) {
	fMax[i] = max(rLogistic.fit$beta[varSel,i])
	fMin[i] = min(rLogistic.fit$beta[varSel,i])
}

quartz()
plot(rLogistic.fit$lambda, fMax, type='n', ylim=c(min(fMin), max(fMax)),xlab=expression(lambda),ylab=expression(beta[j]))
for (i in 1:length(varSel)) {
	lines(rLogistic.fit$lambda, rLogistic.fit$beta[varSel[i],],lwd=2)
}

# Don't know where these came from!
#ixSNP1 = which(colnames(Z) == "rs8192475")
#ixSNP2 = which(colnames(Z) == "rs3885951")

ixSNP = list()
ixSNP[[1]] = which(colnames(Z) == "rs8034191")
ixSNP[[2]] = which(colnames(Z) == "rs1051730")
ixSNP[[3]] = 1081
ixSNP[[4]] = 4973
ixSNP[[5]] = 7932

# Color Brewer RdYlBu 5
colors = list()
colors[[1]] = rgb(215/255, 25/255, 28/255)
colors[[2]] = rgb(253/255, 174/255, 97/255)
colors[[3]] = rgb(255/255, 255/255, 191/255)
colors[[4]] = rgb(171/255, 217/255, 233/255)
colors[[5]] = rgb(44/255, 123/255, 182/255)

for (i in 1:length(ixSNP)) {
	lines(rLogistic.fit$lambda, rLogistic.fit$beta[ixSNP[[i]],], col=colors[[i]], lwd=2)
}

rug(rLogistic.fit$lambda)
title(main=paste("n =",n,"p =",p,"alpha =",alpha))

####------- Test elastic net
library(glmnet)

#alpha = 1
system.time({glmnet.fit = glmnet(Z,Y,family="binomial",alpha=alpha,standardize=T)})

# Fish out the SNPs Chris found
any(names(which(abs(glmnet.fit$beta[,ncol(glmnet.fit$beta)]) > 0)) %in% c("rs8034191"))
any(names(which(abs(glmnet.fit$beta[,ncol(glmnet.fit$beta)]) > 0)) %in% c("rs1051730"))

cutOff = min(which(glmnet.fit$df > 20))
subSample = 1:cutOff

varSel.gn = c()
for (i in subSample) {
	varSel.gn = union(varSel.gn, which(abs(glmnet.fit$beta[,i]) > 0))
}

fMax.gn = c()
fMin.gn = c()
for (i in subSample) {
	fMax.gn[i] = max(glmnet.fit$beta[varSel.gn,i])
	fMin.gn[i] = min(glmnet.fit$beta[varSel.gn,i])
}

quartz()
plot(glmnet.fit$lambda[subSample], fMax.gn, type='n', ylim=c(min(fMin.gn), max(fMax.gn)),xlab=expression(lambda),ylab=expression(beta[j]))
for (i in 1:length(varSel.gn)) {
	lines(glmnet.fit$lambda[subSample], glmnet.fit$beta[varSel.gn[i],subSample],lwd=2)
}

for (i in 1:length(ixSNP)) {
	lines(glmnet.fit$lambda[subSample], glmnet.fit$beta[ixSNP[[i]],subSample], col=colors[[i]], lwd=2)
}
rug(glmnet.fit$lambda[subSample])
title(main=paste("n =",n,"p =",p,"alpha =",alpha))

# BIC curves
# First calculate the degrees of freedom.
activeSet.gn = list()
lrm.fit.gn = list()
dof = c()
BIC.gn = c()
for (i in subSample) {

	activeSet.gn[[i]] = which(abs(glmnet.fit$beta[,i]) > 0)
	sas = length(activeSet.gn[[i]])
	if (length(activeSet.gn[[i]]) > 0) {
		pmatrix = glmnet.fit$lambda[i]*(1-alpha)*diag(1,sas,sas)
		Za = Z[,activeSet.gn[[i]],drop=F]
		lrm.fit.gn[[i]] = lrm(Y ~ Za, penalty.matrix= pmatrix)
		d = svd(Za)$d
		dof[i] = sum(d^2 / (d^2 + (1-alpha)*glmnet.fit$lambda[i]/2))
		L = lrm.fit.gn[[i]]$dev[2]
	} else {
		dof[i] = 0
		linearPredictor = glmnet.fit$a0[i]
		L = -2*sum(Y * linearPredictor - log(1 + exp(linearPredictor)))
	}
	BIC.gn[i] = L + log(n)*dof[[i]]
}

variables.gn = activeSet.gn[[min(which(BIC.gn == min(BIC.gn)))]]

# First calculate the degrees of freedom.
# BIC.L2E part 2
activeSet.L2E = list()
lrm.fit.L2E = list()
dof = c()
BIC.L2E = c()
for (i in 1:ncol(lastHx)) {
	varSel.L2E = which(abs(lastHx[,i]) > 0)
	activeSet.L2E[[i]] = varSel.L2E[-1] - 1
	sas = length(activeSet.L2E[[i]])	
	if (length(activeSet.L2E[[i]]) > 0) {
		Za = Z[,activeSet.L2E[[i]],drop=F]	
		beta0 = lastHx[1,i]
		beta = lastHx[activeSet.L2E[[i]]+1,i]
		system.time(call_logisticL2Eas(Za,Y,0,(1-alpha)*lambda[i],beta0,beta, tol=1e-6))

		rbHx = read.table("betaHx.csv", sep=",")
		beta = t(rbHx[nrow(rbHx),-1,drop=F])
		beta0 = rbHx[nrow(rbHx),1]
		linearPredictor = beta0 + Za %*% beta

		d = svd(Za)$d
		dof[i] = sum(d^2 / (d^2 + (1-alpha)*lambda[i]/2))
		
	} else {
		dof[i] = 0	
		beta0.L2E = lastHx[1,i]
		linearPredictor = beta0.L2E
	}
	L = sum(Y * linearPredictor - log(1 + exp(linearPredictor)))
	BIC.L2E[i] = -2*L + log(n)*dof[[i]]
}

variables.L2E = activeSet.L2E[[min(which(BIC.L2E == min(BIC.L2E)))]]