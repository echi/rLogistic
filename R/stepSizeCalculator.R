stepBound = function(x,w) {	
	return(-0.25 *(1.5*w*x^4 - x^3 - 2*w*x^2 + x + w/2))
}

grad = function(x,w) {
	return(-0.25*(6*w*x^3 - 3*x^2 - 4*w*x + 1))
}

hess = function(x,w) {
	return(matrix(-0.25*(18*w*x^2 - 6*x - 4*w),1,1))
}

stepSizeCalculator = function(w, debugFlag = F, plotFlag = F, wiggle = 0.000001) {
	
	if (w < 0 || w > 1) stop("w must be between 0 and 1")

	start = (-3+sqrt(33))/12
	lower = ifelse(w > 0, max(-1, (1 - sqrt(1+8*w^2))/(6*w)), 0)
	upper = min( 1, (1 + sqrt(1+8*w^2))/(6*w))

	solution = nlminb(start, objective=stepBound, gradient=grad, hessian=hess, lower=lower, upper=upper, w=w)

	if (plotFlag) {
		quartz()
		x = seq(-1,1,0.001)
		y = -1*apply(matrix(x,ncol=1),1,stepBound,w=w)
		yp = -1*apply(matrix(x,ncol=1),1,grad,w=w)
		ypp = -1*apply(matrix(x,ncol=1),1,hess,w=w)
		theRoot = solution$par
		par(mfrow=c(3,1))
		plot(x,y,type='l',ylab="g",main=paste("w = ",w)); abline(h=0,col='black',lty=2); abline(v=theRoot,col='red')
		abline(v = lower, col='blue'); abline(v = upper, col='blue')
		plot(x,yp,type='l',ylab="g'",main=paste("Maximizing Root = ",theRoot)); abline(h=0,col='black',lty=2); abline		(v=theRoot,col='red')
		abline(v = lower, col='blue'); abline(v = upper, col='blue')
		plot(x,ypp,type='l',ylab="g''",main=paste("Optimal K = ",-1*solution$objective)); abline(h=0,col='black',lty=2); 		abline(v=theRoot,col='red')
		abline(v = lower, col='blue'); abline(v = upper, col='blue')
	}

	if (debugFlag) {
		return(list(K = -1*solution$objective+wiggle, solution = solution))	} else {
		return(-1*solution$objective+wiggle)
	}
}

### Manual diagnostics for step size correctness
diagnosticPlots = F

if (diagnosticPlots) {
	
### 1. Derivative on left and right should have opposite signs.
	w = seq(1,0,length.out=1000)
	lower = c(); upper = c(); left = c(); right = c()
	for (i in 1:length(w)) {
		lower[i] = ifelse(w[i] > 0, max(-1, (1 - sqrt(1+8*w[i]^2))/(6*w[i])), 0)
		upper[i] = min( 1, (1 + sqrt(1+8*w[i]^2))/(6*w[i]))
		left[i] = -1*grad(lower[i],w[i])
		right[i] = -1*grad(upper[i], w[i])
	}

	quartz()
	plot(w, left, type='l', col='red', ylim=c(-1,1))
	lines(w, right, type='l', col='blue')
	abline(h=0,col='black',lty=2)
	legend('topright',legend=c("Gradient on left boundary", "Gradient on right boundary"),col=c('red','blue'),lty=1)

	quartz()
	plot(w, lower, type='l', col='red', ylim=c(-1,1), main='Optimal K is between red and blue lines',ylab='K')
	lines(w, upper, type='l', col='blue')
	abline(h=0,col='black',lty=2)

### 2. Plot the function, gradient, and hessian for a couple of w values.

	w = seq(0,1,length.out=5)
	for (i in 1:length(w)) {
		stepSizeCalculator(w[i], plotFlag = T)
	}
}
