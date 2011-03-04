mean.l2e = function(par, X, Y, lambda, alpha) {
  w = par[1]
  theta = matrix(par[-1],ncol=1)
  p.case = plogis(X%*%theta)
  p.cntl = 1 - p.case
  l2e = 0.5*mean( w^2 * (p.case^2 + p.cntl^2) -2 * w * (p.case^Y)*(p.cntl^(1-Y)))
  ridg = 0.5*(1-alpha)*sum(theta[-1]^2)
  lass = alpha*sum(abs(theta[-1]))
  return(l2e + lambda*(ridg+lass))
}
