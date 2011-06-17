#include "util.h"
#include "ENLS.h"
#include "HHSVM.h"

double phiPrime(double x, double t) {
  double out = 0.;
  if (x <= t)
    out = -2.*(1.-t);
  else if (x > t && x <= 1)
    out = -2.*(1.-x);
  return out;
}

void updateWorkingResponse(double *X, double *Y, double beta0, double *beta, double *zeta, double *Xbeta,
                           double t, int n, int p) {
  const char transX = 'N';
  const double alp = 1.;
  const double bet = 0.;
  const int one = 1;
  int i;
  dgemv_(&transX, &n, &p, &alp, X, &n, beta, &one, &bet, Xbeta, &one);
  for (i = n; i--; )
    zeta[i] = Y[i] * phiPrime(Y[i] * (beta0 + Xbeta[i]), t);
}

void HHSVM(double *X, double *Y, double *beta0, double *beta, int *n, int *p,
	   double *lambda, double *alpha, double *t, int *maxiter, double *tol, 
	   double *beta_final) {
  int i, j, iter;
  double *Xbeta = calloc(*n, sizeof(double));
  double *zeta = calloc(*n, sizeof(double));
  double *betaLast = calloc(*p+1, sizeof(double));
  double *Xcentered = calloc((*n)*(*p), sizeof(double));
  double meanZeta, deltaBeta;
  int iIter, maxIter = 3;

  centerColumns(X, Xcentered, *n, *p);

  for (iter = *maxiter; iter--; ) {
    recordCurrentBeta(*beta0, beta, *p, betaLast);
    updateWorkingResponse(Xcentered, Y, *beta0, beta, zeta, Xbeta, *t, *n, *p);
    meanZeta = mean(zeta, *n);
    *beta0 = *beta0 - 0.5*meanZeta;
    for (i = *n; i--; )
      zeta[i] = Xbeta[i] - 0.5*(zeta[i] - meanZeta);
    //    ENLS(Xcentered, zeta, beta, n, p, lambda, alpha, maxiter, tol);
    iIter = ENLS(Xcentered, zeta, beta, n, p, lambda, alpha, &maxIter, tol);
    deltaBeta = twoNorm(*beta0, beta, betaLast, *p);
    Rprintf("oIter %3d: iIters %2d : ||beta - beta_old|| = %7.1e\n", *maxiter - iter, iIter, deltaBeta);

    // Check Outer Loop Convergence                                                                                             
    if (deltaBeta < *tol) {
      break;
    }
  }
  recordCurrentBeta(*beta0, beta, *p, beta_final);

  free(Xbeta);
  free(zeta);
  free(betaLast);
  free(Xcentered);
}
