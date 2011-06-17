#include "util.h"
#include "ENLS.h"
#include "logisticL2E.h"

inline double cplogis(double x) {
  return 0.5 * (1. + tanh(0.5 * x));
}

void updateWorkingResponse(double *X, double *U, double w, double K,
                           double beta0, double *beta, double *zeta, double *Xbeta,
                           int n, int p) {
  const char transX = 'N';
  const double alp = 1.;
  const double bet = 0.;
  const int one = 1;
  double P, D, WZ;
  int i;
  dgemv_(&transX, &n, &p, &alp, X, &n, beta, &one, &bet, Xbeta, &one);
  for (i = n; i--; ) {
    P = cplogis(beta0 + Xbeta[i]);
    D = P + P - 1.;
    WZ = P * (1. - P) * (w * D - U[i]);
    zeta[i] = beta0 + Xbeta[i] - (1. / K) * WZ;
  }
}

void logisticL2E(double *X, int *Y, double *beta0, double *beta, int *n, int *p,
		 double * K, double *w, double *lambda, double *alpha,
		 int *maxiter, double *tol, double *beta_final) {
  int i, j, iter;
  double newLambda, deltaBeta;
  double *Xbeta = calloc(*n, sizeof(double));
  double *U = calloc(*n, sizeof(double));
  double *zeta = calloc(*n, sizeof(double));
  double *betaLast = calloc(*p+1, sizeof(double));
  double *Xcentered = calloc((*n)*(*p), sizeof(double));
  int iIter, maxIter = 3;

  for (i = *n; i--; )

    U[i] = (double) (Y[i] + Y[i] - 1);

  centerColumns(X, Xcentered, *n, *p);

  newLambda = *lambda / ( (*w) * (*K) );

  for (iter = *maxiter; iter--; ) {
    recordCurrentBeta(*beta0, beta, *p, betaLast);
    updateWorkingResponse(Xcentered, U, *w, *K, *beta0, beta, zeta, Xbeta, *n, *p);
    *beta0 = mean(zeta, *n);
    for (i = *n; i--; )
      zeta[i] = zeta[i] - *beta0;
    //    ENLS(Xcentered, zeta, beta, n, p, &newLambda, alpha, maxiter, tol);
    iIter = ENLS(Xcentered, zeta, beta, n, p, &newLambda, alpha, &maxIter, tol);
    deltaBeta = twoNorm(*beta0, beta, betaLast, *p);

    Rprintf("oIter %3d: iIters %2d : ||beta - beta_old|| = %7.1e\n", *maxiter - iter, iIter, deltaBeta);

    // Check Outer Loop Convergence
    if (deltaBeta < *tol) {
      break;
    }
  }
  recordCurrentBeta(*beta0, beta, *p, beta_final);

  free(Xbeta);
  free(U);
  free(zeta);
  free(betaLast);
  free(Xcentered);
}
