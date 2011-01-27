#include <R.h> 
#include <Rmath.h> 
#include <R_ext/BLAS.h>
#include <math.h>
#include <stdio.h>

inline double softThreshold(double x, double lambda) {
  return(sign(x) * fmax2( fabs(x) - lambda, 0.));
}

inline double cplogis(double x) {
  return 0.5 * (1. + tanh(0.5 * x));
}

void writeParametersToFile(FILE *pFile, double beta0, double *beta, int p) {
  int j;
  fprintf(pFile,"%g, ", beta0);
  for (j = 0; j < p - 1; j++)
    fprintf(pFile,"%g, ", beta[j]);
  fprintf(pFile,"%g\n", beta[p-1]);
}

void recordCurrentBeta(double beta0, double *beta, int p, double *betaOld) {
  int j;
  for (j = p; j--; )
    betaOld[j+1] = beta[j];
  betaOld[0] = beta0;
}

int hasConverged(double beta0, double *beta, double *betaOld, int p, double tol) {
  int j;
  double delta = beta0 - betaOld[0];
  double norm = delta * delta;
  for (j = p; j--; ) {
    delta = beta[j] - betaOld[j+1];
    norm += delta * delta;
  }
  if (sqrt(norm) < tol)
    return 1;
  return 0;
}  

void updateIntercept(double *X, double *beta0, double *beta, double *zeta, double *Xbeta, int n, int p) {
  const char transX = 'N';
  const double alp = 1.;
  const double bet = 0.;
  const int one = 1;
  int i;
  *beta0 = 0.;
  dgemv_(&transX, &n, &p, &alp, X, &n, beta, &one, &bet, Xbeta, &one);
  for (i = n; i--; )
    *beta0 += (zeta[i] - Xbeta[i]);
  *beta0 /= (double) n;
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

double innerProduct(double *x, double *y, int n) {
  const int one = 1;
  return ddot_(&n, x, &one, y, &one);
}

void updateBeta(double *X, double beta0, double *beta,
		double *zeta, double *Xbeta, double *partialResidual,
		int n, int p, double lambda, double alpha, double m, double *mS) {
  int i, j;
  const char transX = 'N';
  const double alp = 1.;
  const double bet = 0.;
  const int one = 1;
  double la = lambda * alpha;
  double loma = lambda * (1. - alpha);
  
  for (j = p; j--; ) {
    beta[j] = 0.;
    dgemv_(&transX, &n, &p, &alp, X, &n, beta, &one, &bet, Xbeta, &one);
    for (i = n; i--; )
      partialResidual[i] = zeta[i] - beta0 - Xbeta[i];
    beta[j] = softThreshold(m*innerProduct(&X[n*j], partialResidual, n), la) / (mS[j] + loma);  
  }
}

/* Calculates the logistic L2E loss. X is n-by-p where p is the number of covariates. Note
   Note that the first column of X is _not_ 1.
 */
void logisticL2E(double *X, int *Y, double *beta0, double *beta, int *n, int *p, 
		  double * K, double *w, double *lambda, double *alpha,
		  int *maxiter, double *tol,
		  int *save_flag, char const *filename) {
  int i, j;
  int iterOuter, iterInner;
  double *Xbeta = calloc(*n, sizeof(double));
  double *U = calloc(*n, sizeof(double));
  double *mS = calloc(*p, sizeof(double));
  double *zeta = calloc(*n, sizeof(double));
  double *lastOuterBeta = calloc(*p+1, sizeof(double));
  double *lastInnerBeta = calloc(*p+1, sizeof(double));
  double *partialResidual = calloc(*n, sizeof(double));
  double m;
  FILE *pFile;

  pFile = fopen("betaHx.csv","w");
  writeParametersToFile(pFile, *beta0, beta, *p);

  m = (*w)*(*K)/((double)*n);
  for (j = *p; j--; )
    mS[j] = m*innerProduct(&X[(*n)*j], &X[(*n)*j], *n);

  for (i = *n; i--; )
    U[i] = (double) (Y[i] + Y[i] - 1);
  
  for (iterOuter = *maxiter; iterOuter--; ) {
    recordCurrentBeta(*beta0, beta, *p, lastOuterBeta);

    updateWorkingResponse(X, U, *w, *K, *beta0, beta, zeta, Xbeta, *n, *p);

    for (iterInner = *maxiter; iterInner--; ) {
      recordCurrentBeta(*beta0, beta, *p, lastInnerBeta);

      updateIntercept(X, beta0, beta, zeta, Xbeta, *n, *p);
      updateBeta(X, *beta0, beta, zeta, Xbeta, partialResidual, *n, *p, *lambda, *alpha, m, mS);
      
      // Check Inner Loop Convergence
      if (hasConverged(*beta0, beta, lastInnerBeta, *p, *tol)) {
	Rprintf("Number of inner iterations = %d.\n", *maxiter - iterInner);
	break;
      }
    }
    writeParametersToFile(pFile, *beta0, beta, *p);
    
    // Check Outer Loop Convergence
    if (hasConverged(*beta0, beta, lastOuterBeta, *p, *tol)) {
      Rprintf("Number of outer iterations = %d.\n", *maxiter - iterOuter);
      break;
    }
  }
  fclose(pFile);
  free(Xbeta);
  free(U);
  free(mS);
  free(zeta);
  free(lastOuterBeta);
  free(lastInnerBeta);
  free(partialResidual);
}
