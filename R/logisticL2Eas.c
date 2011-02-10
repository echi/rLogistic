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

void updateBetaInactiveSet(double *X, double beta0, double *beta,
			   double *Xactive, double *betaActive, int sizeActiveSet,
			   double *zeta, double *Xbeta, double *partialResidual,
			   int *inactiveSet, int sizeInactiveSet,
			   int n, int p, double lambda, double alpha, double m, double *mS) {
  int i, j;
  const char transX = 'N';
  const double alp = 1.;
  const double bet = 0.;
  const int one = 1;
  int jIS;
  double la = lambda * alpha;
  double loma = lambda * (1. - alpha);

  dgemv_(&transX, &n, &sizeActiveSet, &alp, Xactive, &n, betaActive, &one, &bet, Xbeta, &one);
  for (i = n; i--; )
    partialResidual[i] = zeta[i] - beta0 - Xbeta[i];
  
  for (j = sizeInactiveSet; j--; ) {
    jIS = inactiveSet[j];
    beta[jIS] = softThreshold(m*innerProduct(&X[n*jIS], partialResidual, n), la) / (mS[jIS] + loma);  
  }
}

/*
int computeAdd2ActiveSet(double *X, double beta0, double *beta,
			 double *zeta, double *Xbeta, double *partialResidual,
			 int n, int sizeInactiveSet, int *inactiveSet,
			 double lambda, double alpha, double m, double *mS) {

  int i, j, iVar, sizeAdd2ActiveSet = 0;
  const char transX = 'N';
  const double alp = 1.;
  const double bet = 0.;
  const int one = 1;
  double la = lambda * alpha;
  double loma = lambda * (1. - alpha);

  dgemv_(&transX, &n, &sizeInactiveSet, &alp, X, &n, beta, &one, &bet, Xbeta, &one);
  for (i = n; i--; )
    partialResidual[i] = zeta[i] - beta0 - Xbeta[i];
  for (j = sizeInactiveSet; j--; ) {
    iVar = inactiveSet[j];
    beta[iVar] = softThreshold(m*innerProduct(&X[n*iVar], partialResidual, n), la) / (mS[iVar] + loma);
    if (fabs(beta[iVar]) > 0.)
      sizeAdd2ActiveSet++;
  }
  return sizeAdd2ActiveSet;
}
*/

void updateActiveAndInactiveSets(int *activeSet, int *sizeActiveSet, int *inactiveSet, int *sizeInactiveSet,
				 double *beta, int p) {
  int j;
  *sizeActiveSet = 0;
  *sizeInactiveSet = 0;
  for (j = 0; j < p; j++) {
    if (fabs(beta[j]) > 0.) {
      activeSet[*sizeActiveSet] = j;
      (*sizeActiveSet)++;
    } else {
      inactiveSet[*sizeInactiveSet] = j;
      (*sizeInactiveSet)++;
    }
  }
}

// updateXactive assumes sizeActiveSet > 0!!!
// updateXactive assumes activeSet indexes from 0.
void updateXactive(double *X, double *Xactive, int nrowsX, int *activeSet, int sizeActiveSet) {
  int j, k;
  for (j = sizeActiveSet; j--; )
    for (k = nrowsX; k--; )
      Xactive[nrowsX*j + k] = X[nrowsX*activeSet[j] + k];
}

// copyBeta2BetaActive assumes sizeActiveSet > 0!!!
void copyBeta2BetaActive(double *beta, double *betaActive, int *activeSet, int sizeActiveSet) {
  int j;
  for (j = sizeActiveSet; j--; )
    betaActive[j] = beta[activeSet[j]];
}

void copyBetaActive2Beta(double *beta, double *betaActive, int *activeSet, int sizeActiveSet) {
  int j;
  for (j = sizeActiveSet; j--; )
    beta[activeSet[j]] = betaActive[j];  
}

/*
void printMatrix(double *X, int nrow, int ncol) {
  int i, j;
  for (j = 0; j < nrow; j++) {
    for (i = 0; i < ncol; i++)
      Rprintf("%g ", X[nrow*i + j]);
    Rprintf("\n");
  }
}

void printMatrixInt(int *X, int nrow, int ncol) {
  int i, j;
  for (j = 0; j < nrow; j++) {
    for (i = 0; i < ncol; i++)
      Rprintf("%d ", X[nrow*i + j]);
    Rprintf("\n");
  }
}
*/

double mean(double *vector, int nElements) {
  int i;
  double sum = 0.;
  for (i = nElements; i--; )
    sum += vector[i];
  return sum/((double)nElements);
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
  double *mSactive = calloc(*p, sizeof(double));
  double *zeta = calloc(*n, sizeof(double));
  double *lastOuterBeta = calloc(*p+1, sizeof(double));
  double *lastInnerBeta = calloc(*p+1, sizeof(double));
  double *partialResidual = calloc(*n, sizeof(double));
  double *Xactive = calloc( (*n) * (*p), sizeof(double));
  double *betaActive = calloc(*p, sizeof(double));
  int *activeSet = calloc(*p, sizeof(int));
  int *inactiveSet = calloc(*p, sizeof(int));
  int sizeActiveSet, sizeInactiveSet, sizeAdd2ActiveSet;
  double m, sum;
  FILE *pFile;

  pFile = fopen("betaHx.csv","w");
  writeParametersToFile(pFile, *beta0, beta, *p);

  m = (*w)*(*K)/((double)*n);
  for (j = *p; j--; )
    mS[j] = m*innerProduct(&X[(*n)*j], &X[(*n)*j], *n);

  for (i = *n; i--; )
    U[i] = (double) (Y[i] + Y[i] - 1);
  
  // Initialize active set variables
  updateActiveAndInactiveSets(activeSet, &sizeActiveSet, inactiveSet, &sizeInactiveSet, beta, *p);

  for (iterOuter = *maxiter; iterOuter--; ) {
    recordCurrentBeta(*beta0, beta, *p, lastOuterBeta);

    updateWorkingResponse(X, U, *w, *K, *beta0, beta, zeta, Xbeta, *n, *p);
    /*
    *beta0 = mean(zeta, *n);
    for (i = *n; i--; )
      zeta[i] = zeta[i] - *beta0;
    */

    for (;;) {
      if (sizeActiveSet > 0) {
	updateXactive(X, Xactive, *n, activeSet, sizeActiveSet);
	copyBeta2BetaActive(beta, betaActive, activeSet, sizeActiveSet);

	for (j = 0; j < sizeActiveSet; j++)
	  mSactive[j] = mS[activeSet[j]];
	
	// Run coordinate descent on active set.
	for (iterInner = *maxiter; iterInner--; ) {
	  recordCurrentBeta(*beta0, betaActive, sizeActiveSet, lastInnerBeta);
	  
	  updateIntercept(Xactive, beta0, betaActive, zeta, Xbeta, *n, sizeActiveSet);
	  updateBeta(Xactive, *beta0, betaActive, zeta, Xbeta, partialResidual, *n, sizeActiveSet, *lambda, *alpha,
		     m, mSactive);
	 
	  // Check Inner Loop Convergence
	  if (hasConverged(*beta0, betaActive, lastInnerBeta, sizeActiveSet, *tol)) {
	    Rprintf("Number of inner iterations = %d.\n", *maxiter - iterInner);
	    break;
	  }
	}
	copyBetaActive2Beta(beta, betaActive, activeSet, sizeActiveSet);
      }
      // Check to see the active set needs to be enlarged.
      /*      sizeAdd2ActiveSet = computeAdd2ActiveSet(X, *beta0, beta, zeta, Xbeta, partialResidual, *n, 
					       sizeInactiveSet, inactiveSet,
					       *lambda, *alpha, m, mS);

      Rprintf("sizeAdd2ActiveSet = %d.\n", sizeAdd2ActiveSet);
      
      for (i = 0; i < sizeActiveSet; i++)
	Rprintf("%d ", activeSet[i]);
      Rprintf("\n");
      */
      //      updateBeta(X, *beta0, beta, zeta, Xbeta, partialResidual, *n, *p, *lambda, *alpha, m, mS);
      updateBetaInactiveSet(X, *beta0, beta, Xactive, betaActive, sizeActiveSet, zeta, Xbeta, partialResidual,
			    inactiveSet, sizeInactiveSet, *n, *p, *lambda, *alpha, m, mS);

      sum = 0.;
      for (j = sizeInactiveSet; j--; )
	sum += fabs(beta[inactiveSet[j]]);
      //      Rprintf("Sum = %g\n", sum);

      updateActiveAndInactiveSets(activeSet, &sizeActiveSet, inactiveSet, &sizeInactiveSet, beta, *p);
     
      if (sum == 0.)
	break;
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
  free(mSactive);
  free(zeta);
  free(Xactive);
  free(betaActive);
  free(activeSet);
  free(inactiveSet);
  free(lastOuterBeta);
  free(lastInnerBeta);
  free(partialResidual);
}
