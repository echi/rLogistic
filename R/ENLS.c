#include <R.h> 
#include <Rmath.h> 
#include <R_ext/BLAS.h>
#include <math.h>
#include <stdio.h>

inline double softThreshold(double x, double lambda) {
  return(sign(x) * fmax2( fabs(x) - lambda, 0.));
}

void copyArray(double *source, double *destination, int nElements) {
  int j;
  for (j = nElements; j--; )
    destination[j] = source[j];
}

int hasConverged(double *currentBeta, double *lastBeta, int p, double tol) {
  int j;
  double delta, norm = 0.;
  for (j = p; j--; ) {
    delta = currentBeta[j] - lastBeta[j];
    norm += delta * delta;
  }
  if (sqrt(norm) < tol)
    return 1;
  return 0;
}  

double innerProduct(double *x, double *y, int n) {
  const int one = 1;
  return ddot_(&n, x, &one, y, &one);
}

// copyBigVectorIntoSmallVector assumes sizeIndexSet > 0!!!                                                                     
void copyBigVectorIntoSmallVector(double *bigVector, double *smallVector,
                                  int *indexSet, int sizeIndexSet) {
  int j;
  for (j = 0; j < sizeIndexSet; j++)
    smallVector[j] = bigVector[indexSet[j]];
}

// copySmallVectorIntoBigVector assumes sizeIndexSet > 0!!!                                                                     
void copySmallVectorIntoBigVector(double *smallVector, double *bigVector,
                                  int *indexSet, int sizeIndexSet) {
  int j;
  for (j = 0; j < sizeIndexSet; j++)
    bigVector[indexSet[j]] = smallVector[j];
}


void updateBeta(double *X, double *Y, double *beta,
		double *Xbeta, double *partialResidual,
		int n, int p, double lambda, double alpha, double *avgVar) {
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
      partialResidual[i] = Y[i] - Xbeta[i];
    beta[j] = softThreshold(innerProduct(&X[n*j], partialResidual, n)/((double)n), la) / (avgVar[j] + loma);  
  }
}

void updateBetaInactiveSet(double *X, double *Y, double *beta,
			   double *Xactive, double *betaActive, int sizeActiveSet,
			   double *Xbeta, double *partialResidual,
			   int *inactiveSet, int sizeInactiveSet,
			   int n, int p, double lambda, double alpha, double *avgVar) {
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
    partialResidual[i] = Y[i] - Xbeta[i];
  
  for (j = sizeInactiveSet; j--; ) {
    jIS = inactiveSet[j];
    beta[jIS] = softThreshold(innerProduct(&X[n*jIS], partialResidual, n)/((double)n), la) / (avgVar[jIS] + loma);  
  }
}

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

void minimizeElasticNetSquaredLoss(double *X, double *Y, double *beta, int n, int p, 
		 double lambda, double alpha, int maxiter, double tol) {
  int i, j, iter;
  double *Xbeta = calloc(n, sizeof(double));
  double *avgVar = calloc(p, sizeof(double));
  double *avgVarActive = calloc(p, sizeof(double));
  double *betaLast = calloc(p, sizeof(double));
  double *partialResidual = calloc(n, sizeof(double));
  double *Xactive = calloc( n * (p), sizeof(double));
  double *betaActive = calloc(p, sizeof(double));
  int *activeSet = calloc(p, sizeof(int));
  int *inactiveSet = calloc(p, sizeof(int));
  int sizeActiveSet, sizeInactiveSet;
  double sum;

  for (j = p; j--; )
    avgVar[j] = innerProduct(&X[n*j], &X[n*j], n)/((double)n);
  
  // Initialize active set variables
  updateActiveAndInactiveSets(activeSet, &sizeActiveSet, inactiveSet, &sizeInactiveSet, beta, p);
  
  for (;;) {
    if (sizeActiveSet > 0) {
      updateXactive(X, Xactive, n, activeSet, sizeActiveSet);
      copyBeta2BetaActive(beta, betaActive, activeSet, sizeActiveSet);
      
      for (j = 0; j < sizeActiveSet; j++)
	avgVarActive[j] = avgVar[activeSet[j]];
	
      // Run coordinate descent on active set.
      for (iter = maxiter; iter--; ) {
	  copyArray(betaActive, betaLast, sizeActiveSet);
	  updateBeta(Xactive, Y, betaActive, Xbeta, partialResidual, n, sizeActiveSet, lambda, alpha, avgVar);
	 
	  // Check Inner Loop Convergence
	  if (hasConverged(betaActive, betaLast, sizeActiveSet, tol)) {
	    Rprintf("Number of inner iterations = %d.\n", maxiter - iter);
	    break;
	  }
	}
	copyBetaActive2Beta(beta, betaActive, activeSet, sizeActiveSet);
      }
      updateBetaInactiveSet(X, Y, beta, Xactive, betaActive, sizeActiveSet, Xbeta, partialResidual,
			    inactiveSet, sizeInactiveSet, n, p, lambda, alpha, avgVar);
      sum = 0.;
      for (j = sizeInactiveSet; j--; )
	sum += fabs(beta[inactiveSet[j]]);

      updateActiveAndInactiveSets(activeSet, &sizeActiveSet, inactiveSet, &sizeInactiveSet, beta, p);
     
      if (sum == 0.)
	break;
  }
  free(Xbeta);
  free(avgVar);
  free(avgVarActive);
  free(betaLast);
  free(partialResidual);
  free(Xactive);
  free(betaActive);
  free(activeSet);
  free(inactiveSet);
}

void screenVariablesSTRONG(double *X, double *Y, int n, int p,
                           double lambda, double alpha,
                           int *keepSet, int *sizeKeepSet) {
  int j;
  const char transX = 'T';
  const double alp = 1.;
  const double bet = 0.;
  const int one = 1;
  double *corr = calloc(p, sizeof(double));
  double *scores = calloc(p, sizeof(double));
  double threshold, lambdaMax = -1.;
  Rprintf("screenVariablesSTRONG\n");
  dgemv_(&transX, &n, &p, &alp, X, &n, Y, &one, &bet, corr, &one);
  for (j = 0; j < p; j++) {
    scores[j] = fabs(corr[j]);
    lambdaMax = fmax2(lambdaMax, scores[j]);
  }
  *sizeKeepSet = 0;
  for (j = 0; j < p; j++) {
    threshold = alpha*(2.*((double)n)*lambda - lambdaMax);
    if (scores[j] >= threshold) {
      keepSet[*sizeKeepSet] = j;
      (*sizeKeepSet)++;
    }
  }
  free(corr);
  free(scores);
}

void screenVariablesSAFE(double *X, double *Y, int n, int p, double lambda, double alpha,
			 int *keepSet, int *sizeKeepSet) {
  int i, j;
  const char transX = 'T';
  const double alp = 1.;
  const double bet = 0.;
  const int one = 1;
  double *corr = calloc(p, sizeof(double));
  double *scores = calloc(p, sizeof(double));
  double varXj, normY, threshold, lambdaMax = -1.;
  Rprintf("screenVariablesSAFE\n");
  dgemv_(&transX, &n, &p, &alp, X, &n, Y, &one, &bet, corr, &one);
  for (j = 0; j < p; j++) {
    scores[j] = fabs(corr[j]);
    lambdaMax = fmax2(lambdaMax, scores[j]);
  }
  normY = sqrt(innerProduct(Y, Y, n));
  *sizeKeepSet = 0;
  for (j = 0; j < p; j++) {
    varXj = innerProduct(&X[n*j], &X[n*j], n);
    /*
    Rprintf("((double)n)*alpha*lambda = %g\n", ((double)n*alpha*lambda));
    Rprintf("normY*sqrt(varXj + ((double)n)*(1.-alpha)*lambda) = %g\n", normY*sqrt(varXj + ((double)n)*(1.-alpha)*lambda));
    Rprintf("Second Term = %g\n", normY*sqrt(varXj + ((double)n)*(1.-alpha)*lambda)*(lambdaMax-((double)n)*lambda)/lambdaMax);
    threshold = ((double)n)*alpha*lambda - (normY*sqrt(varXj + ((double)n)*(1.-alpha)*lambda)*(lambdaMax-((double)n)*lambda)/lambdaMax);
    Rprintf("threshold[%d] = %g\n", j, threshold);
    */
    if (scores[j] >= threshold) {
      keepSet[*sizeKeepSet] = j;
      (*sizeKeepSet)++;
    }
  }
  free(corr);
  free(scores);
}

void ENLS(double *X, double *Y, double *beta, int *n, int *p, 
		 double *lambda, double *alpha, int *maxiter, double *tol) {
  int j;
  double *betaScreened = calloc(*p, sizeof(double));
  int *keepSet = calloc(*p, sizeof(int));
  double *Xscreened;  
  int sizeKeepSet;

  minimizeElasticNetSquaredLoss(X, Y, beta, *n, *p, *lambda, *alpha, *maxiter, *tol);

  /* No screening

  // 1. Center the data

  // 2. Discard variables
  screenVariablesSAFE(X, Y, *n, *p, *lambda, *alpha, keepSet, &sizeKeepSet);
  //  screenVariablesSTRONG(X, Y, *n, *p, *lambda, *alpha, keepSet, &sizeKeepSet);

  Rprintf("Number variables retained: %d\n", sizeKeepSet);

  // CASE I: At least one covariate is kept after screening.
  if (sizeKeepSet > 0) {

    Rprintf("Keeping Variables: ");
    for (j = 0; j < sizeKeepSet; j++)
      Rprintf("%d ", keepSet[j]);
    Rprintf("\n");

    Xscreened = calloc((*n)*sizeKeepSet, sizeof(double));
    updateXactive(X, Xscreened, *n, keepSet, sizeKeepSet);
    copyBigVectorIntoSmallVector(beta, betaScreened, keepSet, sizeKeepSet);
    minimizeElasticNetSquaredLoss(Xscreened, Y, betaScreened, *n, sizeKeepSet, 
				  *lambda, *alpha, *maxiter, *tol);
    for (j = *p; j--; )
      beta[j] = 0.;
    copySmallVectorIntoBigVector(betaScreened, beta, keepSet, sizeKeepSet);
    free(Xscreened);
  }
  // CASE II: No covariates are kept after screening.
  else {
    Rprintf("Keeping no variables.");
    for (j = *p; j--; )
      beta[j] = 0.;
  }

  */

  free(betaScreened);
  free(keepSet);
}

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

void writeParametersToFile(FILE *pFile, double beta0, double *beta, int p) {
  int j;
  fprintf(pFile,"%g, ", beta0);
  for (j = 0; j < p - 1; j++)
    fprintf(pFile,"%g, ", beta[j]);
  fprintf(pFile,"%g\n", beta[p-1]);
}


double mean(double *vector, int nElements) {
  int i;
  double sum = 0.;
  for (i = nElements; i--; )
    sum += vector[i];
  return sum/((double)nElements);
}

void recordCurrentBeta(double beta0, double *beta, int p, double *betaOld) {
  int j;
  for (j = p; j--; )
    betaOld[j+1] = beta[j];
  betaOld[0] = beta0;
}

int hasConverged2(double beta0, double *beta, double *betaOld, int p, double tol) {
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

void logisticL2E(double *X, int *Y, double *beta0, double *beta, int *n, int *p,
		 double * K, double *w, double *lambda, double *alpha,
		 int *maxiter, double *tol,
		 int *save_flag, char const *filename) {
  int i, j, iter;
  double newLambda;
  double *Xbeta = calloc(*n, sizeof(double));
  double *U = calloc(*n, sizeof(double));
  double *zeta = calloc(*n, sizeof(double));
  double *betaLast = calloc(*p+1, sizeof(double));
  FILE *pFile;

  pFile = fopen("betaHx.csv","w");
  writeParametersToFile(pFile, *beta0, beta, *p);

  for (i = *n; i--; )
    U[i] = (double) (Y[i] + Y[i] - 1);

  newLambda = *lambda / ( (*w) * (*K) );

  for (iter = *maxiter; iter--; ) {
    recordCurrentBeta(*beta0, beta, *p, betaLast);
    updateWorkingResponse(X, U, *w, *K, *beta0, beta, zeta, Xbeta, *n, *p);
    *beta0 = mean(zeta, *n);
    for (i = *n; i--; )
      zeta[i] = zeta[i] - *beta0;
    ENLS(X, zeta, beta, n, p, &newLambda, alpha, maxiter, tol);

    writeParametersToFile(pFile, *beta0, beta, *p);

    // Check Outer Loop Convergence                                                                                             
    if (hasConverged2(*beta0, beta, betaLast, *p, *tol)) {
      Rprintf("Number of outer iterations = %d.\n", *maxiter - iter);
      break;
    }
  }
  fclose(pFile);
  free(Xbeta);
  free(U);
  free(zeta);
  free(betaLast);
}
