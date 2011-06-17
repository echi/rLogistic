#include "util.h"
#include "ENLS.h"

inline double softThreshold(double x, double lambda) {
  return(sign(x) * fmax2(fabs(x) - lambda, 0.));
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

int minimizeElasticNetSquaredLoss(double *X, double *Y, double *beta, int n, int p, 
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
  int sizeActiveSet, sizeInactiveSet, inner_iter_sum=0;
  double sum;

  for (j = p; j--; )
    avgVar[j] = innerProduct(&X[n*j], &X[n*j], n)/((double)n);
  
  // Initialize active set variables
  updateActiveAndInactiveSets(activeSet, &sizeActiveSet, inactiveSet, &sizeInactiveSet, beta, p);
  
  for (;;) {
    if (sizeActiveSet > 0) {
      updateXactive(X, Xactive, n, activeSet, sizeActiveSet);
      copyBigVectorIntoSmallVector(beta, betaActive, activeSet, sizeActiveSet);
      
      for (j = 0; j < sizeActiveSet; j++)
	avgVarActive[j] = avgVar[activeSet[j]];
      
      // Run coordinate descent on active set.
      for (iter = maxiter; iter--; ) {
	//	copyArray(betaActive, betaLast, sizeActiveSet);
	updateBeta(Xactive, Y, betaActive, Xbeta, partialResidual, n, sizeActiveSet, lambda, alpha, avgVarActive);

	// Check Inner Loop Convergence
	/*	if (hasConverged(betaActive, betaLast, sizeActiveSet, tol)) {
	  //	  Rprintf("Number of inner iterations = %d.\n", maxiter - iter);
	  break;
	  }*/
      }
      copySmallVectorIntoBigVector(betaActive, beta, activeSet, sizeActiveSet);
      inner_iter_sum += fmin2(maxiter - iter, maxiter);
    }
    updateBetaInactiveSet(X, Y, beta, Xactive, betaActive, sizeActiveSet, Xbeta, partialResidual,
			  inactiveSet, sizeInactiveSet, n, p, lambda, alpha, avgVar);
    sum = 0.;
    for (j = sizeInactiveSet; j--; )
      sum += fabs(beta[inactiveSet[j]]);
    
    // Variables may drop out or come in.
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
  
  return(inner_iter_sum);

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
    //    Rprintf("((double)n)*alpha*lambda = %g\n", ((double)n*alpha*lambda));
    //    Rprintf("normY*sqrt(varXj + ((double)n)*(1.-alpha)*lambda) = %g\n", normY*sqrt(varXj + ((double)n)*(1.-alpha)*lambda));
    //    Rprintf("Second Term = %g\n", normY*sqrt(varXj + ((double)n)*(1.-alpha)*lambda)*(lambdaMax-((double)n)*lambda*alpha)/lambdaMax);
    threshold = ((double)n)*alpha*lambda - (normY*sqrt(varXj + ((double)n)*(1.-alpha)*lambda)*(lambdaMax-((double)n)*lambda*alpha)/lambdaMax);
    //    Rprintf("threshold[%d] = %g\n", j, threshold);
    if (scores[j] >= threshold) {
      keepSet[*sizeKeepSet] = j;
      (*sizeKeepSet)++;
    }
  }
  free(corr);
  free(scores);
}

int ENLS(double *X, double *Y, double *beta, int *n, int *p, 
		 double *lambda, double *alpha, int *maxiter, double *tol) {
  int j;
  double *betaScreened = calloc(*p, sizeof(double));
  int *keepSet = calloc(*p, sizeof(int));
  double *Xscreened;  
  int sizeKeepSet;

  return(minimizeElasticNetSquaredLoss(X, Y, beta, *n, *p, *lambda, *alpha, *maxiter, *tol));

  /*
  // No screening

  // Discard variables
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
