#include <R.h>
#include <Rmath.h>
#include <R_ext/BLAS.h>
#include <math.h>
#include <stdio.h>

double innerProduct(double *x, double *y, int n) {
  const int one = 1;
  return ddot_(&n, x, &one, y, &one);
}

void screenVariablesSAFE(double *X, double *Y, int n, int p,
                         double lambda, double alpha,
                         int *keepSet, int *sizeKeepSet) {
  int j;
  const char transX = 'T';
  const double alp = 1.;
  const double bet = 0.;
  const int one = 1;
  double *corr = calloc(p, sizeof(double));
  double *scores = calloc(p, sizeof(double));
  double varXj, normY, threshold, lambdaMax = -1.;

  dgemv_(&transX, &n, &p, &alp, X, &n, Y, &one, &bet, corr, &one);
  for (j = p; j--; ) {
    scores[j] = fabs(corr[j]);
    lambdaMax = fmax2(lambdaMax, scores[j]);
    /*    if (scores[j] > lambdaMax)
	  lambdaMax = scores[j]; */
  }

  normY = sqrt(innerProduct(Y, Y, n));
  *sizeKeepSet = 0;
  for (j = p; j--; ) {
    varXj = innerProduct(&X[n*j], &X[n*j], n);
    threshold = alpha*lambda - normY*sqrt(varXj + (1.-alpha)*lambda)*(lambdaMax-lambda)/lambdaMax;
    if (scores[j] >= threshold) {
      keepSet[*sizeKeepSet] = j;
      (*sizeKeepSet)++;
    }
  }
  free(corr);
  free(scores);
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

  dgemv_(&transX, &n, &p, &alp, X, &n, Y, &one, &bet, corr, &one);
  for (j = p; j--; ) {
    scores[j] = fabs(corr[j]);
    lambdaMax = fmax2(lambdaMax, scores[j]);
    /*    if (scores[j] > lambdaMax)
	  lambdaMax = scores[j]; */
  }

  *sizeKeepSet = 0;
  for (j = p; j--; ) {
    threshold = alpha*(2.*lambda - lambdaMax);
    if (scores[j] >= threshold) {
      keepSet[*sizeKeepSet] = j;
      (*sizeKeepSet)++;
    }
  }
  free(corr);
  free(scores);
}
  
void selectVariablesSAFE(double *X, double *Y, int *keepSet, int *sizeKeepSet, int *n, int *p,
                         double *lambda, double *alpha) {
  screenVariablesSAFE(X, Y, *n, *p, *lambda, *alpha, keepSet, sizeKeepSet);
}

void selectVariablesSTRONG(double *X, double *Y, int *keepSet, int *sizeKeepSet, int *n, int *p,
                         double *lambda, double *alpha) {
  screenVariablesSTRONG(X, Y, *n, *p, *lambda, *alpha, keepSet, sizeKeepSet);
}
