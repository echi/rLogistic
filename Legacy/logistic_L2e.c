#include <R.h> 
#include <Rmath.h> 
#include <R_ext/BLAS.h>
#include <math.h>
#include <stdio.h>

typedef struct {

  int Nrows;
  int Ncols;
  
  double *data;

} matrix;

inline double softThreshold(double x, double lambda);
inline double cplogis(double x);
void meanL2e(double *loss, double *X, int *Y, double *beta0, double *beta, 
	     int *n, int *p, double *w, double *lambda, double *alpha);
void meanL2eBatch(double *loss, double *X, int *Y, double *beta0, double *beta, 
		  int *n, int *p, double *w, double *lambda, double *alpha, int *m);
void logistic_l2e(double *loss, double *X, int *Y, double *beta0, double *beta, int *n, int *p, 
		  double * K, double *w, double *lambda, double *alpha, int *maxiter, int *check, double *tol,
		  int *save_flag, char const *filename);


matrix *matrix_build(int Nrows, int Ncols){
  matrix *A = (matrix*) calloc(1, sizeof(matrix));
  A->Nrows = Nrows;
  A->Ncols = Ncols;
  A->data = (double*) calloc(Nrows*Ncols, sizeof(double));
  return A;
}

double matrix_get_entry(matrix *A, int row, int col){
  /* column major (compatible with fortran) */
  return A->data[ (row-1) + (col-1)*A->Nrows];
}

inline double softThreshold(double x, double lambda) {
  return(sign(x) * fmax2( fabs(x) - lambda, 0.));
}

//   This function calculates 0.5 times the L2e loss plus an elastic net penalty.
//
//     loss is a scalar double.
//     X is a n by p double.
//     Y is a n-length vector of ints.
//     beta0 is a scalar double.
//     beta is a p-length vector of doubles.
//     w, lambda, and alpha are scalar doubles.

inline double cplogis(double x) {
  return 0.5 * (1. + tanh(0.5 * x));
}

void meanL2e(double *loss, double *X, int *Y, double *beta0, double *beta, 
	     int *n, int *p, double *w, double *lambda, double *alpha) {
  const char transA = 'N';
  const double alp = 1.;
  const double bet = 0.;
  const int one = 1;
  double *u = calloc(*n, sizeof(double));
  double cumsum = 0.;
  double pr, qr, el1 = 0., el2 = 0., penalty;
  int i;

  dgemv_(&transA, n, p, &alp, X, n, beta, &one, &bet, u, &one);

  *loss = 0.;
  for (i = *n; i--; ) {
    //    pr = plogis(*beta0 + u[i], 0.0, 1.0, 1, 0);
    pr = cplogis(*beta0 + u[i]);
    qr = 1. - pr;
    //    cumsum += pr*(*w * pr - ((double) Y[i] + Y[i])) + qr*(*w *qr - ((double) 2 - Y[i] - Y[i]));
    *loss += pr*(*w * pr - ((double) Y[i] + Y[i])) + qr*(*w *qr - ((double) 2 - Y[i] - Y[i]));
  }
  //  cumsum *= 0.5* (*w) / (double) *n;
  *loss *= 0.5* (*w) / (double) *n;

  for (i = *p; i--; ) {
    el1 += fabs(beta[i]);
    el2 += beta[i] * beta[i];
  }

  penalty = *lambda * (*alpha * el1 + 0.5 * (1.-*alpha) * el2);

  *loss += penalty;
  //  *loss = cumsum + penalty;

  free(u);
}

// This code expects the beta matrix to be passed in column major - i.e. each column is a 
// p dimensional vector. So, beta0 is a column vector of length m, and beta is matrix of
// size p by m.
void meanL2eBatch(double *loss, double *X, int *Y, double *beta0, double *beta, 
		  int *n, int *p, double *w, double *lambda, double *alpha, int *m) {
  int i;
  Rprintf("p = %d\n", *p);
  for (i = *m; i--; )
    meanL2e(&(loss[i]), X, Y, &beta0[i], &beta[*p * i], n, p, w, lambda, alpha);
}

/*
void meanL2e2D(double *loss, double *X, double *Y, double *beta0, double *beta1, int *n, int *p,
	       double *w, double *lambda, double *alpha, int *m0, int *m1) {

  int i, j;
  
  for (i = 0; i < *m0; ++i)
    for (j = 0; j < *m1; ++j)
      meanL2e(&(loss[(*m1)*i + j], X, Y, &beta
*/

/* Calculates the logistic L2E loss. X is n by p where p is the number of covariates. Note
   Note that the first column of X is _not_ 1.
 */
void logistic_l2e(double *loss, double *X, int *Y, double *beta0, double *beta, int *n, int *p, 
		  double * K, double *w, double *lambda, double *alpha, int *maxiter, int *check, double *tol,
		  int *save_flag, char const *filename) {

  int i, j, iter = 0;
  const char transX = 'N';
  const double alp = 1.;
  const double bet = 0.;
  const int one = 1;
  double *u = calloc(*n, sizeof(double));
  double *U = calloc(*n, sizeof(double));
  double *mS = calloc(*p, sizeof(double));
  double *WZ = calloc(*n, sizeof(double));
  double P, D;
  double m, curL2eCost, lastCost, mWZ, la = (*lambda)*(*alpha) , loma = (*lambda)*(1.-*alpha);
  double last;
  FILE *pFile;

  //  printf("%c%c\n", filename[0],filename[1]);
  m = (*w) * (*K) / ((double) *n);
  pFile = fopen("mybeta.csv","w");
  //  pFile = fopen(filename,"w");
  for (j = *p; j--; ) {
    fprintf(pFile,"%g, ", beta[j]);
    mS[j] = m * ddot_(n, &X[(*n) * j], &one, &X[(*n) * j], &one);
  }
  fprintf(pFile,"%g, ", *beta0);
  fprintf(pFile,"\n");
  
  for (i = *n; i--; )
    U[i] = (double) (Y[i] + Y[i] - 1);

  if (*check) {
    // Calculate inital loss.
    meanL2e(&lastCost, X, Y, beta0, beta, n, p, w, lambda, alpha);

    // Calculate initial working response.
    dgemv_(&transX, n, p, &alp, X, n, beta, &one, &bet, u, &one);
    for (i = *n; i--; ) {
      //      P = plogis(*beta0 + u[i], 0.0, 1.0, 1, 0);
      P = cplogis(*beta0 + u[i]);
      D = P + P - 1.;
      WZ[i] = P * (1. - P) * (*w * D - U[i]);
    }
    for (iter = *maxiter; iter--; ) {
      // Update intercept
      mWZ = 0.;
      for (i = *n; i--; )
        mWZ += WZ[i];
      *beta0 -= (mWZ / (*K * ((double)*n)));
      // Update working response.
      for (i = *n; i--; ) {
	//        P = plogis(*beta0 + u[i], 0.0, 1.0, 1, 0);
        P = cplogis(*beta0 + u[i]);
        D = P + P - 1.;
        WZ[i] = P * (1. - P) * (*w * D - U[i]);
      }
      // Update beta[j].                                                                                                    
      for (j = *p; j--; ) {
	last = beta[j];
	//	        beta[j] = softThreshold(mS[j] * beta[j]
	//	                      - (*w/((double)*n)) * ddot_(n, &X[(*n) * j], &one, WZ, &one),
	//                                        *lambda* (*alpha)) / (mS[j] + *lambda*(1.-*alpha));
	beta[j] = softThreshold( mS[j]*beta[j] - (*w/((double)*n)) * ddot_(n, &X[(*n)*j], &one, WZ, &one), la) / (mS[j] + loma);
        fprintf(pFile,"%g, ", beta[j]);
        // Update working response.
	for (i = *n; i--; ) {
          u[i] += (beta[j] - last) * X[(*n) * j + i];
	  //          P = plogis(*beta0 + u[i], 0.0, 1.0, 1, 0);
          P = cplogis(*beta0 + u[i]);
          D = P + P - 1.;
          WZ[i] = P * (1. - P) * (*w * D - U[i]);
        }
      }
      fprintf(pFile,"%g, ", *beta0);
      fprintf(pFile,"\n");

      // Check change in cost after an entire cycle through the p + 1 variables.
      meanL2e(&curL2eCost, X, Y, beta0, beta, n, p, w, lambda, alpha);
      if ( fabs(curL2eCost - lastCost) / (fabs(lastCost) + 1.) > *tol) {
	lastCost = curL2eCost;
	//	loss[*maxiter - iter] = curL2eCost;
      }
      else {
	Rprintf("Number of iterations = %d.\n", *maxiter - iter);
	break;
      }
    }
  } else {
    // Calculate initial working response.
    dgemv_(&transX, n, p, &alp, X, n, beta, &one, &bet, u, &one);
    for (i = *n; i--; ) {
      //      P = plogis(*beta0 + u[i], 0.0, 1.0, 1, 0);
      P = cplogis(*beta0 + u[i]);
      D = P + P - 1.;
      WZ[i] = P * (1. - P) * (*w * D - U[i]);
    }
    for (iter = *maxiter; iter--; ) {
      // Update intercept                                                                                                             
      mWZ = 0.;
      for (i = *n; i--; ) 
        mWZ += WZ[i];
      *beta0 -= (mWZ / (*K * ((double)*n)));
      // Update working response.
      for (i = *n; i--; ) {
	//        P = plogis(*beta0 + u[i], 0.0, 1.0, 1, 0);
	P = cplogis(*beta0 + u[i]);
	D = P + P - 1.;
        WZ[i] = P * (1. - P) * (*w * D - U[i]);
      }
      // Update beta[j].
      for (j = *p; j--; ) {
	last = beta[j];
	//	        beta[j] = softThreshold(mS[j] * beta[j] 
	//					- (*w/((double)*n)) * ddot_(n, &X[(*n) * j], &one, WZ, &one), 
	//					*lambda* (*alpha)) / (mS[j] + *lambda*(1.-*alpha));
	beta[j] = softThreshold( mS[j]*beta[j] - (*w/((double)*n)) * ddot_(n, &X[(*n)*j], &one, WZ, &one), la) / (mS[j] + loma);
        fprintf(pFile,"%g, ", beta[j]);
	// Update working response.
	for (i = *n; i--; ) {
	  u[i] += (beta[j] - last) * X[(*n) * j + i];
	  //	  P = plogis(*beta0 + u[i], 0.0, 1.0, 1, 0);
	  P = cplogis(*beta0 + u[i]);
	  D = P + P - 1.;
	  WZ[i] = P * (1. - P) * (*w * D - U[i]);
	}	
      }
      fprintf(pFile,"%g, ", *beta0);
      fprintf(pFile,"\n");
    } 
  } 
  fclose(pFile);

  free(mS);
  free(WZ);
  free(U);
  free(u);

}
