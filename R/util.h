#include <R.h> 
#include <Rmath.h> 
#include <R_ext/BLAS.h>
#include <math.h>
#include <stdio.h>

void copyArray(double *source, double *destination, int nElements);
int hasConverged(double *currentBeta, double *lastBeta, int p, double tol);
double innerProduct(double *x, double *y, int n);
void copyBigVectorIntoSmallVector(double *bigVector, double *smallVector, int *indexSet, int sizeIndexSet);
void copySmallVectorIntoBigVector(double *smallVector, double *bigVector, int *indexSet, int sizeIndexSet);
void recordCurrentBeta(double beta0, double *beta, int p, double *betaOld);
double mean(double *vector, int nElements);
int hasConverged2(double beta0, double *beta, double *betaOld, int p, double tol);
double twoNorm(double beta0, double *beta, double *betaOld, int p);
void centerColumns(double *X, double *C, int n, int p);
