#include "util.h" 

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

double mean(double *vector, int nElements) {
  int i;
  double sum = 0.;
  for (i = nElements; i--; )
    sum += vector[i];
  return sum/((double)nElements);
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

void centerColumns(double *X, double *C, int n, int p) {
  int i, j;
  double mu = 0.;
  
  for (j = p; j--; ) {
    for (i = n; i--; )
      mu += X[i + n*j];
    for (i = n; i--; )
      C[i + n*j] = X[i + n*j] - mu/((double)n);
    mu = 0.;
  }
}
