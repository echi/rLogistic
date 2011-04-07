void updateWorkingResponse(double *X, double *U, double w, double K,
                           double beta0, double *beta, double *zeta, double *Xbeta,
                           int n, int p);
void writeParametersToFile(FILE *pFile, double beta0, double *beta, int p);
void recordCurrentBeta(double beta0, double *beta, int p, double *betaOld);
void logisticL2E(double *X, int *Y, double *beta0, double *beta, int *n, int *p,
		 double * K, double *w, double *lambda, double *alpha,
		 int *maxiter, double *tol,
		 int *save_flag, char const *filename);
