double phiPrime(double x, double t);
void updateWorkingResponse(double *X, double *Y, double beta0, double *beta, double *zeta, double *Xbeta,
                           double t, int n, int p);
void HHSVM(double *X, double *Y, double *beta0, double *beta, int *n, int *p,
	   double *lambda, double *alpha, double *t, int *maxiter, double *tol, 
	   double *beta_final);
