void updateBeta(double *X, double *Y, double *beta,
		double *Xbeta, double *partialResidual,
		int n, int p, double lambda, double alpha, double *avgVar);
void updateBetaInactiveSet(double *X, double *Y, double *beta,
			   double *Xactive, double *betaActive, int sizeActiveSet,
			   double *Xbeta, double *partialResidual,
			   int *inactiveSet, int sizeInactiveSet,
			   int n, int p, double lambda, double alpha, double *avgVar);
void updateActiveAndInactiveSets(int *activeSet, int *sizeActiveSet, int *inactiveSet, int *sizeInactiveSet,
				 double *beta, int p);
void updateXactive(double *X, double *Xactive, int nrowsX, int *activeSet, int sizeActiveSet);
void minimizeElasticNetSquaredLoss(double *X, double *Y, double *beta, int n, int p, 
				   double lambda, double alpha, int maxiter, double tol);
void screenVariablesSTRONG(double *X, double *Y, int n, int p,
                           double lambda, double alpha,
                           int *keepSet, int *sizeKeepSet);
void screenVariablesSAFE(double *X, double *Y, int n, int p, double lambda, double alpha,
			 int *keepSet, int *sizeKeepSet);
void ENLS(double *X, double *Y, double *beta, int *n, int *p, 
	  double *lambda, double *alpha, int *maxiter, double *tol);
