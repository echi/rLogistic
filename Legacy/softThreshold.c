#include <R.h>
#include <Rmath.h>
#include <math.h>

void softThreshold(double *shrunk, double *x, double *lambda, int *n) {

  int i;

  for (i = *n; i--; )
    shrunk[i] = sign(x[i]) * fmax2( fabs(x[i]) - *lambda, 0.);

}
