#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include "invlib.h"

/* test omni2uni */

int main (int argc, char *argv[]) {
  int result;
  long int int_params[10];
  double omniflux,dlogomniflux,uniflux,dloguniflux,real_params[10];

  int_params[0] = 50; /* NA - number of angular gridpoints */
  int_params[1] = -1; /* TEM-1 method */
  int_params[2] = 1; /* 1 = verbose to standard out */
  int_params[3] = 3; /* minimizer, 0=BFGS, 3=NM */
  int_params[4] = 1000; /* maximumn # of iterations */

  real_params[0] = 300.0; /* 300 keV */
  real_params[1] = 40000.0/100.0; /* B/B0 */
  real_params[2] = 6.6; /* Lm */

  omniflux = 1.0E4; /* typical value 1E+4 */
  dlogomniflux = log(2)/2;

  printf("omni2uni_test:inputs = omni=%.5e (dlog=%.5e) @ [%.1f keV, B/B0=%.2f, Lm=%.2f]\n",
	 omniflux,dlogomniflux, real_params[0],real_params[1],real_params[2]);
  result = omni2uni(&omniflux,&dlogomniflux,int_params,real_params,NULL,&uniflux,&dloguniflux);
  printf("omni2uni_test:inputs = omni=%.5e (dlog=%.5e) @ [%.1f keV, B/B0=%.2f, Lm=%.2f]\n",
	 omniflux,dlogomniflux, real_params[0],real_params[1],real_params[2]);
  printf("omni2uni_test:result = %i, uni=%.5e (dlog=%.5e)\n",result,uniflux,dloguniflux);
  return(0);
}
