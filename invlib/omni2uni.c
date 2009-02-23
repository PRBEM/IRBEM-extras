/********************************************************************************************************
OMNI2UNI Section - this section defines omni2uni, wide2uni and the routines they need
********************************************************************************************************/

#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h> /* gsl_finite */
#include "invlib.h" /* defines omni2uni and wide2uni*/
#include "invlib_const.h" /* defines constants, e.g., return codes */
#include "omni2uni.h" /* prototypes for omni2uni and wide2uni */
#include "o2u_TEM1.h" /* routines for TEM1 method */

/**************************
 Other omni2uni and wide2uni routines
***************************/

int wide2uni(const double *wideflux, const double *dlogwideflux,
             const double *PAgrid, const double *H, 
             const long int *int_params, 
             const double *real_params,
             char *outFile,
             const long int *ialpha0,
             double *uniflux, double *dloguniflux) {
  /*
    Uses TEM-1 method to compute wede2uni.
    Arguments are explained in invlib.pdf
   */
  long int NA, verbose, w2u_method;
  long int i;
  int ret; /* return code from method */
  /* input checking */
  if ((!wideflux) || (!dlogwideflux) || (!PAgrid) || (!H) || (!int_params)
      || (!real_params) || (!ialpha0) || (!uniflux) || (!dloguniflux)) {
    return(INVLIB_ERR_NULL);
  }

  if ((*wideflux<=0) || (*dlogwideflux<=0)) {
    return(INVLIB_ERR_DATANAN);
  }

  NA = int_params[0];
  w2u_method = int_params[1];
  verbose = int_params[2];

  if (NA<1) {
    return(INVLIB_ERR_DATAEMPTY);
  }

  /* check for invalid ialpha0 */
  if ((*ialpha0 < 0) || (*ialpha0>=NA)) {
    return(INVLIB_ERR_IALPHA0);
  }

  /* check for invalid inputs in arrays */
  for (i=0; i < NA; i++) {
    if ((PAgrid[i]<0) || (PAgrid[i]>180.0) || (!gsl_finite(H[i]))) {
      return(INVLIB_ERR_DATANAN);
    }
  }
  if ((verbose <0) || (verbose > 3)) {
    return(INVLIB_ERR_VERBOSE);
  } else if ((verbose==3) && (! outFile)) {
    return(INVLIB_ERR_VERBOSE); /* user requested output stream is NULL */
  }

  switch(w2u_method) {
    /* list valid methods here */
  case W2U_METHOD_TEM1:
    ret = wide2uni_TEM1(wideflux,dlogwideflux,PAgrid,H,int_params,real_params,outFile,ialpha0,uniflux,dloguniflux);
    break;
  default:
    return(INVLIB_ERR_INVALIDMETH);
  }


  /* done, success ! */
  return(ret);
}

int omni2uni(const double *omniflux, const double *dlogomniflux,
             const long int *int_params, 
             const double *real_params,
             char *outFile, 
             double *uniflux, double *dloguniflux) {
  /*
    Uses TEM-1 method to compute omni2uni.
    Arguments are explained in invlib.pdf
   */
  /* this routine is just a wrapper for wide2uni */
  int ret; /* return value */
  long int NA; /* number of pitch angle grid points */
  double *PAgrid; /* local pitch angle grid from "small" to 90 */
  double *H; /* local pitch angle grid from "small" to 90 */
  double da;
  long int i;

  if (!int_params) {
    return(INVLIB_ERR_NULL);
  }
  NA = int_params[0];
  PAgrid = malloc(NA*sizeof(double));
  H = malloc(NA*sizeof(double));

  da = 90.0/((double)NA); /* dalpha */ 
  /* alpha runs from da to 90 degrees */
  for (i=0; i <NA; i++) {
    PAgrid[i] = (1.0+i)*da;
    H[i] = sin(PAgrid[i]*M_PI/180.0)*da*M_PI/180.0; /* sin(alpha)dalpha */
  }
  i = NA-1; /* last grid point, i.e., locally-mirroring */
  H[i] /= 2; /* trapezoidal integral */ 
  ret = wide2uni(omniflux,dlogomniflux,PAgrid,H,int_params,real_params,outFile,&i,uniflux,dloguniflux);
  free(PAgrid);
  free(H);
  return(ret);
}

