#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include "invlib_const.h" /* defines constants, e.g., return codes */
#include "ae8_atmocutoff.h" /* defines cutoff pitch angle vs L */

double vampola_sin_n(double Lm) {
  /* PAD exponent in sin^n (alpha) from 
     Vampola's CRRES electron data
     spans L = 3 to 8 in 1/4 L bins.
     Linear interpolate, nearaest extrapolate
     L        n
    3.00    5.380
    3.25    5.078
    3.50    4.669
    3.75    3.916
    4.00    3.095
    4.25    2.494
    4.50    2.151
    4.75    1.998
    5.00    1.899
    5.25    1.942
    5.50    1.974
    5.75    1.939
    6.00    1.970
    6.25    2.136
    6.50    1.775
    6.75    1.438
    7.00    1.254
    7.25    1.194
    7.50    1.046
    7.75    0.989
    8.00    0.852
  */
  const double L0 = 3, Lend = 8, dLinv = 4; /* dL = 1/4 */
  const int N = 21;
  const double n[21] = {5.380,5.078,4.669,3.916,3.095,2.494,
		 2.151,1.998,1.899,1.942,1.974,1.939,
		 1.970,2.136,1.775,1.438,1.254,1.194,
		 1.046,0.989,0.852};
  double c,Li;
  int i;
  if (!gsl_finite(Lm)) {
    /* bad data */
    return(GSL_NAN);
  } else if (Lm <= L0) {
    /* plateau extrapolate */
    return(n[0]);
  } else if (Lm >= Lend) {
    /* plateau extrapolate */
    return(n[N-1]);
  } else {
    /* interpolate */
    i = floor((Lm-L0)*dLinv);
    Li = (L0+((double)i)/dLinv);
    c = (Lm-Li)*dLinv;
    return(n[i]*(1.0-c) + n[i+1]*c);
  }
}


int wide2uni_vampola(const double *wideflux, const double *dlogwideflux,
		     const double *PAgrid, const double *H, 
		     const long int *int_params, 
		     const double *real_params,
		     char *outFile,
		     const long int *ialpha0,
		     double *uniflux, double *dloguniflux) {
  /*
    Uses Vampola's sin^n to compute wede2uni.
    Arguments are same as for wide2uni.
   */
  long int NA, verbose;
  double Lm,BB0,EPAcutoff, PAcutoff;
  double integral=0,n;
  long int i;

  int closeOutFile = 0; /* default: no need to close outFile */
  FILE *outFilePtr = NULL; 

  /* inputs used by all methods already checked: */
  /* inputs already checked for NULLs */
  /* NA, w2u_method, and verbose are already checked for validity */
  NA = int_params[0];
  /* w2u_method = int_params[1]; */
  verbose = int_params[2];
  /* maybe not all methods use these, so check them here: */
  /* minimizer_flag = int_params[3]; */
  /* MaxIter = int_params[4]; */

  /* keV = real_params[0]; */ /* Electron energy, keV */
  BB0 = real_params[1]; /* B/B0, OPQ */
  Lm = real_params[2]; /* McIlwain L, OPQ */

  if ((!gsl_finite(Lm)) || (Lm < 1) || (!gsl_finite(BB0)) || (BB0 < 1)) {
    return(INVLIB_ERR_DATANAN);
  }

  n = vampola_sin_n(Lm); /* get exponent */

  switch(verbose) { /* this switch statement must be consistent with specinv.h */
  case 0:
    outFilePtr = NULL;
    break;
  case 1:
    outFilePtr = stdout;
    break;
  case 2:
    outFilePtr = stderr;
    break;
  case 3:
    if (outFile) {
      outFilePtr = (FILE *)outFile;
    } else {
      return(INVLIB_ERR_VERBOSE); /* can't output to NULL */
    }
    break;
  case 4:
    if (outFile) {
      if (! (outFilePtr = fopen(outFile,"w"))) {
	return(INVLIB_ERR_OUTFILE); /* Couldn't open outFile */
	  }
      closeOutFile = 1;
    } else {
      return(INVLIB_ERR_VERBOSE); /* can't output to NULL */
    }
    break;
  case 5:
    if (outFile) {
      if (! (outFilePtr = fopen(outFile,"a"))) {
	return(INVLIB_ERR_OUTFILE); /* Couldn't open outFile */
	  }
      closeOutFile = 1;
    } else {
      return(INVLIB_ERR_VERBOSE); /* can't output to NULL */
    }
    break;
  default:
    return(INVLIB_ERR_VERBOSE); /* invalid verbose value */
  }

  EPAcutoff = AE8_AtmoCutoff(Lm); /* acute equatorial pitch angle of loss cone */
  PAcutoff = asin(sin(EPAcutoff*M_PI/180.0)*sqrt(BB0))*180.0/M_PI; /* acute local pitch angle */
  if (verbose) {
    fprintf(outFilePtr,"Vampola wide2uni with Lm=%g, BB0 = %g, n=%g\n",Lm,BB0,n);
  }

  /* do numerical integral */
  for (i=0; i < NA; i++) {
    if ((PAgrid[i] > PAcutoff) && (PAgrid[i] < 180-PAcutoff)) {
      integral += pow(sin(PAgrid[i]*M_PI/180.0),n)*H[i]; /* integral sin(alpha)^n*H */
    }
  }

  if ((integral<=0) || (!gsl_finite(integral))) {
    /* any pitch angle grid points left? */
    if (closeOutFile)
      fclose(outFilePtr);
    return(INVLIB_ERR_DATANAN);
  }

  /* get uniflux at alpha0 */
  *uniflux = (*wideflux)/integral*pow(sin(PAgrid[*ialpha0]*M_PI/180.0),n); /* f(a0) = fiso/integral*sin(a0)^n */
  *dloguniflux = *dlogwideflux; /* preserve error */

  /* done, success! */
  if (closeOutFile)
    fclose(outFilePtr);
  return(INVLIB_SUCCESS);
}

