int wide2uni_vampola(const double *wideflux, const double *dlogwideflux,
		  const double *PAgrid, const double *H, 
		  const long int *int_params, 
		  const double *real_params,
		  char *outFile,
		  const long int *ialpha0,
		  double *uniflux, double *dloguniflux);
  /*
    Uses Vampola's sin^n method to compute wede2uni and omni2uni.
    Arguments are same as for wide2uni.
   */
