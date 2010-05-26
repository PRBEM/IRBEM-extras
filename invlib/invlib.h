/*
  Spectral inversion library.
  See invlib.pdf for explanations.
 */


int pc_spec_inv(const double *y, const double *dy, const double *Egrid, const double *H, const double *b,
		const double *mean_log_flux, const double *basis_vectors, const double *basis_variance,
		const long int *int_params, const double *real_params,
		char *outFile, double *flux, double *dlogflux, double *lambda, double *support_data);

int ana_spec_inv(const double *y, const double *dy, const double *Egrid, const double *H, const double *b,
		  const long int *int_params, const double *real_params,
		 char *outFile, double *Eout, double *flux, double *dlogflux, double *lambda, double *support_data);

int ana_spec_inv_multi(const long int Ntimes,
		       const double *y, const double *dy, 
		       const double *Egrid, const double *H0, 
		       const double *dt, const double *b,
		       const long int *int_params, const double *real_params,
		       char *outFile, double *Eout, 
		       double *flux, double *dlogflux, double *lambda, double *support_data, int *result_codes);

int omni2uni(const double *omniflux, const double *dlogomniflux,
             const long int *int_params, 
             const double *real_params,
             char *outFile, 
             double *uniflux, double *dloguniflux);

int wide2uni(const double *wideflux, const double *dlogwideflux,
             const double *PAgrid, const double *H, 
             const long int *int_params, 
             const double *real_params,
             char *outFile,
             const long int *ialpha0,
             double *uniflux, double *dloguniflux);
