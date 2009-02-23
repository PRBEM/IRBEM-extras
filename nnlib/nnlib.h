/* DLL headers for neural network library. See nnlib.pdf. */

double nnlib_fit(const unsigned long int Nt, const unsigned long int Nx, const double *X, 
		 const unsigned long int Nh, 
		 const unsigned long int Ny, const double *Y, 
		 const double *s, const unsigned long int flag, 
		 const unsigned long int MaxIter, 
		 double epsabs,
		 double *theta, 
		 double *xbar, double *ybar, double *sx, double *sy,
		 double *theta_cov);

void nnlib_save_training_set(const char *filename, const unsigned long int Nt, 
			     const unsigned long int Nx, const double *X, 
			     const unsigned long int Ny, const double *Y, 
			     const double *s, const unsigned long int sflag); 

void nnlib_load_training_set(const char *filename, unsigned long int *Nt, unsigned long int *Nx, double *X, 
			     unsigned long int *Ny, double *Y, 
			     double *s, unsigned long int *sflag); 
void nnlib_eval(const unsigned long int Nt,
		const unsigned long int Nx, const double *X, 
		const unsigned long int Nh, const double *theta, 
		const double *xbar, const double *ybar,
		const double *sx, const double *sy,
		const unsigned long int Ny, double *Y, 
		const unsigned long int dY_flag, const double *theta_cov, double *dY);

void nnlib_save_net(const char *filename, const unsigned long int Nx, 
		    const unsigned long int Nh, const unsigned long int Ny, const double *theta, 
		    const double *xbar,const double *ybar, const double *sx, const double *sy,
		    const unsigned long int cov_flag, const double *theta_cov);

void nnlib_load_net(const char *filename, unsigned long int *Nx, 
		    unsigned long int *Nh, unsigned long int *Ny, 
		    double *theta,
		    double *xbar,double *ybar, double *sx, double *sy,
		    unsigned long int *cov_flag, double *theta_cov);

double nnlib_ell(const unsigned long int Nt, const unsigned long int Nx, const double *Z, 
		 const unsigned long int Nh, 
		 const unsigned long int Ny, const double *Y, 
		 const double *xbar, const double *ybar, 
		 const double *sx, const double *sy,
		 const double *s, const unsigned long int sflag,
		 const double *theta,
		 double *grad, double *hess);
