#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "invlib.h" /* defines ana_spec_inv */
#include "invlib_const.h" /* defines constants, e.g., return codes */
#include "invutil.h" /* defines inv_matrix_once */
#include "optim.h" /* optimize and optfunTy */
#include "specinv.h" /* defines macros */

/* penalty function type */
typedef double (*pen_funcTy)(const double lambda,const double y, const void *params, double *grad, double *hess);

/* flux spectrum function type */
typedef double (*fluxFuncTy)(const gsl_vector *q, const double E,const void *params,gsl_vector *grad, gsl_matrix *hess);

/* holds params for likelihood function for a single analytical fit spectrum */
typedef struct {
  long int NY,NE,Nq; /* size of vectors */
  long int verbose; /* verbose setting */
  double *y,*dy; /* counts, relative error for each channel */
  pen_funcTy *pen_funcs; /* penalty function for each channel */
  const gsl_matrix *H; /* measurement matrix */
  const gsl_vector *Egrid, *b; /* energy grid, background counts*/
  fluxFuncTy flux_func; /* analytical flux function being used */
  void *flux_func_params; /* parameters of flux func */
  FILE *outFilePtr; /* output file pointer */
  int closeOutFile; /* does output file pointer need to be closed on exit? */
} asi_ell_paramsTy;


/********************************************************************************************************
ANA_SPEC_INV Section - this section defines ana_spec_inv and the
routines it needs
********************************************************************************************************/

/********************************************************************************************************
penalty functions:
each penalty function takes as input:
 y - observed counts
 lambda - estimated counts
 params - some parameters (not always used)
and returns:
 result: negative natural log of likelihood function without constants (i.e., which don't depend on lambda)
 *grad: gradient with respect to lambda
 *hess: Hessian (second derivative) with respect to lambda
Note:
 if grad is NULL, neither grad nor hess is computed
 if hess is NULL, hess is not computed
********************************************************************************************************/

double nan_pen(const double lambda,const double y, const void *params, double *grad, double *hess) {
  /* dummy penalty function for NaN's */
  /* params not used */
  double pen;
  pen = 0;
  if (grad) {
    *grad = 0;
    if (hess) {
      *hess = 0;
    }
  }
  return(pen);
}

double poiss_pen(const double lambda,const double y, const void *params, double *grad, double *hess) {
/* negative log likelihood of poisson penalty function */
  /* params not used */
  double pen;
  pen = lambda-y*log(lambda);
  if (grad) {
    *grad = 1-y/lambda;
    if (hess) {
      *hess = y/gsl_pow_2(lambda);
    }
  }
  return(pen);
}

double rel_pen(const double lambda,const double y,const void *params, double *grad, double *hess) {
  /* negative log likelihood of gaussian relative error penalty function */
  /* params provides a const double *, with first element giving "sigma" the relative error */
  /* this routine has been validated against matlab counterpart on a test case,
   including the hess calculation*/
  double pen,logy,z,loglam,sigma;
  logy = log(y);
  loglam = log(lambda);
  sigma = *((double *)params);
  z = (logy-loglam)/sigma;
  pen = 0.5*gsl_pow_2(z);
  if (grad) {
    (*grad) = -z/sigma/lambda;
    if (hess) {
      *hess = (1.0+logy-loglam)/gsl_pow_2(sigma*lambda);
    }
  }
  return(pen);
}

/********************************************************************************************************
flux spectrum functions:
each flux function takes as input:
 q - parameters/coefficients of spectral function
 E - energy at which to evaluate
 params - some parameters (not always used)
and returns:
 result: flux at energy E
 grad: gradient with respect to q's
 hess: Hessian (second derivative matrix) with respect to q's
Note:
 if grad is NULL, neither grad nor hess is computed
 if hess is NULL, hess is not computed

********************************************************************************************************/
double flux_exp(const gsl_vector *q, const double E,
		const void *params,
		gsl_vector *grad, gsl_matrix *hess) {
  /* exponential spectrum */
  /* params not used */
  double p1,p2, flux, tmp;
  p1 = gsl_vector_get(q,0);
  p2 = gsl_vector_get(q,1);
  flux = exp(p1+E*p2);
  if (grad) {
    gsl_vector_set(grad,0,flux);
    tmp = flux*E;
    gsl_vector_set(grad,1,tmp);
    if (hess) {
      gsl_matrix_set(hess,0,0,flux);
      gsl_matrix_set(hess,1,0,tmp);
      gsl_matrix_set(hess,0,1,tmp);
      gsl_matrix_set(hess,1,1,tmp*E);
    }
  }
  return(flux);
}

double flux_rm(const gsl_vector *q, const double E,
		const void *params,
		gsl_vector *grad, gsl_matrix *hess) {
  /* relativistic maxwellian spectrum, params = (void *)(&m0c2) */
  /* m0c2 is the rest energy of the particle species, double precision */
  double p1,p2, flux, tmp;
  double m0c2 = *((double *)params);
  p1 = gsl_vector_get(q,0);
  p2 = gsl_vector_get(q,1);
  flux = E*(1+E/2/m0c2)*exp(p1+E*p2);
  if (grad) {
    gsl_vector_set(grad,0,flux);
    tmp = flux*E;
    gsl_vector_set(grad,1,tmp);
    if (hess) {
      gsl_matrix_set(hess,0,0,flux);
      gsl_matrix_set(hess,1,0,tmp);
      gsl_matrix_set(hess,0,1,tmp);
      gsl_matrix_set(hess,1,1,tmp*E);
    }
  }
  return(flux);
}

double flux_rm2(const gsl_vector *q, const double E,
		const void *params,
		gsl_vector *grad, gsl_matrix *hess) {
  /* double relativistic maxwellian spectrum, params = (void *)(&m0c2) */
  /* m0c2 is the rest energy of the particle species, double precision */
  /* grad is a 4-element vector of deriviatives with respect to q */
  /* hess is a 4x4 matrix of second deriviatives with respect to q */
  double p1,p2,p3,p4, flux, flux1,flux2, tmp1,tmp2;
  double m0c2 = *((double *)params);

  /* first population */
  p1 = gsl_vector_get(q,0);
  p2 = gsl_vector_get(q,1);
  flux1 = E*(1+E/2/m0c2)*exp(p1+E*p2);

  /* second population */
  p3 = gsl_vector_get(q,2);
  p4 = gsl_vector_get(q,3);
  flux2 = E*(1+E/2/m0c2)*exp(p3+E*p4);

  /* total */
  flux = flux1+flux2;
  if (grad) {
    /* grad and hess are independent */

    /* first population */
    gsl_vector_set(grad,0,flux1);
    tmp1 = flux1*E;
    gsl_vector_set(grad,1,tmp1);

    /* second population */
    gsl_vector_set(grad,2,flux2);
    tmp2 = flux2*E;
    gsl_vector_set(grad,3,tmp2);

    if (hess) { /* second derivative matrix */
      /* upper-right, and lower-left 2x2 submatrices are zeros */
      gsl_matrix_set_zero(hess);

      /* first population, upper-left */
      gsl_matrix_set(hess,0,0,flux1);
      gsl_matrix_set(hess,1,0,tmp1);
      gsl_matrix_set(hess,0,1,tmp1);
      gsl_matrix_set(hess,1,1,tmp1*E);

      /* second population, lower-right */
      gsl_matrix_set(hess,2,2,flux2);
      gsl_matrix_set(hess,3,2,tmp2);
      gsl_matrix_set(hess,2,3,tmp2);
      gsl_matrix_set(hess,3,3,tmp2*E);
    }
  }
  return(flux);
}

double flux_pl(const gsl_vector *q, const double E,
		const void *params,
		gsl_vector *grad, gsl_matrix *hess) {
  /* power law spectrum */
  /* this routine has been validated against matlab
     counterpart for single test case, including hess */
  double p1,p2, flux, neglogE,tmp;
  neglogE = -log(E);
  p1 = gsl_vector_get(q,0);
  p2 = gsl_vector_get(q,1);
  flux = exp(p1+p2*neglogE);
  if (grad) {
    gsl_vector_set(grad,0,flux);
    tmp = neglogE*flux;
    gsl_vector_set(grad,1,tmp);
    if (hess) {
      gsl_matrix_set(hess,0,0,flux);
      gsl_matrix_set(hess,1,0,tmp);
      gsl_matrix_set(hess,0,1,tmp);
      gsl_matrix_set(hess,1,1,neglogE*tmp);
    }
  }
  return(flux);
}

double flux_ple(const gsl_vector *q, const double E,
		const void *params,
		gsl_vector *grad, gsl_matrix *hess) {
  /* power law spectrum with exponential after E_break */
  double p1,p2, flux,E0,E_break, neglogE,neglogE_break,tmp;
  double *dbl_params = (double *)params;
  E_break = dbl_params[0];
  E0 = dbl_params[1];
  neglogE_break = -log(E_break);
  neglogE = -log(E);
  p1 = gsl_vector_get(q,0);
  p2 = gsl_vector_get(q,1);
  if (E>E_break) {
    flux = exp(p1+p2*neglogE_break-(E-E_break)/E0);
    neglogE = neglogE_break; /* overwrite neglogE for derivatives */
  } else {
    flux = exp(p1+p2*neglogE);
  }
  if (grad) {
    gsl_vector_set(grad,0,flux);
    tmp = neglogE*flux;
    gsl_vector_set(grad,1,tmp);
    if (hess) {
      gsl_matrix_set(hess,0,0,flux);
      gsl_matrix_set(hess,1,0,tmp);
      gsl_matrix_set(hess,0,1,tmp);
      gsl_matrix_set(hess,1,1,neglogE*tmp);
    }
  }
  return(flux);
}

void flux_multi(fluxFuncTy flux_func, 
		const gsl_vector *q, const gsl_vector *Egrid, const void *params,
		gsl_vector *flux,gsl_matrix *grad_matrix) {
  /* evaluate a flux function for an array of energies
     return gradient as grad_matrix [NE x Nq] */
  /* inputs:
     flux_func - one of the flux spectrum functions above (e.g., flux_pl).
     q - parameters/coefficients of spectral function
     Egrid - vector of energies at which to evaluate
     params - some parameters, passed to flux_func
     outputs:
     flux: vector of fluxes at energies in Egrid
     grad_matrix: gradient of flux with respect to q's
     Note:
     if grad_matrix is NULL, no gradient is computed
   */
  long int i;
  double tmp,E;
  gsl_vector_view row; /* use a row view to reference rows of grad_matrix */
  gsl_vector *rowptr=0;
  for (i=1; i <= Egrid->size; i++) {
    if (grad_matrix) {
      row = gsl_matrix_row(grad_matrix,i-1);
      rowptr = &(row.vector);
    } else {
      rowptr = NULL;
    }
    E = gsl_vector_get(Egrid,i-1);
    tmp = (*flux_func)(q,E,params,rowptr,NULL);
    gsl_vector_set(flux,i-1,tmp);
  }
}

void make_lambda(const gsl_vector *q,const asi_ell_paramsTy *asi_ell_params,gsl_vector *lambda) {
  /* sets lambda based on q and info in asi_ell_params: H, b, flux_func, etc */
  gsl_vector *flux = gsl_vector_alloc(asi_ell_params->NE);
  flux_multi(asi_ell_params->flux_func,q,asi_ell_params->Egrid,asi_ell_params->flux_func_params,flux,NULL);
  /* copy background into lambda: lambda = b */
  gsl_vector_memcpy(lambda,asi_ell_params->b);  
  /* lambda = 1.0*H*flux+1.0*b (since lambda=b) */
  gsl_blas_dgemv(CblasNoTrans,1.0,asi_ell_params->H,flux,1.0,lambda);
  gsl_vector_free(flux);
}


/****************************
  asi_ell_combine - the master negative-log-likelihood function
  inputs:
     q - parameters/coefficients of spectral function
     params_void - parameters, an asi_ell_params type, holds all kinds of info for optimization
     grad: gradient with respect to q's
     hess: Hessian (second derivative matrix) with respect to q's
 Note:
    if grad is NULL, neither grad nor hess is computed
    if hess is NULL, hess is not computed
 ****************************/

double asi_ell_combine(const gsl_vector *q, void *params_void, gsl_vector *grad, gsl_matrix *hess) {
  double ell=0; /* neglogp */
  double E,h;
  long int j,i;
  double grad_ell_lambdak,hess_ell_lambdak;
  double *grad_ell_lambdak_ptr=0,*hess_ell_lambdak_ptr=0;
  gsl_vector *lambda=0,*flux=0,*grad_lambda_q=0,*grad_lambda_qj=0;
  gsl_matrix *hess_lambda_q=0,*hess_lambda_qj=0;
  asi_ell_paramsTy *asi_ell_params = (asi_ell_paramsTy *)params_void;

  flux = gsl_vector_alloc(asi_ell_params->NE);
  flux_multi(asi_ell_params->flux_func,q,asi_ell_params->Egrid,asi_ell_params->flux_func_params,flux,NULL);

  lambda = gsl_vector_alloc(asi_ell_params->NY);
  /* copy background into lambda: lambda = b */
  gsl_vector_memcpy(lambda,asi_ell_params->b);  
  /* lambda = 1.0*H*flux+1.0*b (since lambda=b) */
  gsl_blas_dgemv(CblasNoTrans,1.0,asi_ell_params->H,flux,1.0,lambda);

  /* have verified that lambda agrees with matlab result */

  if (grad) { /* prepare space for gradient computation */
    gsl_vector_set_zero(grad);
    grad_lambda_qj = gsl_vector_alloc(q->size);
    grad_ell_lambdak_ptr = &grad_ell_lambdak; /* if grad is NULL, this remains NULL, so won't be computed in ancillary routines */
    if (hess) { /* prepare space for hessian computation */
      gsl_matrix_set_zero(hess);
      grad_lambda_q = gsl_vector_alloc(q->size); /* stores sum over j */
      hess_lambda_q = gsl_matrix_alloc(q->size,q->size);
      hess_lambda_qj = gsl_matrix_alloc(q->size,q->size);
      hess_ell_lambdak_ptr = &hess_ell_lambdak;  /* if hess is NULL, this remains NULL, so won't be computed in ancillary routines */
    }
  }

  for (i=1; i <= asi_ell_params->NY; i++) {
    /* get the ell increment and gradient/hessian as needed */
    /* grad_ell_lambdak, hess_ell_lambdak are scalars */
    ell += (asi_ell_params->pen_funcs[i-1])(gsl_vector_get(lambda,i-1),asi_ell_params->y[i-1],
				    &(asi_ell_params->dy[i-1]), grad_ell_lambdak_ptr, hess_ell_lambdak_ptr);
    if (grad) { /* compute gradient */
      if (hess) {
	/* hess needs to accumulate grad and hess over j, so clear these variables to be accumulators */
	gsl_vector_set_zero(grad_lambda_q); 
	gsl_matrix_set_zero(hess_lambda_q);
      }
      for (j=1; j<=asi_ell_params->NE; j++) {
	h = gsl_matrix_get(asi_ell_params->H,i-1,j-1); /*H_ij*/
	E = gsl_vector_get(asi_ell_params->Egrid,j-1); /* E_j */
	/* now, evaulate flux function at E(j) and get grad/hess as needed */
	asi_ell_params->flux_func(q,E,asi_ell_params->flux_func_params,grad_lambda_qj, hess_lambda_qj);
	/* convert gradients of flux wrt q into gradients of lambda wrt q */
	gsl_vector_scale(grad_lambda_qj,h); /* multiply by H_ij */
	gsl_blas_daxpy(grad_ell_lambdak,grad_lambda_qj,grad); /* grad += grad_ell_lambdak*dlambda_dp */
	if (hess) {
	  gsl_vector_add(grad_lambda_q,grad_lambda_qj); /* accumulate */
	  gsl_matrix_scale(hess_lambda_qj,h); /* multiply by H_ij */
	  gsl_matrix_add(hess_lambda_q,hess_lambda_qj); /* accumulate */
	}
      }
      if (hess) {
	/* convert hessians of flux wrt q into hessians of lambda wrt q */
	/* hess = hess+dlambda_dp*hess_ell_lambda*dlambda_dp' + grad_ell_lambdak*hess_lambda_p; */
	gsl_blas_dger(hess_ell_lambdak,grad_lambda_q,grad_lambda_q,hess); /* hess += dlambda_dp*hess_ell_lambda*dlambda_dp' */
	/* hess += dell_dlambda*hess_lambda_p; */
	gsl_matrix_scale(hess_lambda_q,grad_ell_lambdak); /* hess_lambda_p *= grad_ell_lambdak */
	gsl_matrix_add(hess,hess_lambda_q); /* hess += hess_lambda_p*grad_ell_lambdak */
      }
    }
  }
  /* free memory */
  gsl_vector_free(lambda);
  gsl_vector_free(flux);
  if (grad) {
    gsl_vector_free(grad_lambda_qj);
    if (hess) {
      gsl_vector_free(grad_lambda_q);
      gsl_matrix_free(hess_lambda_q);
      gsl_matrix_free(hess_lambda_qj);
    }
  }
  /* have verified that ell, grad, and hess are right, compared to matlab */
  return(ell);
}

int ana_spec_inv(const double *y, const double *dy, const double *Egrid, const double *H, const double *b,
		 const long int *int_params, const double *real_params,
		 char *outFile, double *Eout, double *flux, double *dlogflux, 
		 double *lambda, double *support_data) {
/****
     see header file for details of inputs/outputs
****/
  long int NY, NE, NEout, fxn_bit_map, minimizer_flag, MaxIter, verbose, dE_mode;
  long int fxn_bit;
  long int i,j; /* loop control vars */
  double m0c2=0,E0=0,E_break=0;
  double dbl_params[10]; /* holder for parameters to flux functions */
  int any_y_positive = 0, Ny_valid=0;
  gsl_vector_view fluxgsl,dlogfluxgsl,grad_matrix_row,lambda_view; 
  gsl_vector *fluxhat=0, *sigmahat=0;
  gsl_vector *q=0; /* parameter vector */
  gsl_vector *grad=0;
  gsl_matrix *hess=0, *covq=0, *grad_matrix=0;
  pen_funcTy *pen_funcs=0;
  asi_ell_paramsTy asi_ell_params;
  optfunTy optfun; 
  gsl_vector *qs[ASI_MAX_POW2+1],*fluxes[ASI_MAX_POW2+1],*sigmas[ASI_MAX_POW2+1],*lambdak[ASI_MAX_POW2+1];
  double ells[ASI_MAX_POW2+1];
  double weights[ASI_MAX_POW2+1],weight_sum,min_ell;
  double sigma_flux_j, flux_j, logflux_j, sigma2_log_flux_j, tmp_dbl;
  gsl_matrix *Htmp=0; /* holds H*dE if necessary */
  /* the following "const" view variables will be declared inline below
  gsl_matrix_const_view Hgsl = gsl_matrix_const_view_array(H,NY,NE); gsl view of H 
  gsl_vector_const_view Egridgsl = gsl_vector_const_view_array(Egrid,NE); gsl view of Egrid
  gsl_vector_const_view Eoutgsl = gsl_vector_const_view_array(Eout,NEout); gsl view of Eout
  gsl_vector_const_view bgsl = gsl_vector_const_view_array(b,NY); gsl view of b
  */

  int clean_return(int return_value) {
    if (asi_ell_params.closeOutFile) {
      fclose(asi_ell_params.outFilePtr);
    }
    return(return_value);
  };

  asi_ell_params.closeOutFile = 0; /* default: no need to close outFile */

  /* test for input NULLs */
  if ((!y) || (!dy) || (!Egrid) || (!H) || (!b) || (!int_params) || (!flux) || (!dlogflux)) {
    return(clean_return(INVLIB_ERR_NULL));
  }

  /* start initialization */
  NY = int_params[0];
  NE = int_params[1];
  NEout = int_params[2];

  /* check for not enough channels */
  if ((NY<=1) || (NE<=1) || (NEout < 1)) {
    return(clean_return(INVLIB_ERR_DATAEMPTY));
  }

  /* check for NaNs, Infs, and inappropriate negative numbers */
  for (i=0; i < NY; i++) {
    if (gsl_finite(y[i])) {
      Ny_valid++;
      if ((! gsl_finite(dy[i])) || (! gsl_finite(b[i]))) {
	return(clean_return(INVLIB_ERR_DATANAN));
	
      }
      if ((y[i]<0) || (dy[i]<0) || (b[i]<0)) {
	return(clean_return(INVLIB_ERR_DATANAN));
      }
      if (y[i]>b[i]) { /* also check for at least one count above background in some channel */
	any_y_positive = 1;
      }
    }
      
    for (j=0; j < NE; j++) {
      if (! gsl_finite(H[NY*j+i])) {
	return(clean_return(INVLIB_ERR_DATANAN));
      }
    }
  }

  /* check for not enough valid channels */
  if (Ny_valid<=1) {
    return(clean_return(INVLIB_ERR_DATAEMPTY));
  }


  /* confirm at least one count in some channel */
  if (! any_y_positive) {
    return(clean_return(INVLIB_ERR_NOCOUNTS));
  }

  /* check for NaNs and Infs or <= in Egrid */
  for (j=0; j < NE; j++) {
    if ((!gsl_finite(Egrid[j])) || (Egrid[j]<=0)) {
      return(clean_return(INVLIB_ERR_DATANAN));
    }
  }

  /* check for NaNs and Infs or <=0 in Eout */
  for (j=0; j < NEout; j++) {
    if (!gsl_finite(Eout[j])) {
      return(clean_return(INVLIB_ERR_DATANAN));
    }
  }

  /* initialize some more parameters */
  fxn_bit_map = int_params[3];
  minimizer_flag = int_params[4];
  MaxIter = int_params[5];
  verbose = int_params[6];
  dE_mode = int_params[7];

  /* check verbose setting */
  asi_ell_params.verbose = verbose;
  switch(verbose) { /* this switch statement must be consistent with specinv.h */
  case 0:
    asi_ell_params.outFilePtr = NULL;
    break;
  case 1:
    asi_ell_params.outFilePtr = stdout;
    break;
  case 2:
    asi_ell_params.outFilePtr = stderr;
    break;
  case 3:
    if (outFile) {
      asi_ell_params.outFilePtr = (FILE *)outFile;
    } else {
      return(clean_return(INVLIB_ERR_VERBOSE)); /* can't output to NULL */
    }
    break;
  case 4:
    if (outFile) {
      if (! (asi_ell_params.outFilePtr = fopen(outFile,"w"))) {
	return(clean_return(INVLIB_ERR_OUTFILE)); /* Couldn't open outFile */
	  }
      asi_ell_params.closeOutFile = 1;
    } else {
      return(clean_return(INVLIB_ERR_VERBOSE)); /* can't output to NULL */
    }
    break;
  case 5:
    if (outFile) {
      if (! (asi_ell_params.outFilePtr = fopen(outFile,"a"))) {
	return(clean_return(INVLIB_ERR_OUTFILE)); /* Couldn't open outFile */
	  }
      asi_ell_params.closeOutFile = 1;
    } else {
      return(clean_return(INVLIB_ERR_VERBOSE)); /* can't output to NULL */
    }
    break;
  default:
    return(clean_return(INVLIB_ERR_VERBOSE)); /* invalid verbose value */
  }

  /* check function bitmap for no set bit */
  if ((fxn_bit_map & ASI_FXN_ALL) == 0) {
    return(clean_return(INVLIB_ERR_NOFXN));
  };

  /* check function bitmap for invalid bit(s) */
  if ((fxn_bit_map & ASI_FXN_ALL) != fxn_bit_map) {
    return(clean_return(INVLIB_ERR_INVALIDFXN));
  }

  /* check for rest mass needed by Relativistic Maxwellian (single or double) */
  if (fxn_bit_map & (ASI_FXN_RM | ASI_FXN_RM2)) {
    if ((!real_params) || (real_params[0]<=0)) {
      return(clean_return(INVLIB_ERR_RME0)); /* Relativistic Maxwellian requiers positive E0=m0c2 */
    }
    m0c2 = real_params[0]; /* store it, if need be */
  }

  /* check for E_break and E0 needed by Power-Law-Exponential */
  if (fxn_bit_map & ASI_FXN_PLE) {
    if ((!real_params) || (real_params[1]<=0) || (real_params[2]<=0)) {
      return(clean_return(INVLIB_ERR_PLE)); /* PLE requires E0 and E_break */
    }
    E_break = real_params[1]; /* store it, if need be */
    E0 = real_params[2]; /* store it, if need be */
  }

  /* check for valid minimizer flag */
  if ((minimizer_flag < 0) || (minimizer_flag > OPTIM_MIN_MAX)) {
    return(clean_return(INVLIB_ERR_INVALIDMIN));
  }

  /* check for valid maximum interations */
  if (MaxIter<=0) {
    return(clean_return(INVLIB_ERR_INVALIDITER));
  }

  if (lambda || support_data) { /* prepare to compute lambda */
    for (i=0; i <= ASI_MAX_POW2; i++) {
      lambdak[i] = gsl_vector_alloc(NY);
    }
  }
  
  /* must declare these here because otherwise compiler treats assignments as violating const-ness */
  gsl_matrix_const_view Hgsl = gsl_matrix_const_view_array(H,NY,NE); /* gsl view of H */
  gsl_vector_const_view Egridgsl = gsl_vector_const_view_array(Egrid,NE);
  gsl_vector_const_view Eoutgsl = gsl_vector_const_view_array(Eout,NEout);
  gsl_vector_const_view bgsl = gsl_vector_const_view_array(b,NY);

  if (dE_mode != ASI_DE_INCLUDED) {
    gsl_vector *dE_tmp = gsl_vector_alloc(NE);
    for (i=1; i < NE-1; i++) {
      gsl_vector_set(dE_tmp,i,(gsl_vector_get(&(Egridgsl.vector),i+1)-gsl_vector_get(&(Egridgsl.vector),i-1))/2);
    }
    gsl_vector_set(dE_tmp,0,(gsl_vector_get(&(Egridgsl.vector),1)-gsl_vector_get(&(Egridgsl.vector),0)));
    gsl_vector_set(dE_tmp,NE-1,(gsl_vector_get(&(Egridgsl.vector),NE-1)-gsl_vector_get(&(Egridgsl.vector),NE-2)));
    if (dE_mode == ASI_DE_TRAPZ) { /* trapz divides the endpoint weights by 2 */
      gsl_vector_set(dE_tmp,0,gsl_vector_get(dE_tmp,0)/2);
      gsl_vector_set(dE_tmp,NE-1,gsl_vector_get(dE_tmp,NE-1)/2);
    }
    /* now create and scale columns of Htmp by dE */
    Htmp = gsl_matrix_alloc(NY,NE);
    gsl_matrix_memcpy(Htmp,&(Hgsl.matrix));
    for (i=0; i < NE; i++) {
      gsl_vector_view Hcol = gsl_matrix_column(Htmp,i);
      gsl_vector_scale(&(Hcol.vector),gsl_vector_get(dE_tmp,i));
    }
    gsl_vector_free(dE_tmp);
  }
  /* initialize other gsl views */
  fluxgsl = gsl_vector_view_array(flux,NEout);
  fluxhat = &(fluxgsl.vector);
  dlogfluxgsl = gsl_vector_view_array(dlogflux,NEout);
  sigmahat = &(dlogfluxgsl.vector);

  /* determine appropriate penalty function */
  pen_funcs = malloc(NY*sizeof(pen_funcTy));
  for (i=0; i < NY; i++) {
    if (! gsl_finite(y[i])) {
      if (verbose) {
	fprintf(asi_ell_params.outFilePtr,"Selecting NaN penalty function for y[%li]=%g\n",i,y[i]);
      }
      pen_funcs[i] = &nan_pen; /* dummy for non-finite y */
    } else if ((dy[i]==0) || (y[i] < 1.0/gsl_pow_2(dy[i]))) {
      if (verbose) {
	fprintf(asi_ell_params.outFilePtr,"Selecting Poisson penalty function for y[%li]=%g\n",i,y[i]);
      }
      pen_funcs[i] = &poiss_pen;
    } else {
      if (verbose) {
	fprintf(asi_ell_params.outFilePtr,"Selecting Gaussian Relative Error penalty function for y[%li]=%g\n",i,y[i]);
      }
      pen_funcs[i] = &rel_pen;
    }
  }

  min_ell = GSL_POSINF; /* initialize to biggest number */

  /* prepare parameters for asi_ell_combine */
  asi_ell_params.NY = NY;
  asi_ell_params.NE = NE;
  asi_ell_params.pen_funcs = pen_funcs;
  if (Htmp) {
    asi_ell_params.H = Htmp;
  } else {
    asi_ell_params.H = &(Hgsl.matrix);
  }
  asi_ell_params.b = &(bgsl.vector);
  asi_ell_params.y = y;
  asi_ell_params.dy = dy;
  asi_ell_params.Egrid = &(Egridgsl.vector);
  /* loop through the analytical functions */
  for (i=0; i <= ASI_MAX_POW2; i++) {
    /* initialize to NULL */
    qs[i] = NULL; /* qs[i]=NULL will later be used to indicate the function isn't used */
    fluxes[i] = NULL;
    sigmas[i] = NULL;
    /* create i'th spectral function bit */
    fxn_bit = (1<<i);
    if (fxn_bit & fxn_bit_map) { /* is this function's bit in the request bitmap? */
      /* initialize for this function */
      switch(fxn_bit) {
      case ASI_FXN_PL:
	/* power law */
	/* allocate, this'll get stored in qs[i] and freed at the end of the function */
	q = gsl_vector_alloc(2); /* two free parameters */
	/* initialize with default params */
	gsl_vector_set(q,0,1);
	gsl_vector_set(q,1,3);
	asi_ell_params.flux_func = flux_pl; /* power law flux spectrum */
	asi_ell_params.flux_func_params = (void *)NULL; /* no constant params */
	if (verbose) {
	  fprintf(asi_ell_params.outFilePtr,"Trying Power-Law\n");
	}
	break;
      case ASI_FXN_EXP:
	/* exponential */
	/* allocate, this'll get stored in qs[i] and freed at the end of the function */
	q = gsl_vector_alloc(2); /* two free parameters */
	/* initialize with default params */
	gsl_vector_set(q,0,1);
	gsl_vector_set(q,1,-1);
	asi_ell_params.flux_func = flux_exp; /* exponential flux spectrum */
	asi_ell_params.flux_func_params = (void *)NULL; /* no constant params */
	if (verbose) {
	  fprintf(asi_ell_params.outFilePtr,"Trying Exponential\n");
	}
	break;
      case ASI_FXN_RM:
	/* relativistic maxwellian */
	/* allocate, this'll get stored in qs[i] and freed at the end of the function */
	q = gsl_vector_alloc(2); /* two free parameters */
	/* initialize with default params */
	gsl_vector_set(q,0,1);
	gsl_vector_set(q,1,-1);
	asi_ell_params.flux_func = flux_rm; /* relativistic Maxwellian */
	asi_ell_params.flux_func_params = (void *)(&m0c2); /* one constant param, E0=m0c2 */
	if (verbose) {
	  fprintf(asi_ell_params.outFilePtr,"Trying Relativistic Maxwellian\n");
	}
	break;
      case ASI_FXN_RM2:
	/* double relativistic maxwellian */
	/* allocate, this'll get stored in qs[i] and freed at the end of the function */
	q = gsl_vector_alloc(4); /* four free parameters */
	/* initialize with default params */
	gsl_vector_set(q,0,1);
	gsl_vector_set(q,1,-1);
	gsl_vector_set(q,2,2);
	gsl_vector_set(q,3,-2);
	asi_ell_params.flux_func = flux_rm2; /* relativistic Maxwellian */
	asi_ell_params.flux_func_params = (void *)(&m0c2); /* one constant param, E0=m0c2 */
	if (verbose) {
	  fprintf(asi_ell_params.outFilePtr,"Trying Double Relativistic Maxwellian\n");
	}
	break;
      case ASI_FXN_PLE:
	/* power law w/ exponential tail*/
	/* allocate, this'll get stored in qs[i] and freed at the end of the function */
	q = gsl_vector_alloc(2); /* two free parameters */
	/* initialize with default params */
	gsl_vector_set(q,0,1);
	gsl_vector_set(q,1,3);
	asi_ell_params.flux_func = flux_ple; /* power law flux spectrum w/ exp tail */
	dbl_params[0] = E_break;
	dbl_params[1] = E0;
	asi_ell_params.flux_func_params = (void *)dbl_params; /* two constant params */
	if (verbose) {
	  fprintf(asi_ell_params.outFilePtr,"Trying Power-Law-Expnential\n");
	}
	break;
      }

      /* store size of q */
      asi_ell_params.Nq = q->size;

      /* now do the optimization */
      optfun.func = &asi_ell_combine;
      optfun.params = (void *)(&asi_ell_params);
      optimize(q,&optfun, minimizer_flag, MaxIter,asi_ell_params.outFilePtr);
      if (verbose) {
	/* report fit coefficients */
	fprintf(asi_ell_params.outFilePtr,"Fit results, q:\n");
	gsl_vector_fprintf(asi_ell_params.outFilePtr,q,"%lf ");
	fprintf(asi_ell_params.outFilePtr,"\n");
      }

      /* store q, flux & errors & covariance */
      qs[i] = q;
      /* allocate space for gradient, hess, and q error covariance */
      grad = gsl_vector_alloc(q->size); 
      hess = gsl_matrix_alloc(q->size,q->size);
      covq = gsl_matrix_alloc(q->size,q->size); /* q error covariance = inv(hess) */
      ells[i] = asi_ell_combine(q,(void*)&asi_ell_params,grad,hess); /* get ell, grad, and hess */
      if (support_data) {
      	/* store fit parameters in support data */
      	support_data[ASI_SD_START(i,NY)] = ells[i];
      	for (j=0; j < q->size; j++) {	
      	  support_data[ASI_SD_START(i,NY)+2+j] = gsl_vector_get(q,j);
      	}
      }
      if (support_data || lambda) {
	/* populate estimated counts */
	make_lambda(q,&asi_ell_params,lambdak[i]);
	if (support_data) {
	  for (j=0; j < NY; j++) {
	    support_data[ASI_SD_START(i,NY)+2+ASI_MAX_NQ+j] = gsl_vector_get(lambdak[i],j);
	  }
	}
      }
      min_ell = gsl_min(ells[i],min_ell); /* update min_ell */
      inv_matrix_once(hess,covq,1);

      gsl_matrix_free(hess); /* don't need this any more */

      /* prepare space for flux and sigma (dlogflux) */
      fluxes[i] = gsl_vector_alloc(NEout);
      sigmas[i] = gsl_vector_alloc(NEout);

      /* will grad_matrix because we're computing a whole spectrum with flux_multi */
      grad_matrix = gsl_matrix_alloc(NEout,q->size);
      /* compute flux spectrum, store gradient wrt q in grad_matrix */
      flux_multi(asi_ell_params.flux_func,q,&(Eoutgsl.vector),asi_ell_params.flux_func_params,fluxes[i],grad_matrix);

      /* loop through j and calculate sigmas[i][j] */
      for (j=1; j <= NEout; j++) {
	gsl_vector_set_zero(grad);
	/* sigma_flux_i = sqrt(flux_grad_p*covq*flux_grad_p'); */
	grad_matrix_row = gsl_matrix_row(grad_matrix,j-1);
	gsl_blas_dgemv(CblasNoTrans,1.0,covq,&(grad_matrix_row.vector),0.0,grad); /* grad = covq*row */
	gsl_blas_ddot(grad,&(grad_matrix_row.vector),&sigma_flux_j);
	sigma_flux_j = sqrt(sigma_flux_j);
	gsl_vector_set(sigmas[i],j-1,sigma_flux_j);
      }

      /* free memory */
      gsl_matrix_free(covq);
      gsl_vector_free(grad);
      gsl_matrix_free(grad_matrix);

    }
  }

  /* compute weights for WEXV combo method
    % subtract min(elli) from elli to avoid exp of very large numbers
  */

  /* compute weights and weight sum for each requested analtyical function */
  weight_sum = 0;
  for (i=0; i <= ASI_MAX_POW2; i++) {
    if (qs[i]) { /* qs[i]=NULL if this function not used */
      /* w = exp(-(elli-min(elli))-Ntheta)'; */
      weights[i] = exp(-(ells[i]-min_ell)-(qs[i]->size));
      weight_sum += weights[i];
    } else {
      weights[i] = 0;
    }
  }

  if (lambda) {
    /* create view of lambda and clear it */
    lambda_view = gsl_vector_view_array(lambda,NY);
    gsl_vector_set_zero(&(lambda_view.vector));
  }


  /* the following code is written to reduce memory alloc/dealloc,
   so array/vector variables get reused */
  /* will acumulate over i "in place" in flux and dlogflux */
  /* for now, fluxhat is actually logflux, sigmahat is x2 */
  gsl_vector_set_zero(fluxhat); 
  gsl_vector_set_zero(sigmahat);
  for (i=0; i <= ASI_MAX_POW2; i++) {
    if (weights[i]>0) {
      weights[i] /= weight_sum; /* normalize */
      if (support_data) {
      	support_data[ASI_SD_START(i,NY)+1] = weights[i];
      }
      for (j=1; j <= NEout; j++) {
	flux_j = gsl_vector_get(fluxes[i],j-1);
	logflux_j = log(flux_j);
	sigma2_log_flux_j = gsl_pow_2(gsl_vector_get(sigmas[i],j-1)/flux_j); /* sigma squared */
	/*
                % weighted mean estimate:
                logflux_hat(iE) = sum(w.*logflux);
                x2 = sum(w.*(logflux.^2+sigma2_log_flux));
                sigma_logflux_hat(iE) = sqrt(x2-logflux_hat(iE)^2);
	*/
	tmp_dbl = gsl_vector_get(fluxhat,j-1);
	gsl_vector_set(fluxhat,j-1,tmp_dbl+weights[i]*logflux_j);
	tmp_dbl = gsl_vector_get(sigmahat,j-1); /* x2 for now */
	gsl_vector_set(sigmahat,j-1,tmp_dbl + weights[i]*(gsl_pow_2(logflux_j)+sigma2_log_flux_j));
      }
      if (lambda) {
      	gsl_blas_daxpy(weights[i],lambdak[i],&(lambda_view.vector)); /* lambda += weights(i)*lambdak(i) */
      }
    }
  }

  /* now logflux -> flux, x2 -> sigma_logflux_hat */
  for (j=1; j <= NEout; j++) {
    logflux_j = gsl_vector_get(fluxhat,j-1);
    /* logflux -> flux */
    gsl_vector_set(fluxhat,j-1,exp(logflux_j));
    /* compute sigma_logflux_hat from x2, logflux_j */
    tmp_dbl = gsl_vector_get(sigmahat,j-1); /* x2 */
    tmp_dbl -= gsl_pow_2(logflux_j); /* x2 - logflux_hat(iE)^2 */
    gsl_vector_set(sigmahat,j-1,sqrt(tmp_dbl)); /* sqrt(x2-logflux_hat(iE)^2) */
  }
  
  /* free memory */
  for (i=0; i <= ASI_MAX_POW2; i++) {
    if (lambda || support_data) {
      gsl_vector_free(lambdak[i]);
    }
    if (qs[i]) {
      gsl_vector_free(qs[i]);
      gsl_vector_free(fluxes[i]);
      gsl_vector_free(sigmas[i]);
    }
  }
  free(pen_funcs);
  if (Htmp) {
    gsl_matrix_free(Htmp);
  }

  /* done, success ! */
  return(clean_return(INVLIB_SUCCESS));
}

int ana_spec_inv_multi(const long int Ntimes,
		       const double *y, const double *dy, 
		       const double *Egrid, const double *H0, 
		       const double *dt, const double *b,
		       const long int *int_params, const double *real_params,
		       char *outFile, double *Eout, 
		       double *flux, double *dlogflux, 
		       double *lambda, double *support_data, int *result_codes) {
  long int i,j,t,NE,NY;
  double *H,*sdptr,*lambdaptr;
  int result_code = INVLIB_SUCCESS;

  /* test for input NULLs */
  if ((!y) || (!dy) || (!Egrid) || (!H0) || (!dt) || (!b) || (!int_params) || (!flux) || (!dlogflux) || (!result_codes)) {
    return(INVLIB_ERR_NULL);
  }

  /* start initialization */
  NY = int_params[0];
  NE = int_params[1];
  if ((Ntimes<1) || (NY<=1) || (NE<=1)) {
    return(INVLIB_ERR_DATAEMPTY);
  }

  H = malloc(NE*NY*sizeof(double));
  for (t=0; t < Ntimes; t++) {
    /* set H = H0*dt */
    for (i = 0; i < NY; i++) {
      for (j=0; j < NE; j++) {
	H[NY*j+i] = H0[NY*j+i]*dt[t];
      }
    }
    if (lambda) {
      lambdaptr = lambda+NY*t;
    } else {
      lambdaptr = NULL;
    }

    if (support_data) {
      sdptr = support_data+ASI_SD_START(1,NY)*(ASI_MAX_POW2+1)*t;
    } else {
      sdptr = NULL;
    }

    result_codes[t] = ana_spec_inv(y+NY*t,dy,Egrid,H,b+NY*t,int_params,real_params,
				   outFile,Eout,flux+NE*t,dlogflux+NE*t,
				   lambdaptr,sdptr);
    if (result_codes[t] != INVLIB_SUCCESS) {
      result_code = result_codes[t]; /* store first error code, to be returned later */
    }
  }
  free(H);
  return(result_code);
}
