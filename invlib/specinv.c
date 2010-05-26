#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "invlib.h" /* defines ana_spec_inv, pc_spec_inv */
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
  int type; /* INV_TYPE_ASI or INV_TYPE_PC */

  long int NY,NE,Nq; /* size of vectors */
  long int verbose; /* verbose setting */
  const double *y,*dy; /* counts, relative error for each channel */
  pen_funcTy *pen_funcs; /* penalty function for each channel */
  const gsl_matrix *H; /* measurement matrix */
  const gsl_vector *Egrid, *b; /* energy grid, background counts*/
  void *flux_func_params; /* parameters of flux func */
  FILE *outFilePtr; /* output file pointer */
  int closeOutFile; /* does output file pointer need to be closed on exit? */

  /* used only in analytical spectral inversion (ana_spec_inv) */
  fluxFuncTy flux_func; /* analytical flux function being used */

  /* used in principal component inversion (pc_spec_inv) */
  const gsl_vector *mean_log_flux, *basis_variance;
  const gsl_matrix *basis_vectors;
} ell_paramsTy;

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

double pc_pen(const gsl_vector *q, const void *params,gsl_vector *grad) {
  /* principal component amplitude penalty function */
  /* hess is zero */
  /* negative log likelihood of principal components penalty function */
  /* pc's are assumed to be gaussian with zero mean and variance given by basis_variance */
  double pen=0;
  ell_paramsTy *ell_params = (ell_paramsTy *)params;
  long int i;
  for (i=0; i < q->size; i++) {
    pen += gsl_pow_2(gsl_vector_get(q,i))/gsl_vector_get(ell_params->basis_variance,i); /* q^2/variance */
  }
  pen /= 2.0; /* Gaussian has 1/2 in exponent */
  if (grad) {
    /* d ell / dq_k = q_k / variance_k */
    gsl_vector_memcpy(grad,q);
    gsl_vector_div(grad,ell_params->basis_variance);
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


void flux_pc(const gsl_vector *q, 
	       const void *params,
	       gsl_vector *flux,
	       gsl_matrix *grad) {
  /* principal component spectrum */
  /* don't need to compute hessian, as it's an analytical function
     of flux and grad: d^2flux_i / dq_k dq_l = flux_i*A_{ik}*A {il} 
     where A is the matrix of basis vectors
  */
  ell_paramsTy *ell_params = (ell_paramsTy *)params;
  long int i,j;
  gsl_vector_view tmp_vector_view;

  /* flux starts out as log flux, and we'll exponentiate at the end */
  gsl_vector_memcpy(flux,ell_params->mean_log_flux); /* flux = mean_log_flux */

  /*  y = \alpha A x + \beta y , where A*x = basis_vectors*q , and  alpha=beta=1 */
  gsl_blas_dgemv(CblasNoTrans, 1.0, ell_params->basis_vectors, q, 1.0, flux);

  /* exponentiate log flux into flux */
  vector_func(flux,flux,&exp);
  if (grad) {
    /* dflux_i / dq_k = flux_i * A_{ik} */
    for (i=0; i < grad->size1; i++) {
      for (j=0; j < grad->size2; j++) {
	gsl_matrix_set(grad,i,j,gsl_vector_get(flux,i)*gsl_matrix_get(ell_params->basis_vectors,i,j));
      }
    }
  }
}


void make_lambda(const gsl_vector *q,const ell_paramsTy *ell_params,gsl_vector *lambda) {
  /* sets lambda based on q and info in ell_params: H, b, flux_func, etc */
  gsl_vector *flux = gsl_vector_alloc(ell_params->NE);
  if (ell_params->type == INV_TYPE_ASI) {
    flux_multi(ell_params->flux_func,q,ell_params->Egrid,ell_params->flux_func_params,flux,NULL);
  } else {
    flux_pc(q,ell_params,flux,NULL);
  }
  /* copy background into lambda: lambda = b */
  gsl_vector_memcpy(lambda,ell_params->b);  
  /* lambda = 1.0*H*flux+1.0*b (since lambda=b) */
  gsl_blas_dgemv(CblasNoTrans,1.0,ell_params->H,flux,1.0,lambda);
  gsl_vector_free(flux);
}

/****************************
  ell_combine - the master negative-log-likelihood function
  inputs:
     q - parameters/coefficients of spectral function
     params_void - parameters, an ell_params type, holds all kinds of info for optimization
     grad: gradient with respect to q's
     hess: Hessian (second derivative matrix) with respect to q's
 Note:
    if grad is NULL, neither grad nor hess is computed
    if hess is NULL, hess is not computed
 ****************************/

double ell_combine(const gsl_vector *q, void *params_void, gsl_vector *grad, gsl_matrix *hess) {
  double ell=0; /* neglogp */
  double E,h;
  long int i,j,k,l;
  double grad_ell_lambdak,hess_ell_lambdak;
  double *grad_ell_lambdak_ptr=0,*hess_ell_lambdak_ptr=0;
  gsl_vector *lambda=0,*flux=0,*grad_lambda_q=0,*grad_lambda_qj=0;
  gsl_matrix *hess_lambda_q=0,*hess_lambda_qj=0,*flux_grad=0;
  ell_paramsTy *ell_params = (ell_paramsTy *)params_void;

  if (grad) { /* prepare space for gradient computation */
    gsl_vector_set_zero(grad);
    if (ell_params->type == INV_TYPE_ASI) {
      flux_grad = gsl_matrix_alloc(ell_params->NE,q->size);
    }
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

  flux = gsl_vector_alloc(ell_params->NE);
  if (ell_params->type == INV_TYPE_ASI) {
    flux_multi(ell_params->flux_func,q,ell_params->Egrid,ell_params->flux_func_params,flux,NULL);
  } else {
    flux_pc(q,ell_params,flux,flux_grad); /* get grad if requested, otherwise just get flux */
  }

  lambda = gsl_vector_alloc(ell_params->NY);
  /* copy background into lambda: lambda = b */
  gsl_vector_memcpy(lambda,ell_params->b);  
  /* lambda = 1.0*H*flux+1.0*b (since lambda=b) */
  gsl_blas_dgemv(CblasNoTrans,1.0,ell_params->H,flux,1.0,lambda);

  /* have verified that lambda agrees with matlab result */

  for (i=1; i <= ell_params->NY; i++) {
    /* get the ell increment and gradient/hessian as needed */
    /* grad_ell_lambdak, hess_ell_lambdak are scalars */
    ell += (ell_params->pen_funcs[i-1])(gsl_vector_get(lambda,i-1),ell_params->y[i-1],
				    &(ell_params->dy[i-1]), grad_ell_lambdak_ptr, hess_ell_lambdak_ptr);
    if (grad) { /* compute gradient */
      if (hess) {
	/* hess needs to accumulate grad and hess over j, so clear these variables to be accumulators */
	gsl_vector_set_zero(grad_lambda_q); 
	gsl_matrix_set_zero(hess_lambda_q);
      }
      for (j=1; j<=ell_params->NE; j++) {
	h = gsl_matrix_get(ell_params->H,i-1,j-1); /*H_ij*/
	/* now, evaulate flux function at E(j) and get grad/hess as needed */
	if (ell_params->type == INV_TYPE_ASI) {
	  E = gsl_vector_get(ell_params->Egrid,j-1); /* E_j */
	  ell_params->flux_func(q,E,ell_params->flux_func_params,grad_lambda_qj, hess_lambda_qj);
	} else {
	  gsl_matrix_get_row(grad_lambda_qj, flux_grad, j); /* dflux(j)/dq for all q */
	  if (hess) {
	    /* d^2flux_j / dq_k dq_l = flux_j*A_{jk}*A_{jl} = dflux_j/dq_k * A_{jl} */
	    for (k=1; k <= q->size; k++) {
	      for (l=1; l <= q->size; l++) {
		gsl_matrix_set(hess_lambda_qj,k,l,gsl_matrix_get(flux_grad,j-1,k-1)*gsl_matrix_get(ell_params->basis_vectors,j-1,l-1));
	      }
	    }
	  }
	}
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

  /* now, add prior penalty for PC inversion */
  if (ell_params->type == INV_TYPE_PC) {
    /* re-use grad_lambda_qj for pen_grad (_qj b/c _q is only initialized for hess) */
    ell+=pc_pen(q,ell_params,grad_lambda_qj);
    if (grad) {
      gsl_vector_add(grad,grad_lambda_qj);
      /* pc_pen contribution to hess is zero */
    }
  };

  /* free memory */
  gsl_vector_free(lambda);
  gsl_vector_free(flux);
  if (grad) {
    gsl_vector_free(grad_lambda_qj);
    if (flux_grad) {
      gsl_matrix_free(flux_grad);
    }
    if (hess) {
      gsl_vector_free(grad_lambda_q);
      gsl_matrix_free(hess_lambda_q);
      gsl_matrix_free(hess_lambda_qj);
    }
  }
  /* have verified that ell, grad, and hess are right, compared to matlab */
  return(ell);
}

int check_common_inputs(const long int *int_params, 
			long int *NYptr, long int *NEptr,
			long int *minimizer_flagptr,
			long int *MaxIterptr,
			long int *verboseptr,
			long int *dE_modeptr,
			const double *y, const double *dy, const double *Egrid, 
			const double *H, const double *b, 
			const double *flux, const double *dlogflux,
			ell_paramsTy *ell_params,const char *outFile) {
  long int i,j; /* loop control vars */
  int any_y_positive = 0, Ny_valid=0;
  long int NY,NE,minimizer_flag,MaxIter,verbose,dE_mode;

  /* test for input NULLs */
  if ((!y) || (!dy) || (!Egrid) || (!H) || (!b) || (!int_params) || (!flux) || (!dlogflux)) {
    return(INVLIB_ERR_NULL);
  }
  NY = int_params[0];
  NE = int_params[1];
  minimizer_flag = int_params[4];
  MaxIter = int_params[5];
  verbose = int_params[6];
  dE_mode = int_params[7];
  *NYptr = NY;
  *NEptr = NE;
  *minimizer_flagptr = minimizer_flag;
  *MaxIterptr = MaxIter;
  *verboseptr = verbose;
  *dE_modeptr = dE_mode;
  /* check for not enough channels */
  if ((NY<=1) || (NE<=1)) {
    return(INVLIB_ERR_DATAEMPTY);
  }

  /* check for NaNs, Infs, and inappropriate negative numbers */
  for (i=0; i < NY; i++) {
    if (gsl_finite(y[i])) {
      Ny_valid++;
      if ((! gsl_finite(dy[i])) || (! gsl_finite(b[i]))) {
	return(INVLIB_ERR_DATANAN);
	
      }
      if ((y[i]<0) || (dy[i]<0) || (b[i]<0)) {
	return(INVLIB_ERR_DATANAN);
      }
      if (y[i]>b[i]) { /* also check for at least one count above background in some channel */
	any_y_positive = 1;
      }
    }
      
    for (j=0; j < NE; j++) {
      if (! gsl_finite(H[NY*j+i])) {
	return(INVLIB_ERR_DATANAN);
      }
    }
  }

  /* check for not enough valid channels */
  if (Ny_valid<=1) {
    return(INVLIB_ERR_DATAEMPTY);
  }


  /* confirm at least one count in some channel */
  if (! any_y_positive) {
    return(INVLIB_ERR_NOCOUNTS);
  }

  /* check for NaNs and Infs or <= in Egrid */
  for (j=0; j < NE; j++) {
    if ((!gsl_finite(Egrid[j])) || (Egrid[j]<=0)) {
      return(INVLIB_ERR_DATANAN);
    }
  }

  /* check for valid minimizer flag */
  if ((minimizer_flag < 0) || (minimizer_flag > OPTIM_MIN_MAX)) {
    return(INVLIB_ERR_INVALIDMIN);
  }

  /* check for valid maximum interations */
  if (MaxIter<=0) {
    return(INVLIB_ERR_INVALIDITER);
  }

  /* check verbose setting */
  ell_params->verbose = verbose;
  switch(verbose) { /* this switch statement must be consistent with specinv.h */
  case 0:
    ell_params->outFilePtr = NULL;
    break;
  case 1:
    ell_params->outFilePtr = stdout;
    break;
  case 2:
    ell_params->outFilePtr = stderr;
    break;
  case 3:
    if (outFile) {
      ell_params->outFilePtr = (FILE *)outFile;
    } else {
      return(INVLIB_ERR_VERBOSE); /* can't output to NULL */
    }
    break;
  case 4:
    if (outFile) {
      if (! (ell_params->outFilePtr = fopen(outFile,"w"))) {
	return(INVLIB_ERR_OUTFILE); /* Couldn't open outFile */
	  }
      ell_params->closeOutFile = 1;
    } else {
      return(INVLIB_ERR_VERBOSE); /* can't output to NULL */
    }
    break;
  case 5:
    if (outFile) {
      if (! (ell_params->outFilePtr = fopen(outFile,"a"))) {
	return(INVLIB_ERR_OUTFILE); /* Couldn't open outFile */
	  }
      ell_params->closeOutFile = 1;
    } else {
      return(INVLIB_ERR_VERBOSE); /* can't output to NULL */
    }
    break;
  default:
    return(INVLIB_ERR_VERBOSE); /* invalid verbose value */
  }
  return(INVLIB_SUCCESS);
}

pen_funcTy *setup_pen_funcs(const long int NY, const long int NE, const long int dE_mode,
			    const double *y, const double *dy, 
			    const gsl_matrix *Hgsl, const gsl_vector *Egridgsl, const gsl_vector *bgsl,
			    ell_paramsTy *ell_params,gsl_matrix **Htmpptr) {
  pen_funcTy *pen_funcs=0;
  long int i;

  pen_funcs = malloc(NY*sizeof(pen_funcTy));

  if (dE_mode != SPECINV_DE_INCLUDED) {
    gsl_vector *dE_tmp = gsl_vector_alloc(NE);
    for (i=1; i < NE-1; i++) {
      gsl_vector_set(dE_tmp,i,(gsl_vector_get(Egridgsl,i+1)-gsl_vector_get(Egridgsl,i-1))/2);
    }
    gsl_vector_set(dE_tmp,0,(gsl_vector_get(Egridgsl,1)-gsl_vector_get(Egridgsl,0)));
    gsl_vector_set(dE_tmp,NE-1,(gsl_vector_get(Egridgsl,NE-1)-gsl_vector_get(Egridgsl,NE-2)));
    if (dE_mode == SPECINV_DE_TRAPZ) { /* trapz divides the endpoint weights by 2 */
      gsl_vector_set(dE_tmp,0,gsl_vector_get(dE_tmp,0)/2);
      gsl_vector_set(dE_tmp,NE-1,gsl_vector_get(dE_tmp,NE-1)/2);
    }
    /* now create and scale columns of Htmp by dE */
    (*Htmpptr) = gsl_matrix_alloc(NY,NE);
    gsl_matrix_memcpy((*Htmpptr),Hgsl);
    for (i=0; i < NE; i++) {
      gsl_vector_view Hcol = gsl_matrix_column((*Htmpptr),i);
      gsl_vector_scale(&(Hcol.vector),gsl_vector_get(dE_tmp,i));
    }
    gsl_vector_free(dE_tmp);
  }

  /* determine appropriate penalty function */
  for (i=0; i < NY; i++) {
    if (! gsl_finite(y[i])) {
      if (ell_params->verbose) {
	fprintf(ell_params->outFilePtr,"Selecting NaN penalty function for y[%li]=%g\n",i,y[i]);
      }
      pen_funcs[i] = &nan_pen; /* dummy for non-finite y */
    } else if ((dy[i]==0) || (y[i] < 1.0/gsl_pow_2(dy[i]))) {
      if (ell_params->verbose) {
	fprintf(ell_params->outFilePtr,"Selecting Poisson penalty function for y[%li]=%g\n",i,y[i]);
      }
      pen_funcs[i] = &poiss_pen;
    } else {
      if (ell_params->verbose) {
	fprintf(ell_params->outFilePtr,"Selecting Gaussian Relative Error penalty function for y[%li]=%g\n",i,y[i]);
      }
      pen_funcs[i] = &rel_pen;
    }
  }
  ell_params->NY = NY;
  ell_params->NE = NE;
  ell_params->pen_funcs = pen_funcs;
  if (*Htmpptr) {
    ell_params->H = *Htmpptr;
  } else {
    ell_params->H = Hgsl;
  }
  ell_params->b = bgsl;
  ell_params->y = y;
  ell_params->dy = dy;
  ell_params->Egrid = Egridgsl;

  return(pen_funcs);
}

int ana_spec_inv(const double *y, const double *dy, const double *Egrid, const double *H, const double *b,
		 const long int *int_params, const double *real_params,
		 char *outFile, double *Eout, double *flux, double *dlogflux, 
		 double *lambda, double *support_data) {
/****
     see header file for details of inputs/outputs
****/
  int result;
  long int NY, NE, NEout, fxn_bit_map, minimizer_flag, MaxIter, verbose, dE_mode;
  long int fxn_bit;
  long int i,j; /* loop control vars */
  double m0c2=0,E0=0,E_break=0;
  double dbl_params[10]; /* holder for parameters to flux functions */
  gsl_vector_view fluxgsl,dlogfluxgsl,grad_matrix_row,lambda_view; 
  gsl_vector *fluxhat=0, *sigmahat=0;
  gsl_vector *q=0; /* parameter vector */
  gsl_vector *grad=0;
  gsl_matrix *hess=0, *covq=0, *grad_matrix=0;
  pen_funcTy *pen_funcs=0;
  ell_paramsTy ell_params;
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
    if (ell_params.closeOutFile) {
      fclose(ell_params.outFilePtr);
    }
    return(return_value);
  };

  ell_params.closeOutFile = 0; /* default: no need to close outFile */
  ell_params.type = INV_TYPE_ASI; /* set ana-spec-inv as the type */

  /* check/initialize common inputs/vars */
  if ((result=check_common_inputs(int_params,
				  &NY,&NE,&minimizer_flag,&MaxIter,&verbose,&dE_mode,
				  y,dy,Egrid,H,b,flux,dlogflux,&ell_params,outFile
				  )) != INVLIB_SUCCESS) {
    return(clean_return(result));
  }

  /* check/initialize asi specific inputs/vars */
  NEout = int_params[2];

  /* check for not enough channels */
  if (NEout < 1) {
    return(clean_return(INVLIB_ERR_DATAEMPTY));
  }

  /* check for NaNs and Infs or <=0 in Eout */
  for (j=0; j < NEout; j++) {
    if (!gsl_finite(Eout[j])) {
      return(clean_return(INVLIB_ERR_DATANAN));
    }
  }

  /* initialize some more parameters */
  fxn_bit_map = int_params[3];

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

  if (lambda || support_data) { /* prepare to compute lambda */
    for (i=0; i <= ASI_MAX_POW2; i++) {
      lambdak[i] = gsl_vector_alloc(NY);
    }
  }
  
  /* must declare these here because otherwise compiler treats assignments as violating const-ness */
  gsl_matrix_const_view Hgsl = gsl_matrix_const_view_array(H,NY,NE); /* gsl view of H */
  gsl_vector_const_view Egridgsl = gsl_vector_const_view_array(Egrid,NE);
  gsl_vector_const_view bgsl = gsl_vector_const_view_array(b,NY);
  gsl_vector_const_view Eoutgsl = gsl_vector_const_view_array(Eout,NEout);

  pen_funcs = setup_pen_funcs(NY,NE,dE_mode,y,dy,&(Hgsl.matrix),&(Egridgsl.vector),&(bgsl.vector),&ell_params,&Htmp);
  
  /* initialize other gsl views */
  fluxgsl = gsl_vector_view_array(flux,NEout);
  fluxhat = &(fluxgsl.vector);
  dlogfluxgsl = gsl_vector_view_array(dlogflux,NEout);
  sigmahat = &(dlogfluxgsl.vector);

  min_ell = GSL_POSINF; /* initialize to biggest number */

  /* prepare parameters for ell_combine */
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
	ell_params.flux_func = flux_pl; /* power law flux spectrum */
	ell_params.flux_func_params = (void *)NULL; /* no constant params */
	if (verbose) {
	  fprintf(ell_params.outFilePtr,"Trying Power-Law\n");
	}
	break;
      case ASI_FXN_EXP:
	/* exponential */
	/* allocate, this'll get stored in qs[i] and freed at the end of the function */
	q = gsl_vector_alloc(2); /* two free parameters */
	/* initialize with default params */
	gsl_vector_set(q,0,1);
	gsl_vector_set(q,1,-1);
	ell_params.flux_func = flux_exp; /* exponential flux spectrum */
	ell_params.flux_func_params = (void *)NULL; /* no constant params */
	if (verbose) {
	  fprintf(ell_params.outFilePtr,"Trying Exponential\n");
	}
	break;
      case ASI_FXN_RM:
	/* relativistic maxwellian */
	/* allocate, this'll get stored in qs[i] and freed at the end of the function */
	q = gsl_vector_alloc(2); /* two free parameters */
	/* initialize with default params */
	gsl_vector_set(q,0,1);
	gsl_vector_set(q,1,-1);
	ell_params.flux_func = flux_rm; /* relativistic Maxwellian */
	ell_params.flux_func_params = (void *)(&m0c2); /* one constant param, E0=m0c2 */
	if (verbose) {
	  fprintf(ell_params.outFilePtr,"Trying Relativistic Maxwellian\n");
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
	ell_params.flux_func = flux_rm2; /* relativistic Maxwellian */
	ell_params.flux_func_params = (void *)(&m0c2); /* one constant param, E0=m0c2 */
	if (verbose) {
	  fprintf(ell_params.outFilePtr,"Trying Double Relativistic Maxwellian\n");
	}
	break;
      case ASI_FXN_PLE:
	/* power law w/ exponential tail*/
	/* allocate, this'll get stored in qs[i] and freed at the end of the function */
	q = gsl_vector_alloc(2); /* two free parameters */
	/* initialize with default params */
	gsl_vector_set(q,0,1);
	gsl_vector_set(q,1,3);
	ell_params.flux_func = flux_ple; /* power law flux spectrum w/ exp tail */
	dbl_params[0] = E_break;
	dbl_params[1] = E0;
	ell_params.flux_func_params = (void *)dbl_params; /* two constant params */
	if (verbose) {
	  fprintf(ell_params.outFilePtr,"Trying Power-Law-Expnential\n");
	}
	break;
      }

      /* store size of q */
      ell_params.Nq = q->size;

      /* now do the optimization */
      optfun.func = &ell_combine;
      optfun.params = (void *)(&ell_params);
      optimize(q,&optfun, minimizer_flag, MaxIter,ell_params.outFilePtr);
      if (verbose) {
	/* report fit coefficients */
	fprintf(ell_params.outFilePtr,"Fit results, q:\n");
	gsl_vector_fprintf(ell_params.outFilePtr,q,"%lf ");
	fprintf(ell_params.outFilePtr,"\n");
      }

      /* store q, flux & errors & covariance */
      qs[i] = q;
      /* allocate space for gradient, hess, and q error covariance */
      grad = gsl_vector_alloc(q->size); 
      hess = gsl_matrix_alloc(q->size,q->size);
      covq = gsl_matrix_alloc(q->size,q->size); /* q error covariance = inv(hess) */
      ells[i] = ell_combine(q,(void*)&ell_params,grad,hess); /* get ell, grad, and hess */
      if (support_data) {
      	/* store fit parameters in support data */
      	support_data[ASI_SD_START(i,NY)] = ells[i];
      	for (j=0; j < q->size; j++) {	
      	  support_data[ASI_SD_START(i,NY)+2+j] = gsl_vector_get(q,j);
      	}
      }
      if (support_data || lambda) {
	/* populate estimated counts */
	make_lambda(q,&ell_params,lambdak[i]);
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
      flux_multi(ell_params.flux_func,q,&(Eoutgsl.vector),ell_params.flux_func_params,fluxes[i],grad_matrix);

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



int pc_spec_inv(const double *y, const double *dy, const double *Egrid, const double *H, const double *b,
		const double *mean_log_flux, const double *basis_vectors, const double *basis_variance,
		const long int *int_params, const double *real_params,
		char *outFile, double *flux, double *dlogflux, double *lambda, double *support_data) {

  int result;
  long int NY, NE, minimizer_flag, MaxIter, verbose, dE_mode;
  long int Nbases,Nactive_bases;
  long int i,j;
  ell_paramsTy ell_params;
  pen_funcTy *pen_funcs=0;
  optfunTy optfun; 
  gsl_vector *grad,*q;
  gsl_vector_view fluxgsl;
  gsl_matrix *hess,*covq,*grad_matrix;
  double ell;
  gsl_matrix *Htmp=0; /* holds H*dE if necessary */
  /* the following "const" view variables will be declared inline below
  gsl_matrix_const_view Hgsl = gsl_matrix_const_view_array(H,NY,NE); gsl view of H 
  gsl_vector_const_view Egridgsl = gsl_vector_const_view_array(Egrid,NE); gsl view of Egrid
  gsl_vector_const_view bgsl = gsl_vector_const_view_array(b,NY); gsl view of b
  */


  int clean_return(int return_value) {
    if (ell_params.closeOutFile) {
      fclose(ell_params.outFilePtr);
    }
    return(return_value);
  };

  ell_params.closeOutFile = 0; /* default: no need to close outFile */
  ell_params.type = INV_TYPE_PC; /* set principal components as the type */

  /* check/initialize common inputs/vars */
  if ((result=check_common_inputs(int_params,
				  &NY,&NE,&minimizer_flag,&MaxIter,&verbose,&dE_mode,
				  y,dy,Egrid,H,b,flux,dlogflux,&ell_params,outFile
				  )) != INVLIB_SUCCESS) {
    return(clean_return(result));
  }

  /* notes: Eout = Egrid */
  Nbases = int_params[2];
  Nactive_bases = int_params[3];

  if ((Nbases < 1) || (Nactive_bases < 1)) {
    return(clean_return(INVLIB_ERR_DATAEMPTY));
  }

  /* check for NaNs and Infs or <=0 in PC definitions */
  for (j=0; j < Nactive_bases; j++) {
    if ((!gsl_finite(basis_variance[j])) || (basis_variance[j]<=0)) {
	return(clean_return(INVLIB_ERR_DATANAN));
    }
    for (i=0; i < NE; i++) {
      if (!gsl_finite(basis_vectors[NE*j+i])) {
	return(clean_return(INVLIB_ERR_DATANAN));
      }
    }
  }
  for (i=0; i < NE; i++) {
    if (!gsl_finite(mean_log_flux[i])) {
	return(clean_return(INVLIB_ERR_DATANAN));
    }
  }

  /* must declare these here because otherwise compiler treats assignments as violating const-ness */
  gsl_matrix_const_view Hgsl = gsl_matrix_const_view_array(H,NY,NE); /* gsl view of H */
  gsl_vector_const_view Egridgsl = gsl_vector_const_view_array(Egrid,NE);
  gsl_vector_const_view bgsl = gsl_vector_const_view_array(b,NY);

  pen_funcs = setup_pen_funcs(NY,NE,dE_mode,y,dy,&(Hgsl.matrix),&(Egridgsl.vector),&(bgsl.vector),&ell_params,&Htmp);

  /* will have to do some fancy gsl view manipulation to get from basis_vectors to ell_params.basis_vectors, etc */

  gsl_vector_const_view mean_log_flux_view = gsl_vector_const_view_array(mean_log_flux,NE); /* view of original array */
  gsl_matrix_const_view basis_vectors_full_view = gsl_matrix_const_view_array(basis_vectors,NE,Nbases); /* view of original array */
  gsl_vector_const_view basis_variance_full_view = gsl_vector_const_view_array(basis_variance,Nbases); /* view of original array */
  gsl_matrix_const_view basis_vectors_active_view = gsl_matrix_const_submatrix(&(basis_vectors_full_view.matrix),0,0,NE,Nactive_bases); /* view of active columns of matrix */
  gsl_vector_const_view basis_variance_active_view = gsl_vector_const_subvector(&(basis_variance_full_view.vector),0,Nactive_bases); /* view of active part of array */


  q = gsl_vector_calloc(Nactive_bases); /* initialize q to 0 */
  /* store size of q */
  ell_params.Nq = q->size;
  ell_params.mean_log_flux = &(mean_log_flux_view.vector);
  ell_params.basis_vectors = &(basis_vectors_active_view.matrix);
  ell_params.basis_variance = &(basis_variance_active_view.vector);

  /* now do the optimization */
  optfun.func = &ell_combine;
  optfun.params = (void *)(&ell_params);
  optimize(q,&optfun, minimizer_flag, MaxIter,ell_params.outFilePtr);
  if (verbose) {
    /* report fit coefficients */
    fprintf(ell_params.outFilePtr,"Fit results, q:\n");
    gsl_vector_fprintf(ell_params.outFilePtr,q,"%lf ");
    fprintf(ell_params.outFilePtr,"\n");
  }

  /* now do equivalent to code following line 927 */
  /* call ell_combine w/ grad and hess */
  grad = gsl_vector_alloc(q->size); /* dell/dq */
  hess = gsl_matrix_alloc(q->size,q->size); /* d^2ell/dq^2 */
  covq = gsl_matrix_alloc(q->size,q->size); /* q error covariance = inv(hess) */
  ell = ell_combine(q,(void*)&ell_params,grad,hess); /* get ell, grad, and hess */

  /* store lambda in lambda */
  if (lambda) {
    gsl_vector_view lambda_view = gsl_vector_view_array(lambda,NY);
    make_lambda(q,&ell_params,&(lambda_view.vector));
  }
  if (support_data) {
    /* store ell, q in support_data */
  }

  /* compute dlogflux = sqrt(diag(dq/dlogflux * d^2ell/dq^2 * dq/dlogflux)) */
  inv_matrix_once(hess,covq,1);
  gsl_matrix_free(hess);
  /* store flux in flux, get grad_matrix */
  /* will get grad_matrix because we're computing a whole spectrum with flux_pc */
  grad_matrix = gsl_matrix_alloc(NE,q->size);
  fluxgsl = gsl_vector_view_array(flux,NE);
  flux_pc(q,&ell_params,&(fluxgsl.vector),grad_matrix);

  if (support_data) { /* ell, q[0],q[1],... */
    support_data[0] = ell;
    for (j=0; j < q->size; j++) {
      support_data[j+1] = gsl_vector_get(q,j);
    }
  }
  
  /* loop through j and calculate sigmas[i][j] */
  for (j=1; j <= NE; j++) {
    gsl_vector_set_zero(grad);
    /* sigma_flux_i = sqrt(flux_grad_p*covq*flux_grad_p'); */
    gsl_vector_view grad_matrix_row = gsl_matrix_row(grad_matrix,j-1);
    gsl_blas_dgemv(CblasNoTrans,1.0,covq,&(grad_matrix_row.vector),0.0,grad); /* grad = covq*row */
    double sigma_flux_j=0;
    gsl_blas_ddot(grad,&(grad_matrix_row.vector),&sigma_flux_j);
    sigma_flux_j = sqrt(sigma_flux_j);
    dlogflux[j-1] = sigma_flux_j;
  }

  free(pen_funcs);
  gsl_vector_free(q);
  gsl_vector_free(grad);
  gsl_matrix_free(covq);
  gsl_matrix_free(grad_matrix);
}

