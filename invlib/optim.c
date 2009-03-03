/********************************************************************************************************
OPTIMIZE routines - this section defines optimize and the routines
it needs
********************************************************************************************************/


#include <stdlib.h>
#include <stdio.h>
#include <gsl/gsl_multimin.h>
#include <gsl/gsl_multiroots.h>
#include <gsl/gsl_errno.h> /* strerror */
#include <gsl/gsl_blas.h> /* nrm2 */
#include "optim.h" /* includes gsl vector and matrix headers */

/*************************************************************
Wrapper functions used by optimize and gsl multimin routines
*************************************************************/

double optim_f( const gsl_vector *q, void *params) {
  /* just evaluate ell_combine at q */
  optfunTy *optfun = params;
  double neglogp = optfun->func(q,optfun->params,NULL,NULL);
  return(neglogp);
}

void optim_df(const gsl_vector *q, void *params, gsl_vector *df) {
  /* just evaluate gradient of ell_combine wrt q */
  optfunTy *optfun = params;
  optfun->func(q,optfun->params,df,NULL);
}

void optim_fdf( const gsl_vector *q, void *params, double *f, gsl_vector *df ) {
  /* evaluate ell_combine and get gradient  wrt q */
  optfunTy *optfun = params;
  *f = optfun->func(q,optfun->params,df,NULL);
}

double optimize(gsl_vector *q, optfunTy *optfun, const long int minimizer_flag, const long int MaxIter, FILE *fid) {
  /* see optim.h for helps */
  const gsl_multimin_fdfminimizer_type *Tfdf;
  const gsl_multimin_fminimizer_type *Tf;
  gsl_multimin_fminimizer *s_f;
  gsl_multimin_fdfminimizer *s_fdf;
  gsl_vector *ss=0; /* step size for simplex */
  gsl_multimin_function my_f;
  gsl_multimin_function_fdf my_fdf;
  int status;
  double size;
  size_t iter = 0;
  double fval;
  unsigned long int soOften = 1;  /* this sets the print status frequency initially */

  /* Nelder Mead Simples (minimizer_flag == OPTIM_MIN_NM) is
     a special case because it does not require the gradient,
     but does, instead, require an initial guess at the
     step size.
   */

  /* select the minimizer, either Tfdf or Tf */
  switch(minimizer_flag) {
  case OPTIM_MIN_FR:
    Tfdf = gsl_multimin_fdfminimizer_conjugate_fr;
    break;
  case OPTIM_MIN_PR:
    Tfdf = gsl_multimin_fdfminimizer_conjugate_pr;
    break;
  case OPTIM_MIN_NM:
    Tf = gsl_multimin_fminimizer_nmsimplex;
    break;
  case OPTIM_MIN_BFGS:
  default:
    Tfdf = gsl_multimin_fdfminimizer_vector_bfgs2;
    break;
  }

  /* set up the multimin function my_f or my_fdf, and set the minimizer */
  if (OPTIM_MIN_ISFDF(minimizer_flag)) {
    s_fdf = gsl_multimin_fdfminimizer_alloc (Tfdf, q->size);
    my_fdf.n = q->size;
    my_fdf.f = &optim_f;
    my_fdf.df = &optim_df;
    my_fdf.fdf = &optim_fdf;
    my_fdf.params = (void *)optfun;
    gsl_multimin_fdfminimizer_set (s_fdf, &my_fdf, q, 1e-6, 1e-6);
    if (fid) {
      fprintf(fid,"%s invoked with minimizer %s, MaxIter=%li\n",__func__,
	      gsl_multimin_fdfminimizer_name(s_fdf),MaxIter);
    }
  } else {
    ss = gsl_vector_alloc(q->size); /* initial step size for simplex */
    gsl_vector_set_all(ss,1.0);
    s_f = gsl_multimin_fminimizer_alloc (Tf, q->size);
    my_f.n = q->size;
    my_f.f = &optim_f;
    my_f.params = (void *)optfun;
    gsl_multimin_fminimizer_set (s_f, &my_f, q, ss);
    if (fid) {
      fprintf(fid,"%s invoked with minimizer %s, MaxIter=%li\n",__func__,
	      gsl_multimin_fminimizer_name(s_f),MaxIter);
    }
  }

  /* pretty standard minimizer loop from the GSL examples */
      
  do
    {
      iter++;
      /* do one iteration */
      if (OPTIM_MIN_ISFDF(minimizer_flag)) {
	status = gsl_multimin_fdfminimizer_iterate (s_fdf); /* iterate */
	fval = gsl_multimin_fdfminimizer_minimum(s_fdf);  /* get current function value */
      } else {
	status = gsl_multimin_fminimizer_iterate (s_f); /* iterate */
	fval = s_f->fval; /* get current function value */
      }
      
      if (status) { /* maybe we're done, but probably not for a happy reason */
	break;
      }
      
      /* now test stopping criteria */
      if (OPTIM_MIN_ISFDF(minimizer_flag)) {
	status = gsl_multimin_test_gradient (s_fdf->gradient, 1e-8); /* test gradient against stop */
      } else {
	size = gsl_multimin_fminimizer_size(s_f); /* get simplex size */
	status = gsl_multimin_test_size(size,1e-6); /* test simplex size against stop */
      }
      
      if (fid && ((iter % soOften) == 0)) {
	/* update when to print status -- only when first digit of iter changes */
	if (iter >= soOften*10) {
	  soOften = (soOften < 1000) ? soOften*10 : 1000;
	}
	/* print status */
	fprintf(fid,"%s: %lu/%lu: %.14e\n",__func__,(unsigned long int)iter,(unsigned long int)MaxIter,fval);
      }      
    }
  while ((status == GSL_CONTINUE) && (iter < MaxIter)); /* stop on success or run out of iterations */
  
  if (fid) {
    fprintf(fid,"%s: completed after %lu/%lu: %s\n",__func__,(unsigned long int)iter,(unsigned long int)MaxIter,gsl_strerror(status));
  }

  /* copy solution into q and fval */
  if (OPTIM_MIN_ISFDF(minimizer_flag)) {
    gsl_vector_memcpy(q,s_fdf->x);
    fval = gsl_multimin_fdfminimizer_minimum(s_fdf);
  } else {
    gsl_vector_memcpy(q,s_f->x);
    fval = s_f->fval;
  }

  /* free memory */
  if (OPTIM_MIN_ISFDF(minimizer_flag)) {
    gsl_multimin_fdfminimizer_free (s_fdf);
  } else {
    gsl_multimin_fminimizer_free (s_f);
    free(ss);
  }
  return(fval);
}

/*************************************************************
Wrapper functions used fzero and gsl multiroots routines
*************************************************************/

int fzero_f( const gsl_vector *q, void *params, gsl_vector *f) {
  /* evaluate gradient (f) of ell at q */
  optfunTy *optfun = params;
  optfun->func(q,optfun->params,f,NULL);
  return(GSL_SUCCESS);
}

int fzero_df(const gsl_vector *q, void *params, gsl_matrix *df) {
  /* evaluate hessian (df) of ell at q */
  optfunTy *optfun = params;
  gsl_vector *f = gsl_vector_alloc(q->size); /* dummy variable becaues hess requires grad */
  optfun->func(q,optfun->params,f,df);
  gsl_vector_free(f);
  return(GSL_SUCCESS);
}

int fzero_fdf( const gsl_vector *q, void *params, gsl_vector *f, gsl_matrix *df ) {
  /* evaluate gradient (f) and hessian (df) of ell at q */
  optfunTy *optfun = params;
  optfun->func(q,optfun->params,f,df);
  return(GSL_SUCCESS);
}

double fzero(gsl_vector *q, optfunTy *optfun, const long int minimizer_flag, const long int MaxIter, FILE *fid) {
  const gsl_multiroot_fsolver_type *T_f;
  const gsl_multiroot_fdfsolver_type *T_fdf;
  gsl_multiroot_fsolver *s_f;
  gsl_multiroot_fdfsolver *s_fdf;
  int status;
  size_t iter = 0;
  gsl_multiroot_function my_f;
  gsl_multiroot_function_fdf my_fdf;
  unsigned long int soOften = 1;  /* this sets the print status frequency initially */
  double tmp;

  /* f and fdf solvers have sligthly different syntax */

  switch(minimizer_flag) {
  case OPTIM_FZ_HYBRID:
    T_f = gsl_multiroot_fsolver_hybrid;
  case OPTIM_FZ_DNEWTON:
    T_f = gsl_multiroot_fsolver_dnewton;
  case OPTIM_FZ_BROYDEN:
    T_f = gsl_multiroot_fsolver_broyden;
  case OPTIM_FZ_HYBRIDJ:
    T_fdf = gsl_multiroot_fdfsolver_hybridj;
  case OPTIM_FZ_NEWTON:
    T_fdf = gsl_multiroot_fdfsolver_newton;
  case OPTIM_FZ_GNEWTON:
    T_fdf = gsl_multiroot_fdfsolver_gnewton;
  case OPTIM_FZ_HYBRIDSJ:
  default:
    T_fdf = gsl_multiroot_fdfsolver_hybridsj;
  }

  if (OPTIM_FZ_ISFDF(minimizer_flag)) {
    s_fdf = gsl_multiroot_fdfsolver_alloc (T_fdf, q->size);
    my_fdf.n = q->size;
    my_fdf.f = &fzero_f;
    my_fdf.df = &fzero_df;
    my_fdf.fdf = &fzero_fdf;
    my_fdf.params = (void *)optfun;
    gsl_multiroot_fdfsolver_set (s_fdf, &my_fdf, q);
    if (fid) {
      fprintf(fid,"%s Invoked with solver %s, MaxIter=%li\n",__func__,
	      gsl_multiroot_fdfsolver_name(s_fdf),MaxIter);
    }
  } else {
    s_f = gsl_multiroot_fsolver_alloc (T_f, q->size);
    my_f.n = q->size;
    my_f.f = &fzero_f;
    my_f.params = (void *)optfun;
    gsl_multiroot_fsolver_set (s_f, &my_f, q);
    if (fid) {
      fprintf(fid,"%s Invoked with solver %s, MaxIter=%li\n",__func__,
	      gsl_multiroot_fsolver_name(s_f),MaxIter);
    }
  }

  do
    {
      iter++;
      if (OPTIM_FZ_ISFDF(minimizer_flag)) {
	status = gsl_multiroot_fdfsolver_iterate (s_fdf);
      } else {
	status = gsl_multiroot_fsolver_iterate (s_f);
      }

      if (status) { /* check if solver is stuck */
	break;
      }

      if (OPTIM_FZ_ISFDF(minimizer_flag)) {
	status = gsl_multiroot_test_residual (s_fdf->f, 1e-7);
      } else {
	status = gsl_multiroot_test_residual (s_f->f, 1e-7);
      }

      if (fid && ((iter % soOften) == 0)) {
	/* update when to print status -- only when first digit of iter changes */
	if (iter >= soOften*10) {
	  soOften = (soOften < 1000) ? soOften*10 : 1000;
	}
	/* print status */
	if (OPTIM_FZ_ISFDF(minimizer_flag)) {
	  tmp = gsl_blas_dnrm2(s_fdf->f);
	} else {
	  tmp = gsl_blas_dnrm2(s_f->f);
	}
	fprintf(fid,"%s: %lu/%lu: %lg:\n",__func__,(unsigned long int)iter,(unsigned long int)MaxIter,tmp);
      }      

    }
  while (status == GSL_CONTINUE && iter < MaxIter);

  if (fid) {
    fprintf(fid,"%s: completed after %lu/%lu: %s\n",__func__,(unsigned long int)iter,(unsigned long int)MaxIter,gsl_strerror(status));
  }

  /* copy solution into q */
  if (OPTIM_FZ_ISFDF(minimizer_flag)) {
    gsl_vector_memcpy(q,s_fdf->x);
  } else {
    gsl_vector_memcpy(q,s_f->x);
  }
  tmp = optfun->func(q,optfun->params,NULL,NULL); /* evaluate ell */

  /* free memory */
  if (OPTIM_FZ_ISFDF(minimizer_flag)) {
    gsl_multiroot_fdfsolver_free (s_fdf);
  } else {
    gsl_multiroot_fsolver_free (s_f);
  }
  return(tmp);
}
