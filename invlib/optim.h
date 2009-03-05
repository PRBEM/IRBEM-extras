/* header for generic optimize function */
#ifndef OPTIM_H

#include <stdio.h> /* FILE * */
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>

typedef struct {
  /* by having the function itself be a parameter, we can make the function optimize truly abstract */
  double (*func)(const gsl_vector *q, void *params, gsl_vector *grad, gsl_matrix *hess); /* function that computes fval, grad, hess */
  void *params; /* params argument for func */
} optfunTy;

double optimize(gsl_vector *q, optfunTy *optfun, const long int minimizer_flag, const long int MaxIter, FILE *fid);
  /* wrapper for GSL multimin routines
     computes q that minimizes optfun->func, with optfun->params 
     optfun->func returns fval, and grad, hess on request when grad, hess args not NULL
     inputs:
     q - initial guess at coefficients
     optfun - provides func to minimize and its params
     minimizer_flag controls which minimizer is used, as passed to ana_spec_inv
     MaxIter - maximum number of iterations for minimizer
     fid - NULL for no output, otherwise various messages are printed to fid with fprintf
     returns fval = ell = neglogp 
  */

double fzero(gsl_vector *q, optfunTy *optfun, const long int minimizer_flag, const long int MaxIter, FILE *fid);
  /* wrapper for GSL multiroot routines
     computes q that solves grad(optfun->func)=0
     optfun->func returns fval, and grad, hess on request when grad, hess args not NULL
     inputs:
     q - initial guess at coefficients
     optfun - provides func to minimize and its params
     minimizer_flag controls which solver is used, as passed to ana_spec_inv
     MaxIter - maximum number of iterations for solver
     fid - NULL for no output, otherwise various messages are printed to fid with fprintf
     returns fval = ell = neglogp 
  */

/* macros for minimizer */
/* fdf solvers */
#define OPTIM_MIN_BFGS  (0)
#define OPTIM_MIN_FR    (1)
#define OPTIM_MIN_PR    (2)
/* f solvers */
#define OPTIM_MIN_NM    (3)
/* other macros */
#define OPTIM_MIN_MAX   (3)
#define OPTIM_MIN_ISFDF(x) (((x)>=0) && ((x)<OPTIM_MIN_NM))

/* fdf solvers */
#define OPTIM_FZ_HYBRIDSJ (0)
#define OPTIM_FZ_HYBRIDJ  (1)
#define OPTIM_FZ_NEWTON   (2)
#define OPTIM_FZ_GNEWTON  (3)
/* f solvers */
#define OPTIM_FZ_HYBRIDS  (4)
#define OPTIM_FZ_HYBRID   (5)
#define OPTIM_FZ_DNEWTON  (6)
#define OPTIM_FZ_BROYDEN  (7)
/* other macros */
#define OPTIM_FZ_MAX      (7)
#define OPTIM_FZ_ISFDF(x) (((x)>=0) && (x<OPTIM_FZ_HYBRIDS))

#define OPTIM_H 1
#endif
