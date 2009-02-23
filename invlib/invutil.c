/* utility functions for inv lib */

#include "invutil.h" /* includes gsl matrix, defines inv_matrix_once */
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>

/* matrix inverse functions */
void inv2x2_once(const gsl_matrix *A, gsl_matrix *invA) {
  /* quick inverse of 2 x 2 matrix */
  double a11,a12,a21,a22, det;
  a11 = gsl_matrix_get(A,0,0);
  a12 = gsl_matrix_get(A,0,1);
  a21 = gsl_matrix_get(A,1,0);
  a22 = gsl_matrix_get(A,1,1);
  det = a11*a22-a12*a21;
  gsl_matrix_set(invA,0,0,a22/det);
  gsl_matrix_set(invA,0,1,-a12/det);
  gsl_matrix_set(invA,1,0,-a21/det);
  gsl_matrix_set(invA,1,1,a11/det);
}

void inv_matrix_once(const gsl_matrix *A, gsl_matrix *invA, int symmetric) {
  /* this routine has been verified by multilpying invA*A */
  gsl_matrix *LU;
  gsl_permutation *p;
  int signum;
  
  if (A->size1==2) {
    inv2x2_once(A,invA);
  } else {
    LU = gsl_matrix_alloc(A->size1,A->size2);
    p = gsl_permutation_alloc(A->size1);
    gsl_matrix_memcpy(LU,A);
    gsl_linalg_LU_decomp(LU,p,&signum);
    gsl_linalg_LU_invert(LU,p,invA);
    if (symmetric) { /* force symmetry by taking invA = (invA'+invA)/2 */
      gsl_matrix_transpose_memcpy(LU,invA);
      gsl_matrix_add(invA,LU);
      gsl_matrix_scale(invA,0.5);
    }
    gsl_matrix_free(LU);
    gsl_permutation_free(p);
  }
}

gsl_vector *vector_func(const gsl_vector *x, gsl_vector *funcx, double (*func)(double)) {
  gsl_vector *out;
  double tmp;
  unsigned long int i;
  if (funcx) {
    out = funcx;
  } else {
    out = gsl_vector_alloc(x->size);
  }
  for (i=0; i < x->size; i++) {
    tmp = gsl_vector_get(x,i);
    tmp = (*func)(tmp);
    gsl_vector_set(out,i,tmp);
  }
  return(out);
}

gsl_matrix *matrix_func(const gsl_matrix *x, gsl_matrix *funcx, double (*func)(double)) {
  gsl_matrix *out;
  double tmp;
  unsigned long int i,j;
  if (funcx) {
    out = funcx;
  } else {
    out = gsl_matrix_alloc(x->size1,x->size2);
  }
  for (i=0; i < x->size1; i++) {
    for (j=0; j < x->size2; j++) {
      tmp = gsl_matrix_get(x,i,j);
      tmp = (*func)(tmp);
      gsl_matrix_set(out,i,j,tmp);
    }
  }
  return(out);
}
