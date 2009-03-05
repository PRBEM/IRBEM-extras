/* utility functions for invlib */

#ifndef INVUTIL_H

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>

/* invert matrix, if symmetric, foces invA to be symmetric */
void inv_matrix_once(const gsl_matrix *A, gsl_matrix *invA, int symmetric);

/* these functions apply func to all elements of a vector or matrix */
/* if funcx is NULL, then a new fector is created and returned, must be deallocated */
/* if funcx is not NULL, then the result is stored in funcx */
/* can operate in place if funcx = x */
gsl_vector *vector_func(const gsl_vector *x, gsl_vector *funcx, double (*func)(double));
gsl_matrix *matrix_func(const gsl_matrix *x, gsl_matrix *funcx, double (*func)(double));

#define INVUTIL_H 1
#endif
