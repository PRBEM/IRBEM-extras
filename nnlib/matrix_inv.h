/* matrix inverse for 2x2 and 3x3 matrices */

#ifndef MATRIX_INV_H

#include <stdint.h> 
#include "subscripts.h"

void inv3x3(uint32_t Nt, double *A);
  /* inverse in place of an array of 3 x 3 matrices */

void inv2x2(uint32_t Nt, double *A);
  /* inverse in place of an array of 2 x 2 matrices */

#define MATRIX_INV_H 1
#endif
