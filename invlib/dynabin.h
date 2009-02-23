#ifndef DYNABIN_H

#define DYNABIN_H 1

typedef struct {
  long int Ntree,Ndims;

  /* btree definitions */
  unsigned char *c;
  long int *left, *right;
  double *Xc;

} dynabin;

dynabin dynabin_read(const char *filename);
/* read dynabin from file, allocates memory */

void dynabin_free(dynabin db);
/* free memory allocated by dynabin_read */

long int dynabin_lookup(const dynabin db, const double *x); 
/* returns 1-based bin index, or 0 if x is unbinable */
/* x is an array of db.Ndims coordinate values */

#endif
