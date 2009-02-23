#include <stdio.h>
#include "dynabin.h"

long int dynabin_lookup(const dynabin db, const double *x) {
  /* fortran call might look something like this: 
     FUNCTION DYNABIN_LOOKUP(x,Ndims,Ntree,c,Xc,left,right)
     REAL*8 x(Ndims),Xc(Ntree)
     INTEGER*1 c(Ntree)
     INTEGER*4 left(Ntree),right(Ntree)
   */
  long int i=1;
  /* search the btree until a negative leaf index is found */
  while (i>0) {
    if (x[db.c[i-1]-1] <= db.Xc[i-1]) {
      i = db.left[i-1];
    } else {
      i = db.right[i-1];
    }
  }
  /* negative leaf index is the negative of the bin index */
  return(-i);
}

dynabin dynabin_read(const char *filename) {
  /* only needs Ndims, Ntree, c, Xc, left, right */
  /* all other info in file is discarded */
  FILE *filep;
  dynabin db;
  long int i,Nbins,c,left,right;
  double Xc;
  if (!(filep = fopen(filename,"rt"))) {
    perror("dynabin_read:Unable to open file");
    fprintf(stderr,"%s: file name: %s\n",__func__,filename);
    exit(-1);
  }
  fscanf(filep,"%*[NDIMS:]%li\n",&(db.Ndims));
  fscanf(filep,"%*[BINS:]%li\n",&Nbins);
  for (i=1; i <= Nbins; i++) {
    fscanf(filep,"%*[^\n]\n"); /* read line */
    
  }
  fscanf(filep,"%*[BTREE:]%li\n",&(db.Ntree));
  db.c = malloc(db.Ntree*sizeof(unsigned char));
  db.left = malloc(db.Ntree*sizeof(long int));
  db.right = malloc(db.Ntree*sizeof(long int));
  db.Xc = malloc(db.Ntree*sizeof(double));
  for (i=1; i <= db.Ntree; i++) {
    /* i,c,Xc,left,right,parent */
    fscanf(filep,"%*li,%li,%lg,%li,%li,%*li\n",&c,&Xc,&left,&right); /* read line */
    db.c[i-1] = c; /* loss of precission is OK */
    db.Xc[i-1] = Xc;
    db.left[i-1] = left;
    db.right[i-1] = right;
  }
  fclose(filep);
  return(db);
}

void dynabin_free(dynabin db) {
  free(db.c);
  free(db.left);
  free(db.right);
  free(db.Xc);
}


#ifdef DYNABIN_TEST
int main (int argc, char *argv[]) {
  FILE *filep;
  dynabin db;
  long int i,i2,Nbins,errors=0;
  unsigned char c;
  double *mid;
  char filename[] = "test.dynabin";

  db = dynabin_read(filename);
  mid = malloc(db.Ndims*sizeof(double));

  if (!(filep = fopen(filename,"rt"))) {
    perror("main:Unable to open file");
    fprintf(stderr,"%s: file name: %s\n",__func__,filename);
    exit(-1);
  }
  fscanf(filep,"%*[NDIMS:]%*li\n");
  fscanf(filep,"%*[BINS:]%li\n",&Nbins);
  for (i=1; i <= Nbins; i++) {
    /* i,EdgeMin(:),min(:),mid(:),max(:),EdgeMax(:),N,Uspan,parent */
    fscanf(filep,"%*li,"); /* i */
    /* EdgeMin,min */
    for (c=1; c <= db.Ndims*2; c++) {
      fscanf(filep,"%*lg,");
    }
    /* mid */
    for (c=1; c <= db.Ndims; c++) {
      fscanf(filep,"%lg,",&(mid[c-1]));
      printf("%g,",mid[c-1]);
    }
    i2 = dynabin_lookup(db,mid);
    if (i==i2) {
      printf("%li=%li\n",i,i2);
    } else {
      printf("%li!=%li **********\n",i,i2);
      errors++;
    }
    fscanf(filep,"%*[^\n]\n"); /* rest of line */
  }
  fclose(filep);

  free(mid);
  dynabin_free(db);
  printf("%li Errors\n",errors);
  return(1);
}
#endif
