#include <stdlib.h>
#include <stdio.h>
#include "invlib.h"

/* test ana_spec_inv */

int main (int argc, char *argv[]) {
  FILE *fid;
  int result;
  long int i,j,NC,NE,NEout,int_params[10];
  double *c,*dc,*Egrid,*H,*b,*Eout,*flux,*dlogflux,real_params[10];
  char infile[] = "specinv_test.in1";
  char outfile[] = "specinv_test.out1";

  /* read c, dc, Egrid, H, b from ascii file */
  if (!(fid = fopen(infile,"r"))) {
    printf("Unable to load infile %s\n",infile);
    exit(-1);
  }
  
  fscanf(fid,"%li",&NC);
  fscanf(fid,"%li",&NE);
  fscanf(fid,"%li",&NEout);

  c = malloc(NC*sizeof(double));
  dc = malloc(NC*sizeof(double));
  Egrid = malloc(NE*sizeof(double));
  H = malloc(NC*NE*sizeof(double));
  b = malloc(NC*sizeof(double));
  Eout = malloc(NEout*sizeof(double));
  flux = malloc(NEout*sizeof(double));
  dlogflux = malloc(NEout*sizeof(double));

  for (i=0; i < NC; i++) {
    fscanf(fid,"%lf",&(c[i]));
  }

  for (i=0; i < NC; i++) {
    fscanf(fid,"%lf",&(dc[i]));
  }

  for (i=0; i < NE; i++) {
    fscanf(fid,"%lf",&(Egrid[i]));
  }

  for (i=0; i < NE*NC; i++) {
    fscanf(fid,"%lf",&(H[i]));
  }

  for (i=0; i < NC; i++) {
    fscanf(fid,"%lf",&(b[i]));
  }

  for (i=0; i < NEout; i++) {
    fscanf(fid,"%lf",&(Eout[i]));
  }

  fclose(fid);

  printf("specinv_test: NC=%li, NE=%li, NEout=%li\n",NC,NE,NEout);

  int_params[0] = NC;
  int_params[1] = NE;
  int_params[2] = NEout;
  int_params[3] = 1+2; /* analtyical functions bitmap*/
  int_params[4] = 0; /* minimizer, 0=BFGS, 3=NM */
  int_params[5] = 1000; /* maximumn # of iterations */
  int_params[6] = 1; /* 1 = verbose to standard out */
  int_params[7] = 0; /* 0 - dE already included in H */
  int_params[8] = 0; /* reserved */
  int_params[9] = 0; /* reserved */

  real_params[0] = 0.511; /* electron rest energy, MeV */
  real_params[1] = 100; /* E_break for proton PLE, MeV */
  real_params[2] = 345; /* E0 for proton PLE, MeV */
  real_params[3] = 0; /* reserved */
  real_params[4] = 0; /* reserved */
  real_params[5] = 0; /* reserved */
  real_params[6] = 0; /* reserved */
  real_params[7] = 0; /* reserved */
  real_params[8] = 0; /* reserved */
  real_params[9] = 0; /* reserved */

  result = ana_spec_inv(c,dc,Egrid,H,b,int_params,real_params,NULL,Eout,flux,dlogflux,NULL,NULL);
  printf("specinv_test: ana_spec_inv result= %i\n",result);

  /*
  for (j=0; j < NEout; j++) {
    printf("E=%lg\tflux=%lg\tdlogflux=%lg\n",Eout[j],flux[j],dlogflux[j]);
  }
  */

  if (!(fid = fopen(outfile,"w"))) {
    printf("specinv_test: Unable to load infile %s\n",outfile);
    exit(-1);
  }
  
  for (j=0; j < NEout; j++) {
//     fprintf(fid,"%lg,%lg,%lg\n",Eout[j],flux[j],dlogflux[j]);
    fprintf(fid,"%6.4f,%6.4e,%6.4e\n",Eout[j],flux[j],dlogflux[j]);
  }
  fclose(fid);

  free(c);
  free(dc);
  free(Egrid);
  free(H);
  free(b);
  free(Eout);
  free(flux);
  free(dlogflux);


  return(0);
}
