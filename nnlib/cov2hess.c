/*MINGW/MSYS:
gcc -I/usr/local/include -L/usr/local/lib -o cov2hess.exe nn.c nnio.c matrix_inv.c cov2hess.c -lgsl -lgslcblas -static
*/
/* test theta cov calculation */

#include <math.h>
#include <stdlib.h>
#include <string.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include "nn.h"

#define DEFAULT_SETFILE "logflux_test.set"

void print_usage() {
  fprintf(stderr,"cov2hess <netfile> <outfile> <setfile>\n");
  fprintf(stderr,"Evaluates <netfile> on <setfile> and saves hess to <outfile>\n");
  fprintf(stderr,"<setfile> defaults to %s\n",DEFAULT_SETFILE);
  fprintf(stderr,"<outfile> defaults to <netfile>\n");
  fprintf(stderr,"<netfile> is a required input\n");
  fprintf(stderr,"\n\n");
}

int main (int argc, char *argv[]) {
  unsigned long int i,j;
  double fval;
  net_type net;
  set_type set;
  char setname[1024] = DEFAULT_SETFILE;
  char netname[1024];
  char outnetname[1024];
  double *sinv;
  gsl_vector *grad;
  gsl_matrix *hess;
  double *hess_dbl;
  long int Ntheta;
  FILE *filep;
  double *Y, *dY;

  if (argc < 1) {
    fprintf(stderr,"Need at least name of net to read\n");
    print_usage();
    exit(-1);
  }

  strcpy(netname,argv[1]);
  if (argc>=3) {
    strcpy(outnetname,argv[2]);
  } else {
    strcpy(outnetname,netname);
  }

  if (argc>=4) {
    strcpy(setname,argv[3]);
  }

  printf("Converting %s to %s (hess) using %s\n",netname,outnetname,setname);

  printf("Loading net: %s\n",netname);
  net = load_net(netname);
  printf("%s: Nx=%li, Nh=%li, Ny=%li,\ntheta[0]=%lg,cov_flag=%li\n\n",
	 netname,net.Nx,net.Nh,net.Ny,net.theta[0],net.cov_flag);

  if ((net.cov_flag & 3) == 2) {
    printf("%s already has hess\n",netname);
    if (! strcmp(netname,outnetname)) {
      printf("Saving %s to %s\n",netname,outnetname);
      save_net(outnetname,net);
    }
  } else {
    
    printf("Loading set: %s\n",setname);
    set = load_set(setname);
    printf("%s: Nt=%li, Nx=%li, Ny=%li,\nX[0]=%lg,Y[0]=%lg,s[0]=%lg,sflag=%li\n\n",
	   setname,set.Nt,set.Nx,set.Ny,set.X[0],set.Y[0],set.s[0],set.sflag);
    
    if ((net.Nx != set.Nx) || (net.Ny != set.Ny)) {
      fprintf(stderr,"Net and Set have different Nx or Ny\nAborting.\n");
      exit(-1);
    }

    sinv = make_sinv(set.Nt,set.Ny,set.s,set.sflag);
    
    Ntheta = theta_Ntheta(net.Nx,net.Nh,net.Ny);
    grad = gsl_vector_alloc(Ntheta);
    hess = gsl_matrix_alloc(Ntheta,Ntheta);
    
    fval = nn_fit_eval(set.Nt,set.Nx,set.X,net.Nh,set.Ny,set.Y,sinv,set.sflag,net.theta,grad,hess);
    hess_dbl = malloc(sizeof(double)*Ntheta*Ntheta);
    for (i=1; i <= Ntheta; i++) {
      for (j=1; j <= Ntheta; j++) {
	hess_dbl[ss2(Ntheta,Ntheta,i,j)] = gsl_matrix_get(hess,i-1,j-1);
      }
    }
    printf("Saving net: %s\n",outnetname);
    nnlib_save_net(outnetname,net.Nx,net.Nh,net.Ny,net.theta,2,hess_dbl);
    free(sinv);
    free(hess_dbl);
    free(grad);
    gsl_matrix_free(hess);
    free_set(&set);
  }
  free_net(&net);
  printf("Memory Freed\n Done.");
  return(0);
}
