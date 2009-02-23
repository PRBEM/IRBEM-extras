/*MINGW/MSYS:
gcc -I/usr/local/include -L/usr/local/lib -o cov_test.exe nn.c nnio.c matrix_inv.c cov_test.c -lgsl -lgslcblas -static
*/
/* test theta cov calculation */

#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include "nn.h"

int main () {
  unsigned long int i,j;
  double fval;
  net_type net;
  set_type set;
  char setname[1024] = "logflux_test.set";
  char netname[1024] = "logflux_cov.net";
  char outnetname[1024] = "logflux_hess.net";
  double *sinv;
  gsl_vector *grad;
  gsl_matrix *hess;
  double *hess_dbl;
  long int Ntheta;
  FILE *filep;
  double *Y, *dY;

  printf("Starting nntest\n");
  set = load_set(setname);
  printf("%s: Nt=%li, Nx=%li, Ny=%li,\nX[0]=%lg,Y[0]=%lg,s[0]=%lg,sflag=%li\n\n",
	 setname,set.Nt,set.Nx,set.Ny,set.X[0],set.Y[0],set.s[0],set.sflag);

  if (1) {
    sinv = make_sinv(set.Nt,set.Ny,set.s,set.sflag);

    net = load_net(netname);
    printf("%s: Nx=%li, Nh=%li, Ny=%li,\ntheta[0]=%lg,cov_flag=%li\n\n",
	   netname,net.Nx,net.Nh,net.Ny,net.theta[0],net.cov_flag);
    
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
    nnlib_save_net(outnetname,net.Nx,net.Nh,net.Ny,net.theta,2,hess_dbl);
    free(hess_dbl);

    filep = fopen("test.grad","wb");
    gsl_vector_fwrite(filep,grad);
    fclose(filep);
    
    filep = fopen("test.hess","wb");
    gsl_matrix_fwrite(filep,hess);
    fclose(filep);
    free(sinv);
    free(grad);
    gsl_matrix_free(hess);
  } else {
    net = load_net(outnetname);
  }

  Y = malloc(sizeof(double)*set.Nt*set.Ny);
  dY = malloc(sizeof(double)*set.Nt*set.Ny*set.Ny);
  nnlib_eval(set.Nt,set.Nx,set.X,net.Nh,net.theta,set.Ny,Y,4+2,net.theta_cov,dY);
  free(Y);
  free(dY);

  free_net(&net);
  free_set(&set);

  printf("Memory Freed\n");
  return(0);
}
