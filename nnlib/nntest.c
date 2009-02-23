/*MINGW/MSYS:
gcc -I/usr/local/include -L/usr/local/lib -o nntest.exe nn.c nnio.c matrix_inv.c nntest.c -lgsl -lgslcblas -static
*/

#include <math.h>
#include <stdlib.h>
#include <gsl/gsl_math.h>
#include "nn.h"

void test_grad_hess(unsigned long int Nt, unsigned long int Nx, unsigned long int Nh, unsigned long int Ny, 
		    double *X, double *Y, double *theta, double *s,unsigned long int sflag) {
  /* evaulate net and its grad, hess */
  double *sinv=0;
  double theta_tmp,dw,fval2,gi,gi2,dg;
  unsigned long int Ntheta,i,j;
  double hi,hi2,dh;
  gsl_vector *grad2=0;
  double fval;
  gsl_vector *grad=0;
  gsl_matrix *hess=0;
  char buf[255],buf2[255];

  Ntheta = theta_Ntheta(Nx,Nh,Ny);
  sinv = make_sinv(Nt,Ny,s,sflag);

  grad = gsl_vector_alloc(Ntheta);
  hess = gsl_matrix_alloc(Ntheta,Ntheta);
  fval = nn_fit_eval(Nt,Nx,X,Nh,Ny,Y,sinv,sflag,theta,grad,hess);
  printf("fval: %lg\n\n",fval);
  printf("gradient test:\n");
  for (i=1; i <= Ntheta; i++) { 
    gi = gsl_vector_get(grad,i-1);
    theta_tmp = theta[i-1];
    dw = gsl_max(fabs(theta_tmp)/1e8,1e-10);
    theta[i-1] = theta_tmp+dw;
    fval2 = nn_fit_eval(Nt,Nx,X,Nh,Ny,Y,sinv,sflag,theta,0,0);
    theta[i-1] = theta_tmp;
    gi2 = (fval2-fval)/dw;
    dg = fabs(gi-gi2)/(fabs(gi)+fabs(gi2))*2.0;
    if (0 || ((dg>1e-3) && (gsl_max(fabs(gi),fabs(gi2))>=1e-8))) {
      theta_id(Nx,Nh,Ny,i,buf);
      printf("%li;%s](%lg,%lg) %lg-%lg = %lg\n",i,buf,theta_tmp,dw,gi,gi2,dg);
    }
  }
  printf("\n");

  printf("hessian test:\n");
  grad2 = gsl_vector_alloc(Ntheta);
  for (i=1; i <= Ntheta; i++) {
    theta_id(Nx,Nh,Ny,i,buf);
    theta_tmp = theta[i-1];
    dw = gsl_max(fabs(theta_tmp)/1e6,1e-8);
    theta[i-1] = theta_tmp+dw;
    fval2 = nn_fit_eval(Nt,Nx,X,Nh,Ny,Y,sinv,sflag,theta,grad2,0);
    theta[i-1] = theta_tmp;
    for (j=1; j <= Ntheta; j++) {
      hi = gsl_matrix_get(hess,i-1,j-1);
      hi2 = (gsl_vector_get(grad2,j-1)-gsl_vector_get(grad,j-1))/dw;
      dh = fabs(hi-hi2)/(fabs(hi)+fabs(hi2))*2.0;
      if (0 || ((dh>1e-4) && (gsl_max(fabs(hi),fabs(hi2))>1e-8))) {
	theta_id(Nx,Nh,Ny,j,buf2);
	printf("%li,%li;%s,%s](%lg,%lg) %lg-%lg = %lg\n",i,j,buf,buf2,theta_tmp,dw,hi,hi2,dh);
      }
    }
  }
  printf("\n");
  
  printf("\n Done with tests\n");

  if (sinv) {
    free(sinv);
  }

  if (grad) {
    gsl_vector_free(grad);
  };

  if (grad2) {
    gsl_vector_free(grad2);
  };

  if (hess) {
    gsl_matrix_free(hess);
  }
}


int main () {
  double fval;
  net_type net;
  set_type set;
  char setname[] = "c:\\temp2\\test1.set";
  char netname[] = "c:\\temp2\\test.net";
  char outnetname[] = "c:\\temp2\\test1.net";

  printf("Starting nntest\n");
  set = load_set(setname);
  printf("%s: Nt=%li, Nx=%li, Ny=%li,\nX[0]=%lg,Y[0]=%lg,s[0]=%lg,sflag=%li\n\n",
	 setname,set.Nt,set.Nx,set.Ny,set.X[0],set.Y[0],set.s[0],set.sflag);

  net = alloc_net(set.Nx,20,set.Ny,1);

  /*
  net = load_net(netname);
  printf("%s: Nx=%li, Nh=%li, Ny=%li,\ntheta[0]=%lg,cov_flag=%li\n\n",
	 netname,net.Nx,net.Nh,net.Ny,net.theta[0],net.cov_flag);
  */

  /* test_grad_hess(Nt,Nx,Nh,Ny,X,Y,theta,s,sflag);  */
  fval = nnlib_fit(set.Nt,set.Nx,set.X,net.Nh,set.Ny,set.Y,set.s,set.sflag | 4 | 8 | NN_OPT_BFGS,100000,net.theta,net.theta_cov);
  save_net(outnetname,net);

  free_net(&net);
  free_set(&set);

  printf("Memory Freed\n");
  return(0);
}
