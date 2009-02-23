/* neural network codes
   train and evaluate neural networks
 */

#ifndef NN_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include "nnlib.h"
#include "subscripts.h"
#include "matrix_inv.h"

/* theta ordering is not intuitive, but built to match Matlab */
#define theta_w(Nx,Nh,Ny,j,i) (((i)-1)*(Nx) + (j)-1)
#define theta_v(Nx,Nh,Ny,i,k) ((Nh)*(Nx) + ((k)-1)*(Nh)+(i)-1) 
#define theta_w0(Nx,Nh,Ny,i) ((Nh)*(Nx) + (Ny)*(Nh) + (i)-1)
#define theta_v0(Nx,Nh,Ny,k) ((Nh)*(Nx) + (Ny)*(Nh) + (Nh) + (k)-1)
#define theta_Ntheta(Nx,Nh,Ny) ((Ny) + (Ny)*(Nh) + (Nh) + (Nh)*(Nx))

#define NN_OPT_BFGS (0)
#define NN_OPT_FR   (32)
#define NN_OPT_PR   (64)
#define NN_OPT_NM   (96)
#define NN_OPT_BITS   (96)

typedef struct {
  unsigned long int Nx, Nh, Ny, cov_flag;
  double *theta,*theta_cov;
  double *xbar,*ybar,*sx,*sy; /* sample mean, std of training set */
} net_type;

typedef struct {
  unsigned long int Nt, Nx, Ny, sflag;
  double *X,*Y,*s;
} set_type;

unsigned long int get_Ns(const unsigned long int sflag,const unsigned long int Nt,const unsigned long int Ny);

double *make_sinv(const unsigned long int Nt, const unsigned long int Ny, const double *s, const unsigned long int flag);

double nn_fit_eval(const unsigned long int Nt, const unsigned long int Nx, const double *Z, 
		   const unsigned long int Nh, 
		   const unsigned long int Ny, const double *Y, 
		   const double *ybar, const double *sy,
		   const double *sinv, const unsigned long int flag,
		   const double *theta,
		   gsl_vector *grad, gsl_matrix *hess);
  /* flag has same meaning as in nnlib_fit, sinv has same size */


void theta_id(const unsigned long int Nx, const unsigned long int Nh, const unsigned long int Ny, const unsigned long int n, char *buf);
/* returns string description of theta index n: e.g. w(7,21) (1-based indexing)
   buf must be big enough to hold string (100 chars should more than do it */

net_type init_net();
net_type alloc_net(unsigned long int Nx,unsigned long int Nh,unsigned long int Ny,unsigned long int cov_flag);
void free_net(net_type *net);
net_type load_net(const char *filename);
void save_net(const char *filename,net_type net);

set_type init_set();
set_type alloc_set();
void free_set(set_type *set);
set_type load_set(const char *filename);
void save_set(const char *filename,set_type set);

/* mpi routines */
int nn_mpi_init(int *argc,char **argv[]);
set_type nn_mpi_dispatch(set_type set, net_type net);
void nn_mpi_slave();
void nn_mpi_cleanup();

#define NN_H 1
#endif
