/*
see Makefile for compile switches
see function print_usage below for syntax

to enable parallelization by MPI use -DUSEMPI (important for nn.c)
to enable parallelization by OpenMP use -DUSEOMP (important for nn.c)

*/


#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "nn.h"

void print_usage(){
  fprintf(stderr,"nntrain.exe <setfile> <Nh> <netfile> ... \n");
  fprintf(stderr,"nntrain.x <setfile> <Nh> <netfile> ... \n");
  fprintf(stderr,"trains a net with <Nh> hidden nodes on <setfile>, saves to <netfile>\n");
  fprintf(stderr,"(nntrain.mpix) Uses MPI to divide the job among multiple processors\n");
  fprintf(stderr,"(nntrain.ompx) Uses OpenMP to divide the job among multiple processors\n");
  fprintf(stderr,"nntrain ... -bfgs  (train with vector bfgs)\n");
  fprintf(stderr,"... -fr    (train with congugate fr)\n");
  fprintf(stderr,"... -pr    (train with congugate pr)\n");
  fprintf(stderr,"... -nm    (train with Nelder-Mead Simplex)\n");
  fprintf(stderr,"... -mI <maxIter> (set maximum iterations, to be done each cycle, default is 50000)\n");
  fprintf(stderr,"... -# <maxIter> (depricated, same as mI)\n");
  fprintf(stderr,"... -c <NumCycles> (set number of cycles, default is 3)\n");
  fprintf(stderr,"... -nm0   (pre-train with 1000 iterations of N-M before each cycle)\n");
  fprintf(stderr,"... -nm0 <nm0Iter> (set number of iterations of N-M for nm0, default is 1000)\n");
  fprintf(stderr,"... -e <epsabs> (set fit tolerance, default is 1e-4)\n");
  fprintf(stderr,"... -val <setfile> (compute out-of-sample penalty on <setfile>)\n");
  fprintf(stderr,"\n\n");
}

void cleanup(int exit_code){
  /* cleans up, including sending end op code 0 to mpi children */
  /* if exit_code !=0, exits */
  nn_mpi_cleanup(); /* kill slaves as needed */
  if(exit_code) {
    exit(exit_code);
  }
}

int main (int argc, char *argv[]) {
  unsigned long int sflag;
  unsigned long int i, nm0Iter=1000, maxIter=50000,NumCycles=3;
  double fval;
  net_type net;
  set_type set,subset;
  char setname[1024];
  char outnetname[1024];
  char valname[1024];
  double epsabs = 1e-4;
  int do_val = 0;
  int pre_train_nm = 0;
  unsigned long int Nh=10;
  unsigned long int nn_opt_method=0;


  if (nn_mpi_init(&argc,&argv) != 0) { /* initialize mpi, divert for slave */
    nn_mpi_slave(); /* prefab nn_mpi_slave function */
    return(0);
  } /* else master */


  /* read arguments */
  printf("%s: reading %d args\n",__func__,argc); /* debug */
  if (argc < 4) {
    if (!((argc == 2) && (! strcmp(argv[1],"--help")))) {
      fprintf(stderr,"%s: too few arguments\n",__func__);
    }
    print_usage();
    cleanup(-1);
  }
  strcpy(setname,argv[1]);
  Nh = atol(argv[2]);
  strcpy(outnetname,argv[3]);
  for (i=4; i < argc; i++) {
    if (argv[i][0] == '-') {
      if (! strcmp(argv[i],"-bfgs")) {
	nn_opt_method = NN_OPT_BFGS;
      } else if (! strcmp(argv[i],"-fr")) {
	nn_opt_method = NN_OPT_FR;
      } else if (! strcmp(argv[i],"-pr")) {
	nn_opt_method = NN_OPT_PR;
      } else if (! strcmp(argv[i],"-nm0")) {
	pre_train_nm = 1;
	if ((i+1 < argc) && (argv[i+1][0] != '-')) {
	  /* a non-switch arg follows -nm0, so it's nm0Iter */
	  i++; 
	  nm0Iter = atol(argv[i]);
	}
      } else if (! strcmp(argv[i],"-nm")) {
	nn_opt_method = NN_OPT_NM;
      } else if ((! strcmp(argv[i],"-mI")) || (! strcmp(argv[i],"-#"))) {
	i++;
	maxIter = atol(argv[i]);
      } else if (! strcmp(argv[i],"-c")) {
	i++;
	NumCycles = atol(argv[i]);
      } else if (! strcmp(argv[i],"-e")) {
	i++;
	epsabs = atof(argv[i]);
      } else if (! strcmp(argv[i],"-val")) {
	i++;
	do_val = 1;
	strcpy(valname,argv[i]);
      } else {
	switch(argv[i][1]) {
	  /* other options here */
	default:
	  fprintf(stderr,"%s: unknown argument '%s'",__func__,argv[i]);
	  print_usage();
	  cleanup(-1);
	  break;
	} /* end switch */
      } /* end else */
    } /* end if */
  } /* end for */
  
  set = load_set(setname);
  printf("%s: set info: Nt=%li, Nx=%li, Ny=%li,\nX[0]=%lg,Y[0]=%lg,s[0]=%lg,sflag=%li\n\n",
	 setname,set.Nt,set.Nx,set.Ny,set.X[0],set.Y[0],set.s[0],set.sflag);
  net = alloc_net(set.Nx,Nh,set.Ny,2); /* alloc but use theta hess instead of cov */
  
  subset = nn_mpi_dispatch(set,net); /* dispatch subsets, create virtual subset for master (mpi) */

  printf("%s: Starting Training\n",__func__); /* debug */
  
  for (i=1; i <= NumCycles; i++) {
    
    if (i==1) {
      sflag = subset.sflag | 8; /* random initial theta (no 16), 8 is for verbose output */
    } else {
      sflag = subset.sflag | 8 | 16; /* start at theta, don't compute hess (no 128), verbose */
    }
    
    if (pre_train_nm) {
      fval = nnlib_fit(subset.Nt,subset.Nx,subset.X,net.Nh,subset.Ny,subset.Y,subset.s,sflag | NN_OPT_NM,
		       nm0Iter,epsabs,net.theta,net.xbar,net.ybar,net.sx,net.sy,net.theta_cov);
      sflag |= 16; /* start at theta hereafter */
    };
    
    if (i==NumCycles) {
      sflag |= 128; /* compute theta hess */
    }
    fval = nnlib_fit(subset.Nt,subset.Nx,subset.X,net.Nh,subset.Ny,subset.Y,subset.s,sflag | nn_opt_method,
		     maxIter,epsabs,net.theta,net.xbar,net.ybar,net.sx,net.sy,net.theta_cov);
    printf("After cycle %lu/%lu: fval=%lg\n",i,NumCycles,fval);
    
    printf("%s: Nx=%li, Nh=%li, Ny=%li\n",outnetname,net.Nx,net.Nh,net.Ny);
    save_net(outnetname,net);
  }
  
  free_set(&set);
  
  if (do_val) {
    
    set = load_set(valname);
    if ((set.Nx != net.Nx) || (set.Ny != net.Ny)) {
      fprintf(stderr,"ERRROR: Cannot do out-of-sample test, %s and %s do not have same Nx, Ny\n",setname,valname);
    } else {
      fval = nnlib_ell(set.Nt,set.Nx,set.X,net.Nh,net.Ny,set.Y,
		       net.xbar,net.ybar,net.sx,net.sy,set.s,set.sflag,net.theta,NULL,NULL);
      printf("Out of sample (%s): fval=%lg\n",valname,fval);
    }
    free_set(&set);
    free_net(&net);
  }
  
  printf("Memory Freed\n");
  
  cleanup(0);
  
  return(0);
}
