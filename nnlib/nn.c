/* 
   main code for neural network library
   see Makefile for compile switches
   to enable parallelization by MPI use -DUSEMPI (important for nn.c)
*/


#include <stdlib.h>
#include <time.h>
#include <string.h>
#ifdef USEMPI
#include <mpi.h>
#endif
#include <gsl/gsl_math.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_permutation.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_multimin.h>
#include "nn.h"
#include "matrix_inv.h"

#ifdef USEMPI

/* define op codes for MPI */
#define NN_MPI_OP_EXIT (-1)
#define NN_MPI_OP_INIT (0)
#define NN_MPI_OP_ELL  (1)
#define NN_MPI_OP_GRAD (2)
#define NN_MPI_OP_HESS (3)

/* define tag codes for MPI */
#define NN_MPI_TAG_OP     (0)
#define NN_MPI_TAG_ELL    (1)
#define NN_MPI_TAG_GRAD   (2)
#define NN_MPI_TAG_HESS   (3)
#define NN_MPI_TAG_NT     (10)
#define NN_MPI_TAG_NX     (11)
#define NN_MPI_TAG_NY     (12)
#define NN_MPI_TAG_SFLAG  (13)
#define NN_MPI_TAG_X      (21)
#define NN_MPI_TAG_Y      (22)
#define NN_MPI_TAG_S      (23)
#define NN_MPI_TAG_XBAR   (31)
#define NN_MPI_TAG_YBAR   (32)
#define NN_MPI_TAG_SX     (41)
#define NN_MPI_TAG_SY     (42)
#define NN_MPI_TAG_NH     (100)
#define NN_MPI_TAG_THETA  (101)

#endif

#ifdef USEOMP
/* use OpenMP */
#include <omp.h>
#endif

double ofunc(double u) { /* o(u) */
  const double min_u = -500.0;
  if (u > min_u) { /* avoid overflow, i.e., inf */
    return(exp(-u));
  } else {
    return(exp(-min_u));
  }
}
double gfunc(double o) { /* g(o), o = o(u) */
  return(1.0/(1.0+o));
}

unsigned long int get_Ns(unsigned long int sflag,unsigned long int Nt,unsigned long int Ny) {
  unsigned long int Ns;
  switch(sflag & 3) {
  case 0: /* scalar */
  default:
    Ns = 1;
    break;
  case 1: /* 1 x Ny */
    Ns = Ny;
    break;
  case 2: /* Nt x Ny */
    Ns = Ny*Nt;
    break;
  case 3: /* Nt x Ny x Ny */
    Ns = Ny*Ny*Nt;
    break;
  }
  return(Ns);
}

unsigned long int get_Nst(unsigned long int sflag,unsigned long int Ny) {
  /* returns number of s values per time step */
  unsigned long int Nst;
  switch(sflag & 3) {
  case 0: /* scalar */
  case 1: /* 1 x Ny */
  default:
    Nst = 0;
    break;
  case 2: /* Nt x Ny */
    Nst = Ny;
    break;
  case 3: /* Nt x Ny x Ny */
    Nst = Ny*Ny;
    break;
  }
  return(Nst);
}

void stats_Z(const unsigned long int Nt,
	     const unsigned long int Nx, const double *X, 
	     const unsigned long int Ny, const double *Y, 
	     double *xbar, double *ybar,
	     double *sx, double *sy, double *Z,
	     const unsigned long int flag) {
  /* sets xbar, ybar, sx, sy, and (if not NULL), Z */
  /* requires xbar, ybar, sx, sy allocated */
  /* flag == 0 set xbar, sx, ybar, sy, else use provided values */

  unsigned long int t,j,k;
  unsigned long int nx[Nx], ny[Ny];
  double tmp;


  if (flag == 0) {
    /* init means and stds to zero */
    
    for (j=0; j < Nx; j++) {
      nx[j] = 0;
      xbar[j] = 0;
      sx[j] = 0;
    }
    
    for (j=0; j < Ny; j++) {
      ny[j] = 0;
      ybar[j] = 0;
      sy[j] = 0;
    }
    
    /* aggregate sum of x in xbar, sum of x^2 in sx, same for y */
    
    for (t=1; t <= Nt; t++) {
      for (j=1; j <= Nx; j++) {
	tmp = X[ss2(Nt,Nx,t,j)]; 
	if (gsl_finite(tmp)) {
	  k = ss1(Nx,j);
	  nx[k]++;
	  xbar[k] += tmp;
	  sx[k] += gsl_pow_2(tmp);
	}
      }
      for (j=1; j <= Ny; j++) {
	tmp = Y[ss2(Nt,Ny,t,j)]; 
	if (gsl_finite(tmp)) {
	  k = ss1(Ny,j);
	  ny[k]++;
	  ybar[k] += tmp;
	  sy[k] += gsl_pow_2(tmp);
	}
      }
    }
    
    /* normalize xbar and sx, same for y */
    
    for (j=1; j <= Nx; j++) {
      k = ss1(Nx,j);
      xbar[k] /= nx[k]; /* xbar = mean(x) */
      sx[k] /= nx[k]; /* sx = mean(x^2) */
      sx[k] -= gsl_pow_2(xbar[k]); /* sx = mean(x^2)-mean(x)^2 = var(x) */
      sx[k] = sqrt(sx[k]); /* sx = std(x) */
    }
    
    for (j=1; j <= Ny; j++) {
      k = ss1(Ny,j);
      ybar[k] /= ny[k]; /* ybar = mean(y) */
      sy[k] /= ny[k]; /* sy = mean(y^2) */
      sy[k] -= gsl_pow_2(ybar[k]); /* sy = mean(y^2)-mean(y)^2 = var(y) */
      sy[k] = sqrt(sy[k]); /* sy = std(y) */
    }
    
  }

  if (Z) {
    /* compute Z = (X-xbar)/s(x) */

    for (t=1; t <= Nt; t++) {
      for (j=1; j <= Nx; j++) {
	k = ss2(Nt,Nx,t,j);
	Z[k] = (X[k]-xbar[ss1(Nx,j)])/sx[ss1(Nx,j)];
      }
    }
  }

}

void make_gt(const unsigned long int Nt,const unsigned long int Nx, const double *Z, 
	     const unsigned long int Nh, const double *theta, 
	     const unsigned long int Ny, const unsigned long int t,
	     double *gt, double *dgt, double *ddgt) {
  /*
    computs g(t,:) for X(t,:) and derivatives (if points not NULL) 
    gt , dgt, ddgt are arrays of Nh doubles
  */

  unsigned long int i,j;
  double oi,ui;
  for (i=1; i <= Nh; i++) {
    ui = theta[theta_w0(Nx,Nh,Ny,i)]; /* w0i */
    for (j=1; j <= Nx; j++) {
      ui += theta[theta_w(Nx,Nh,Ny,j,i)]*Z[ss2(Nt,Nx,t,j)]; /* wji*Xtj */
    }
    oi = ofunc(ui);
    gt[ss1(Nh,i)] = gfunc(oi);

    if (dgt) {
      dgt[ss1(Nh,i)] = gsl_pow_2(gt[ss1(Nh,i)])*oi; /* dg/du */
      if (ddgt) {
	ddgt[ss1(Nh,i)] = (oi-1.0)*gt[ss1(Nh,i)]*dgt[ss1(Nh,i)]; /* d^2g/du^2 */
      }
    }
  }
}

void make_ht(const unsigned long int Nt,const unsigned long int Nx, const double *Z, 
	     const unsigned long int Nh, const double *theta, 
	     const unsigned long int Ny, const unsigned long int t, 
	     const double *ybar,
	     const double *sy,
	     double *ht, gsl_matrix *dht, double *ddht) {
  /* 
     computes h(t,:) for X(t,:)  [ h(t,:) = Y(t,:) ]
     ht is an array of Ny doubles
     dht is a matrix of Ntheta x Ny doubles
     ddht is a tensor of Ntheta x Ntheta x Ny doubles
   */

  unsigned long int i,k,ip,jp,jpp,Ntheta,ind1,ind2,ind3;
  double hk,temp,temp2,temp3, syk;
  double *gt=0, *dgt=0, *ddgt=0;

  gt = calloc(Nh,sizeof(double));
  if (dht) {
    gsl_matrix_set_zero(dht);
    dgt = calloc(Nh,sizeof(double));
  }

  if (ddht) {
    Ntheta = theta_Ntheta(Nx,Nh,Ny);
    memset(ddht,0,sizeof(double)*Ntheta*Ntheta*Ny);
    ddgt = calloc(Nh,sizeof(double));
  }

  make_gt(Nt,Nx,Z,Nh,theta,Ny,t,gt,dgt,ddgt);
  for (k=1; k <= Ny; k++) {
    hk = theta[theta_v0(Nx,Nh,Ny,k)]; /*v0k*/
    for (i=1; i <= Nh; i++) {
      hk += theta[theta_v(Nx,Nh,Ny,i,k)]*gt[ss1(Nh,i)]; /*v(i,k)*/
    }
    syk = sy[ss1(Ny,k)];
    hk *= syk; /* scale by sy[k] */
    hk += ybar[ss1(Ny,k)]; /* add ybar[k] */

    ht[ss1(Ny,k)] = hk;

    if (dht) {
      /* compute dhk/dtheta */
      /* note: kp = k, due to delta function*/

      /* d/dv0 */
      gsl_matrix_set(dht,theta_v0(Nx,Nh,Ny,k),ss1(Ny,k),syk); /* v0(kp) */

      for (ip = 1; ip <= Nh; ip++) {
	/* d/dv */
	gsl_matrix_set(dht,theta_v(Nx,Nh,Ny,ip,k),ss1(Ny,k),syk*gt[ss1(Nh,ip)]); /* v(ip,kp) */
	
	/* d/dw0 */
	temp = syk*theta[theta_v(Nx,Nh,Ny,ip,k)]*dgt[ss1(Nh,ip)];
	gsl_matrix_set(dht,theta_w0(Nx,Nh,Ny,ip),ss1(Ny,k),temp); /* w0(ip) */
	
	/* d/dw */
	for (jp = 1; jp <= Nx; jp++) {
	  gsl_matrix_set(dht,theta_w(Nx,Nh,Ny,jp,ip),ss1(Ny,k),temp*Z[ss2(Nt,Nx,t,jp)]); /* w(jp,ip) */
	} /* end for jp */
      } /* end for ip */

      if (ddht) {
	/* hessians involving /dv0/d* or /dv/dv are all zero */

	for (ip = 1; ip <= Nh; ip++) {
	  /* note: ipp = ip, kpp = kp = k, due to delta functions */
	  /* d^2hk/dv/dw0 */
	  ind1 = 1+theta_v(Nx,Nh,Ny,ip,k); /* v(ip,kp) */
	  ind2 = 1+theta_w0(Nx,Nh,Ny,ip); /* w0(ip) */
	  temp = syk*dgt[ss1(Nh,ip)];
	  ddht[ss3(Ntheta,Ntheta,Ny,ind1,ind2,k)] = temp;
	  ddht[ss3(Ntheta,Ntheta,Ny,ind2,ind1,k)] = temp;
	  /* d^hk/dv/dw */
	  for (jpp = 1; jpp <= Nx; jpp++) {
	    ind1 = 1+theta_v(Nx,Nh,Ny,ip,k); /* v(ip,kp) */
	    ind2 = 1+theta_w(Nx,Nh,Ny,jpp,ip); /* w(jpp,ipp) */
	    temp2 = syk*dgt[ss1(Nh,ip)]*Z[ss2(Nt,Nx,t,jpp)];
	    ddht[ss3(Ntheta,Ntheta,Ny,ind1,ind2,k)] = temp2;
	    ddht[ss3(Ntheta,Ntheta,Ny,ind2,ind1,k)] = temp2;
	  } /* end for jpp */

	  /* d^2hk/dw0/dw0 */
	  ind1 = 1+theta_w0(Nx,Nh,Ny,ip);
	  temp = syk*theta[theta_v(Nx,Nh,Ny,ip,k)]*ddgt[ss1(Nh,ip)];
	  ddht[ss3(Ntheta,Ntheta,Ny,ind1,ind1,k)] = temp;
	  for (jpp = 1; jpp <= Nx; jpp++) {
	    /* d^2hk/dw0/dw */
	    ind2 = 1+theta_w(Nx,Nh,Ny,jpp,ip);
	    temp2 = temp*Z[ss2(Nt,Nx,t,jpp)];
	    ddht[ss3(Ntheta,Ntheta,Ny,ind1,ind2,k)] = temp2;
	    ddht[ss3(Ntheta,Ntheta,Ny,ind2,ind1,k)] = temp2;
	    for (jp = 1; jp <= Nx; jp++) {
	      /* d^2hk/dw/dw */
	      ind3 = 1+theta_w(Nx,Nh,Ny,jp,ip);
	      temp3 = temp2*Z[ss2(Nt,Nx,t,jp)];
	      ddht[ss3(Ntheta,Ntheta,Ny,ind3,ind2,k)] = temp3;
	      ddht[ss3(Ntheta,Ntheta,Ny,ind2,ind3,k)] = temp3;
	    } /* end for jp */
	  } /* end for jpp */
	} /* end for ip */
      } /* end if ddht */
    } /* end if dht */
  } /* end for k */

  free(gt);

  if (dgt) {
    free(dgt);
  }
  if (ddgt) {
    free(ddgt);
  }

} /* end make_ht */

void nnlib_eval(const unsigned long int Nt,
		const unsigned long int Nx, const double *X, 
		const unsigned long int Nh, const double *theta, 
		const double *xbar, const double *ybar,
		const double *sx, const double *sy,
		const unsigned long int Ny, double *Y, 
		const unsigned long int dY_flag, const double *theta_cov, double *dY) {
/* 
   dY_flag & 3 = 0: don't compute dY (theta_cov ignored)
   dY_flag & 3 = 1: for Nt x Ny dY (standard errors)
   dY_flag & 3 = 2: for Nt x Ny x Ny dY (covariance matrix)
   dy_flag & 4 = 4 theta_cov is really theta_hess
   OUTPUTS: Y, dY
 */

  unsigned long int t,i,j,k,Ntheta;
  double *ht,dtemp;
  gsl_matrix *dht=0;
  gsl_matrix *cov_theta_gsl=0, *cov_theta_times_dy=0, *cov_y=0;
  gsl_permutation *hess_theta_perm=0;
  gsl_vector_view lu_column, dht_column;
  int signum=0; /* used by LU decomp */
  double *Z;

  Z = malloc(Nt*Nx*sizeof(double));

  stats_Z(Nt,Nx,X,Ny,Y,xbar,ybar,sx,sy,Z,1);  /* just compute Z */

  ht = calloc(Ny,sizeof(double)); /* ht holds h(t,:) for x(t,:) */

  if ((dY_flag & 3) > 0) {
    Ntheta = theta_Ntheta(Nx,Nh,Ny);
    dht = gsl_matrix_calloc(Ntheta,Ny); /* dht holds dh/dtheta for x(t,:) */
    cov_theta_times_dy = gsl_matrix_calloc(Ntheta,Ny); /* temporary matrix: holds cov(theta)*dh */
    cov_y = gsl_matrix_alloc(Ny,Ny);
    cov_theta_gsl = gsl_matrix_alloc(Ntheta,Ntheta);
    for (i=1; i <= Ntheta; i++) {
      for (j=1; j <= Ntheta; j++) {
	gsl_matrix_set(cov_theta_gsl,i-1,j-1,theta_cov[ss2(Ntheta,Ntheta,i,j)]);
      }
    }

    if ((dY_flag & 4) == 4) {
      hess_theta_perm = gsl_permutation_alloc(Ntheta);
      gsl_linalg_LU_decomp(cov_theta_gsl,hess_theta_perm,&signum);
      /* now cov_theta_gsl holds LU decomp of hess_theta */
    }
  }

  for (t=1; t <= Nt; t++) {

    make_ht(Nt,Nx,Z,Nh,theta,Ny,t,ybar,sy,ht,dht,0); /* get h(t,:) and its grad with theta */


    /* copy h(t,:) into Y(t,:) */
    for (k=1; k <= Ny; k++) {
      Y[ss2(Nt,Ny,t,k)] = ht[ss1(Ny,k)];
    }

    if ((dY_flag & 3) > 0) {
      
      /* compute cov(yt) = dht'*cov_theta*dht */

      if ((dY_flag & 4) == 0) {
      
	/*  cov_theta_times_dy = cov_theta*dht */
	/* cov_theta is symmetric */
	/* gsl_blas_dsymm(CblasLeft,CblasUpper,1.0,cov_theta_gsl,dht,0.0,cov_theta_times_dy); */
	gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,cov_theta_gsl,dht,0.0,cov_theta_times_dy);
      } else /* dY_flag & 4 == 4, cov_theta_gsl is really hess_theta */ {
	/* cov_theta_times_dy = cov_theta*dht */
	/* hess_theta*cov_theta_times_dy = dht */
	/* already have LU decomp of hess_theta, so use LU solver */
	for (k=1; k <= Ny; k++) { /* do one column at a time */
	  /* access columns by view */
	  dht_column = gsl_matrix_column(dht,k-1); /* column ought not be modified */
	  lu_column = gsl_matrix_column(cov_theta_times_dy,k-1); /* column will be modified */
	  gsl_linalg_LU_solve(cov_theta_gsl,hess_theta_perm,&(dht_column.vector),&(lu_column.vector));
	}
      }
	
      /* cov_y = dy' * cov_theta_times_dy */
      gsl_blas_dgemm(CblasTrans,CblasNoTrans,1.0,dht,cov_theta_times_dy,0.0,cov_y);
      
      if ((dY_flag & 3) ==1) {
	/* compute standard erros from diagonal of covariance matrix */
	for (k=1; k <= Ny; k++) {
	  dY[ss2(Nt,Ny,t,k)] = sqrt(gsl_matrix_get(cov_y,k-1,k-1));
	}
      } else /* dY_flag & 3 == 2 */ {
	/* copy covariance matrix to dY */
	for (i=1; i <= Ny; i++) {
	  dY[ss3(Nt,Ny,Ny,t,i,i)] = gsl_matrix_get(cov_y,i-1,i-1);
	  for (j=1; j < i; j++) {
	    /* force symmetry */
	    dtemp = (gsl_matrix_get(cov_y,i-1,j-1) + gsl_matrix_get(cov_y,j-1,i-1))/2.0; 
	    dY[ss3(Nt,Ny,Ny,t,i,j)] = dtemp; 
	    dY[ss3(Nt,Ny,Ny,t,j,i)] = dtemp; 
	  }
	}
      }
    }
  }
  
  free(ht);
  free(Z);
  if ((dY_flag & 3) >0) {
      gsl_matrix_free(dht);
      gsl_matrix_free(cov_theta_gsl);
      gsl_matrix_free(cov_theta_times_dy);
      gsl_matrix_free(cov_y);
      if ((dY_flag & 4) == 4) {
	gsl_permutation_free(hess_theta_perm);
      }
  }

}

double nn_fit_eval_noreg(const unsigned long int Nt, const unsigned long int Nx, const double *Z, 
			 const unsigned long int Nh, 
			 const unsigned long int Ny, const double *Y, 
			 const double *ybar, const double *sy,
			 const double *sinv, const unsigned long int sflag,
			 const double *theta,
			 gsl_vector *grad, gsl_matrix *hess) {
  /* computes ell, grad, hess without regularization */
  /* returns ell ~ -log(p) */
  /* sflag has same meaning as in nnlib_fit, sinv has same size */

  double ell=0; /* -log(L) */
  double temp,si; /* scratch for single value from sinv */
  double *et,*ht,*ddht=0;
  gsl_matrix *dht=0;
  unsigned long int t,k,m,np,npp,Ntheta;

  ht = calloc(Ny,sizeof(double)); /* ht holds h(t,:) for x(t,:) */
  et = calloc(Ny,sizeof(double)); /* et holds h(t,:)-Y(t,:) for x(t,:) */

  Ntheta = theta_Ntheta(Nx,Nh,Ny);

  if (grad) {
    gsl_vector_set_zero(grad);
    dht = gsl_matrix_alloc(Ntheta,Ny);
    if (hess) {
      gsl_matrix_set_zero(hess);
      ddht = calloc(Ntheta*Ntheta*Ny,sizeof(double));
    }
  }

  for (t=1; t <= Nt; t++) {
    /* compute ht and deriv, hess */
    make_ht(Nt,Nx,Z,Nh,theta,Ny,t,ybar,sy,ht,dht,ddht);

    /* compute et */
    for (k=1; k <= Ny; k++) {
      temp = Y[ss2(Nt,Ny,t,k)];
      if (gsl_finite(temp)) {
	et[ss1(Ny,k)] = ht[ss1(Ny,k)]-temp;
      } else {
	et[ss1(Ny,k)] = 0; /* no error for non-finite y */
      }
    } /* end for k */

    for (k=1; k <= Ny; k++) {
      for (m=1; m <= Ny; m++) {
	switch(sflag & 3) {
	case 3: /* sinv: Nt x Ny x Ny */
	  si = sinv[ss3(Nt,Ny,Ny,t,k,m)];
	  break;
	case 2: /* sinv: Nt x Ny */
	  if (k==m) {
	    si = sinv[ss2(Nt,Ny,t,k)];
	  } else {
	    si = 0;
	  } 
	  break;
	case 1: /* sinv: 1 x Ny */
	  if (k==m) {
	    si = sinv[ss1(Ny,k)];
	  } else {
	    si = 0;
	  } 
	  break;
	case 0: /* sinv: scalar */
	default:
	  si = sinv[0];
	} /* end switch */
	if (si) {
	  /* increment ell */
	  ell += si*et[ss1(Ny,k)]*et[ss1(Ny,m)];
	  
	  if (grad) {
	    /* compute grad */
	    for (np=1; np <= Ntheta; np++) {
	      temp =  si * et[ss1(Ny,k)] * gsl_matrix_get(dht,np-1,m-1);
	      gsl_vector_set(grad,np-1,gsl_vector_get(grad,np-1)+temp);

	      if (hess) {
		/* compute hess */
		for (npp=1; npp <= np; npp++) {
		  temp =  si * ( et[ss1(Ny,k)] * ddht[ss3(Ntheta,Ntheta,Ny,np,npp,m)]
				 + gsl_matrix_get(dht,np-1,k-1) * gsl_matrix_get(dht,npp-1,m-1));
		  temp += gsl_matrix_get(hess,np-1,npp-1);
		  gsl_matrix_set(hess,np-1,npp-1,temp);
		  gsl_matrix_set(hess,npp-1,np-1,temp);
		} /* end for npp */
	      } /* end for np */
	    } /* end if hess */
	  } /* end if grad */
	} /* end if si */
      } /* end for m */
    } /* end for k */
  } /* end for t */
  
  if (dht) {
    gsl_matrix_free(dht);
  }
  if (ddht) {
    free(ddht);
  }
  free(et);
  free(ht);

  ell /= 2.0;
  return(ell);

}


double nn_fit_eval(const unsigned long int Nt, const unsigned long int Nx, const double *Z, 
		   const unsigned long int Nh, 
		   const unsigned long int Ny, const double *Y, 
		   const double *ybar, const double *sy,
		   const double *sinv, const unsigned long int sflag,
		   const double *theta,
		   gsl_vector *grad, gsl_matrix *hess) {

  double ell=0; /* -log(L) */
  unsigned long int np, Ntheta;
  gsl_vector_view vecview; /* used as thetaview and, in mpi mode, gradview */

  Ntheta = theta_Ntheta(Nx,Nh,Ny);

#ifdef USEMPI

  /* this block of code relies on pre-existing slave threads to share computation of
     ell, grad, hess on pre-shared subsets */
  /* this code should only execute in the master (mpi_rank=0) */
  int mpi_size=1, mpi_slave, mpi_op_code;
  double *mpi_buff=0;
  MPI_Status status;
  gsl_matrix_view hessview;

  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);	/* get number of processes */  

  if (mpi_size>1) {
    if (grad) {
      if (hess) {
	mpi_op_code = NN_MPI_OP_HESS; /* ell, grad, hess */
      } else {
	mpi_op_code = NN_MPI_OP_GRAD; /* ell, grad */
      }
    } else {
      mpi_op_code = NN_MPI_OP_ELL; /* ell only */
    }
    
    for (mpi_slave = 1; mpi_slave < mpi_size; mpi_slave++) {
      /* send op code & theta to slaves - MPI */
      MPI_Send(&mpi_op_code,1,MPI_INT,mpi_slave,NN_MPI_TAG_OP,MPI_COMM_WORLD);
      /* send theta */
      MPI_Send(theta,Ntheta,MPI_DOUBLE,mpi_slave,NN_MPI_TAG_THETA,MPI_COMM_WORLD);
    }
  }

  /* master: do my own calculation, while others are churning away */
  ell += nn_fit_eval_noreg(Nt,Nx,Z,Nh,Ny,Y,ybar,sy,sinv,sflag,theta,grad,hess);

  if (mpi_size>1) {
    /* receive ell, grad, hess into mpi_buff */
    mpi_buff = malloc(Ntheta*Ntheta*sizeof(double)); /* max space needed Ntheta x Ntheta for hess */
    vecview = gsl_vector_view_array(mpi_buff,Ntheta); /* used as view of gradient */
    hessview = gsl_matrix_view_array(mpi_buff,Ntheta,Ntheta);
    for (mpi_slave = 1; mpi_slave < mpi_size; mpi_slave++) {
      MPI_Recv(mpi_buff,1,MPI_DOUBLE,mpi_slave,NN_MPI_TAG_ELL,MPI_COMM_WORLD,&status);
      ell += mpi_buff[0];
      
      if (grad) {
	MPI_Recv(mpi_buff,Ntheta,MPI_DOUBLE,mpi_slave,NN_MPI_TAG_GRAD,MPI_COMM_WORLD,&status);
	gsl_vector_add(grad,&(vecview.vector));
	if (hess) {
	  MPI_Recv(mpi_buff,Ntheta*Ntheta,MPI_DOUBLE,mpi_slave,NN_MPI_TAG_HESS,MPI_COMM_WORLD,&status);
	  gsl_matrix_add(hess,&(hessview.matrix));
	}
      }
    }
    free(mpi_buff);
  }
  /* end of MPI block */
#elif defined(USEOMP)
  /* OpenMP version */
  /* this code is done before parallelization */
  unsigned long int Nst; /* number of pointer steps/time in sinv */

  Nst = get_Nst(sflag,Ny);
  
  if (grad) {
    gsl_vector_set_zero(grad); /* grad is the accumulator, initialize it */
    if (hess) {
      gsl_matrix_set_zero(hess); /* hess is the accumulator, initialize it */
    }
  }
  
#pragma omp parallel 
  { /* start omp parallelization, variables local to this loop are private, all others shared */
    int num_threads = 1,tid=0;
    unsigned long int Ntsub; /* Number of time steps per sub */
    
    double ell_tmp = 0;
    unsigned long int myNt;
    gsl_vector *grad_tmp=0;
    gsl_matrix *hess_tmp=0;
    unsigned long int t0;
    
    tid = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    Ntsub = Nt/num_threads;

    if (tid == num_threads-1) {
      myNt = Nt-Ntsub*(num_threads-1);
    } else {
      myNt = Ntsub;
    }
    
    /* allocate temporary variables  */
    if (grad) {
      grad_tmp = gsl_vector_alloc(Ntheta);
      if (hess) {
	hess_tmp = gsl_matrix_alloc(Ntheta,Ntheta);
      }
    }
    
    t0 = Ntsub*tid; /* t offset to subset in Z, Y and sinv */
    /* execute */
    ell_tmp = nn_fit_eval_noreg(myNt,Nx,Z+t0*Nx,Nh,Ny,Y+t0*Ny,ybar,sy,sinv+t0*Nst,sflag,theta,grad_tmp,hess_tmp);
    /* sum  */
#pragma omp critical      
    { /* only one process can do this block at a time */
      ell += ell_tmp;
      if (grad) {
	gsl_vector_add(grad,grad_tmp);
	if (hess) {
	  gsl_matrix_add(hess,hess_tmp);
	}
      }
    } /* end of critical */
    /* free temporary variables */
    if (grad_tmp) {
      gsl_vector_free(grad_tmp);
      if (hess_tmp) {
	gsl_matrix_free(hess_tmp);
      }
    }
  } /* end of omp parallel */
  /* end of OMP block */
#else
  /* single-processor/serial version */
  ell += nn_fit_eval_noreg(Nt,Nx,Z,Nh,Ny,Y,ybar,sy,sinv,sflag,theta,grad,hess);
#endif



  /* now do regularization criterion */
  for (np = 1; np <= Ntheta; np++) {
    ell += gsl_pow_2(theta[ss1(Ntheta,np)])/2;
  }

  if (grad) {
    vecview = gsl_vector_view_array(theta,Ntheta); /* temporarily holds gsl view of theta */
    gsl_vector_add(grad,&(vecview.vector)); /* add theta to grad */
    if (hess) {
      for (np = 1; np <= Ntheta; np++) {
	/* add identy to hess */
	gsl_matrix_set(hess,np-1,np-1,1.0+gsl_matrix_get(hess,np-1,np-1));
      }
    }
  }
  
  return(ell);
}

double nnlib_ell(const unsigned long int Nt, const unsigned long int Nx, const double *X, 
			 const unsigned long int Nh, 
			 const unsigned long int Ny, const double *Y, 
			 const double *xbar, const double *ybar, 
			 const double *sx, const double *sy,
			 const double *s, const unsigned long int sflag,
			 const double *theta,
			 double *grad, double *hess){
  double ell;
  unsigned long int Ntheta;
  gsl_vector_view grad_gsl;
  gsl_vector *gradptr;
  gsl_matrix_view hess_gsl;
  gsl_matrix *hessptr;
  double *sinv,*Z;

  Ntheta = theta_Ntheta(Nx,Nh,Ny);

  if (grad) {
    grad_gsl = gsl_vector_view_array(grad,Ntheta);
    gradptr = &(grad_gsl.vector);
  } else {
    gradptr = NULL;
  }
  if (hess) {
    hess_gsl = gsl_matrix_view_array(hess,Ntheta,Ntheta);
    hessptr = &(hess_gsl.matrix);
  } else {
    hessptr = NULL;
  }

  sinv = make_sinv(Nt,Ny,s,sflag);
  
  /* compute Z = (X-xbar)/s(x) */
  Z = malloc(Nt*Nx*sizeof(double));
  stats_Z(Nt,Nx,X,Ny,Y,xbar,ybar,sx,sy,Z,1);  /* just compute Z */

  ell = nn_fit_eval(Nt, Nx,Z,Nh,Ny,Y,ybar,sy,sinv,sflag,theta,gradptr,hessptr);
  free(Z);
  free(sinv);
  return(ell);
}


double *make_sinv(const unsigned long int Nt, const unsigned long int Ny, const double *s, const unsigned long int sflag) {
  double *sinv;
  gsl_matrix *stemp, *sinvtemp;
  gsl_permutation *perm;
  int signum=0;
  unsigned long int t,i,j;
  switch(sflag & 3) {
  case 3: /* sinv: Nt x Ny x Ny */
    sinv = malloc(Nt*Ny*Ny*sizeof(double));
    if (Ny==1) {
      for (i=1; i <= Nt; i++) {
	sinv[i-1] = 1.0/s[i-1];
      }
    } else if (Ny==2) {
      memcpy(sinv,s,Nt*Ny*Ny*sizeof(double));
      inv2x2(Nt,sinv);
    } else if (Ny==3) {
      memcpy(sinv,s,Nt*Ny*Ny*sizeof(double));
      inv3x3(Nt,sinv);
    } else {
      /* the really hard case, Ny>3 */
      /* use GSL's LU decomp machinery */
      stemp = gsl_matrix_alloc(Ny,Ny);
      perm = gsl_permutation_alloc(Ny);
      sinvtemp = gsl_matrix_alloc(Ny,Ny);
      
      for (t=1; t <= Nt; t++) {

	/* copy s(t,:,:) to stemp */
	for (i=1; i <= Ny; i++) {
	  for (j=1; j <= Ny; j++) {
	    gsl_matrix_set(stemp,i-1,j-1,s[ss3(Nt,Ny,Ny,t,i,j)]);
	  }
	}

	/* invert stemp */
	gsl_linalg_LU_decomp(stemp,perm,&signum);
	gsl_linalg_LU_invert(stemp,perm,sinvtemp);

	/* copy stemp to s(t,:,:) */
	for (i=1; i <= Ny; i++) {
	  for (j=1; j <= Ny; j++) {
	    sinv[ss3(Nt,Ny,Ny,t,i,j)] = gsl_matrix_get(sinvtemp,i-1,j-1);
	  }
	}

      }
      gsl_matrix_free(stemp);
      gsl_permutation_free(perm);
      gsl_matrix_free(sinvtemp);
    } /* end else */
    break;
  case 2: /* sinv: Nt x Ny */
    sinv = malloc(Nt*Ny*sizeof(double));
    for (i=1; i <= Nt*Ny; i++) {
      sinv[i-1] = 1.0/s[i-1];
    }
    break;
  case 1: /* sinv: 1 x Ny */
    sinv = malloc(Ny*sizeof(double));
    for (i=1; i <= Ny; i++) {
      sinv[i-1] = 1.0/s[i-1];
    }
    break;
  case 0: /* sinv: scalar */
  default:
    sinv = malloc(sizeof(double));
    sinv[0] = 1.0/s[0];
  } /* end switch */
  return(sinv);
}

typedef struct {
  unsigned long int Nt,Nx,Nh,Ny,sflag;
  const double *Z,*Y,*sinv;
  const double *ybar,*sy;
} nn_opt_params_type;


double nn_opt_f(const gsl_vector *x,void *params) {
  double fval;
  const nn_opt_params_type *par;
  double *theta;
  unsigned long int i;

  theta = malloc((x->size)*sizeof(double));
  for (i=1; i<= x->size; i++) {
    theta[i-1] = gsl_vector_get(x,i-1);
  }

  par = (const nn_opt_params_type *)params;
  fval = nn_fit_eval(par->Nt,par->Nx,par->Z,par->Nh,par->Ny,par->Y,par->ybar,par->sy,par->sinv,par->sflag,theta,0,0);

  free(theta);

  return(fval);
}

void nn_opt_df(const gsl_vector * x, void * params, gsl_vector * g) {
  const nn_opt_params_type *par;
  double *theta;
  unsigned long int i;

  theta = malloc((x->size)*sizeof(double));
  for (i=1; i<= x->size; i++) {
    theta[i-1] = gsl_vector_get(x,i-1);
  }

  par = (const nn_opt_params_type *)params;
  nn_fit_eval(par->Nt,par->Nx,par->Z,par->Nh,par->Ny,par->Y,par->ybar,par->sy,par->sinv,par->sflag,theta,g,0);

  free(theta);
}

void nn_opt_fdf(const gsl_vector * x, void * params, double * f, gsl_vector * g) {
  const nn_opt_params_type *par;
  double *theta;
  unsigned long int i;

  theta = malloc((x->size)*sizeof(double));
  for (i=1; i<= x->size; i++) {
    theta[i-1] = gsl_vector_get(x,i-1);
  }

  par = (const nn_opt_params_type *)params;
  *f = nn_fit_eval(par->Nt,par->Nx,par->Z,par->Nh,par->Ny,par->Y,par->ybar,par->sy,par->sinv,par->sflag,theta,g,0);

  free(theta);
}



double nnlib_fit(const unsigned long int Nt, const unsigned long int Nx, const double *X, 
		 const unsigned long int Nh, 
		 const unsigned long int Ny, const double *Y, 
		 const double *s, const unsigned long int flag, 
		 const unsigned long int MaxIter, 
		 double epsabs,
		 double *theta, 
		 double *xbar, double *ybar, double *sx, double *sy,
		 double *theta_cov) {
  /* sets theta, xbar, ybar, sx, sy, and optionally theta_cov */
/* 
   flag & 3 = 0: for 1x1 s
   flag & 3 = 1: for 1 x Ny s
   flag & 3 = 2: for Nt x Ny s
   flag & 3 = 3: for Nt x Ny x Ny s
   flag & 8 = 0: quiet
   flag & 8 = 8: verbose
   flag & 16 = 0: initialize at random point
   flag & 16 = 16: start at initial theta
   flag & 96 = 0: use vector_bfgs minimizer
   flag & 96 = 32: use conjugate_fr minimizer
   flag & 96 = 64: use conjugate_pr minimizer
   flag & 96 = 96: use nelder-mead simplex
   flag & 132 = 0: don't compute theta_cov or theta_hess
   flag & 132 = 4: compute theta_cov
   flag & 132 = 128: compute theta_hess
   flag & 256 = 256: don't overwrite xbar, ybar, sx, sy with values from X

   OUTPUTS: theta, theta_cov
   returns ell = neglogp
 */

  unsigned long int i,j,Ntheta;
  double *Z;
  double fval;
  double *sinv=0;
  gsl_permutation *perm;
  int signum=0;
  gsl_vector *grad;
  gsl_matrix *temp,*cov_theta_gsl;
  nn_opt_params_type par;
  size_t iter = 0;
  int status;
  double size; /* size of simplex */
  double dtemp; /* temporary double */
  const gsl_multimin_fdfminimizer_type *Tfdf;
  const gsl_multimin_fminimizer_type *Tf;
  gsl_multimin_fminimizer *s_f;
  gsl_multimin_fdfminimizer *s_fdf;
  gsl_vector *x; /* x is theta */
  gsl_vector *ss=0; /* step size for simplex */
  gsl_multimin_function my_f;
  gsl_multimin_function_fdf my_fdf;
  unsigned long int soOften = 1;  /* this sets the print status frequency initially */
  double initial_step_size = 1e-4;
  double linesearch_tol = epsabs;
 
  if (flag & 8) {
    printf("%s:Making sinv\n",__func__);
  }
  sinv = make_sinv(Nt,Ny,s,flag);

  if (flag & 8) {
    printf("%s:Making Z\n",__func__);
  }

  Z = malloc(Nt*Nx*sizeof(double));

  stats_Z(Nt,Nx,X,Ny,Y,xbar,ybar,sx,sy,Z,flag & 256); /* initialize xbar, ybar, sx, sy, & Z */

  par.Nt = Nt;
  par.Nx = Nx;
  par.Z = Z;
  par.Nh = Nh;
  par.Ny = Ny;
  par.Y = Y;
  par.ybar = ybar;
  par.sy = sy;
  par.sinv = sinv;
  par.sflag = flag;
  
  Ntheta = theta_Ntheta(Nx,Nh,Ny);

  x = gsl_vector_alloc(Ntheta);
  if (flag & 16) {
    for (i=1; i <= Ntheta; i++) {
      gsl_vector_set(x,i-1,theta[i-1]);
    }
  } else {
    /* Start at random point, in [-0.5 0.5]^n*/
    if (flag & 8) {
      printf("%s: Initializing weights (theta) - random\n",__func__);
    }
    srand( time(NULL) ); /* seed random number generator */
    for (i=1; i <= Ntheta; i++) {
      gsl_vector_set(x,i-1,((double)rand()/RAND_MAX)-0.5);
    }
  }

  switch(flag & NN_OPT_BITS) {
  case NN_OPT_FR:
    Tfdf = gsl_multimin_fdfminimizer_conjugate_fr;
    break;
  case NN_OPT_PR:
    Tfdf = gsl_multimin_fdfminimizer_conjugate_pr;
    break;
  case NN_OPT_NM:
    Tf = gsl_multimin_fminimizer_nmsimplex;
    break;
  default:
  case NN_OPT_BFGS:
#ifdef OLDGSL
    Tfdf = gsl_multimin_fdfminimizer_vector_bfgs; /* old version of GSL has no bfgs2 */
#else
    Tfdf = gsl_multimin_fdfminimizer_vector_bfgs2;
#endif
    break;
  }

  if ((flag & NN_OPT_BITS) == NN_OPT_NM) {
    ss = gsl_vector_alloc(Ntheta);
    gsl_vector_set_all(ss,1.0);
    s_f = gsl_multimin_fminimizer_alloc (Tf, Ntheta);
    my_f.n = Ntheta;
    my_f.f = &nn_opt_f;
    my_f.params = &par;
    gsl_multimin_fminimizer_set (s_f, &my_f, x, ss);
    if (flag & 8) {
      printf("%s: Minimizing with: %s\n",__func__,gsl_multimin_fminimizer_name(s_f));
    }
  } else {
    my_fdf.n = Ntheta;
    my_fdf.f = &nn_opt_f;
    my_fdf.df = &nn_opt_df;
    my_fdf.fdf = &nn_opt_fdf;
    my_fdf.params = &par;
    s_fdf = gsl_multimin_fdfminimizer_alloc (Tfdf, Ntheta);
    gsl_multimin_fdfminimizer_set (s_fdf, &my_fdf, x, initial_step_size, linesearch_tol);
    if (flag & 8) {
      printf("%s:Minimizing with: %s\n",__func__,gsl_multimin_fdfminimizer_name(s_fdf));
    }
  }
  
  
  if (flag & 8) {
    printf("%s:Starting loop\n",__func__);
  }
  do
    {
      iter++;
      if ((flag & NN_OPT_BITS) == NN_OPT_NM) {
	status = gsl_multimin_fminimizer_iterate (s_f);
      } else {
	status = gsl_multimin_fdfminimizer_iterate (s_fdf);
      }
      
      if (status) {
	break;
      }
      
      if ((flag & NN_OPT_BITS) == NN_OPT_NM) {
	size = gsl_multimin_fminimizer_size(s_f);
	status = gsl_multimin_test_size(size,epsabs);
      } else {
	status = gsl_multimin_test_gradient (s_fdf->gradient, epsabs);
      }

      if (flag & 8) {
	if (iter >= soOften*10) {
	  soOften = (soOften < 1000) ? soOften*10 : 1000;
	}
	if ((iter % soOften) == 0) {
	  if ((flag & NN_OPT_BITS) == NN_OPT_NM) {
	    printf("%s:%lu %.14e\n",__func__,(unsigned long int)iter,s_f->fval);
	  } else {
	    printf("%s:%lu %.14e\n",__func__,(unsigned long int)iter,gsl_multimin_fdfminimizer_minimum(s_fdf));
	  }
	}      
      }

    }
  while ((status == GSL_CONTINUE) && (iter < MaxIter));

  if (flag & 8) {
    printf("%s: completed after %lu/%lu: %s\n",__func__,(unsigned long int)iter,(unsigned long int)MaxIter,gsl_strerror(status));
    if ((flag & NN_OPT_BITS) == NN_OPT_NM) {
      printf("%s:size = %le, epsabs = %le\n",__func__,gsl_multimin_fminimizer_size(s_f),epsabs);
    } else {
      printf("%s:gradient = %le, epsabs = %le\n",__func__,gsl_blas_dnrm2(s_fdf->gradient),epsabs);
    }
  }

  /* copy s_fdf->x into theta */
  if ((flag & NN_OPT_BITS) == NN_OPT_NM) {
    gsl_vector_memcpy(x,s_f->x);
    fval = s_f->fval;
  } else {
    gsl_vector_memcpy(x,s_fdf->x);
    fval = gsl_multimin_fdfminimizer_minimum(s_fdf);
  }
  for (i=1; i <= Ntheta; i++) {
    theta[i-1] = gsl_vector_get(x,i-1);
  }
  

  if (flag & 8) {
    printf("%s:Final penalty function value: %lg, after %lu iterations\n",__func__,fval,(unsigned long int)iter);
  }

  if ((flag & NN_OPT_BITS) == NN_OPT_NM) {
    gsl_multimin_fminimizer_free (s_f);
    free(ss);
  } else {
    gsl_multimin_fdfminimizer_free (s_fdf);
  }
  gsl_vector_free (x);

  if (flag & 132) {
    /* compute theta_cov or theta_hess */

    if (flag & 8) {
      printf("%s:Computing theta_cov/hess\n",__func__);
    }

    grad = gsl_vector_alloc(Ntheta);
    cov_theta_gsl = gsl_matrix_alloc(Ntheta,Ntheta);

    if ((flag & 132) == 4) { /* theta cov */
      temp = gsl_matrix_alloc(Ntheta,Ntheta);
      nn_fit_eval(Nt,Nx,Z,Nh,Ny,Y,ybar,sy,sinv,flag,theta,grad,temp);
      perm = gsl_permutation_alloc(Ntheta);
      /* use GSL's LU decomp machinery */
      gsl_linalg_LU_decomp(temp,perm,&signum);
      gsl_linalg_LU_invert(temp,perm,cov_theta_gsl);
      gsl_matrix_free(temp);
      gsl_permutation_free(perm);
    } else { /* theta hess */
      nn_fit_eval(Nt,Nx,Z,Nh,Ny,Y,ybar,sy,sinv,flag,theta,grad,cov_theta_gsl);
    }


    /* copy result into theta_cov */
    for (i=1; i <= Ntheta; i++) {
      theta_cov[ss2(Ntheta,Ntheta,i,i)] = gsl_matrix_get(cov_theta_gsl,i-1,i-1);
      for (j=1; j < i; j++) {
	/* force symmetry */
	dtemp = (gsl_matrix_get(cov_theta_gsl,i-1,j-1) + gsl_matrix_get(cov_theta_gsl,j-1,i-1))/2.0;
	theta_cov[ss2(Ntheta,Ntheta,i,j)] = dtemp;
	theta_cov[ss2(Ntheta,Ntheta,j,i)] = dtemp;
      }
    }
  
    gsl_vector_free(grad);
    gsl_matrix_free(cov_theta_gsl);
    
    }

  if (Z) {
    free(Z);
  }
  if (sinv) {
    free(sinv);
  }
  if (flag & 8) {
    printf("%s:Done fitting\n",__func__);
  }

  return(fval);
}


void theta_id(const unsigned long int Nx, const unsigned long int Nh, const unsigned long int Ny, const unsigned long int n, char *buf) {
  unsigned long int Nw, Nv,ntemp,i,j;
  Nw = Nx*Nh;
  Nv = Ny*Nh;

  if (n <= Nw) {
    ntemp = n-1;
    i = 1+(ntemp % Nx);
    j = 1+(ntemp / Nx);
    sprintf(buf,"w(%lu,%lu)",i,j);
  } else if (n-Nw <= Nv) {
    ntemp = n-1;
    i = 1+(ntemp % Nh);
    j = 1+(ntemp / Nh);
    sprintf(buf,"v(%lu,%lu)",i,j);
  } else if (n-Nw-Nv <= Nh) {
    sprintf(buf,"w0(%lu)",n-Nw-Nv);
  } else if (n-Nw-Nv-Nh <= Ny) {
    sprintf(buf,"v0(%lu)",n-Nw-Nv-Nh);
  } else {
    sprintf(buf,"%s: ERROR %lu>%lu\n",__func__,n,theta_Ntheta(Nx,Nh,Ny));
  }
}

net_type init_net() {
  net_type net;
  net.Nx = 0;
  net.Nh = 0;
  net.Ny = 0;
  net.cov_flag = 0;
  net.theta = 0;
  net.theta_cov = 0;
  net.xbar = 0;
  net.ybar = 0;
  net.sx = 0;
  net.sy = 0;
  return(net);
}

void free_net(net_type *net) {
  if (net->theta) {
    free(net->theta);
  }
  if (net->theta_cov) {
    free(net->theta_cov);
  }
  if (net->xbar) {
    free(net->xbar);
  }
  if (net->ybar) {
    free(net->ybar);
  }
  if (net->sy) {
    free(net->sx);
  }
  if (net->sy) {
    free(net->sy);
  }
  *net = init_net();
}

net_type alloc_net(unsigned long int Nx,unsigned long int Nh,unsigned long int Ny,unsigned long int cov_flag) {
  net_type net;
  unsigned long int Ntheta;

  net = init_net();
  net.Nx = Nx;
  net.Nh = Nh;
  net.Ny = Ny;
  net.cov_flag = cov_flag;
  Ntheta = theta_Ntheta(Nx,Nh,Ny);

  net.theta = malloc(Ntheta*sizeof(double));
  net.xbar = malloc(Nx*sizeof(double));
  net.ybar = malloc(Ny*sizeof(double));
  net.sx = malloc(Nx*sizeof(double));
  net.sy = malloc(Ny*sizeof(double));
  if (cov_flag) {
    net.theta_cov = malloc(Ntheta*Ntheta*sizeof(double));
  }
  return(net);
}

net_type load_net(const char *filename) {
  net_type net;
  /* load a net */
  net.cov_flag = 99999;
  nnlib_load_net(filename,&(net.Nx),&(net.Nh),&(net.Ny),net.theta,net.xbar,net.ybar,net.sx,net.sy,&(net.cov_flag),net.theta_cov);
  net = alloc_net(net.Nx,net.Nh,net.Ny,net.cov_flag);
  nnlib_load_net(filename,&(net.Nx),&(net.Nh),&(net.Ny),net.theta,net.xbar,net.ybar,net.sx,net.sy,&(net.cov_flag),net.theta_cov);
  return(net);
}

void save_net(const char *filename,const net_type net) {
  nnlib_save_net(filename,net.Nx,net.Nh,net.Ny,net.theta,net.xbar,net.ybar,net.sx,net.sy,net.cov_flag,net.theta_cov);
}

set_type init_set() {
  set_type set;
  set.Nt = 0;
  set.Nx = 0;
  set.Ny = 0;
  set.sflag = 0;
  set.X = 0;
  set.Y = 0;
  set.s = 0;
  return(set);
}

void free_set(set_type *set) {
  if (set->X) {
    free(set->X);
  }
  if (set->Y) {
    free(set->Y);
  }
  if (set->s) {
    free(set->s);
  }
  *set = init_set();
}

set_type alloc_set(unsigned long int Nt,unsigned long int Nx,unsigned long int Ny,unsigned long int sflag) {
  set_type set;
  unsigned long int Ns;
  set = init_set();
  set.Nt = Nt;
  set.Nx = Nx;
  set.Ny = Ny;
  set.sflag = sflag;

  set.X = malloc(Nt*Nx*sizeof(double));
  set.Y = malloc(Nt*Ny*sizeof(double));

  Ns = get_Ns(set.sflag,set.Nt,set.Ny);
  set.s = malloc(Ns*sizeof(double));
  return(set);
}

set_type load_set(const char *filename) {
  /* load a training set */
  set_type set;
  set = init_set();
  set.sflag = 99999;
  nnlib_load_training_set(filename,&(set.Nt),&(set.Nx),set.X,&(set.Ny),set.Y,set.s,&(set.sflag));
  set = alloc_set(set.Nt,set.Nx,set.Ny,set.sflag);
  nnlib_load_training_set(filename,&(set.Nt),&(set.Nx),set.X,&(set.Ny),set.Y,set.s,&(set.sflag));
  return(set);
}

void save_set(const char *filename,const set_type set) {
  nnlib_save_training_set(filename,set.Nt,set.Nx,set.X,set.Ny,set.Y,set.s,set.sflag);
}

void nn_mpi_slave() {
#ifdef USEMPI

  int mpi_rank=-1, mpi_op_code, initialized=0;
  net_type net;
  set_type subset;
  MPI_Status status;
  double *gradbuf=0,*hessbuf=0, *Z=0, ell=0, *sinv=0;
  gsl_vector_view grad_view;
  gsl_matrix_view hess_view;
  gsl_vector *gradptr;
  gsl_matrix *hessptr;
  unsigned long int Ns,Ntheta;
  /* local sub-function, operates on local nn_mpi_slave's variables, doesn't work on coto-c1's mpicc */
    void nn_mpi_slave_free(void) {
      if (initialized) {
        free(Z);
        free(sinv);
	free(gradbuf);
	free(hessbuf);
	free_net(&net);
	free_set(&subset);
      }
      initialized = 0;
    } 

  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);	/* get current process id */  
  
  /* wait & do calculations on trial nets */
  
  do {

    /* wait/receive mpi op code */
    MPI_Recv(&mpi_op_code,1,MPI_INT,0,NN_MPI_TAG_OP,MPI_COMM_WORLD,&status);
    
    if (mpi_op_code == NN_MPI_OP_INIT) {
      
      printf("%s[%d]: op code %i received - initializing subset, net\n",__func__,mpi_rank,mpi_op_code);

      nn_mpi_slave_free(); /* free memory if needed */
      
      /* receive training subset */
      MPI_Recv(&(subset.Nt),1,MPI_UNSIGNED_LONG,0,NN_MPI_TAG_NT,MPI_COMM_WORLD,&status);
      MPI_Recv(&(subset.Nx),1,MPI_UNSIGNED_LONG,0,NN_MPI_TAG_NX,MPI_COMM_WORLD,&status);
      MPI_Recv(&(subset.Ny),1,MPI_UNSIGNED_LONG,0,NN_MPI_TAG_NY,MPI_COMM_WORLD,&status);
      MPI_Recv(&(subset.sflag),1,MPI_UNSIGNED_LONG,0,NN_MPI_TAG_SFLAG,MPI_COMM_WORLD,&status);
      printf("%s[%d]: Slave receiving Nt=%li,Nx=%li,Ny=%li,sflag=%li\n",__func__,mpi_rank,subset.Nt,subset.Nx,subset.Ny,subset.sflag);
      subset = alloc_set(subset.Nt,subset.Nx,subset.Ny,subset.sflag);
      MPI_Recv(subset.X,subset.Nt*subset.Nx,MPI_DOUBLE,0,NN_MPI_TAG_X,MPI_COMM_WORLD,&status);
      MPI_Recv(subset.Y,subset.Nt*subset.Ny,MPI_DOUBLE,0,NN_MPI_TAG_Y,MPI_COMM_WORLD,&status);
      Ns = get_Ns(subset.sflag,subset.Nt,subset.Ny);
      MPI_Recv(subset.s,Ns,MPI_DOUBLE,0,NN_MPI_TAG_S,MPI_COMM_WORLD,&status);
      printf("%s[%d]: Slave received X, Y, s\n",__func__,mpi_rank);
      
      /* receive nH */
      MPI_Recv(&(net.Nh),1,MPI_UNSIGNED_LONG,0,NN_MPI_TAG_NH,MPI_COMM_WORLD,&status);
      /* alloc net */
      net = alloc_net(subset.Nx,net.Nh,subset.Ny,0); /* alloc net, don't need space for cov */
      /* receive xbar--sy */
      MPI_Recv(net.xbar,net.Nx,MPI_DOUBLE,0,NN_MPI_TAG_XBAR,MPI_COMM_WORLD,&status);
      MPI_Recv(net.ybar,net.Ny,MPI_DOUBLE,0,NN_MPI_TAG_YBAR,MPI_COMM_WORLD,&status);
      MPI_Recv(net.sx,net.Nx,MPI_DOUBLE,0,NN_MPI_TAG_SX,MPI_COMM_WORLD,&status);
      MPI_Recv(net.sy,net.Ny,MPI_DOUBLE,0,NN_MPI_TAG_SY,MPI_COMM_WORLD,&status);
      
      printf("%s[%d]: Slave received xbar,ybar,sx,sy\n",__func__,mpi_rank);
      Ntheta = theta_Ntheta(net.Nx,net.Nh,net.Ny);
      
      /* compute Z = (X-xbar)/s(x) */
      Z = malloc(subset.Nt*subset.Nx*sizeof(double));
      stats_Z(subset.Nt,subset.Nx,subset.X,subset.Ny,subset.Y,net.xbar,net.ybar,net.sx,net.sy,Z,1);  /* just compute Z */
      
      /* invert s */
      sinv = make_sinv(subset.Nt,subset.Ny,subset.s,subset.sflag);
      
      /* allocate space for grad and hess (do it this way to ensure minimum contiguous blocks) */
      gradbuf = malloc(Ntheta*sizeof(double));
      hessbuf = malloc(Ntheta*Ntheta*sizeof(double));
      
      /* create gsl views */
      grad_view = gsl_vector_view_array(gradbuf,Ntheta);
      hess_view = gsl_matrix_view_array(hessbuf,Ntheta,Ntheta);
      
      printf("%s[%d]: Slave received subset and allocated net and buffers\n",__func__,mpi_rank);
      initialized = 1;
      
    } else if ((mpi_op_code >= NN_MPI_OP_ELL) && (mpi_op_code <= NN_MPI_OP_HESS)) {
      /* receive xbar, ybar, sx, sy, theta */
      MPI_Recv(net.theta,Ntheta,MPI_DOUBLE,0,NN_MPI_TAG_THETA,MPI_COMM_WORLD,&status);
      /* check mpi_op_code. Set NULL for un-needed results from nn_fit_eval_noreg */
      if (mpi_op_code>=NN_MPI_OP_GRAD) {
	gradptr = &(grad_view.vector);
      } else {
	gradptr = NULL;
      }
      if (mpi_op_code>=NN_MPI_OP_HESS) {
	hessptr = &(hess_view.matrix);
      } else {
	hessptr = NULL;
      }
      /* evaluate net on subset */
      ell = nn_fit_eval_noreg(subset.Nt,subset.Nx,Z,net.Nh,subset.Ny,subset.Y,net.ybar,net.sy,sinv,subset.sflag,net.theta,gradptr,hessptr);
      
      /* send ell, grad, hess */
      MPI_Send(&ell,1,MPI_DOUBLE,0,NN_MPI_TAG_ELL,MPI_COMM_WORLD); /* send ell */
      if (mpi_op_code>=NN_MPI_OP_GRAD) { /* send grad */
	MPI_Send(gradbuf,Ntheta,MPI_DOUBLE,0,NN_MPI_TAG_GRAD,MPI_COMM_WORLD);
	if (mpi_op_code>=NN_MPI_OP_HESS) { /* send hess */
	  MPI_Send(hessbuf,Ntheta*Ntheta,MPI_DOUBLE,0,NN_MPI_TAG_HESS,MPI_COMM_WORLD);
	}
      }
    } /* else ignore or exit */
    
  } while (mpi_op_code != NN_MPI_OP_EXIT); /* end slave wait/process loop */
  
  /* done, master instructed slave to end */
  
  printf("%s[%d]: exiting \n",__func__,mpi_rank);
  nn_mpi_slave_free(); /* free memory if needed */
  MPI_Finalize();

#endif
}

int nn_mpi_init(int *argc,char **argv[]) {
  int mpi_rank=0;
#ifdef USEMPI
  MPI_Init (argc, argv);	/* starts MPI */  
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);	/* get current process id */  
  printf("%s: MPI initialized, rank=%d\n",__func__,mpi_rank);
#endif
  return(mpi_rank);
}

set_type nn_mpi_dispatch(set_type set, net_type net) {
  /* create slaves if mpi enambled, return subset owned by this process, computes net.xbar-sy */
  set_type subset;

#ifdef USEMPI

  /* pre-compute xbar, sx, ybar, sy from whole set, needed for mpi case */
  stats_Z(set.Nt,set.Nx,set.X,set.Ny,set.Y,net.xbar,net.ybar,net.sx,net.sy,NULL,0); 
  set.sflag |= 256; /* indicates net has precomputed xbar,ybar,sx,sy, will be returned in subset.sflag  */

  int mpi_slave,mpi_size=1, mpi_op_code=NN_MPI_OP_EXIT;
  unsigned long int Ns,subsetNt;
  
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);	/* get number of processes */  

  if (mpi_size>1) {
    subset = set;
    subsetNt = set.Nt/mpi_size;
    for (mpi_slave=1; mpi_slave < mpi_size; mpi_slave++) {
      subset.X += set.Nx*subsetNt;
      subset.Y += set.Ny*subsetNt;

      subset.s += get_Nst(set.sflag,set.Ny)*subsetNt;
      
      if (mpi_slave==mpi_size-1) { /* special case for last one, limit Nt to what's left */
	subset.Nt = set.Nt-subsetNt*mpi_slave;
      } else {
	subset.Nt = subsetNt;
      }
      
      /* send init op code */
      mpi_op_code = NN_MPI_OP_INIT;
      MPI_Send(&mpi_op_code,1,MPI_INT,mpi_slave,NN_MPI_TAG_OP,MPI_COMM_WORLD);
      
      /* send subset */
      MPI_Send(&(subset.Nt),1,MPI_UNSIGNED_LONG,mpi_slave,NN_MPI_TAG_NT,MPI_COMM_WORLD);
      MPI_Send(&(subset.Nx),1,MPI_UNSIGNED_LONG,mpi_slave,NN_MPI_TAG_NX,MPI_COMM_WORLD);
      MPI_Send(&(subset.Ny),1,MPI_UNSIGNED_LONG,mpi_slave,NN_MPI_TAG_NY,MPI_COMM_WORLD);
      MPI_Send(&(subset.sflag),1,MPI_UNSIGNED_LONG,mpi_slave,NN_MPI_TAG_SFLAG,MPI_COMM_WORLD);
      MPI_Send(subset.X,subset.Nt*subset.Nx,MPI_DOUBLE,mpi_slave,NN_MPI_TAG_X,MPI_COMM_WORLD);
      MPI_Send(subset.Y,subset.Nt*subset.Ny,MPI_DOUBLE,mpi_slave,NN_MPI_TAG_Y,MPI_COMM_WORLD);
      Ns = get_Ns(subset.sflag,subset.Nt,subset.Ny);
      MPI_Send(subset.s,Ns,MPI_DOUBLE,mpi_slave,NN_MPI_TAG_S,MPI_COMM_WORLD);
      /* send net */
      MPI_Send(&(net.Nh),1,MPI_UNSIGNED_LONG,mpi_slave,NN_MPI_TAG_NH,MPI_COMM_WORLD);
      MPI_Send(net.xbar,net.Nx,MPI_DOUBLE,mpi_slave,NN_MPI_TAG_XBAR,MPI_COMM_WORLD);
      MPI_Send(net.ybar,net.Ny,MPI_DOUBLE,mpi_slave,NN_MPI_TAG_YBAR,MPI_COMM_WORLD);
      MPI_Send(net.sx,net.Nx,MPI_DOUBLE,mpi_slave,NN_MPI_TAG_SX,MPI_COMM_WORLD);
      MPI_Send(net.sy,net.Ny,MPI_DOUBLE,mpi_slave,NN_MPI_TAG_SY,MPI_COMM_WORLD);
      
    }
    /* now point subset at set */
    subset = set;
    subset.Nt = subsetNt;
  }
#else
  /* non-MPI case */
  /* now point subset at set */
  subset = set;
#endif
  return(subset);
}


void nn_mpi_cleanup() {
#ifdef USEMPI
  int mpi_slave, mpi_size=1, mpi_op_code;

  MPI_Comm_size (MPI_COMM_WORLD, &mpi_size);	/* get number of processes */  
  printf("%s: cleaning up %d nodes\n",__func__,mpi_size);

  if (mpi_size>1) {
    /* send op code exit to slaves, i.e., slaves die */
    mpi_op_code = NN_MPI_OP_EXIT;
    for (mpi_slave=1; mpi_slave < mpi_size; mpi_slave++) {
      MPI_Send(&mpi_op_code,1,MPI_INT,mpi_slave,NN_MPI_TAG_OP,MPI_COMM_WORLD);
    }
  }
  
  printf("%s: finalizing\n",__func__);
  MPI_Finalize();  
#endif
}
