#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "kdtree.h"

#ifdef USEOMP
/* use OpenMP */
#include <omp.h>
#endif

/* 
kdtree_build build k-d tree for fast nearest neighbors.  This routine
is meant to be called as a dll from Matlab--that's why no structures.
*/

long int kdtree_build_idl( int argc, void *argv[]) {
  /* expects all arguments passed by reference */
  return(kdtree_build((double *)argv[0], /* *X */
		      *(unsigned long int *)argv[1], /* Nx */
		      *(unsigned long int *)argv[2], /* Nc */
		      *(int *)argv[3], /* flags */
		      (unsigned long int *)argv[4], /* *root */
		      (unsigned short int *)argv[5], /* *c */
		      (unsigned long int *)argv[6], /* *parent */
		      (unsigned long int *)argv[7], /* *left */
		      (unsigned long int *)argv[8])); /* *right */
}
int kdtree_build(const double *X, 
		 const unsigned long int Nx,
		 const unsigned long int Nc,
		 const int flags,
		 unsigned long int *root,
		 unsigned short int *c,
		 unsigned long int *parent,
		 unsigned long int *left,
		 unsigned long int *right) {

  const short int XmajorFlag = (flags & 1);
  
  inline unsigned long int ss2(const unsigned long int Nrows, 
				      const unsigned long int Ncols,
				      const unsigned long int i,
				      const unsigned long int j) { 
    /* 1-based 2-D subscript to 0-based 1-D subscript */
    /* XmajorFlag - 0 for ColumnMajor, 1 for RowMajor, 
       (controls majority of all input and output matrices) */
    return(XmajorFlag?(Ncols*(i-1)+(j-1)):(Nrows*(j-1)+(i-1)));
  }

  unsigned long int *new_range(const unsigned long int i1, const unsigned long int i2, unsigned long int *NI) {
    unsigned long int i,*I = malloc((i2-i1+1)*sizeof(unsigned long int));
    for (i=i1; i <= i2; i++) {
      I[i-i1] = i;
    }
    if (NI) {
      *NI = (i2-i1+1);
    }
    return(I);
  }

  unsigned long int *new_reindex(const unsigned long int *I, const unsigned long int *i1, unsigned long int Ni1) {
    /* return I(i1), where i1 is a set of Ni1 1-based indexes */
    unsigned long int i,*I1 = malloc(Ni1*sizeof(unsigned long int));
    for (i=0; i < Ni1; i++) {
      I1[i] = I[i1[i]-1];
    }
    return(I1);
  }

  unsigned long int kdtree_split(const unsigned long int *I, const unsigned long int NI,const short int c_in,const unsigned long int iparent) {

    int XIcompare(const  void *i1, const void *i2) {
      double x1 = X[ss2(Nx,Nc,*((unsigned long int*)i1),c_in)];
      double x2 = X[ss2(Nx,Nc,*((unsigned long int*)i2),c_in)];
      if (x1==x2) {
	return(0);
      } else if (x1 < x2) {
	return(-1);
      } else {
	return(+1);
      }
    }

    unsigned long int *Isorted = malloc(NI*sizeof(unsigned long int));

    /* get Isorted that sorts X(Isorted,c) */
    memcpy(Isorted,I,NI*sizeof(unsigned long int));
    qsort(Isorted,NI,sizeof(unsigned long int),&XIcompare);

    unsigned long int jmid = ceil((double)NI/2.0);
    unsigned long int leaf_index = Isorted[jmid-1];

    if (iparent==0) {
      *root = leaf_index;
    }
    c[leaf_index-1] = c_in;
    parent[leaf_index-1] = iparent;

    /* prepare to call split */
    short int c_out = 1+(c_in % Nc);

    if (jmid>1) {
      unsigned long int NI1,*i1 = new_range(1,jmid-1,&NI1);
      unsigned long int *I1 = new_reindex(Isorted,i1,NI1);
      left[leaf_index-1] = kdtree_split(I1,NI1,c_out,leaf_index);
      free(i1);
      free(I1);
    }

    if (jmid < NI) {
      unsigned long int NI2,*i2 = new_range(jmid+1,NI,&NI2);
      unsigned long int *I2 = new_reindex(Isorted,i2,NI2);
      right[leaf_index-1] = kdtree_split(I2,NI2,c_out,leaf_index);
      free(i2);
      free(I2);
    }
    free(Isorted);
    return(leaf_index);
  }
  
  *root = 0;
  /* assume c, parent, left, right are allocated already */
  /* set to zero */
  memset(c,0,Nx*sizeof(*c));
  memset(parent,0,Nx*sizeof(*parent));
  memset(left,0,Nx*sizeof(*left));
  memset(right,0,Nx*sizeof(*right));

  unsigned long int NI,*I = new_range(1,Nx,&NI);

  kdtree_split(I,NI,1,0);

  free(I);

  return(1); /* return value presently not used */
}

int kdtree_kNN_range(const unsigned long int i1,
		     const unsigned long int i2,
		     const double *X, 
		     const unsigned long int Nx,
		     const unsigned long int Nc,
		     const int flags,
		     const unsigned long int root,
		     const unsigned short int *c,
		     const unsigned long int *parent,
		     const unsigned long int *left,
		     const unsigned long int *right,
		     const double *X0,
		     const unsigned long int NX0,
		     const unsigned long int k,
		     const double *DistScale,
		     unsigned long int *index,
		     double *R2) {

  /* i1, i2 1-based indices for which X0's to do in this call (helps parallilization) */

  const short int XmajorFlag = (flags & 1);
  const double double_inf = +1./0.0;  /* positive infinity */

    inline unsigned long int ss2(const unsigned long int Nrows, 
				 const unsigned long int Ncols,
				 const unsigned long int i,
				 const unsigned long int j) { 
      /* 1-based 2-D subscript to 0-based 1-D subscript */
      /* XmajorFlag - 0 for ColumnMajor, 1 for RowMajor, 
	 (controls majority of all input and output matrices) */
      return(XmajorFlag?(Ncols*(i-1)+(j-1)):(Nrows*(j-1)+(i-1)));
    }
    
    double dist_fun1(unsigned long int iX, unsigned long int iX0,unsigned long int j) {
      /* signed, scaled distance between X(iX,:) and X0(iX0,:) in the j coordinate */
      double dist = (X[ss2(Nx,Nc,iX,j)]-X0[ss2(NX0,Nc,iX0,j)]);
      if (DistScale) {
	dist *= DistScale[j-1];
      }
      return(dist);
    }
    
    double dist_fun_sumsq(unsigned long int iX, unsigned long int iX0) {
      /* sum squared, scaled distance between X(iX,:) and X0(iX0,:) */
      unsigned long int j;
      double dist = 0;
      for (j=1; j <= Nc; j++) {
	double tmp = dist_fun1(iX,iX0,j);
	dist += tmp*tmp;
      }
      return(dist);
    }

    double dist_fun_maxsq(unsigned long int iX, unsigned long int iX0) {
      /* maximum squared, scaled distance between X(iX,:) and X0(iX0,:) */
      unsigned long int j;
      double dist = 0;
      for (j=1; j <= Nc; j++) {
	double tmp = dist_fun1(iX,iX0,j);
	tmp *= tmp; /* tmp^2 */
	if (tmp>dist) {
	  dist = tmp;
	}
      }
      return(dist);
    }

    /* flags & 2 == 0 for pythagorean distance, 2 for max absolute difference */
  double (*dist_fun)(unsigned long int, unsigned long int);
  if (flags & 2) {
    dist_fun = &dist_fun_maxsq;
  } else {
    dist_fun = &dist_fun_sumsq;
  }

    
    void check_node(unsigned long int i,
		    unsigned long int iX0,
		    short int check_parent,
		    short int check_left,
		    short int check_right) {
      /* check node and recursively check its
	 parent, left child, and/or right child
	 i is an index into the tree variables (X, parent, right, left, c, etc)
	 iX0 is an index into the center points matrix (X0)
      */
      if (i<=0) {
	return;
      }
      
      double dxc = dist_fun1(i,iX0,c[i-1]); /* signed distance along c */
      double dxc2 = dxc*dxc; /* squared distance along c */
      unsigned long int ssR2; /* subscript for iR2 */
      if (dxc2 <= R2[ssR2 = ss2(NX0,k,iX0,k)]) { /* check against farthest nearest neighbor */
	/* this node and its children could be nearest neighbors */
	double r2i = (*dist_fun)(i,iX0); /* full Nc-d squared distance */
	if (r2i <= R2[ssR2]) { /* this is an active nearest neighbor */
	  unsigned long int iR2 = k,ssR2m1; /* ssR2m1 - subscript for iR2-1 */
	  /* walk inward from farthest nearest neighbor,
	     copying distances larger than r2i outward,
	     until we find a spot for this node
	  */
	  while ((iR2 > 1) && (r2i <= R2[ssR2m1 = ss2(NX0,k,iX0,iR2-1)])) {
	    R2[ssR2] = R2[ssR2m1];
	    index[ssR2] = index[ssR2m1];
	    iR2--;
	    ssR2 = ssR2m1;
	  } /* end while */
	  /* now we have a spot for this node, store it */
	  R2[ssR2] = r2i;
	  index[ssR2] = i;
	} /* end if r2i */
	
	/* ok, check both children if allowed */
	if (check_right) {
	  check_node(right[i-1],iX0,0,1,1);
	}
	if (check_left) {
	  check_node(left[i-1],iX0,0,1,1);
	}
      } else {
	/* this node is too far away, but its kids might be closer */
	if (dxc<0) { /* only kids to the right can have smaller dxc2, and be in range */
	  if (check_right) {
	    check_node(right[i-1],iX0,0,1,1);
	  }
	} else { /* only kids to the left can have smaller dxc2, and be in range */
	  if (check_left) {
	    check_node(left[i-1],iX0,0,1,1);
	  }
	} /* end if dxc<0 */
      } /*end if dxc2 */
      if (check_parent && (parent[i-1]>0)) {
	/* check parent, but don't recheck this node */
	short int isright = (right[parent[i-1]-1] == i); /* is this a right child of parent? */
	check_node(parent[i-1],iX0,1,isright,!isright);
      }
      
    }

    unsigned long int iX0;
    for (iX0 = i1; iX0 <= i2; iX0++) {
      unsigned long int i,ilast;
      /* prepare distance array w/ infs */
      for (i=1; i <= k; i++) {
      	R2[ss2(NX0,k,iX0,i)] = double_inf;
      }
      /* descend to closest leaf node */
      i = root;
      while (i>0) {
      	ilast = i;
      	if (X0[ss2(NX0,Nc,iX0,c[i-1])] > X[ss2(Nx,Nc,i,c[i-1])]) {
      	  i = right[i-1];
      	} else {
      	  i = left[i-1];
      	}
      } /* end while */
      /* ilast now points to nearest leaf node (may have one branch in away direction */
      check_node(ilast,iX0,1,1,1); /* check node, its parent and any children */
    } /* end of for loop */
    
  return(1); /* return value presently not used */
}

	       
long int kdtree_kNN_idl( int argc, void *argv[]) {
  /* expects all arguments passed by reference */
  return(kdtree_kNN(argv[0], /* *X */
		    *(unsigned long int*)argv[1], /* Nx */
		    *(unsigned long int*)argv[2], /* Nc */
		    *(int*)argv[3], /* flags */
		    *(unsigned long int*)argv[4], /* root */
		    (unsigned short int*)argv[5], /* *c */
		    (unsigned long int*)argv[6], /* *parent */
		    (unsigned long int*)argv[7], /* *left */
		    (unsigned long int*)argv[8], /* *right */
		    (double*)argv[9], /* *X0 */
		    *(unsigned long int*)argv[10], /* NX0 */
		    *(unsigned long int*)argv[11], /* k */
		    (double*)argv[12], /* *DistScale */
		    (unsigned long int*)argv[13], /* *index */
		    (double*)argv[14])); /* *R2 */
}

int kdtree_kNN(const double *X, 
	       const unsigned long int Nx,
	       const unsigned long int Nc,
	       const int flags,
	       const unsigned long int root,
	       const unsigned short int *c,
	       const unsigned long int *parent,
	       const unsigned long int *left,
	       const unsigned long int *right,
	       const double *X0,
	       const unsigned long int NX0,
	       const unsigned long int k,
	       const double *DistScale,
	       unsigned long int *index,
	       double *R2) {
  short int freeR2=0;
  int result = 0;

  if (!R2) {
    /* allocate distance array */
    R2 = malloc(NX0*k*sizeof(double));
    freeR2 = 1; /* will need to free later */
  }

#ifdef USEOMP
#pragma omp parallel
  { /* start parallel */
    int num_threads = 1,tid=0;
    /* i1, i2 1-based indices for which X0's to do in the sub-call */
    /* Ntsub Number of time steps per sub */
    unsigned long int i1,i2,Ntsub;

    tid = omp_get_thread_num();
    num_threads = omp_get_num_threads();
    Ntsub = NX0/num_threads;

    i1 = 1+tid*Ntsub;
    if (tid == num_threads-1) {
      i2 = NX0;
    } else {
      i2 = i1+Ntsub-1;
    }
    
    kdtree_kNN_range(i1,i2,X,Nx,Nc,flags,root,c,parent,left,right,X0,NX0,k,DistScale,index,R2);
    result = 1;
  } /* end parallel */
#else
  result = kdtree_kNN_range(1,NX0,X,Nx,Nc,flags,root,c,parent,left,right,X0,NX0,k,DistScale,index,R2);
#endif

  /* free memory */
  if (freeR2) {
    free(R2);
  }

  return(result);
};


