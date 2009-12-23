#include <octave/oct.h>
#include <octave/dMatrix.h>
#include <octave/dNDArray.h>
#include <octave/uint64NDArray.h>
#include <octave/uint32NDArray.h>
#include <octave/uint16NDArray.h>
#include <octave/int64NDArray.h>
#include <octave/int32NDArray.h>
#include <octave/int16NDArray.h>
#include <string.h>
#ifdef __cplusplus
extern "C"
{
#endif
#include "kdtree.h"
#ifdef __cplusplus
}  /* end extern "C" */
#endif

DEFUN_DLD (kdtree_build_oct,args,nargout,"kdtree_buld help, see kdtree.m")
{

/*

Inputs:
0 X - matrix of sample points (Nx x Nc) (double *)
1 Nx - number of rows of X (unsigned long int)
2 Nc - number of columns of X (unsigned long int)
3 flags : (int)
  flags & 1 == 0 for ColumnMajor, 1 for RowMajor, 
  (controls majority of all input and output matrices)

Outputs:
0 root - root node of kdtree (scalar) (unsigned long int *)
1 c - c (column index) array of kdtree (Nx x 1) (unsigned short int *)
2 parent - parent index array of kdtree (Nx x 1) (unsigned long int *)
3 left - parent index array of kdtree (Nx x 1) (unsigned long int *)
4 right - parent index array of kdtree (Nx x 1) (unsigned long int *)
 (root, c, parent, left, and right are 1-based indices into the tree)
 (assume c, parent, left, right are allocated externally)

return value (1) presently not used
*/
  /* inputs */
  unsigned long int Nx, Nc;
  int flags;
  double *X;
  
  /* outputs */
  octave_value_list retval;
  unsigned long int *root, *parent, *left, *right;
  unsigned short int *c;
  
  /* internals */
  int result = 0;
  
  /* typecast inputs */
  
  /* when getting fortran_vec, must first copy
     input argument into *NDArray of right type.
     Otherwise, memory fortran_vec points to is freed after assignment
  */

  /* doubles - always 64-bit */
  NDArray X_array = args(0).array_value();
  X = X_array.fortran_vec(); /* (double *) */

#ifdef LINUX64
  // /* 64-bit long int */
  Nx = (unsigned long int)(args(1).uint64_array_value().elem(0)); /* (unsigned long int) */
  Nc = (unsigned long int)(args(2).uint64_array_value().elem(0)); /* (unsigned long int) */
  flags = (int)(args(3).int32_array_value().elem(0)); /* (int) */

#else
  // /* 32-bit long int */
  Nx = (unsigned long int)(args(1).uint32_array_value().elem(0)); /* (unsigned long int) */
  Nc = (unsigned long int)(args(2).uint32_array_value().elem(0)); /* (unsigned long int) */
  flags = (long int)(args(3).int32_array_value().elem(0)); /* (int) */ /* int and long int same on win32 */
#endif

  /* define outputs */
  dim_vector dv1(2),dvNx(2);
  dv1(0) = 1;  dv1(1) = 1;
  dvNx(0) = Nx;  dvNx(1) = 1;

#ifdef LINUX64
  // /* 64-bit long int */
  uint64NDArray root_array(dv1); /* 1 x 1 */
  uint64NDArray parent_array(dvNx); /* Nx x 1 */
  uint64NDArray left_array(dvNx); /* Nx x 1 */
  uint64NDArray right_array(dvNx); /* Nx x 1 */

#else
  // /* 32-bit long int */
  uint32NDArray root_array(dv1); /* 1x1 */
  uint32NDArray parent_array(dvNx); /* Nx x 1 */
  uint32NDArray left_array(dvNx); /* Nx x 1 */
  uint32NDArray right_array(dvNx); /* Nx x 1 */
  
#endif

  root = (unsigned long int *)(root_array.fortran_vec()); /* (unsigned long int*) */
  parent = (unsigned long int *)(parent_array.fortran_vec()); /* (unsigned long int *) */
  left = (unsigned long int *)(left_array.fortran_vec()); /* (unsigned long int *) */
  right = (unsigned long int *)(right_array.fortran_vec()); /* (unsigned long int *) */

  /* short int, always 16-bit */
  uint16NDArray c_array(dvNx);
  c = (unsigned short int *)(c_array.fortran_vec()); /* (unsigned short int *) */
  
  result = kdtree_build(X,Nx,Nc,flags,root,c,parent,left,right);

  retval(0) = root_array;
  retval(1) = c_array;
  retval(2) = parent_array;
  retval(3) = left_array;
  retval(4) = right_array;
  return retval;

}

DEFUN_DLD (kdtree_multikNN_oct,args,nargout,"kdtree_kNN help, see kdtree.m")
{

/*

Inputs:
0  const double *X - matrix of sample points (Nx x Nc)
1  const unsigned long int Nx - number of rows of X
2  const unsigned long int Nc - number of columns of X
3  const int flags :
   flags & 1 == 0 for ColumnMajor, 1 for RowMajor, 
   (controls majority of all input and output matrices)
4  const unsigned long int root - root node of kdtree
5  const unsigned short int *c - c (column index) array of kdtree (Nx x 1)
6  const unsigned long int *parent - parent index array of kdtree (Nx x 1)
7  const unsigned long int *left - parent index array of kdtree (Nx x 1)
8  const unsigned long int *right - parent index array of kdtree (Nx x 1)
  (root, c, parent, left, and right are 1-based indices into the tree)
9  const double *X0 - matrix of centers for Nearest Neighbors (NX0 x Nc)
10 const unsigned long int NX0 - number of centers
11 const unsigned long int k - number of nearest neighbrost to find at each center
12 const double *DistScale - scale multiplier for each dimension (1 x Nc)
  (distance in c'th dimension is scaled by DistScale[c],
   provide NULL to omit)

Outputs:
13 unsigned long int *index - matrix of nearest neighbor indices (NX0 x k)
14 double *R2 - distance to nearest neighbors (NX0 x k)
 (provide NULL to omit)

 return value (1) presently not used
*/
  /* inputs */
  unsigned long int Nx, Nc, root, NX0, k, *parent, *left, *right;
  int flags;
  unsigned short int *c;
  double *X, *X0, *DistScale=0;
  
  /* outputs */
  octave_value_list retval;
  unsigned long int *index;
  double *R2;
  
  /* internals */
  int result = 0;
  
  /* typecast inputs */
  
  /* when getting fortran_vec, must first copy
     input argument into *NDArray of right type.
     Otherwise, memory fortran_vec points to is freed after assignment
  */
  
  /* doubles - always 64-bit */
  NDArray X_array = args(0).array_value();
  X = X_array.fortran_vec(); /* (double *) */
  NDArray X0_array = args(9).array_value();
  X0 = X0_array.fortran_vec(); /* (double *) */
  NDArray DistScale_array = args(12).array_value();
  if (DistScale_array.nelem()) {
    DistScale = DistScale_array.fortran_vec(); /* (double *) */
  } else {
    DistScale = NULL;
  }
  
  /* dummy size for index array */
  dim_vector index_dv(2);
  index_dv(0) = 1;
  index_dv(1) = 1;

#ifdef LINUX64
  // /* 64-bit long int */
  Nx = (unsigned long int)(args(1).uint64_array_value().elem(0)); /* (unsigned long int) */
  Nc = (unsigned long int)(args(2).uint64_array_value().elem(0)); /* (unsigned long int) */
  root = (unsigned long int)(args(4).uint64_array_value().elem(0)); /* (unsigned long int) */

  uint64NDArray parent_array = args(6).uint64_array_value();
  uint64NDArray left_array = args(7).uint64_array_value();
  uint64NDArray right_array = args(8).uint64_array_value();

  NX0 = (unsigned long int)(args(10).uint64_array_value().elem(0)); /* (unsigned long int) */
  k = (unsigned long int)(args(11).uint64_array_value().elem(0)); /* (unsigned long int) */
  flags = (int)(args(3).int32_array_value().elem(0)); /* (int) */

  uint64NDArray index_array(index_dv);
#else
  // /* 32-bit long int */
  Nx = (unsigned long int)(args(1).uint32_array_value().elem(0)); /* (unsigned long int) */
  Nc = (unsigned long int)(args(2).uint32_array_value().elem(0)); /* (unsigned long int) */
  root = (unsigned long int)(args(4).uint32_array_value().elem(0)); /* (unsigned long int) */
  
  uint32NDArray parent_array = args(6).uint32_array_value();
  uint32NDArray left_array = args(7).uint32_array_value();
  uint32NDArray right_array = args(8).uint32_array_value();
  
  NX0 = (unsigned long int)(args(10).uint32_array_value().elem(0)); /* (unsigned long int) */
  k = (unsigned long int)(args(11).uint32_array_value().elem(0)); /* (unsigned long int) */
  flags = (long int)(args(3).int32_array_value().elem(0)); /* (int) */ /* int and long int same on win32 */

  uint32NDArray index_array(index_dv);
#endif


  parent = (unsigned long int *)(parent_array.fortran_vec()); /* (unsigned long int *) */
  left = (unsigned long int *)(left_array.fortran_vec()); /* (unsigned long int *) */
  right = (unsigned long int *)(right_array.fortran_vec()); /* (unsigned long int *) */

  index_dv(0) = NX0;
  index_dv(1) = k;
  index_array.resize(index_dv);
  
  /* short int, always 16-bit */
  uint16NDArray c_array = args(5).uint16_array_value();
  c = (unsigned short int *)(c_array.fortran_vec()); /* (unsigned short int *) */
  
  /* prepare outputs */

  index = (unsigned long int *)(index_array.fortran_vec());  /* (unsigned long int *) */
  
  Matrix R2_array(1,1); /* trivial memory allocation */
  if (nargout == 2) {
    R2_array.resize(NX0,k);
    R2 = R2_array.fortran_vec(); /* (double *) */
  } else {
    R2 = NULL;
  }

  /* print debug info
  octave_stdout << "Nx " << Nx << "\n";
  octave_stdout << "Nc " << Nc << "\n";
  octave_stdout << "flags " << flags << "\n";
  octave_stdout << "root " << root << "\n";
  octave_stdout << "X[0],X[1] " << X[0] << "," << X[1] << "\n";
  octave_stdout << "c[0],c[1] " << c[0] << "," << c[1] << "\n";
  octave_stdout << "parent[0],parent[1] " << parent[0] << "," << parent[1] << "\n";
  octave_stdout << "left[0],left[1] " << left[0] << "," << left[1] << "\n";
  octave_stdout << "right[0],right[1] " << right[0] << "," << right[1] << "\n";
  octave_stdout << "X0[0],X0[1] " << X0[0] << "," << X0[1] << "\n";
  octave_stdout << "NX0 " << NX0 << "\n";
  octave_stdout << "k " << k << "\n";
  if (DistScale) {
    octave_stdout << "DistScale[0],DistSCale[1] " << DistScale[0] << "," << DistScale[1] << "\n";
  } else {
    octave_stdout << "DistScale is NULL\n";
  }
  octave_stdout << "index[0],index[1] " << index[0] << "," << index[1] << "\n";
  if (R2) {
    octave_stdout << "R2[0],R2[1] " << R2[0] << "," << R2[1] << "\n";
  } else {
    octave_stdout << "R2 is NULL \n";
  }

  */
  
  result = kdtree_kNN(X,Nx,Nc,flags,root,c,parent,left,right,X0,NX0,k,DistScale,index,R2);

  retval(0) = index_array;
  if (nargout == 2) {
    retval(1) = R2_array;
  }
  
  return retval;

}


/*

cygwin/octave3.0 on win32
float 4
double 8
long double 12
char 1
short int 2
int 4
long int 4
long long int 8
 */
