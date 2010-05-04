/* 
kdtree_build build k-d tree for fast nearest neighbors.  This routine
is meant to be called as a dll from Matlab--that's why no structures.
*/

int kdtree_build(const double *X, 
		 const unsigned long int Nx,
		 const unsigned long int Nc,
		 const int flags,
		 unsigned long int *root,
		 unsigned short int *c,
		 unsigned long int *parent,
		 unsigned long int *left,
		 unsigned long int *right);
/*

Inputs:
X - matrix of sample points (Nx x Nc)
Nx - number of rows of X
Nc - number of columns of X
flags :
  flags & 1 == 0 for ColumnMajor, 1 for RowMajor, 
  (controls majority of all input and output matrices)

Outputs:
root - root node of kdtree (scalar)
c - c (column index) array of kdtree (Nx x 1)
parent - parent index array of kdtree (Nx x 1)
left - parent index array of kdtree (Nx x 1)
right - parent index array of kdtree (Nx x 1)
 (root, c, parent, left, and right are 1-based indices into the tree)
 (assume c, parent, left, right are allocated externally)

return value (1) presently not used
*/


/* 
kdtree_kNN find k nearest neighbors from a pre-constructed tree
use Matlab kdtree.m to build the tree. This routine is
meant to be called as a dll from Matlab--that's why no structures
are used and that's why there's no routine to build
the kdtree.
*/

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
	       double *R2);
/*

Inputs:
X - matrix of sample points (Nx x Nc)
Nx - number of rows of X
Nc - number of columns of X
flags :
  flags & 1 == 0 for ColumnMajor, 1 for RowMajor, 
  (controls majority of all input and output matrices)
  flags & 2 == 0 for pythagorean distance, 2 for max absolute difference,
root - root node of kdtree
c - c (column index) array of kdtree (Nx x 1)
parent - parent index array of kdtree (Nx x 1)
left - parent index array of kdtree (Nx x 1)
right - parent index array of kdtree (Nx x 1)
 (root, c, parent, left, and right are 1-based indices into the tree)
X0 - matrix of centers for Nearest Neighbors (NX0 x Nc)
NX0 - number of centers
k - number of nearest neighbrost to find at each center
DistScale - scale multiplier for each dimension (1 x Nc)
  (distance in c'th dimension is scaled by DistScale[c],
   provide NULL to omit)

Outputs:
index - matrix of nearest neighbor indices (NX0 x k)
R2 - squared distance to nearest neighbors (NX0 x k)
 (provide NULL to omit)

return value (1) presently not used
*/


/* IDL wrappers */
/* form: return_type example(int argc, void *argv[]) */
/* expects all arguments passed by reference. In IDL
   this requires an extra "VALUE" keyword to be set
   to a byte array of zeros with length argc. Strings
   should be passed as null-terminated arrays of bytes
   (this can be accomplished by converting them in IDL
   or by setting the corresponding entry in VALUE to 1).
 */

long int kdtree_build_idl( int argc, void *argv[]);
/* IDL call:
   call_extern,'kdtree.so','kdtree_build_idl_', \
   X,Nx,Nc,flags,root,c,parent,left,right, \
   value=bytearry(9))
 */
long int kdtree_kNN_idl( int argc, void *argv[]);
/* IDL call:
   call_extern,'kdtree.so','kdtree_kNN_idl_', \
   X,Nx,Nc,flags,root,c,parent,left,right, \
   X0, NX0,k,DistScale,index,R2, \
   value=bytearry(15))
 */
