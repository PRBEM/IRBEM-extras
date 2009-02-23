#include "matrix_inv.h"

void inv3x3(uint32_t Nt, double *A) {
  /* inverse in place of an array of 3 x 3 matrices */
  uint32_t t,i,j;
  uint32_t it,ib,jl,jr; /* itop,ibot,ileft, jleft, jright */
  double a,b,c,d,det2x2,det3x3,B[9];

  for (t=1; t <= Nt; t++) {
    /* copy A(t,:,:) to B */
    for (i=1; i <= 3; i++) {
      for (j=1; j <= 3; j++) {
	B[ss2(3,3,i,j)] = A[ss3(Nt,3,3,t,i,j)];
      }
    }

    /* now compute determinants in 2x2 submatrices */
    for (i=1; i <= 3; i++) {
      for (j=1; j <= 3; j++) {
	it = 1+(j==1); /* 2, 1, 1 */
	ib = 3-(j==3); /* 3, 3, 2 */
	jl = 1+((i+(j==2)) % 3); /* complicated */
	jr = 1+((i+(j!=2)) % 3); /* complicated */
	a = B[ss2(3,3,it,jl)];
	b = B[ss2(3,3,it,jr)];
	c = B[ss2(3,3,ib,jl)];
	d = B[ss2(3,3,ib,jr)];
	det2x2 = a*d-b*c;
	A[ss3(Nt,3,3,t,i,j)] = det2x2;
      }
    }

    /* compute 3x3 determinant */
    det3x3 = A[ss3(Nt,3,3,t,1,1)]*B[ss2(3,3,1,1)]
      + A[ss3(Nt,3,3,t,2,1)]*B[ss2(3,3,1,2)]
      + A[ss3(Nt,3,3,t,3,1)]*B[ss2(3,3,1,3)];
      

    /* now divide through by det3x3 */
    for (i=1; i <= 3; i++) {
      for (j=1; j <= 3; j++) {
	A[ss3(Nt,3,3,t,i,j)] /= det3x3; 
      }
    }

  }
}

void inv2x2(uint32_t Nt, double *A) {
  /* inverse in place of an array of 2 x 2 matrices */
  uint32_t t;
  double a11,a12,a21,a22, det;
  for (t=1; t <= Nt; t++) {
    a11 = A[ss3(Nt,2,2,t,1,1)];
    a12 = A[ss3(Nt,2,2,t,1,2)];
    a21 = A[ss3(Nt,2,2,t,2,1)];
    a22 = A[ss3(Nt,2,2,t,2,2)];
    det = a11*a22-a12*a21;
    A[ss3(Nt,2,2,t,1,1)] = a22/det;
    A[ss3(Nt,2,2,t,1,2)] = -a12/det;
    A[ss3(Nt,2,2,t,2,1)] = -a21/det;
    A[ss3(Nt,2,2,t,2,2)] = a11/det;
  }
}

/* test code
int main() {
  double A[4] = {0.655741,0.0357117,0.849129,0.933993};
  double B[9] = {0.765517,0.489764,0.709365,0.7952,0.445586,0.754687,0.186873,0.646313,0.276025};
  int i,j;
  inv2x2(1,A);
  for (i=1; i <= 2; i++) {
    for (j=1; j<= 2; j++) {
      printf("invA[%i,%i]=%g\n",i,j,A[ss2(2,2,i,j)]);
    }
  }
  inv3x3(1,B);
  for (i=1; i <= 3; i++) {
    for (j=1; j<= 3; j++) {
      printf("invB[%i,%i]=%g\n",i,j,B[ss2(3,3,i,j)]);
    }
  }
}
*/
