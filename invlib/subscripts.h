#ifndef SUBSCRIPTS_H

/* subscripts vector, matrix, tensor */
#define ss1(N,i) ((i)-1)
#define ss2(Nr,Nc,i,j) (((i)-1)*(Nc)+(j)-1)
#define ss3(N1,N2,N3,i,j,k) ((((i)-1)*(N2)+(j)-1)*(N3)+(k)-1)

#define SUBSCRIPTS_H 1
#endif
