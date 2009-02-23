/* provides load/save routines for
   neural networks and training sets */

#include <stdlib.h>
#include <stdio.h>
#include <stdint.h>
#include "nn.h"

void nnlib_save_training_set(const char *filename, const unsigned long int Nt, 
			     const unsigned long int Nx, const double *X, 
			     const unsigned long int Ny, const double *Y, 
			     const double *s, const unsigned long int flag) {

  FILE *filep;
  size_t Ns;
  uint32_t temp;

  if (!(filep = fopen(filename,"wb"))) {
    perror("main:Unable to create file");
    fprintf(stderr,"%s: file name: %s\n",__func__,filename);
    exit(-1);
  }

  temp = Nt;
  fwrite(&temp,sizeof(temp),1,filep);
  temp = Nx;
  fwrite(&temp,sizeof(temp),1,filep);
  temp = Ny;
  fwrite(&temp,sizeof(temp),1,filep);
  temp = flag & 3;
  fwrite(&temp,sizeof(temp),1,filep);
  fwrite(X,sizeof(double),Nt*Nx,filep);
  fwrite(Y,sizeof(double),Nt*Ny,filep);
  Ns = get_Ns(flag,Nt,Ny);
  fwrite(s,sizeof(double),Ns,filep);
  fclose(filep);
} 

void nnlib_load_training_set(const char *filename, unsigned long int *Nt, 
			     unsigned long int *Nx, double *X, 
			     unsigned long int *Ny, double *Y, 
			     double *s, unsigned long int *flag) {
/* !!!!! X and Y must have enough space to accommodate net 
Call with *flag = 99999 to just set sizes (Nt, Nx, Ny and flag).
Then allocate memory and call again.
*/

  FILE *filep;
  size_t Ns;
  unsigned long int inflag;
  uint32_t temp;

  if (!(filep = fopen(filename,"rb"))) {
    perror("main:Unable to open file");
    fprintf(stderr,"%s: file name: %s\n",__func__,filename);
    exit(-1);
  }

  inflag = *flag;
  fread(&temp, sizeof(uint32_t), 1, filep);
  *Nt = temp;
  fread(&temp, sizeof(uint32_t), 1, filep);
  *Nx = temp;
  fread(&temp, sizeof(uint32_t), 1, filep);
  *Ny = temp;
  fread(&temp, sizeof(uint32_t), 1, filep);
  *flag = temp;
  if (inflag != 99999) {
    fread(X,sizeof(double),(*Nt)*(*Nx),filep);
    fread(Y,sizeof(double),(*Nt)*(*Ny),filep);
    Ns = get_Ns(*flag,*Nt,*Ny);
    fread(s,sizeof(double),Ns,filep);
  }
  fclose(filep);
}


void nnlib_save_net(const char *filename, const unsigned long int Nx, 
		    const unsigned long int Nh, const unsigned long int Ny, const double *theta, 
		    const double *xbar,const double *ybar, const double *sx, const double *sy,
		    const unsigned long int cov_flag, const double *theta_cov) {
  FILE *filep;
  uint32_t temp,Ntheta;

  if (!(filep = fopen(filename,"wb"))) {
    perror("main:Unable to create file");
    fprintf(stderr,"%s: file name: %s\n",__func__,filename);
    exit(-1);
  }

  Ntheta = theta_Ntheta(Nx,Nh,Ny);

  /* file starts with a 2-uint32 header
     first uint32 is zero (spot held by Nx in old format)
     second uint32 is "type" for future reverse compatibility
   */
  temp = 0; /* indicate "new" file format */
  fwrite(&temp,sizeof(temp),1,filep);
  temp = 1; /* indicate specific type of new file format */
  fwrite(&temp,sizeof(temp),1,filep);
  temp = Nx;
  fwrite(&temp,sizeof(temp),1,filep);
  temp = Nh;
  fwrite(&temp,sizeof(temp),1,filep);
  temp = Ny;
  fwrite(&temp,sizeof(temp),1,filep);
  if (theta_cov) {
    temp = cov_flag & 3;
    } else {
      temp = 0; /* override user request because theta_cov is empty */
    }
  fwrite(&temp,sizeof(temp),1,filep);
  fwrite(theta,sizeof(double),Ntheta,filep);
  fwrite(xbar,sizeof(double),Nx,filep);
  fwrite(ybar,sizeof(double),Ny,filep);
  fwrite(sx,sizeof(double),Nx,filep);
  fwrite(sy,sizeof(double),Ny,filep);
  if (cov_flag && theta_cov) {
    fwrite(theta_cov,sizeof(double),Ntheta*Ntheta,filep);
  }
  fclose(filep);
}

void nnlib_load_net(const char *filename, unsigned long int *Nx, 
		    unsigned long int *Nh, unsigned long int *Ny, double *theta, 
		    double *xbar,double *ybar, double *sx, double *sy,
		    unsigned long int *cov_flag, double *theta_cov) {

  FILE *filep;
  unsigned long int inflag,filetype, j;
  uint32_t temp,Ntheta;

  if (! filename) {
    fprintf(stderr,"%s: filename is NULL\n",__func__);
    return;
  }

  if (! Nx) {
    fprintf(stderr,"%s: Nx is NULL\n",__func__);
    return;
  }

  if (! Nh) {
    fprintf(stderr,"%s: Nh is NULL\n",__func__);
    return;
  }

  if (! Ny) {
    fprintf(stderr,"%s: Ny is NULL\n",__func__);
    return;
  }

  if (! cov_flag) {
    fprintf(stderr,"%s: cov_flag is NULL\n",__func__);
    return;
  }

  inflag = *cov_flag;
 
  if (!(filep = fopen(filename,"rb"))) {
    perror("main:Unable to open file");
    fprintf(stderr,"%s: file name: %s\n",__func__,filename);
    exit(-1);
  }

  /* file starts with Nx or 2-uint32 header
     If first uint32 is zero then second uint32 defines type
     Otherwise just use Nx and "old" format, provide dummy for variables in "new" format
   */
  fread(&temp, sizeof(uint32_t), 1, filep);
  if (temp > 0) {
    filetype = 0;
    /* temp is Nx */
  } else {
    /* new format, read second int32 */
    fread(&temp, sizeof(uint32_t), 1, filep);
    filetype = temp;
    fread(&temp, sizeof(uint32_t), 1, filep); /* read Nx */
  }
  *Nx = temp;
  fread(&temp, sizeof(uint32_t), 1, filep);
  *Nh = temp;
  fread(&temp, sizeof(uint32_t), 1, filep);
  *Ny = temp;
  fread(&temp, sizeof(uint32_t), 1, filep);
  *cov_flag = temp;
  if (inflag != 99999) {
    if (! theta) {
      fprintf(stderr,"%s: theta is NULL for cov_flag==%li, trying to read theta\n",__func__,inflag);
      fclose(filep);
      return;
    }
    Ntheta = theta_Ntheta(*Nx,*Nh,*Ny);
    fread(theta,sizeof(double),Ntheta,filep);

    if (filetype == 0) {
      /* xbar, ybar, sx, sy dummy values */
      for (j=0; j < *Nx; j++ ) {
	xbar[j] = 0;
	sx[j] = 1;
      }
      for (j=0; j < *Ny; j++ ) {
	ybar[j] = 0;
	sy[j] = 1;
      }
    } else {
      /* xbar, ybar, sx, sy stored in file */
      fread(xbar,sizeof(double),*Nx,filep);
      fread(ybar,sizeof(double),*Ny,filep);
      fread(sx,sizeof(double),*Nx,filep);
      fread(sy,sizeof(double),*Ny,filep);
    }

    if (cov_flag) {
      if (! theta_cov) {
	fprintf(stderr,"%s: theta_cov is NULL for cov_flag==%li, trying to read theta_cov\n",__func__,inflag);
	fclose(filep);
	return;
      }
      fread(theta_cov,sizeof(double),Ntheta*Ntheta,filep);
    }
  }
  fclose(filep);
}

