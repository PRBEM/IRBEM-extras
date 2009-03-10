#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_blas.h>
#include "invlib_const.h" /* defines constants, e.g., return codes */
#include "invutil.h" /* defines inv_matrix_once, vector_func */
#include "nnlib.h" /* used by TEM1 routines */
#include "optim.h" /* optimize, fzero and optfunTy */
#include "ae8_atmocutoff.h" /* AE8_AtmoCutoff() */

/* TEM1 uses AE8 atmospheric cutoff model */
#define TEM1_AtmoCutoff(Lm) AE8_AtmoCutoff(Lm)

/* holds params for likelihood function for a single wide2omni fit */
typedef struct {
  double y; /* wide-angle flux */
  gsl_vector *H; /* measurement weights for y */
  gsl_vector *mu; /* a priori value of log flux, in log-normal */
  gsl_matrix *Sigma_inv; /* inverse of flux vs angle covariance matrix */
  FILE *outFilePtr; /* output file pointer */
  int closeOutFile; /* does output file pointer need to be closed on exit? */
} w2u_TEM1_paramsTy;



/**************************
 TEM-1 routines 
***************************/

double TEM1_alpha_rho(double Lm, double keV, double EPA1deg, double EPA2deg) {
  /* compute TEM1 Gaussian covariance (e.g., for log-normal assumption)
     between EPA1 and EPA2 (equatorial pitch angles, degrees)
     at Lm and energy keV
   */
  double rho,absdA,Amid;
  /* TEM1 covariance works in degrees */
  absdA = fabs(EPA2deg-EPA1deg);
  Amid = (EPA1deg+EPA2deg)/2.0;
  /* this is the TEM-1 full covariance function 
     it gives the log of the covariance: 
-2.38045422354778740000e-004 x absdL*logEmid*logEmid
-1.94237761524617050000e-003 x absdlogE*logEmid*Lmid
-4.57775467646251600000e-010 x absdA*Emid*Emid
6.92048589352182450000e-009 x absdE*absdL*Emid
-2.77718729137489650000e-002 x absdlogE*logEmid
1.27959165482814540000e-008 x absdE*absdA*Amid
-7.25627663580207180000e-005 x absdA*Lmid*Lmid
-9.88350656479867160000e-003 x absdL*absdL*Lmid
1.08291924325609880000e-006 x absdA*absdL*Emid
-1.00412194407433640000e-008 x absdL*Emid*Emid
-2.07259507044728790000e-001 x absdL
  */
  /* We only need terms without absdE, absdlogE or absdL */
  rho = -4.57775467646251600000e-010 * absdA*keV*keV
    -7.25627663580207180000e-005 * absdA*Lm*Lm;
  rho = exp(rho);
  return(rho);
}

/* TEM1 neural network */
/* Neural network coefficients from TEM1_logflux.net */

const unsigned long int TEM1_Nx=3,TEM1_Nh=24,TEM1_Ny=2;
const double TEM1_xbar = 0, TEM1_ybar = 0, TEM1_sx = 1, TEM1_sy = 1;
const double TEM1_theta[] = {
-2.32233269714936560000e+000,-5.29141132661746380000e+000,-2.86449825922063940000e+000,-1.28394731147018230000e+000,
-8.44929290787063390000e-004,-1.19299967375961050000e+000,-3.42719877762713440000e+000,-2.79518133365517630000e+000,
-1.36504990698968650000e+001,-1.89049697096156560000e+001,-2.58884679322468040000e+000,-5.64695338815610890000e+000,
9.90034062324596230000e-001,-1.54961721651657340000e-002,-2.31138349818970430000e+000,-4.36362575121701380000e+000,
-1.47590330705035710000e+000,-2.03628784794812080000e+000,1.59556914835272500000e+000,-4.81676512106896380000e-002,
3.01329518088655310000e+000,4.71293641111188060000e+000,2.38223259056419190000e+000,4.38221731611645190000e+000,
1.42812030348646780000e+000,-2.10391851980329920000e+000,-4.37566204275786670000e+000,-3.11148607442081020000e+000,
3.41734416144306370000e+000,7.69088043241130890000e-001,7.23782244585444600000e-001,-5.08214396430565160000e+000,
-1.02357006988981460000e+001,-2.60691165297326100000e+000,-2.82115339954197930000e+000,-9.64978620839669430000e+000,
7.14419742385879890000e+000,-5.02735737843328440000e-002,3.49410894435860400000e+000,-3.89350342809173180000e+000,
-5.98440795540902130000e-001,2.47096686812677800000e-001,-1.07157822420712520000e-001,8.00933812054480750000e+000,
-6.49993870079576210000e+000,1.57738793711490090000e+001,5.02375291938038430000e+000,4.43699581506083970000e+000,
4.54574091554282500000e+000,-2.03471473639842460000e+000,-3.87043509559156630000e+000,-8.46720545155019040000e-001,
1.26756732392653920000e-002,1.12610485972773460000e+000,-2.81849219178339490000e-001,-4.46585719519653250000e+000,
4.78079201619137530000e+000,7.20241813527515440000e+000,-3.73539650354132750000e+000,-1.53253593031524620000e+000,
-4.83266386164988670000e+000,-2.69215304726191770000e+000,-5.94876012846670130000e+000,-3.26116785408658890000e+000,
-1.03878478056284110000e+001,2.20335032441057950000e-001,-1.19236656733315520000e+000,1.20543876550083220000e-002,
-5.49538456138832300000e+000,-6.77226851602575790000e-001,3.49531777060749890000e-001,3.38712764710896510000e+000,
-3.71605752757315950000e+000,1.08409514551553240000e+001,-2.68720444900039590000e+000,-3.05290853171085220000e+000,
5.39950916817901040000e+000,4.73137741535640900000e+000,-3.28456178643768880000e+000,-2.15107899564235930000e+000,
-9.33989025600738470000e+000,3.31662436821627440000e+000,4.95024974495926170000e+000,-2.20792212499064800000e+000,
-4.49705337613929770000e+000,8.65299280590953050000e+000,-3.26951707315755780000e+000,1.86358902474869480000e+000,
-2.01426084787260870000e+000,1.65016822294933550000e+001,-7.51547435319420360000e-001,2.01675038169942230000e+000,
-5.76538289032190130000e+000,-5.91442553914867110000e+000,7.98309992875715050000e+000,-2.43302340800836390000e+000,
-2.67114664397515390000e+000,9.10303624353602990000e+000,-1.46115147698390890000e+000,-5.26413427741111480000e+000,
8.11166525653914760000e+000,4.23998771800293990000e+000,-5.37378873456432780000e+000,1.95824154806895010000e+000,
-2.02102196375255170000e-001,2.04869779344129730000e+000,9.37067508011091380000e+000,-3.47025641687998300000e+000,
1.12350585841447700000e+000,1.02066013127507220000e+001,-6.85744943329824870000e+000,-3.36975325045859670000e+000,
-1.22661519270283680000e+000,1.75601053561052250000e+001,-5.37623299682625700000e+000,1.18450964443872840000e+000,
-1.04892755291358660000e+000,1.06211512121456760000e+000,5.63624612250501760000e+000,-3.43904235535493540000e+000,
-1.27872691170037990000e+000,1.61852708083980980000e+001,1.05418561368941090000e-001,-2.30010019670539290000e+000,
3.94689719401651340000e+000,2.36533125884515630000e+000,5.63734786514848360000e+000,2.46281582342676720000e+000,
-5.92228945119838370000e+000,-2.70009498741793230000e+000,-2.04548381986317460000e+000,4.71304256547800640000e+000,
8.94998852004973240000e+000,8.71106721821782060000e+000,2.71335197749214570000e+000,-3.97793392793331440000e+000,
5.50837217098035840000e+000,9.57779167146826320000e-001,4.92851977521583700000e+000,-6.26527150241898760000e+000,
5.35286235718743160000e+000,7.19619696074851880000e-001,1.87104306304927700000e+001,3.77853177159461500000e+000,
-6.81655938670052210000e+000,-2.07386404446184350000e+000,
  };

void TEM1_net_eval(const double keV, const double EPAdeg, const double Lm, double y[TEM1_Ny] ) {
  /* evaluate TEM1 logflux neural network for one energy (keV),
     equatorial pitch angle (EPAdeg, degrees), and Lm (McIlwain L) */

  /* this routine has been validated against matlab */

  double x[TEM1_Nx];
  x[0] = log(keV);
  x[1] = EPAdeg;
  x[2] = Lm;
  nnlib_eval(1,TEM1_Nx,x,TEM1_Nh,TEM1_theta, &TEM1_xbar, &TEM1_ybar, &TEM1_sx, &TEM1_sy,TEM1_Ny,y,0,NULL,NULL);
}

void TEM1_mu_sigma(const double keV, const double EPAdeg, const double Lm, double *mu, double *sigma) {
  /* get parameters of log-normal distribution for
     TEM1 (initialized with init_TEM1_net)
     energy in keV
     EPAdeg - equatorial pitch angle in degrees
     Lm - McIlwain L in OPQ
     returns:
     mu = mean of lognormal
     sigma = standard deviation of lognormal
   */

  /* this routine has been validated against matlab */

  double y[2]; /* stores nn output */
  TEM1_net_eval(keV,EPAdeg,Lm,y);
  *mu = y[0];
  *sigma = (y[1]-(*mu))/1.6448536270; /* norminv(0.95) = 1.64... */
  *sigma = GSL_MAX(log(2)/2,(*sigma)); /* force at least factor of 2 (95%) natural variability */
}

double w2u_TEM1_ell_combine(const gsl_vector *x, void *params_void, gsl_vector *grad, gsl_matrix *hess) {
  /* 
     compute TEM1 combined -log(p) penalty function
     at x = log(flux), with parameters params_void (actually a w2u_TEM1_paramsTy *),
     Note: The first x->size-1 elements of x are log fluxes, the last element is a LaGrange multiplier
     if grad is NULL grad and hess aren't calculated.
     if hess is NULL hess isn't calculated.
*/
  double fval=0;
  w2u_TEM1_paramsTy *params = (w2u_TEM1_paramsTy*)params_void;
  gsl_vector *flux, *dx, *siginv_dx, *Hf, *grad_logflux;
  gsl_vector_view hess_diag,grad_logflux_view;
  gsl_matrix_view hess_logflux_view;
  double lambda, dx_siginv_dx, lagrange, tmp, dloglamy;
  long int i;
  /* consts logflux_view and logflux will be declared in-line below */

  lagrange = gsl_vector_get(x,x->size-1); /* get lagrange multiplier */

  gsl_vector_const_view logflux_view = gsl_vector_const_subvector(x,0,x->size-1);
  const gsl_vector *logflux = &(logflux_view.vector);

  flux = gsl_vector_alloc(x->size-1);
  Hf = gsl_vector_alloc(flux->size);
  dx = gsl_vector_alloc(flux->size);
  siginv_dx = gsl_vector_alloc(flux->size);

  gsl_vector_memcpy(dx,logflux);
  gsl_vector_sub(dx,params->mu); /* dx = x-mu */
  /* undo log */
  vector_func(logflux,flux,&exp);
  gsl_vector_memcpy(Hf,params->H);
  gsl_vector_mul(Hf,flux); /* Hf(i) = H(i)*f(i) */
  gsl_blas_ddot(params->H,flux,&lambda); /* lambda = H'*flux */
  gsl_vector_set_zero(siginv_dx);
  gsl_blas_dgemv(CblasNoTrans,1.0,params->Sigma_inv,dx,1.0,siginv_dx); /* siginv_dx = Sigma_inv*dx */
  gsl_blas_ddot(dx,siginv_dx,&dx_siginv_dx); /* dx_siginv_dx = dx'*Signa_inv*dx */
  dloglamy = (log(lambda)-log(params->y));
  fval = lagrange*dloglamy + dx_siginv_dx/2.0;
  /*
  printf("lagrange=%lg, lam=%lg, y=%lg, lam-y=%lg, part1=%lg, part2=%lg\n",
	 lagrange, lambda,params->y,dloglamy,
	 lagrange*dloglamy,dx_siginv_dx/2.0);
  */

  if (grad) {
    /* grad = lagrange/lambda*Hf + siginv_dx */

    /* set up & set the logflux part, using views */
    grad_logflux_view = gsl_vector_subvector(grad,0,grad->size-1);
    grad_logflux = &(grad_logflux_view.vector);

    /* start with lagrange/lambda*Hf part */
    gsl_vector_memcpy(grad_logflux,Hf);
    gsl_vector_scale(grad_logflux,lagrange/lambda);

    /* add siginv_dx part */
    gsl_vector_add(grad_logflux,siginv_dx);

    /* d/dlagrange part */
    gsl_vector_set(grad,grad->size-1,dloglamy);

    if (hess) {
      gsl_matrix_set_zero(hess); /* clear hess */

      /* logflux part of hess = -lagrange*(Hf)(Hf)'/lambda^2 + lagrange*diag(Hf)/lambda + Sigma_inv */
      /* prepare to operate only on logflux part */
      hess_logflux_view = gsl_matrix_submatrix(hess,0,0,grad->size-1,grad->size-1);

      /* first handle diag part due to Hf */
      hess_diag = gsl_matrix_diagonal(&(hess_logflux_view.matrix)); /* access diagonal of hess as vector view */
      gsl_vector_memcpy(&(hess_diag.vector),Hf); /* set Hf along diagonal */
      gsl_vector_scale(&(hess_diag.vector),lagrange/lambda); /* mult diagonal by lagrange/lambda */

      /* now add -lagrange*(Hf)(Hf)'/lambda^2 part */
      gsl_blas_dger(-lagrange/gsl_pow_2(lambda),Hf,Hf,&(hess_logflux_view.matrix));

      /* next, add part due to Sigma_inv*/
      gsl_matrix_add(&(hess_logflux_view.matrix),params->Sigma_inv); /* Sigma_inv part */

      /* now set outside rows/columns d^ell/dlagrange/dx = Hf/lambda*/
      for (i=0; i < Hf->size; i++) {
	tmp = gsl_vector_get(Hf,i)/lambda;
	gsl_matrix_set(hess,hess->size1-1,i,tmp);
	gsl_matrix_set(hess,i,hess->size2-1,tmp);
      }

    }
  }

  /* free memory */
  gsl_vector_free(flux);
  gsl_vector_free(Hf);
  gsl_vector_free(dx);
  gsl_vector_free(siginv_dx);

  return(fval); /* return neglogp */
}


int wide2uni_TEM1(const double *wideflux, const double *dlogwideflux,
             const double *PAgrid, const double *H, 
             const long int *int_params, 
             const double *real_params,
             char *outFile,
             const long int *ialpha0,
             double *uniflux, double *dloguniflux) {
  /*
    Uses TEM-1 method to compute wede2uni.
    Arguments are same as for wide2uni.
   */
  long int NA, minimizer_flag, MaxIter, verbose;
  long int ifirst=0; /* first PA grid point out of northern loss cone */
  long int ilast=0; /* first PA grid point out of southern loss cone */
  long int NA2; /* number of PA points outside loss cone (reduced vector) */
  long int i2alpha0; /* index of ialpha0 in reduced vector */
  double Lm,BB0,keV,EPAcutoff, PAcutoff;
  double mui,si,sj,tmp,EPAdegi,EPAdegj,rho;
  gsl_vector *mu=0, *sigma=0, *EPAdeg=0,*flux0=0;
  gsl_matrix *Sigma=0, *Sigma_inv=0;
  w2u_TEM1_paramsTy w2u_TEM1_params;
  long int i,j;
  gsl_vector *x=0, *grad=0, *logflux=0;
  gsl_vector_view logflux_view;
  gsl_matrix *covx=0, *hess=0;
  optfunTy optfun;
  /* Hgsl will be declared below */

  int clean_return(int return_value) {
    if (w2u_TEM1_params.closeOutFile) {
      fclose(w2u_TEM1_params.outFilePtr);
    }
    return(return_value);
  }

  w2u_TEM1_params.closeOutFile = 0; /* default: no need to close outFile */

  /* inputs used by all methods already checked: */
  /* inputs already checked for NULLs */
  /* NA, w2u_method, and verbose are already checked for validity */
  NA = int_params[0];
  /* w2u_method = int_params[1]; */
  verbose = int_params[2];
  /* maybe not all methods use these, so check them here: */
  minimizer_flag = int_params[3];
  MaxIter = int_params[4];


  /* check for valid minimizer flag */
  if ((minimizer_flag < 0) || (minimizer_flag > OPTIM_MIN_MAX)) {
    return(clean_return(INVLIB_ERR_INVALIDMIN));
  }

  /* check for valid maximum interations */
  if (MaxIter<=0) {
    return(clean_return(INVLIB_ERR_INVALIDITER));
  }

  keV = real_params[0]; /* Electron energy, keV */
  BB0 = real_params[1]; /* B/B0, OPQ */
  Lm = real_params[2]; /* McIlwain L, OPQ */

  if ((!gsl_finite(Lm)) || (Lm < 1) || (!gsl_finite(BB0)) || (BB0 < 1) || (!gsl_finite(keV)) || (keV < 0)) {
    return(clean_return(INVLIB_ERR_DATANAN));
  }

  EPAcutoff = TEM1_AtmoCutoff(Lm); /* acute equatorial pitch angle of loss cone */
  PAcutoff = asin(sin(EPAcutoff*M_PI/180.0)*sqrt(BB0))*180.0/M_PI; /* acute local pitch angle */
  if (PAgrid[(*ialpha0)] < PAcutoff) {
    return(clean_return(INVLIB_ERR_IALPHA0));
  }
  ifirst=0;
  while ((ifirst < NA) && (PAgrid[ifirst] < PAcutoff)) {
    ifirst++;
  }
  if (ifirst == NA) {
    return(clean_return(INVLIB_ERR_IALPHA0)); 
  }

  ilast=NA-1;
  while ((ilast >= 0) && (PAgrid[ilast] > 180-PAcutoff)) {
    ilast--;
  }
  if (ilast < 0) {
    return(clean_return(INVLIB_ERR_IALPHA0)); 
  }

  NA2 = ilast-ifirst+1;
  /* any pitch angle grid points left? */
  if (NA2 <=0) {
    return(clean_return(INVLIB_ERR_DATANAN));
  }

  /* calculate ialpha0 in new reduced vector */
  i2alpha0 = *ialpha0-ifirst;

  /* build vectors */
  EPAdeg = gsl_vector_alloc(NA2);
  mu = gsl_vector_alloc(NA2);
  sigma = gsl_vector_alloc(NA2);
  for (i=0; i < NA2; i++) {
    EPAdegi = asin(sin(PAgrid[ifirst+i]*M_PI/180.0)/sqrt(BB0))*180.0/M_PI;
    gsl_vector_set(EPAdeg,i,EPAdegi);
    TEM1_mu_sigma(keV,EPAdegi,Lm,&mui,&si);
    gsl_vector_set(mu,i,mui);
    gsl_vector_set(sigma,i,si);
  }
  /* declare this here because it's a const view */
  gsl_vector_const_view Hgsl = gsl_vector_const_view_array(H+ifirst,NA2); /* subset of H starting at ifirst */

  /* build Sigma */
  Sigma = gsl_matrix_alloc(NA2,NA2);
  for (i=0; i < NA2; i++) {
    EPAdegi = gsl_vector_get(EPAdeg,i);
    si = gsl_vector_get(sigma,i);
    gsl_matrix_set(Sigma,i,i,si*si);
    for (j=(i+1); j < NA2; j++ ){
      EPAdegj = gsl_vector_get(EPAdeg,j);
      sj = gsl_vector_get(sigma,j);
      rho = TEM1_alpha_rho(Lm,keV,EPAdegi,EPAdegj);
      tmp = si*sj*rho;
      gsl_matrix_set(Sigma,i,j,tmp);
      gsl_matrix_set(Sigma,j,i,tmp);
    }
  }
  /* won't need these any more, free them */
  gsl_vector_free(EPAdeg);
  gsl_vector_free(sigma);

  /* build sigma inv */
  Sigma_inv = gsl_matrix_alloc(NA2,NA2);
  inv_matrix_once(Sigma,Sigma_inv,1);
  gsl_matrix_free(Sigma); /* done with this */


  /* prepare to call optimizer */
  x = gsl_vector_alloc(NA2+1); /* logflux and 1 lagrange multiplier */
  gsl_vector_set(x,NA2,1); /* initialize the lagrange multiplier to 1 */
  /* create view of just logflux part */
  logflux_view = gsl_vector_subvector(x,0,NA2); 
  logflux = &(logflux_view.vector);

  w2u_TEM1_params.y = (*wideflux);
  w2u_TEM1_params.H = &(Hgsl.vector); /* will not preserve constness */
  w2u_TEM1_params.mu = mu;
  w2u_TEM1_params.Sigma_inv = Sigma_inv;

  /* construct good initial guess: wide/(H'*exp(mu)) */
  flux0 = vector_func(mu,NULL,&exp); /* exp(mu) */
  gsl_blas_ddot(&(Hgsl.vector),flux0,&tmp);
  gsl_vector_scale(flux0,(*wideflux)/tmp);
  vector_func(flux0,logflux,&log); /* undo log */
  gsl_vector_free(flux0);


  switch(verbose) { /* this switch statement must be consistent with specinv.h */
  case 0:
    w2u_TEM1_params.outFilePtr = NULL;
    break;
  case 1:
    w2u_TEM1_params.outFilePtr = stdout;
    break;
  case 2:
    w2u_TEM1_params.outFilePtr = stderr;
    break;
  case 3:
    if (outFile) {
      w2u_TEM1_params.outFilePtr = (FILE *)outFile;
    } else {
      return(clean_return(INVLIB_ERR_VERBOSE)); /* can't output to NULL */
    }
    break;
  case 4:
    if (outFile) {
      if (! (w2u_TEM1_params.outFilePtr = fopen(outFile,"w"))) {
	return(clean_return(INVLIB_ERR_OUTFILE)); /* Couldn't open outFile */
	  }
      w2u_TEM1_params.closeOutFile = 1;
    } else {
      return(clean_return(INVLIB_ERR_VERBOSE)); /* can't output to NULL */
    }
    break;
  case 5:
    if (outFile) {
      if (! (w2u_TEM1_params.outFilePtr = fopen(outFile,"a"))) {
	return(clean_return(INVLIB_ERR_OUTFILE)); /* Couldn't open outFile */
	  }
      w2u_TEM1_params.closeOutFile = 1;
    } else {
      return(clean_return(INVLIB_ERR_VERBOSE)); /* can't output to NULL */
    }
    break;
  default:
    return(clean_return(INVLIB_ERR_VERBOSE)); /* invalid verbose value */
  }

  optfun.func = &w2u_TEM1_ell_combine;
  optfun.params = (void *)(&w2u_TEM1_params);
  fzero(x,&optfun, minimizer_flag, MaxIter,w2u_TEM1_params.outFilePtr);

  /* get uniflux */
  *uniflux = exp(gsl_vector_get(logflux,i2alpha0));

  /* now get error bar */
  grad = gsl_vector_alloc(NA2+1);
  hess = gsl_matrix_alloc(NA2+1,NA2+1);
  optfun.func(x,optfun.params,grad,hess); /* need hess */

  /* free memory */
  gsl_vector_free(x);
  gsl_vector_free(mu);
  gsl_matrix_free(Sigma_inv);
  gsl_vector_free(grad);

  covx = gsl_matrix_alloc(NA2+1,NA2+1);
  inv_matrix_once(hess,covx,1); /* covx = inv(hess) */
  gsl_matrix_free(hess);
  *dloguniflux = sqrt(gsl_pow_2(*dlogwideflux) + gsl_matrix_get(covx,i2alpha0,i2alpha0));
  gsl_matrix_free(covx);
  
  /* done, success! */
  return(clean_return(INVLIB_SUCCESS));
}
