//
//  TS07.cpp
//  UBJDevelopment
//
//  Created by KYUNGGUK MIN on 7/8/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#include "TS07.h"
#include "Geopack.h"
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace UBK {
    using namespace std;

    /*! @abstract TS07 context definition.
     @discussion TS07 model requires model coefficient, PARAM_A, that are specific time of interest. In addition, it requires lengthy shielding parameters. Since the model is still in development and these coefficients can be possibly changed, Context abstract class is provided, instead of hard coding the coefficients. Although it's an abstract class, it's struct definition with default constructor and destructor, and public member variables for easy used and for easy adoptation to C. User can subclass to implement his/her own coefficient loading routines, or just instanciate Context variable and assign the coefficients.
     @note All zeroth indexes of coefficient C array are dummy. Writing should start from index 1 at every rank.
     */
    struct TS07Context {
        // MODEL PARAMETER FOR A SPECIFIC TIME MOMENT
        //double const* PARAM_A;// PARAM_A[102]; // Zeroth index is dummy. Coefficients should be filled starting index 1.
        // SHIELDING FIELD PARAMETERS
        double TSS_TSS[81][6];// TSS_TSS[81][6]; Zeroth indexes are dummy, i.e. TSS_TSS[0][*] and TSS_TSS[*][0]. Coefficients should be filled starting index 1.
        double TSO_TSO[81][6][5];// TSO_TSO[81][6][5]; Zeroth indexes are dummy, i.e. TSO_TSO[*][0][0], TSO_TSO[0][*][0] and TSO_TSO[0][0][*]. Coefficients should be filled starting index 1.
        double TSE_TSE[81][6][5];// TSE_TSE[81][6][5]; Zeroth indexes are dummy, i.e. TSE_TSE[*][0][0], TSE_TSE[0][*][0] and TSE_TSE[0][0][*]. Coefficients should be filled starting index 1.

        TS07Context();
    };
    static TS07Context const ts07shielding;

    ////////////////////////////////////////////////////////////////////////////////
    // Implementation
    ////////////////////////////////////////////////////////////////////////////////
    TS07::TS07 (Geopack const* geopack, double const parmod[]) : TSExternalField(geopack, 1, parmod)
    {
        if (NULL == parmod) {
            throw invalid_argument("PARMOD must be 102-element array (PARMOD[0] is Pdyn and and PARMOD[1:101] are the fitting coefficients).");
        }
        this->setParmod(parmod);
    }

    /*! Entry function to TS07 external model. The required arguments are exactly same as the previous models except Context variable which was unavoidable. This is the only visible function of the model.
     In the original fortran code, PARMOD is dummy argument and PDYN, only required parameter from PARMOD, is passed via COMMON BLOCK variable. In this version, however, the first element of PARMOD array is used to pass PDYN. Since the first element of PARMOD array in the previous models are PDYN, this way makes much more sense and causes little confusion.
     @note Unlike the original fortran code, GEOPACK_08 is not part of the code. User should calculate dipole tilt angle as well as main field somewhere else in his/her own code.
     */
    static void TS07D (int const IOPT,double const PARMOD[]/*zero index based 10 element C array*/,double const PS,
                       double const X,double const Y,double const Z,
                       double *BX,double *BY,double *BZ,
                       TS07Context const *ctx);
    void TS07::getFieldInGSW_atPoint(Point *bOut, const Point ptgsw) const
    {
        TS07D(this->iopt(), this->parmod(), this->geopack()->psi(), ptgsw.x, ptgsw.y, ptgsw.z, &bOut->x, &bOut->y, &bOut->z, &ts07shielding);
    }


    ////////////////////////////////////////////////////////////////
    // TS07 Model
    ////////////////////////////////////////////////////////////////
#pragma mark - TS07 Model
    /*-----------------------------Model---------------------*/
    /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     P r o g r a m   TS07D

     c----------------------------------------------------------------------
     c
     REAL*8 A,PDYN,TSS,TSO,TSE

     CHARACTER*80 NAMETSS(5),NAMETSO(5,4),NAMETSE(5,4)

     PARAMETER (NTOT=101)

     COMMON /GEOPACK1/ AAA(10),SPS,CPS,BBB(3),PSI,CCC(18)
     COMMON /PARAM/ A(NTOT)
     COMMON /INPUT/  PDYN

     COMMON /TSS/ TSS(80,5)
     COMMON /TSO/ TSO(80,5,4)
     COMMON /TSE/ TSE(80,5,4)

     DIMENSION PARMOD(10)
     DIMENSION CXY(101,101),CJY(101,101),XX(101),ZZ(101)

     DATA NAMETSS/                                  !  SHIELDING FIELD PARAMETERS:
     *'../TSG_DYN_PAR/tailamebhr1.par',         !   MAKE SURE TO MODIFY THE PATH TO THE ACTUAL
     *'../TSG_DYN_PAR/tailamebhr2.par',         !   STORAGE FOLDER, IF NECESSARY
     *'../TSG_DYN_PAR/tailamebhr3.par',
     *'../TSG_DYN_PAR/tailamebhr4.par',
     *'../TSG_DYN_PAR/tailamebhr5.par'/
     DATA NAMETSO/
     *'../TSG_DYN_PAR/tailamhr_o_11.par',
     *'../TSG_DYN_PAR/tailamhr_o_21.par',
     *'../TSG_DYN_PAR/tailamhr_o_31.par',
     *'../TSG_DYN_PAR/tailamhr_o_41.par',
     *'../TSG_DYN_PAR/tailamhr_o_51.par',
     *'../TSG_DYN_PAR/tailamhr_o_12.par',
     *'../TSG_DYN_PAR/tailamhr_o_22.par',
     *'../TSG_DYN_PAR/tailamhr_o_32.par',
     *'../TSG_DYN_PAR/tailamhr_o_42.par',
     *'../TSG_DYN_PAR/tailamhr_o_52.par',
     *'../TSG_DYN_PAR/tailamhr_o_13.par',
     *'../TSG_DYN_PAR/tailamhr_o_23.par',
     *'../TSG_DYN_PAR/tailamhr_o_33.par',
     *'../TSG_DYN_PAR/tailamhr_o_43.par',
     *'../TSG_DYN_PAR/tailamhr_o_53.par',
     *'../TSG_DYN_PAR/tailamhr_o_14.par',
     *'../TSG_DYN_PAR/tailamhr_o_24.par',
     *'../TSG_DYN_PAR/tailamhr_o_34.par',
     *'../TSG_DYN_PAR/tailamhr_o_44.par',
     *'../TSG_DYN_PAR/tailamhr_o_54.par'/
     DATA NAMETSE/
     *'../TSG_DYN_PAR/tailamhr_e_11.par',
     *'../TSG_DYN_PAR/tailamhr_e_21.par',
     *'../TSG_DYN_PAR/tailamhr_e_31.par',
     *'../TSG_DYN_PAR/tailamhr_e_41.par',
     *'../TSG_DYN_PAR/tailamhr_e_51.par',
     *'../TSG_DYN_PAR/tailamhr_e_12.par',
     *'../TSG_DYN_PAR/tailamhr_e_22.par',
     *'../TSG_DYN_PAR/tailamhr_e_32.par',
     *'../TSG_DYN_PAR/tailamhr_e_42.par',
     *'../TSG_DYN_PAR/tailamhr_e_52.par',
     *'../TSG_DYN_PAR/tailamhr_e_13.par',
     *'../TSG_DYN_PAR/tailamhr_e_23.par',
     *'../TSG_DYN_PAR/tailamhr_e_33.par',
     *'../TSG_DYN_PAR/tailamhr_e_43.par',
     *'../TSG_DYN_PAR/tailamhr_e_53.par',
     *'../TSG_DYN_PAR/tailamhr_e_14.par',
     *'../TSG_DYN_PAR/tailamhr_e_24.par',
     *'../TSG_DYN_PAR/tailamhr_e_34.par',
     *'../TSG_DYN_PAR/tailamhr_e_44.par',
     *'../TSG_DYN_PAR/tailamhr_e_54.par'/
     C
     DO 1001 IREAD=1,5
     OPEN (UNIT=1,FILE=NAMETSS(IREAD))
     READ (1,200) (TSS(KK,IREAD),KK=1,80)
     200  FORMAT(G17.10)
     1001 CLOSE(1)

     DO 1002 IREAD=1,5
     DO 1003 KREAD=1,4
     OPEN (UNIT=1,FILE=NAMETSO(IREAD,KREAD))
     READ (1,200) (TSO(KK,IREAD,KREAD),KK=1,80)
     1003   CONTINUE
     1002 CLOSE(1)

     DO 1004 IREAD=1,5
     DO 1005 KREAD=1,4
     OPEN (UNIT=1,FILE=NAMETSE(IREAD,KREAD))
     READ (1,200) (TSE(KK,IREAD,KREAD),KK=1,80)
     1005   CONTINUE
     1004 CLOSE(1)
     C
     PRINT *, '   SHIELDING COEFFICIENTS HAS BEEN READ INTO RAM'
     c
     OPEN (UNIT=1,FILE='../MR98_2/Best_par_so_far.par')  !  MODEL PARAMETER FILE FOR
     READ (1,100) (A(I),I=1,NTOT)                             !  A SPECIFIC TIME MOMENT
     100  FORMAT(G15.6)                                            !  MAKE SURE TO MODIFY THE PATH
     CLOSE(1)                                                 !  IF NECESSARY
     C
     C
     PRINT *, '   ENTER PDYN (IN NANOPASCALS)'
     READ *, PDYN

     c -------------------------------------------

     PRINT *,'  ENTER DATE & UT:  YEAR, DAYOFYEAR, HOUR, MIN, SEC'
     READ *, IYEAR,IDAY,IHOUR,MIN,ISEC

     c ------------------------------------------------
     C
     VXGSE=-400.  !  GSE COMPONENTS OF SOLAR WIND VELOCITY VECTOR; THIS PARTICULAR CHOICE
     VYGSE=   0.  !   IS MADE IN ORDER TO MAKE THE GSW SYSTEM IDENTICAL TO THE STANDARD GSM
     VZGSE=   0.

     call RECALC_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,VXGSE,VYGSE,VZGSE) ! CALCULATES TILT ANGLE AND
     C                                                                  UPDATES MAIN FIELD COEFFICIENTS
     C
     PRINT *, '  Enter S/C position: XGSM,YGSM,ZGSM '
     READ *, XGSM,YGSM,ZGSM

     CALL IGRF_GSW_08 (XGSM,YGSM,ZGSM,HXGSM,HYGSM,HZGSM)  !  HERE HX,HY,HZ ARE IN THE STANDARD GSM
     C                                                                  DUE TO THE ABOVE ASSUMED VY=VZ=0.
     CALL EXTMODEL (1,PARMOD,PSI,XGSM,YGSM,ZGSM,BXGSM,BYGSM,BZGSM)

     PRINT *,'     Main field:'
     PRINT *,' HXGSM,HYGSM,HZGSM=',HXGSM,HYGSM,HZGSM
     PRINT *,' '
     PRINT *,'   External field:'
     PRINT *,' BXGSM,BYGSM,BZGSM=',BXGSM,BYGSM,BZGSM
     PRINT *,' '
     PRINT *,'  Geodipole tilt (degrees):', PSI*57.29578

     PAUSE
     c
     END
     */

    /*! Internal working space. Do not touch.
     */
    struct TS_EXTERNAL_COMMON_ {
        // COMMON /TAIL/
        double TAIL_D;
        // COMMON /BIRKPAR/
        double BIRKPAR_XKAPPA1,BIRKPAR_XKAPPA2;
        // COMMON /G/
        double G_G,G_TW;
        // COMMON /RH0/
        double RH0_RH0;
        // COMMON /DPHI_B_RHO0/
        double DPHI_B_RHO0_DPHI,DPHI_B_RHO0_B,DPHI_B_RHO0_RHO_0,DPHI_B_RHO0_XKAPPA;
        // COMMON /MODENUM/
        int MODENUM_M;
        // COMMON /DTHETA/
        double DTHETA_DTHETA;
    };

    // ===================================================================
    //     DOUBLE PRECISION FUNCTION bessj0(x)
    static double bessj0(double const x) { return j0(x);
        //IMPLICIT REAL*8 (A-H,O-Z)
        //SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
        //*s5,s6
        static double const p1=1.e0,p2=-.1098628627e-2,p3=.2734510407e-4,p4=-.2073370639e-5,p5=.2093887211e-6;
        static double const q1=-.1562499995e-1,q2=.1430488765e-3,q3=-.6911147651e-5,q4=.7621095161e-6,q5=-.934945152e-7;

        static double const r1=57568490574.e0,r2=-13362590354.e0,r3=651619640.7e0,r4=-11214424.18e0,r5=77392.33017e0,r6=-184.9052456e0;
        static double const s1=57568490411.e0,s2=1029532985.e0,s3=9494680.718e0,s4=59272.64853e0,s5=267.8532712e0,s6=1.e0;

        if(fabs(x) <  8.e0) {//if(Dabs(x).lt.8.D0)then
            double y=x*x;
            return (r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) / (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))));
        } else {//else
            double ax=fabs(x);
            double z=8.e0/ax;
            double y=z*z;
            double xx=ax-.785398164e0;
            return sqrt(.636619772/ax) * (cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))));
        }//endif
        //return
    }//END
    // (C) Copr. 1986-92 Numerical Recipes Software .)+1YX39'=K+1.

    //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
    //    DOUBLE PRECISION FUNCTION bessj1(x)
    static double bessj1(double const x) { return j1(x);
        //IMPLICIT REAL*8 (A-H,O-Z)
        //SAVE p1,p2,p3,p4,p5,q1,q2,q3,q4,q5,r1,r2,r3,r4,r5,r6,s1,s2,s3,s4,
        //*s5,s6
        static double const r1=72362614232.e0,r2=-7895059235.e0,r3=242396853.1e0,r4=-2972611.439e0,r5=15704.48260e0,r6=-30.16036606e0;
        static double const s1=144725228442.e0,s2=2300535178.e0,s3=18583304.74e0,s4=99447.43394e0,s5=376.9991397e0,s6=1.e0;

        static double const p1=1.e0,p2=.183105e-2,p3=-.3516396496e-4,p4=.2457520174e-5,p5=-.240337019e-6;
        static double const q1=.04687499995e0,q2=-.2002690873e-3,q3=.8449199096e-5,q4=-.88228987e-6,q5=.105787412e-6;

        if(fabs(x) < 8.e0) {//if(Dabs(x).lt.8.D0)then
            double y=x*x;
            return x*(r1+y*(r2+y*(r3+y*(r4+y*(r5+y*r6))))) / (s1+y*(s2+y*(s3+y*(s4+y*(s5+y*s6)))));
        } else {//else
            double ax=fabs(x);
            double z=8.e0/ax;
            double y=z*z;
            double xx=ax-2.356194491e0;
            return sqrt(.636619772/ax) * (cos(xx)*(p1+y*(p2+y*(p3+y*(p4+y*p5))))-z*sin(xx)*(q1+y*(q2+y*(q3+y*(q4+y*q5)))))*(x>=0. ? 1.:-1.);//Dsign(1.D0,x)
        }//endif
        //return
    }//END
    //  (C) Copr. 1986-92 Numerical Recipes Software .)+1YX39'=K+1.

    // ===================================================================
    static double bessj(int const n,double const x) { return jn(n, x);
        //IMPLICIT REAL*8 (A-H,O-Z)
        static int const IACC=40;
        static double const BIGNO=1.e10,BIGNI=1.e-10;//PARAMETER (IACC=40,BIGNO=1.D10,BIGNI=1.D-10)
        //U    USES bessj0,bessj1
        assert(n>=2);//if(n.lt.2)pause 'bad argument n in bessj'

        double bessjRet;
        double ax=fabs(x);
        if(ax == 0.e0) {//if(ax.eq.0.D0)then
            bessjRet=0.;
        } else if(ax > (double)(n)) {//else if(ax.gt.Dfloat(n))then
            double tox=2.e0/ax;
            double bjm=bessj0(ax);
            double bj=bessj1(ax);
            for (int j=1; j<=n-1; j++) {//do 11 j=1,n-1
                double bjp=j*tox*bj-bjm;
                bjm=bj;
                bj=bjp;
            }//11      continue
            bessjRet=bj;
        } else {//else
            double tox=2.e0/ax;
            int m=2*((n+(int)(sqrt((double)(IACC*n))))/2);
            bessjRet=0.e0;
            int jsum=0;
            double sum=0.e0;
            double bjp=0.e0;
            double bj=1.e0;
            for (int j=m; j>=1; j--) {//do 12 j=m,1,-1
                double bjm=j*tox*bj-bjp;
                bjp=bj;
                bj=bjm;
                if(fabs(bj) > BIGNO) {//if(Dabs(bj).gt.BIGNO)then
                    bj=bj*BIGNI;
                    bjp=bjp*BIGNI;
                    bessjRet=bessjRet*BIGNI;
                    sum=sum*BIGNI;
                }//endif
                if(jsum != 0)sum=sum+bj;
                jsum=1-jsum;
                if(j == n)bessjRet=bjp;
            }//12      continue
            sum=2.e0*sum-bj;
            bessjRet=bessjRet/sum;
        }//endif
        if(x < 0.e0 && n%2 == 1)bessjRet=-bessjRet;
        return bessjRet;
    }//END
    //  (C) Copr. 1986-92 Numerical Recipes Software .)+1YX39'=K+1.

    //=====================================================================================
    static double R_S(double const A[],double const R,double const THETA) {
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION A(31)

        return R+A[2]/R+A[3]*R/sqrt(R*R+A[11]*A[11])+A[4]*R/(R*R+A[12]*A[12])
        +(A[5]+A[6]/R+A[7]*R/sqrt(R*R+A[13]*A[13])+A[8]*R/(R*R+A[14]*A[14]))*
        cos(THETA)
        +(A[9]*R/sqrt(R*R+A[15]*A[15])+A[10]*R/pow((R*R+A[16]*A[16]), 2))
        *cos(2.e0*THETA);

    }//END

    //-----------------------------------------------------------------------------
    static double THETA_S(double const A[],double const R,double const THETA) {
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION A(31)

        return THETA+(A[17]+A[18]/R+A[19]/(R*R)
                      +A[20]*R/sqrt(R*R+A[27]*A[27]))*sin(THETA)
        +(A[21]+A[22]*R/sqrt(R*R+A[28]*A[28])
          +A[23]*R/(R*R+A[29]*A[29]))*sin(2.e0*THETA)
        +(A[24]+A[25]/R+A[26]*R/(R*R+A[30]*A[30]))*sin(3.e0*THETA);

        //RETURN
    }//END

    /*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     c
     SUBROUTINE FIALCOS(R,THETA,PHI,BTHETA,BPHI,N,THETA0,DT)
     C
     C  CONICAL MODEL OF BIRKELAND CURRENT FIELD; BASED ON THE OLD S/R FIALCO (OF 1990-91)
     C  SEE THE OLD NOTEBOOK 1985-86-88, NOTE OF MARCH 5, BUT HERE BOTH INPUT AND OUTPUT ARE IN SPHERICAL CDS.

     C  BTN, AND BPN ARE THE ARRAYS OF BTHETA AND BPHI (BTN(i), BPN(i) CORRESPOND TO i-th MODE).
     C   ONLY FIRST  N  MODE AMPLITUDES ARE COMPUTED (N<=10).
     C    THETA0 IS THE ANGULAR HALF-WIDTH OF THE CONE, DT IS THE ANGULAR H.-W. OF THE CURRENT LAYER

     C   NOTE:  BR=0  (BECAUSE ONLY RADIAL CURRENTS ARE PRESENT IN THIS MODEL)
     */
    static void FIALCOS(double const R,double const THETA,double const PHI,
                        double *BTHETA,double *BPHI,
                        int const N,double const THETA0,double const DT) {
        //IMPLICIT REAL*8 (A-H,O-Z)
        double BTN[11],BPN[11],CCOS[11],SSIN[11];//DIMENSION  BTN(10),BPN(10),CCOS(10),SSIN(10)

        double SINTE=sin(THETA);
        double RO=R*SINTE;
        double COSTE=cos(THETA);
        double SINFI=sin(PHI);
        double COSFI=cos(PHI);
        double TG=SINTE/(1.e0+COSTE);//   !        TAN(THETA/2)
        double CTG=SINTE/(1.e0-COSTE);//  !        CTG(THETA/2)


        double TETANP=THETA0+DT;
        double TETANM=THETA0-DT;
        double TGP,TGM,TGM2,TGP2;
        if(THETA >= TETANM) {//IF(THETA.LT.TETANM) GOTO 1
            TGP=tan(TETANP*0.5e0);
            TGM=tan(TETANM*0.5e0);
            TGM2=TGM*TGM;
            TGP2=TGP*TGP;
        }//1   CONTINUE

        double COSM1=1.e0;
        double SINM1=0.e0;
        double TM=1.e0;
        double TGM2M=1.e0;
        double TGP2M=1.e0;

        double T,DTT,DTT0,FC,FC1,TGM2M1,TG21;
        for (int M=1; M<=N; M++) {//DO 2 M=1,N
            TM=TM*TG;
            CCOS[M]=COSM1*COSFI-SINM1*SINFI;
            SSIN[M]=SINM1*COSFI+COSM1*SINFI;
            COSM1=CCOS[M];
            SINM1=SSIN[M];
            if(THETA < TETANM) {//IF(THETA.LT.TETANM) THEN
                T=TM;
                DTT=0.5e0*M*TM*(TG+CTG);
                DTT0=0.e0;
            } else if (THETA < TETANP) {//ELSE IF(THETA.LT.TETANP) THEN
                TGM2M=TGM2M*TGM2;
                FC=1.e0/(TGP-TGM);
                FC1=1.e0/(2*M+1);
                TGM2M1=TGM2M*TGM;
                TG21=1.e0+TG*TG;
                T=FC*(TM*(TGP-TG)+FC1*(TM*TG-TGM2M1/TM));
                DTT=0.5e0*M*FC*TG21*(TM/TG*(TGP-TG)-FC1*(TM-TGM2M1/(TM*TG)));
                DTT0=0.5e0*FC*((TGP+TGM)*(TM*TG-FC1*(TM*TG-TGM2M1/TM))+
                               TM*(1.e0-TGP*TGM)-(1.e0+TGM2)*TGM2M/TM);
            } else {//ELSE
                TGP2M=TGP2M*TGP2;
                TGM2M=TGM2M*TGM2;
                FC=1.e0/(TGP-TGM);
                FC1=1.e0/(2*M+1);
                T=FC*FC1*(TGP2M*TGP-TGM2M*TGM)/TM;
                DTT=-T*M*0.5e0*(TG+CTG);
            }//ENDIF

            BTN[M]=M*T*CCOS[M]/RO;
            BPN[M]=-DTT*SSIN[M]/R;   // 2
        }

        *BTHETA=BTN[N] *800.;
        *BPHI  =BPN[N] *800.;

        return;
    }//END

    /*
     C-------------------------------------------------------------------------
     C
     SUBROUTINE ONE_CONE(A,X,Y,Z,BX,BY,BZ)
     c
     c  RETURNS FIELD COMPONENTS FOR A DEFORMED CONICAL CURRENT SYSTEM, FITTED TO A BIOSAVART FIELD
     c    BY SIM_14.FOR.  HERE ONLY THE NORTHERN CONE IS TAKEN INTO ACCOUNT.
     */
    static void ONE_CONE(double const A[],double const X,double const Y,double const Z,
                         double *BX,double *BY,double *BZ,
                         struct TS_EXTERNAL_COMMON_ const* TS_EXTERNAL_COMMON) {
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION A(31)

        double const* DTHETA=&TS_EXTERNAL_COMMON->DTHETA_DTHETA;//COMMON /DTHETA/ DTHETA
        int const* M=&TS_EXTERNAL_COMMON->MODENUM_M;//COMMON /MODENUM/ M

        static double const DR=1.e-6,DT=1.e-6;//  !   JUST FOR NUMERICAL DIFFERENTIATION

        double THETA0=A[31];

        double RHO2=X*X+Y*Y;
        double RHO=sqrt(RHO2);
        double R=sqrt(RHO2+Z*Z);
        double THETA=atan2(RHO,Z);
        double PHI=atan2(Y,X);
        /*
         C   MAKE THE DEFORMATION OF COORDINATES:
         */
        double RS=R_S(A,R,THETA);
        double THETAS=THETA_S(A,R,THETA);
        double PHIS=PHI;
        /*
         C   CALCULATE FIELD COMPONENTS AT THE NEW POSITION (ASTERISKED):
         */
        double BTAST,BFAST;
        FIALCOS (RS,THETAS,PHIS,&BTAST,&BFAST,*M,THETA0,*DTHETA);//    !   MODE #M
        /*
         C   NOW TRANSFORM B{R,T,F}_AST BY THE DEFORMATION TENSOR:
         C
         C      FIRST OF ALL, FIND THE DERIVATIVES:
         */
        double DRSDR=(R_S(A,R+DR,THETA)-R_S(A,R-DR,THETA))/(2.e0*DR);
        double DRSDT=(R_S(A,R,THETA+DT)-R_S(A,R,THETA-DT))/(2.e0*DT);
        double DTSDR=(THETA_S(A,R+DR,THETA)-THETA_S(A,R-DR,THETA))/(2.e0*DR);
        double DTSDT=(THETA_S(A,R,THETA+DT)-THETA_S(A,R,THETA-DT))/(2.e0*DT);
        double STSST=sin(THETAS)/sin(THETA);
        double RSR=RS/R;

        double BR     =-RSR/R*STSST*BTAST*DRSDT;//                 !   NB#6, P.43    BRAST DOES NOT ENTER HERE
        double BTHETA = RSR*STSST*BTAST*DRSDR;//                  !               (SINCE IT IS ZERO IN OUR CASE)
        double BPHI   = RSR*BFAST*(DRSDR*DTSDT-DRSDT*DTSDR);//

        double S=RHO/R;
        double C=Z/R;
        double SF=Y/RHO;
        double CF=X/RHO;

        double BE=BR*S+BTHETA*C;

        *BX=A[1]*(BE*CF-BPHI*SF);
        *BY=A[1]*(BE*SF+BPHI*CF);
        *BZ=A[1]*(BR*C-BTHETA*S);

        return;
    }//END

    /*
     C=========================================================================
     c
     SUBROUTINE TWOCONSS (A,X,Y,Z,BX,BY,BZ)
     C
     C   DIFFERS FROM TWOCONES:  THIS S/R IS FOR THE "SYMMETRIC" MODE OF BIRKELAND CURRENTS IN THAT
     C                           HERE THE FIELD IS ROTATED BY 90 DEGS FOR M=1 AND BY 45 DEGS FOR M=2
     C
     C   ADDS FIELDS FROM TWO CONES (NORTHERN AND SOUTHERN), WITH A PROPER SYMMETRY OF THE CURRENT AND FIELD,
     C     CORRESPONDING TO THE REGION 1 BIRKELAND CURRENTS. (SEE NB #6, P.58).
     */
    static void TWOCONSS (double const A[],
                          double const X,double const Y,double const Z,
                          double *BX,double *BY,double *BZ,
                          struct TS_EXTERNAL_COMMON_ const* TS_EXTERNAL_COMMON) {
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION A(31)
        int const* M=&TS_EXTERNAL_COMMON->MODENUM_M;//COMMON /MODENUM/ M

        static double const HSQR2=0.707106781e0;

        double XAS,YAS;
        if (*M == 1) {//   !   ROTATION BY 90 DEGS
            XAS = Y;
            YAS =-X;
        } else {//ELSE               !   ROTATION BY 45 DEGS
            XAS = (X+Y)*HSQR2;
            YAS = (Y-X)*HSQR2;
        }//ENDIF

        double BXN,BYN,BZN;
        ONE_CONE (A,XAS,YAS,Z,&BXN,&BYN,&BZN,TS_EXTERNAL_COMMON);
        double BXS,BYS,BZS;
        ONE_CONE (A,XAS,-YAS,-Z,&BXS,&BYS,&BZS,TS_EXTERNAL_COMMON);

        double BXAS=BXN-BXS;
        double BYAS=BYN+BYS;
        *BZ=BZN+BZS;

        if (*M == 1) {//   !   ROTATION BY 90 DEGS
            *BX =-BYAS;
            *BY = BXAS;
        } else {//ELSE
            *BX=(BXAS-BYAS)*HSQR2;
            *BY=(BXAS+BYAS)*HSQR2;
        }//ENDIF

        return;
    }//END

    /*
     C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     C
     SUBROUTINE BIRSH_SY (A,PS,X_SC,X,Y,Z,BX,BY,BZ)
     C
     C   THIS S/R IS QUITE SIMILAR TO BIRK_SHL, BUT IT IS FOR THE SYMMETRIC MODE OF BIRKELAND CURRENT FIELD
     C     AND FOR THAT REASON THE FIELD COMPONENTS HAVE A DIFFERENT KIND OF SYMMETRY WITH RESPECT TO Y_gsm
     */
    static void BIRSH_SY (double const A[],double const PS,double const X_SC,
                          double const X,double const Y,double const Z,
                          double *BX,double *BY,double *BZ) {
        //IMPLICIT  REAL * 8  (A - H, O - Z)
        //DIMENSION A(86)

        double CPS=cos(PS);
        double SPS=sin(PS);

        double S3PS=2.e0*CPS;

        double PST1=PS*A[85];
        double PST2=PS*A[86];

        double ST1=sin(PST1);
        double CT1=cos(PST1);
        double ST2=sin(PST2);
        double CT2=cos(PST2);

        double X1=X*CT1-Z*ST1;
        double Z1=X*ST1+Z*CT1;
        double X2=X*CT2-Z*ST2;
        double Z2=X*ST2+Z*CT2;

        int L=0;
        double GX=0.e0;
        double GY=0.e0;
        double GZ=0.e0;

        double P,Q,CYPI,CYQI,SYPI,SYQI,R,S,SZRK,CZSK,CZRK,SZSK,SQPR,SQQS,EPR,EQS,FX,FY,FZ,HX,HY,HZ,HXR,HZR;
        for (int M=1; M<=2; M++) {//DO 1 M=1,2     !    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
            //                          AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
            for (int I=1; I<=3; I++) {//DO 2 I=1,3
                P=A[72+I];
                Q=A[78+I];
                CYPI=cos(Y/P);
                CYQI=cos(Y/Q);
                SYPI=sin(Y/P);
                SYQI=sin(Y/Q);

                for (int K=1; K<=3; K++) {//DO 3 K=1,3
                    R=A[75+K];
                    S=A[81+K];
                    SZRK=sin(Z1/R);
                    CZSK=cos(Z2/S);
                    CZRK=cos(Z1/R);
                    SZSK=sin(Z2/S);
                    SQPR=sqrt(1.e0/(P*P)+1.e0/(R*R));
                    SQQS=sqrt(1.e0/(Q*Q)+1.e0/(S*S));
                    EPR=exp(X1*SQPR);
                    EQS=exp(X2*SQQS);

                    for (int N=1; N<=2; N++) {//DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
                        //                                AND N=2 IS FOR THE SECOND ONE

                        for (int NN=1; NN<=2; NN++) {//DO 5 NN=1,2 !   NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
                            //                                         TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE

                            if (M == 1) {
                                FX= SQPR*EPR*SYPI*SZRK;
                                FY=EPR*CYPI*SZRK/P;
                                FZ= EPR*SYPI*CZRK/R;
                                if (N == 1) {
                                    if (NN == 1) {
                                        HX=FX;
                                        HY=FY;
                                        HZ=FZ;
                                    } else {//ELSE
                                        HX=FX*X_SC;
                                        HY=FY*X_SC;
                                        HZ=FZ*X_SC;
                                    }//ENDIF
                                } else {//ELSE
                                    if (NN == 1) {
                                        HX=FX*CPS;
                                        HY=FY*CPS;
                                        HZ=FZ*CPS;
                                    } else {//ELSE
                                        HX=FX*CPS*X_SC;
                                        HY=FY*CPS*X_SC;
                                        HZ=FZ*CPS*X_SC;
                                    }//ENDIF
                                }//ENDIF

                            } else {//ELSE                            !   M.EQ.2
                                FX= SPS*SQQS*EQS*SYQI*CZSK;
                                FY= SPS/Q*EQS*CYQI*CZSK;
                                FZ=-SPS/S*EQS*SYQI*SZSK;
                                if (N == 1) {
                                    if (NN == 1) {
                                        HX=FX;
                                        HY=FY;
                                        HZ=FZ;
                                    } else {//ELSE
                                        HX=FX*X_SC;
                                        HY=FY*X_SC;
                                        HZ=FZ*X_SC;
                                    }//ENDIF
                                } else {//ELSE
                                    if (NN == 1) {
                                        HX=FX*S3PS;
                                        HY=FY*S3PS;
                                        HZ=FZ*S3PS;
                                    } else {//ELSE
                                        HX=FX*S3PS*X_SC;
                                        HY=FY*S3PS*X_SC;
                                        HZ=FZ*S3PS*X_SC;
                                    }//ENDIF
                                }//ENDIF
                            }//ENDIF
                            L=L+1;

                            if (M == 1) {
                                HXR=HX*CT1+HZ*ST1;
                                HZR=-HX*ST1+HZ*CT1;
                            } else {//ELSE
                                HXR=HX*CT2+HZ*ST2;
                                HZR=-HX*ST2+HZ*CT2;
                            }//ENDIF

                            GX=GX+HXR*A[L];
                            GY=GY+HY *A[L];
                            GZ=GZ+HZR*A[L];  // 5
                        }

                    }//4   CONTINUE
                }//3   CONTINUE
            }//2   CONTINUE
        }//1   CONTINUE

        *BX=GX;
        *BY=GY;
        *BZ=GZ;

        return;
    }//END

    /*
     c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     c
     SUBROUTINE BIR1N2SY (NUMB,MODE,PS,X,Y,Z,BX,BY,BZ)        !   SEE NB# 6, P.60 and NB#7, P.35-...
     C
     C   THIS CODE IS VERY SIMILAR TO BIRK_1N2, BUT IT IS FOR THE "SYMMETRICAL" MODE, IN WHICH J_parallel
     C     IS A SYMMETRIC (EVEN) FUNCTION OF Ygsm
     C
     C  CALCULATES COMPONENTS  OF REGION 1/2 FIELD IN SPHERICAL COORDS.  DERIVED FROM THE S/R DIPDEF2C (WHICH
     C    DOES THE SAME JOB, BUT INPUT/OUTPUT THERE WAS IN SPHERICAL COORDS, WHILE HERE WE USE CARTESIAN ONES)
     C
     C   INPUT:  NUMB=1 (2) FOR REGION 1 (2) CURRENTS
     C           MODE=1 YIELDS SIMPLE SINUSOIDAL MLT VARIATION, WITH MAXIMUM CURRENT AT DAWN/DUSK MERIDIAN
     C     WHILE MODE=2 YIELDS THE SECOND HARMONIC.
     C
     */
    static void BIR1N2SY (int const NUMB,int const MODE,
                          double const PS,double const X,double const Y,double const Z,
                          double *BX,double *BY,double *BZ,
                          struct TS_EXTERNAL_COMMON_ *TS_EXTERNAL_COMMON) {//        !   SEE NB# 6, P.60 and NB#7, P.35-...
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION A11(31),A12(31),A21(31),A22(31)
        int *M=&TS_EXTERNAL_COMMON->MODENUM_M;//COMMON /MODENUM/ M
        double *DTHETA=&TS_EXTERNAL_COMMON->DTHETA_DTHETA;//COMMON /DTHETA/ DTHETA

        double *DPHI=&TS_EXTERNAL_COMMON->DPHI_B_RHO0_DPHI,*B=&TS_EXTERNAL_COMMON->DPHI_B_RHO0_B,*RHO_0=&TS_EXTERNAL_COMMON->DPHI_B_RHO0_RHO_0,*XKAPPA=&TS_EXTERNAL_COMMON->DPHI_B_RHO0_XKAPPA;//COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! THESE PARAMETERS CONTROL DAY-NIGHT ASYMMETRY OF F.A.C., AS FOLLOWS:
        /*
         C  (1) DPHI:   HALF-DIFFERENCE (IN RADIANS) BETWEEN DAY AND NIGHT LATITUDE OF FAC OVAL AT IONOSPHERIC ALTITUDE;
         C              TYPICAL VALUE: 0.06
         C  (2) B:      AN ASYMMETRY FACTOR AT HIGH-ALTITUDES;  FOR B=0, THE ONLY ASYMMETRY IS THAT FROM DPHI
         C              TYPICAL VALUES: 0.35-0.70
         C  (3) RHO_0:  A FIXED PARAMETER, DEFINING THE DISTANCE RHO, AT WHICH THE LATITUDE SHIFT GRADUALLY SATURATES AND
         C              STOPS INCREASING
         C              ITS VALUE WAS ASSUMED FIXED, EQUAL TO 7.0.
         C  (4) XKAPPA: AN OVERALL SCALING FACTOR, WHICH CAN BE USED FOR CHANGING THE SIZE OF THE F.A.C. OVAL
         */
        static double const BETA=0.9e0,RH=10.e0,EPS=3.e0;// ! parameters of the tilt-dependent deformation of the untilted F.A.C. field

        static double const A11[32] = {NAN/*dummy*/,.1618068350,-.1797957553,2.999642482,-.9322708978,
            -.6811059760,.2099057262,-8.358815746,-14.86033550,.3838362986,
            -16.30945494,4.537022847,2.685836007,27.97833029,6.330871059,
            1.876532361,18.95619213,.9651528100,.4217195118,-.08957770020,
            -1.823555887,.7457045438,-.5785916524,-1.010200918,.01112389357,
            .09572927448,-.3599292276,8.713700514,.9763932955,3.834602998,
            2.492118385,.7113544659};
        static double const A12[32] = {NAN/*dummy*/,.7058026940,-.2845938535,5.715471266,-2.472820880,
            -.7738802408,.3478293930,-11.37653694,-38.64768867,.6932927651,
            -212.4017288,4.944204937,3.071270411,33.05882281,7.387533799,
            2.366769108,79.22572682,.6154290178,.5592050551,-.1796585105,
            -1.654932210,.7309108776,-.4926292779,-1.130266095,-.009613974555,
            .1484586169,-.2215347198,7.883592948,.02768251655,2.950280953,
            1.212634762,.5567714182};
        static double const A21[32] = {NAN/*dummy*/,.1278764024,-.2320034273,1.805623266,-32.37241440,
            -.9931490648,.3175085630,-2.492465814,-16.21600096,.2695393416,
            -6.752691265,3.971794901,14.54477563,41.10158386,7.912889730,
            1.258297372,9.583547721,1.014141963,.5104134759,-.1790430468,
            -1.756358428,.7561986717,-.6775248254,-.04014016420,.01446794851,
            .1200521731,-.2203584559,4.508963850,.8221623576,1.779933730,
            1.102649543,.8867880020};
        static double const A22[32] = {NAN/*dummy*/,.4036015198,-.3302974212,2.827730930,-45.44405830,
            -1.611103927,.4927112073,-.003258457559,-49.59014949,.3796217108,
            -233.7884098,4.312666980,18.05051709,28.95320323,11.09948019,
            .7471649558,67.10246193,.5667096597,.6468519751,-.1560665317,
            -1.460805289,.7719653528,-.6658988668,.2515179349E-05,
            .02426021891,.1195003324,-.2625739255,4.377172556,.2421190547,
            2.503482679,1.071587299,.7247997430};

        *B=0.5;
        *RHO_0=7.0;

        *M=MODE;
        if (NUMB == 1) {//IF (NUMB.EQ.1) THEN
            *DPHI=0.055e0;
            *DTHETA=0.06e0;
        }//ENDIF

        if (NUMB == 2) {//IF (NUMB.EQ.2) THEN
            *DPHI=0.030e0;
            *DTHETA=0.09e0;
        }//ENDIF

        double Xsc=X**XKAPPA;
        double Ysc=Y**XKAPPA;
        double Zsc=Z**XKAPPA;
        double RHO=sqrt(Xsc*Xsc+Zsc*Zsc);

        double Rsc=sqrt(Xsc*Xsc+Ysc*Ysc+Zsc*Zsc);//                                 !  SCALED
        double RHO2=*RHO_0**RHO_0;
        double PHI;
        if (Xsc == 0.e0 && Zsc == 0.e0) {//IF (Xsc.EQ.0.D0.AND.Zsc.EQ.0.D0) THEN
            PHI=0.e0;
        } else {//ELSE
            PHI=atan2(-Zsc,Xsc);//  !  FROM CARTESIAN TO CYLINDRICAL (RHO,PHI,Y)
        }//ENDIF

        double SPHIC=sin(PHI);
        double CPHIC=cos(PHI);//  !  "C" means "CYLINDRICAL", TO DISTINGUISH FROM SPHERICAL PHI

        double BRACK=*DPHI+*B*RHO2/(RHO2+1.e0)*(RHO*RHO-1.e0)/(RHO2+RHO*RHO);
        double R1RH=(Rsc-1.e0)/RH;
        double PSIAS=BETA*PS/pow((1.e0+pow(R1RH, EPS)), (1.e0/EPS));

        double PHIS=PHI-BRACK*sin(PHI) -PSIAS;
        double DPHISPHI=1.e0-BRACK*cos(PHI);
        double DPHISRHO=-2.e0**B*RHO2*RHO/pow((RHO2+RHO*RHO), 2)*sin(PHI)+BETA*PS*pow(R1RH, (EPS-1.e0))*RHO/(RH*Rsc*pow((1.e0+pow(R1RH, EPS)), (1.e0/EPS+1.e0)));
        double DPHISDY= BETA*PS*pow(R1RH, (EPS-1.e0))*Ysc/(RH*Rsc*pow((1.e0+pow(R1RH, EPS)), (1.e0/EPS+1.e0)));

        double SPHICS=sin(PHIS);
        double CPHICS=cos(PHIS);

        double XS= RHO*CPHICS;
        double ZS=-RHO*SPHICS;

        double BXS,BYAS,BZS;
        if (NUMB == 1) {
            if (MODE == 1) TWOCONSS (A11,XS,Ysc,ZS,&BXS,&BYAS,&BZS,TS_EXTERNAL_COMMON);
            if (MODE == 2) TWOCONSS (A12,XS,Ysc,ZS,&BXS,&BYAS,&BZS,TS_EXTERNAL_COMMON);
        } else {//ELSE
            if (MODE == 1) TWOCONSS (A21,XS,Ysc,ZS,&BXS,&BYAS,&BZS,TS_EXTERNAL_COMMON);
            if (MODE == 2) TWOCONSS (A22,XS,Ysc,ZS,&BXS,&BYAS,&BZS,TS_EXTERNAL_COMMON);
        }//ENDIF

        double BRHOAS=BXS*CPHICS-BZS*SPHICS;
        double BPHIAS=-BXS*SPHICS-BZS*CPHICS;

        double BRHO_S=BRHOAS*DPHISPHI                             **XKAPPA;//        ! SCALING
        double BPHI_S=(BPHIAS-RHO*(BYAS*DPHISDY+BRHOAS*DPHISRHO)) **XKAPPA;
        double BY_S=BYAS*DPHISPHI                                 **XKAPPA;

        *BX=BRHO_S*CPHIC-BPHI_S*SPHIC;
        *BY=BY_S;
        *BZ=-BRHO_S*SPHIC-BPHI_S*CPHIC;

        return;
    }//END

    /*
     C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     C
     SUBROUTINE BIRTOTSY (PS,X,Y,Z,BX11,BY11,BZ11,BX12,BY12,BZ12,
     BX21,BY21,BZ21,BX22,BY22,BZ22)
     C
     C   THIS S/R IS ALMOST IDENTICAL TO BIRK_TOT, BUT IT IS FOR THE SYMMETRIC MODE, IN WHICH
     C     J_parallel IS AN EVEN FUNCTION OF Ygsm.
     C
     C
     C      IOPBS -  BIRKELAND FIELD MODE FLAG:
     C         IOPBS=0 - ALL COMPONENTS
     C         IOPBS=1 - REGION 1, MODES 1 & 2 (SYMMETRIC !)
     C         IOPBS=2 - REGION 2, MODES 1 & 2 (SYMMETRIC !)
     */
    static void BIRTOTSY (double const PS,double const X,double const Y,double const Z,
                          double *BX11,double *BY11,double *BZ11,double *BX12,double *BY12,double *BZ12,
                          double *BX21,double *BY21,double *BZ21,double *BX22,double *BY22,double *BZ22,
                          struct TS_EXTERNAL_COMMON_ *TS_EXTERNAL_COMMON) {
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION SH11(86),SH12(86),SH21(86),SH22(86)
        double const* XKAPPA1=&TS_EXTERNAL_COMMON->BIRKPAR_XKAPPA1,*XKAPPA2=&TS_EXTERNAL_COMMON->BIRKPAR_XKAPPA2;//COMMON /BIRKPAR/ XKAPPA1,XKAPPA2   !  INPUT PARAMETERS, SPECIFIED FROM A MAIN PROGRAM
        //                                            (JOINT WITH  BIRK_TOT  FOR THE ANTISYMMETRICAL MODE)

        double /*&DPHI=ctx.DPHI_B_RHO0_DPHI,&B=ctx.DPHI_B_RHO0_B,&RHO_0=ctx.DPHI_B_RHO0_RHO_0,*/
        *XKAPPA=&TS_EXTERNAL_COMMON->DPHI_B_RHO0_XKAPPA;//COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! PARAMETERS, CONTROL DAY-NIGHT ASYMMETRY OF F.A.C.

        static double const SH11[87] = {NAN/*dummy*/,4956703.683,-26922641.21,-11383659.85,29604361.65,
            -38919785.97,70230899.72,34993479.24,-90409215.02,30448713.69,
            -48360257.19,-35556751.23,57136283.60,-8013815.613,30784907.86,
            13501620.50,-35121638.52,50297295.45,-84200377.18,-46946852.58,
            107526898.8,-39003263.47,59465850.17,47264335.10,-68892388.73,
            3375901.533,-9181255.754,-4494667.217,10812618.51,-17351920.97,
            27016083.00,18150032.11,-33186882.96,13340198.63,-19779685.30,
            -17891788.15,21625767.23,16135.32442,133094.0241,-13845.61859,
            -79159.98442,432.1215298,-85438.10368,1735.386707,41891.71284,
            18158.14923,-105465.8135,-11685.73823,62297.34252,-10811.08476,
            -87631.38186,9217.499261,52079.94529,-68.29127454,56023.02269,
            -1246.029857,-27436.42793,-11972.61726,69607.08725,7702.743803,
            -41114.36810,12.08269108,-21.30967022,-9.100782462,18.26855933,
            -7.000685929,26.22390883,6.392164144,-21.99351743,2.294204157,
            -16.10023369,-1.344314750,9.342121230,148.5493329,99.79912328,
            70.78093196,35.23177574,47.45346891,58.44877918,139.8135237,
            91.96485261,6.983488815,9.055554871,19.80484284,2.860045019,
            .8213262337E-01,-.7962186676E-05};

        static double const SH12[87] = {NAN/*dummy*/,-1210748.720,-52324903.95,-14158413.33,19426123.60,
            6808641.947,-5138390.983,-1118600.499,-4675055.459,2059671.506,
            -1373488.052,-114704.4353,-1435920.472,1438451.655,61199067.17,
            16549301.39,-22802423.47,-7814550.995,5986478.728,1299443.190,
            5352371.724,-2994351.520,1898553.337,203158.3658,2270182.134,
            -618083.3112,-25950806.16,-7013783.326,9698792.575,3253693.134,
            -2528478.464,-546323.4095,-2217735.237,1495336.589,-914647.4222,
            -114374.1054,-1200441.634,-507068.4700,1163189.975,998411.8381,
            -861919.3631,5252210.872,-11668550.16,-4113899.385,6972900.950,
            -2546104.076,7704014.310,2273077.192,-5134603.198,256205.7901,
            -589970.8086,-503821.0170,437612.8956,-2648640.128,5887640.735,
            2074286.234,-3519291.144,1283847.104,-3885817.147,-1145936.942,
            2589753.651,-408.7788403,1234.054185,739.8541716,-965.8068853,
            3691.383679,-8628.635819,-2855.844091,5268.500178,-1774.372703,
            5515.010707,1556.089289,-3665.434660,204.8672197,110.7748799,
            87.36036207,5.522491330,31.06364270,73.57632579,281.5331360,
            140.3461448,17.07537768,6.729732641,4.100970449,2.780422877,
            .8742978101E-01,-.1028562327E-04};

        static double const SH21[87] = {NAN/*dummy*/,-67763516.61,-49565522.84,10123356.08,51805446.10,
            -51607711.68,164360662.1,-4662006.024,-191297217.6,-7204547.103,
            30372354.93,-750371.9365,-36564457.17,61114395.65,45702536.50,
            -9228894.939,-47893708.68,47290934.33,-149155112.0,4226520.638,
            173588334.5,7998505.443,-33150962.72,832493.2094,39892545.84,
            -11303915.16,-8901327.398,1751557.110,9382865.820,-9054707.868,
            27918664.50,-788741.7146,-32481294.42,-2264443.753,9022346.503,
            -233526.0185,-10856269.53,-244450.8850,1908295.272,185445.1967,
            -1074202.863,41827.75224,-241553.7626,-20199.12580,123235.6084,
            199501.4614,-1936498.464,-178857.4074,1044724.507,121044.9917,
            -946479.9247,-91808.28803,532742.7569,-20742.28628,120633.2193,
            10018.49534,-61599.11035,-98709.58977,959095.1770,88500.43489,
            -517471.5287,-81.56122911,816.2472344,55.30711710,-454.5368824,
            25.74693810,-202.5007350,-7.369350794,104.9429812,58.14049362,
            -685.5919355,-51.71345683,374.0125033,247.9296982,159.2471769,
            102.3151816,15.81062488,34.99767599,133.0832773,219.6475201,
            107.9582783,10.00264684,7.718306072,25.22866153,5.013583103,
            .8407754233E-01,-.9613356793E-05};

        static double const SH22[87] = {NAN/*dummy*/,-43404887.31,8896854.538,-8077731.036,-10247813.65,
            6346729.086,-9416801.212,-1921670.268,7805483.928,2299301.127,
            4856980.170,-1253936.462,-4695042.690,54305735.91,-11158768.10,
            10051771.85,12837129.47,-6380785.836,12387093.50,1687850.192,
            -10492039.47,-5777044.862,-6916507.424,2855974.911,7027302.490,
            -26176628.93,5387959.610,-4827069.106,-6193036.589,2511954.143,
            -6205105.083,-553187.2984,5341386.847,3823736.361,3669209.068,
            -1841641.700,-3842906.796,281561.7220,-5013124.630,379824.5943,
            2436137.901,-76337.55394,548518.2676,42134.28632,-281711.3841,
            -365514.8666,-2583093.138,-232355.8377,1104026.712,-131536.3445,
            2320169.882,-174967.6603,-1127251.881,35539.82827,-256132.9284,
            -19620.06116,131598.7965,169033.6708,1194443.500,107320.3699,
            -510672.0036,1211.177843,-17278.19863,1140.037733,8347.612951,
            -303.8408243,2405.771304,174.0634046,-1248.722950,-1231.229565,
            -8666.932647,-754.0488385,3736.878824,227.2102611,115.9154291,
            94.34364830,3.625357304,64.03192907,109.0743468,241.4844439,
            107.7583478,22.36222385,6.282634037,27.79399216,2.270602235,
            .8708605901E-01,-.1256706895E-04};

        *XKAPPA=*XKAPPA1;//        !  FORWARDED IN BIR1N2SY
        double X_SC=*XKAPPA1-1.1e0;//    !  FORWARDED IN BIRSH_SY
        double FX11,FY11,FZ11,HX11,HY11,HZ11;
        BIR1N2SY (1,1,PS,X,Y,Z,&FX11,&FY11,&FZ11,TS_EXTERNAL_COMMON);//           !  REGION 1, MODE 1
        BIRSH_SY (SH11,PS,X_SC,X,Y,Z,&HX11,&HY11,&HZ11);

        *BX11=FX11+HX11;
        *BY11=FY11+HY11;
        *BZ11=FZ11+HZ11;

        double FX12,FY12,FZ12,HX12,HY12,HZ12;
        BIR1N2SY (1,2,PS,X,Y,Z,&FX12,&FY12,&FZ12,TS_EXTERNAL_COMMON);//           !  REGION 1, MODE 2
        BIRSH_SY (SH12,PS,X_SC,X,Y,Z,&HX12,&HY12,&HZ12);

        *BX12=FX12+HX12;
        *BY12=FY12+HY12;
        *BZ12=FZ12+HZ12;

        *XKAPPA=*XKAPPA2;//        !  FORWARDED IN BIR1N2SY
        X_SC=*XKAPPA2-1.0e0;//    !  FORWARDED IN BIRSH_SY

        double FX21,FY21,FZ21,HX21,HY21,HZ21;
        BIR1N2SY (2,1,PS,X,Y,Z,&FX21,&FY21,&FZ21,TS_EXTERNAL_COMMON);//           !  REGION 2, MODE 1
        BIRSH_SY (SH21,PS,X_SC,X,Y,Z,&HX21,&HY21,&HZ21);

        *BX21=FX21+HX21;
        *BY21=FY21+HY21;
        *BZ21=FZ21+HZ21;

        double FX22,FY22,FZ22,HX22,HY22,HZ22;
        BIR1N2SY (2,2,PS,X,Y,Z,&FX22,&FY22,&FZ22,TS_EXTERNAL_COMMON);//           !  REGION 2, MODE 2
        BIRSH_SY (SH22,PS,X_SC,X,Y,Z,&HX22,&HY22,&HZ22);

        *BX22=FX22+HX22;
        *BY22=FY22+HY22;
        *BZ22=FZ22+HZ22;

        return;
    }//END

    /*
     C=========================================================================
     c
     SUBROUTINE TWOCONES (A,X,Y,Z,BX,BY,BZ)
     C
     C   ADDS FIELDS FROM TWO CONES (NORTHERN AND SOUTHERN), WITH A PROPER SYMMETRY OF THE CURRENT AND FIELD,
     C     CORRESPONDING TO THE REGION 1 BIRKELAND CURRENTS. (SEE NB #6, P.58).
     */
    static void TWOCONES (double const A[],double const X,double const Y,double const Z,
                          double *BX,double *BY,double *BZ,
                          struct TS_EXTERNAL_COMMON_ const* TS_EXTERNAL_COMMON) {
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION A(31)
        double BXN,BYN,BZN,BXS,BYS,BZS;
        ONE_CONE (A,X,Y,Z,&BXN,&BYN,&BZN,TS_EXTERNAL_COMMON);
        ONE_CONE (A,X,-Y,-Z,&BXS,&BYS,&BZS,TS_EXTERNAL_COMMON);
        *BX=BXN-BXS;
        *BY=BYN+BYS;
        *BZ=BZN+BZS;

        return;
    }//END

    /*
     c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     c
     c
     SUBROUTINE BIRK_1N2 (NUMB,MODE,PS,X,Y,Z,BX,BY,BZ)        !   SEE NB# 6, P.60
     C
     C  CALCULATES COMPONENTS  OF REGION 1/2 FIELD IN SPHERICAL COORDS.  DERIVED FROM THE S/R DIPDEF2C (WHICH
     C    DOES THE SAME JOB, BUT INPUT/OUTPUT THERE WAS IN SPHERICAL COORDS, WHILE HERE WE USE CARTESIAN ONES)
     C
     C   INPUT:  NUMB=1 (2) FOR REGION 1 (2) CURRENTS
     C           MODE=1 YIELDS SIMPLE SINUSOIDAL MLT VARIATION, WITH MAXIMUM CURRENT AT DAWN/DUSK MERIDIAN
     C     WHILE MODE=2 YIELDS THE SECOND HARMONIC.
     C
     */
    static void BIRK_1N2 (int const NUMB,int const MODE,
                          double const PS,double const X,double const Y,double const Z,
                          double *BX,double *BY,double *BZ,
                          struct TS_EXTERNAL_COMMON_ *TS_EXTERNAL_COMMON) {//        !   SEE NB# 6, P.60
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION A11(31),A12(31),A21(31),A22(31)
        int *M=&TS_EXTERNAL_COMMON->MODENUM_M;//COMMON /MODENUM/ M
        double *DTHETA=&TS_EXTERNAL_COMMON->DTHETA_DTHETA;//COMMON /DTHETA/ DTHETA

        double *DPHI=&TS_EXTERNAL_COMMON->DPHI_B_RHO0_DPHI,*B=&TS_EXTERNAL_COMMON->DPHI_B_RHO0_B,*RHO_0=&TS_EXTERNAL_COMMON->DPHI_B_RHO0_RHO_0,*XKAPPA=&TS_EXTERNAL_COMMON->DPHI_B_RHO0_XKAPPA;//COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! THESE PARAMETERS CONTROL DAY-NIGHT ASYMMETRY OF F.A.C., AS FOLLOWS:
        /*
         C  (1) DPHI:   HALF-DIFFERENCE (IN RADIANS) BETWEEN DAY AND NIGHT LATITUDE OF FAC OVAL AT IONOSPHERIC ALTITUDE;
         C              TYPICAL VALUE: 0.06
         C  (2) B:      AN ASYMMETRY FACTOR AT HIGH-ALTITUDES;  FOR B=0, THE ONLY ASYMMETRY IS THAT FROM DPHI
         C              TYPICAL VALUES: 0.35-0.70
         C  (3) RHO_0:  A FIXED PARAMETER, DEFINING THE DISTANCE RHO, AT WHICH THE LATITUDE SHIFT GRADUALLY SATURATES AND
         C              STOPS INCREASING
         C              ITS VALUE WAS ASSUMED FIXED, EQUAL TO 7.0.
         C  (4) XKAPPA: AN OVERALL SCALING FACTOR, WHICH CAN BE USED FOR CHANGING THE SIZE OF THE F.A.C. OVAL
         */
        static double const BETA=0.9e0,RH=10.e0,EPS=3.e0;// ! parameters of the tilt-dependent deformation of the untilted F.A.C. field

        static double const A11[32] = {NAN/*dummy*/,.1618068350,-.1797957553,2.999642482,-.9322708978,
            -.6811059760,.2099057262,-8.358815746,-14.86033550,.3838362986,
            -16.30945494,4.537022847,2.685836007,27.97833029,6.330871059,
            1.876532361,18.95619213,.9651528100,.4217195118,-.08957770020,
            -1.823555887,.7457045438,-.5785916524,-1.010200918,.01112389357,
            .09572927448,-.3599292276,8.713700514,.9763932955,3.834602998,
            2.492118385,.7113544659};
        static double const A12[32] = {NAN/*dummy*/,.7058026940,-.2845938535,5.715471266,-2.472820880,
            -.7738802408,.3478293930,-11.37653694,-38.64768867,.6932927651,
            -212.4017288,4.944204937,3.071270411,33.05882281,7.387533799,
            2.366769108,79.22572682,.6154290178,.5592050551,-.1796585105,
            -1.654932210,.7309108776,-.4926292779,-1.130266095,-.009613974555,
            .1484586169,-.2215347198,7.883592948,.02768251655,2.950280953,
            1.212634762,.5567714182};
        static double const A21[32] = {NAN/*dummy*/,.1278764024,-.2320034273,1.805623266,-32.37241440,
            -.9931490648,.3175085630,-2.492465814,-16.21600096,.2695393416,
            -6.752691265,3.971794901,14.54477563,41.10158386,7.912889730,
            1.258297372,9.583547721,1.014141963,.5104134759,-.1790430468,
            -1.756358428,.7561986717,-.6775248254,-.04014016420,.01446794851,
            .1200521731,-.2203584559,4.508963850,.8221623576,1.779933730,
            1.102649543,.8867880020};
        static double const A22[32] = {NAN/*dummy*/,.4036015198,-.3302974212,2.827730930,-45.44405830,
            -1.611103927,.4927112073,-.003258457559,-49.59014949,.3796217108,
            -233.7884098,4.312666980,18.05051709,28.95320323,11.09948019,
            .7471649558,67.10246193,.5667096597,.6468519751,-.1560665317,
            -1.460805289,.7719653528,-.6658988668,.2515179349E-05,
            .02426021891,.1195003324,-.2625739255,4.377172556,.2421190547,
            2.503482679,1.071587299,.7247997430};

        *B=0.5;
        *RHO_0=7.0;

        *M=MODE;
        if (NUMB == 1) {//IF (NUMB.EQ.1) THEN
            *DPHI=0.055e0;
            *DTHETA=0.06e0;
        }//ENDIF

        if (NUMB == 2) {//IF (NUMB.EQ.2) THEN
            *DPHI=0.030e0;
            *DTHETA=0.09e0;
        }//ENDIF

        double Xsc=X**XKAPPA;
        double Ysc=Y**XKAPPA;
        double Zsc=Z**XKAPPA;
        double RHO=sqrt(Xsc*Xsc+Zsc*Zsc);

        double Rsc=sqrt(Xsc*Xsc+Ysc*Ysc+Zsc*Zsc);//                                 !  SCALED
        double RHO2=*RHO_0**RHO_0;
        double PHI;
        if (Xsc == 0.e0 && Zsc == 0.e0) {//IF (Xsc.EQ.0.D0.AND.Zsc.EQ.0.D0) THEN
            PHI=0.e0;
        }else {//ELSE
            PHI=atan2(-Zsc,Xsc);//  !  FROM CARTESIAN TO CYLINDRICAL (RHO,PHI,Y)
        }//ENDIF

        double SPHIC=sin(PHI);
        double CPHIC=cos(PHI);//  !  "C" means "CYLINDRICAL", TO DISTINGUISH FROM SPHERICAL PHI

        double BRACK=*DPHI+*B*RHO2/(RHO2+1.e0)*(RHO*RHO-1.e0)/(RHO2+RHO*RHO);
        double R1RH=(Rsc-1.e0)/RH;
        double PSIAS=BETA*PS/pow((1.e0+pow(R1RH, EPS)), (1.e0/EPS));

        double PHIS=PHI-BRACK*sin(PHI) -PSIAS;
        double DPHISPHI=1.e0-BRACK*cos(PHI);
        double DPHISRHO=-2.e0**B*RHO2*RHO/pow((RHO2+RHO*RHO), 2) *sin(PHI)+BETA*PS*pow(R1RH, (EPS-1.e0))*RHO/(RH*Rsc*pow((1.e0+pow(R1RH, EPS)), (1.e0/EPS+1.e0)));
        double DPHISDY= BETA*PS*pow(R1RH, (EPS-1.e0))*Ysc/(RH*Rsc*pow((1.e0+pow(R1RH, EPS)), (1.e0/EPS+1.e0)));

        double SPHICS=sin(PHIS);
        double CPHICS=cos(PHIS);

        double XS= RHO*CPHICS;
        double ZS=-RHO*SPHICS;

        double BXS,BYAS,BZS;
        if (NUMB == 1) {//IF (NUMB.EQ.1) THEN
            if (MODE == 1) TWOCONES (A11,XS,Ysc,ZS,&BXS,&BYAS,&BZS,TS_EXTERNAL_COMMON);
            if (MODE == 2) TWOCONES (A12,XS,Ysc,ZS,&BXS,&BYAS,&BZS,TS_EXTERNAL_COMMON);
        } else {//ELSE
            if (MODE == 1) TWOCONES (A21,XS,Ysc,ZS,&BXS,&BYAS,&BZS,TS_EXTERNAL_COMMON);
            if (MODE == 2) TWOCONES (A22,XS,Ysc,ZS,&BXS,&BYAS,&BZS,TS_EXTERNAL_COMMON);
        }//ENDIF

        double BRHOAS=BXS*CPHICS-BZS*SPHICS;
        double BPHIAS=-BXS*SPHICS-BZS*CPHICS;

        double BRHO_S=BRHOAS*DPHISPHI                             **XKAPPA;//        ! SCALING
        double BPHI_S=(BPHIAS-RHO*(BYAS*DPHISDY+BRHOAS*DPHISRHO)) **XKAPPA;
        double BY_S=BYAS*DPHISPHI                                 **XKAPPA;

        *BX=BRHO_S*CPHIC-BPHI_S*SPHIC;
        *BY=BY_S;
        *BZ=-BRHO_S*SPHIC-BPHI_S*CPHIC;

        return;
    }//END

    /*
     C-------------------------------------------------------------------------
     C
     */
    static void BIRK_SHL (double const A[],double const PS,double const X_SC,
                          double const X,double const Y,double const Z,
                          double *BX,double *BY,double *BZ) {

        //IMPLICIT  REAL * 8  (A - H, O - Z)
        //DIMENSION A(86)

        double CPS=cos(PS);
        double SPS=sin(PS);

        double S3PS=2.e0*CPS;

        double PST1=PS*A[85];
        double PST2=PS*A[86];

        double ST1=sin(PST1);
        double CT1=cos(PST1);
        double ST2=sin(PST2);
        double CT2=cos(PST2);

        double X1=X*CT1-Z*ST1;
        double Z1=X*ST1+Z*CT1;
        double X2=X*CT2-Z*ST2;
        double Z2=X*ST2+Z*CT2;

        int L=0;
        double GX=0.e0;
        double GY=0.e0;
        double GZ=0.e0;

        double P,Q,CYPI,CYQI,SYPI,SYQI,R,S,SZRK,CZSK,CZRK,SZSK,SQPR,SQQS,EPR,EQS,FX,FY,FZ,HX,HY,HZ,HXR,HZR;
        for (int M=1; M<=2; M++) {//DO 1 M=1,2     !    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
            //                          AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
            for (int I=1; I<=3; I++) {//DO 2 I=1,3
                P=A[72+I];
                Q=A[78+I];
                CYPI=cos(Y/P);
                CYQI=cos(Y/Q);
                SYPI=sin(Y/P);
                SYQI=sin(Y/Q);

                for (int K=1; K<=3; K++) {//DO 3 K=1,3
                    R=A[75+K];
                    S=A[81+K];
                    SZRK=sin(Z1/R);
                    CZSK=cos(Z2/S);
                    CZRK=cos(Z1/R);
                    SZSK=sin(Z2/S);
                    SQPR=sqrt(1.e0/(P*P)+1.e0/(R*R));
                    SQQS=sqrt(1.e0/(Q*Q)+1.e0/(S*S));
                    EPR=exp(X1*SQPR);
                    EQS=exp(X2*SQQS);

                    for (int N=1; N<=2; N++) {//DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
                        //                                AND N=2 IS FOR THE SECOND ONE

                        for (int NN=1; NN<=2; NN++) {//DO 5 NN=1,2 !   NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
                            //                                         TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE

                            if (M == 1) {//IF (M.EQ.1) THEN
                                FX=-SQPR*EPR*CYPI*SZRK;
                                FY=EPR*SYPI*SZRK/P;
                                FZ=-EPR*CYPI*CZRK/R;
                                if (N == 1) {//IF (N.EQ.1) THEN
                                    if (NN == 1) {//IF (NN.EQ.1) THEN
                                        HX=FX;
                                        HY=FY;
                                        HZ=FZ;
                                    } else {//ELSE
                                        HX=FX*X_SC;
                                        HY=FY*X_SC;
                                        HZ=FZ*X_SC;
                                    }//ENDIF
                                } else {//ELSE
                                    if (NN == 1) {//IF (NN.EQ.1) THEN
                                        HX=FX*CPS;
                                        HY=FY*CPS;
                                        HZ=FZ*CPS;
                                    } else {//ELSE
                                        HX=FX*CPS*X_SC;
                                        HY=FY*CPS*X_SC;
                                        HZ=FZ*CPS*X_SC;
                                    }//ENDIF
                                }//ENDIF

                            } else {//ELSE                            !   M.EQ.2
                                FX=-SPS*SQQS*EQS*CYQI*CZSK;
                                FY=SPS/Q*EQS*SYQI*CZSK;
                                FZ=SPS/S*EQS*CYQI*SZSK;
                                if (N == 1) {//IF (N.EQ.1) THEN
                                    if (NN == 1) {//IF (NN.EQ.1) THEN
                                        HX=FX;
                                        HY=FY;
                                        HZ=FZ;
                                    } else {//ELSE
                                        HX=FX*X_SC;
                                        HY=FY*X_SC;
                                        HZ=FZ*X_SC;
                                    }//ENDIF
                                } else {//ELSE
                                    if (NN == 1) {//IF (NN.EQ.1) THEN
                                        HX=FX*S3PS;
                                        HY=FY*S3PS;
                                        HZ=FZ*S3PS;
                                    } else {//ELSE
                                        HX=FX*S3PS*X_SC;
                                        HY=FY*S3PS*X_SC;
                                        HZ=FZ*S3PS*X_SC;
                                    }//ENDIF
                                }//ENDIF
                            }//ENDIF
                            L=L+1;

                            if (M == 1) {//IF (M.EQ.1) THEN
                                HXR=HX*CT1+HZ*ST1;
                                HZR=-HX*ST1+HZ*CT1;
                            } else {//ELSE
                                HXR=HX*CT2+HZ*ST2;
                                HZR=-HX*ST2+HZ*CT2;
                            }//ENDIF

                            GX=GX+HXR*A[L];
                            GY=GY+HY *A[L];
                            GZ=GZ+HZR*A[L]; // 5
                        }

                    }//4   CONTINUE
                }//3   CONTINUE
            }//2   CONTINUE
        }//1   CONTINUE

        *BX=GX;
        *BY=GY;
        *BZ=GZ;

        return;
    }//END

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    static void BIRK_TOT (double const PS,double const X,double const Y,double const Z,
                          double *BX11,double *BY11,double *BZ11,double *BX12,double *BY12,double *BZ12,
                          double *BX21,double *BY21,double *BZ21,double *BX22,double *BY22,double *BZ22,
                          struct TS_EXTERNAL_COMMON_ *TS_EXTERNAL_COMMON) {

        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION SH11(86),SH12(86),SH21(86),SH22(86)
        double const* XKAPPA1=&TS_EXTERNAL_COMMON->BIRKPAR_XKAPPA1,*XKAPPA2=&TS_EXTERNAL_COMMON->BIRKPAR_XKAPPA2;//COMMON /BIRKPAR/ XKAPPA1,XKAPPA2   !  INPUT PARAMETERS, SPECIFIED FROM A MAIN PROGRAM
        double /*&DPHI=ctx.DPHI_B_RHO0_DPHI,&B=ctx.DPHI_B_RHO0_B,&RHO_0=ctx.DPHI_B_RHO0_RHO_0,*/
        *XKAPPA=&TS_EXTERNAL_COMMON->DPHI_B_RHO0_XKAPPA;//COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! PARAMETERS, CONTROL DAY-NIGHT ASYMMETRY OF F.A.C.

        static double const SH11[87] = {NAN/*dummy*/,46488.84663,-15541.95244,-23210.09824,-32625.03856,
            -109894.4551,-71415.32808,58168.94612,55564.87578,-22890.60626,
            -6056.763968,5091.368100,239.7001538,-13899.49253,4648.016991,
            6971.310672,9699.351891,32633.34599,21028.48811,-17395.96190,
            -16461.11037,7447.621471,2528.844345,-1934.094784,-588.3108359,
            -32588.88216,10894.11453,16238.25044,22925.60557,77251.11274,
            50375.97787,-40763.78048,-39088.60660,15546.53559,3559.617561,
            -3187.730438,309.1487975,88.22153914,-243.0721938,-63.63543051,
            191.1109142,69.94451996,-187.9539415,-49.89923833,104.0902848,
            -120.2459738,253.5572433,89.25456949,-205.6516252,-44.93654156,
            124.7026309,32.53005523,-98.85321751,-36.51904756,98.88241690,
            24.88493459,-55.04058524,61.14493565,-128.4224895,-45.35023460,
            105.0548704,-43.66748755,119.3284161,31.38442798,-92.87946767,
            -33.52716686,89.98992001,25.87341323,-48.86305045,59.69362881,
            -126.5353789,-44.39474251,101.5196856,59.41537992,41.18892281,
            80.86101200,3.066809418,7.893523804,30.56212082,10.36861082,
            8.222335945,19.97575641,2.050148531,4.992657093,2.300564232,
            .2256245602,-.05841594319};

        static double const SH12[87] = {NAN/*dummy*/,210260.4816,-1443587.401,-1468919.281,281939.2993,
            -1131124.839,729331.7943,2573541.307,304616.7457,468887.5847,
            181554.7517,-1300722.650,-257012.8601,645888.8041,-2048126.412,
            -2529093.041,571093.7972,-2115508.353,1122035.951,4489168.802,
            75234.22743,823905.6909,147926.6121,-2276322.876,-155528.5992,
            -858076.2979,3474422.388,3986279.931,-834613.9747,3250625.781,
            -1818680.377,-7040468.986,-414359.6073,-1295117.666,-346320.6487,
            3565527.409,430091.9496,-.1565573462,7.377619826,.4115646037,
            -6.146078880,3.808028815,-.5232034932,1.454841807,-12.32274869,
            -4.466974237,-2.941184626,-.6172620658,12.64613490,1.494922012,
            -21.35489898,-1.652256960,16.81799898,-1.404079922,-24.09369677,
            -10.99900839,45.94237820,2.248579894,31.91234041,7.575026816,
            -45.80833339,-1.507664976,14.60016998,1.348516288,-11.05980247,
            -5.402866968,31.69094514,12.28261196,-37.55354174,4.155626879,
            -33.70159657,-8.437907434,36.22672602,145.0262164,70.73187036,
            85.51110098,21.47490989,24.34554406,31.34405345,4.655207476,
            5.747889264,7.802304187,1.844169801,4.867254550,2.941393119,
            .1379899178,.06607020029};

        static double const SH21[87] = {NAN/*dummy*/,162294.6224,503885.1125,-27057.67122,-531450.1339,
            84747.05678,-237142.1712,84133.61490,259530.0402,69196.05160,
            -189093.5264,-19278.55134,195724.5034,-263082.6367,-818899.6923,
            43061.10073,863506.6932,-139707.9428,389984.8850,-135167.5555,
            -426286.9206,-109504.0387,295258.3531,30415.07087,-305502.9405,
            100785.3400,315010.9567,-15999.50673,-332052.2548,54964.34639,
            -152808.3750,51024.67566,166720.0603,40389.67945,-106257.7272,
            -11126.14442,109876.2047,2.978695024,558.6019011,2.685592939,
            -338.0004730,-81.99724090,-444.1102659,89.44617716,212.0849592,
            -32.58562625,-982.7336105,-35.10860935,567.8931751,-1.917212423,
            -260.2023543,-1.023821735,157.5533477,23.00200055,232.0603673,
            -36.79100036,-111.9110936,18.05429984,447.0481000,15.10187415,
            -258.7297813,-1.032340149,-298.6402478,-1.676201415,180.5856487,
            64.52313024,209.0160857,-53.85574010,-98.52164290,14.35891214,
            536.7666279,20.09318806,-309.7349530,58.54144539,67.45226850,
            97.92374406,4.752449760,10.46824379,32.91856110,12.05124381,
            9.962933904,15.91258637,1.804233877,6.578149088,2.515223491,
            .1930034238,-.02261109942};

        static double const SH22[87] = {NAN/*dummy*/,-131287.8986,-631927.6885,-318797.4173,616785.8782,
            -50027.36189,863099.9833,47680.20240,-1053367.944,-501120.3811,
            -174400.9476,222328.6873,333551.7374,-389338.7841,-1995527.467,
            -982971.3024,1960434.268,297239.7137,2676525.168,-147113.4775,
            -3358059.979,-2106979.191,-462827.1322,1017607.960,1039018.475,
            520266.9296,2627427.473,1301981.763,-2577171.706,-238071.9956,
            -3539781.111,94628.16420,4411304.724,2598205.733,637504.9351,
            -1234794.298,-1372562.403,-2.646186796,-31.10055575,2.295799273,
            19.20203279,30.01931202,-302.1028550,-14.78310655,162.1561899,
            .4943938056,176.8089129,-.2444921680,-100.6148929,9.172262228,
            137.4303440,-8.451613443,-84.20684224,-167.3354083,1321.830393,
            76.89928813,-705.7586223,18.28186732,-770.1665162,-9.084224422,
            436.3368157,-6.374255638,-107.2730177,6.080451222,65.53843753,
            143.2872994,-1028.009017,-64.22739330,547.8536586,-20.58928632,
            597.3893669,10.17964133,-337.7800252,159.3532209,76.34445954,
            84.74398828,12.76722651,27.63870691,32.69873634,5.145153451,
            6.310949163,6.996159733,1.971629939,4.436299219,2.904964304,
            .1486276863,.06859991529};

        /* ====   LEAST SQUARES FITTING ONLY:
         C       BX11=0.D0
         C       BY11=0.D0
         C       BZ11=0.D0
         C       BX12=0.D0
         C       BY12=0.D0
         C       BZ12=0.D0
         C       BX21=0.D0
         C       BY21=0.D0
         C       BZ21=0.D0
         C       BX22=0.D0
         C       BY22=0.D0
         C       BZ22=0.D0
         C===================================*/

        *XKAPPA=*XKAPPA1;//        !  FORWARDED IN BIRK_1N2
        double X_SC=*XKAPPA1-1.1e0;//    !  FORWARDED IN BIRK_SHL
        double FX11,FY11,FZ11,HX11,HY11,HZ11;
        BIRK_1N2 (1,1,PS,X,Y,Z,&FX11,&FY11,&FZ11,TS_EXTERNAL_COMMON);//           !  REGION 1, MODE 1
        BIRK_SHL (SH11,PS,X_SC,X,Y,Z,&HX11,&HY11,&HZ11);
        *BX11=FX11+HX11;
        *BY11=FY11+HY11;
        *BZ11=FZ11+HZ11;
        double FX12,FY12,FZ12,HX12,HY12,HZ12;
        BIRK_1N2 (1,2,PS,X,Y,Z,&FX12,&FY12,&FZ12,TS_EXTERNAL_COMMON);//           !  REGION 1, MODE 2
        BIRK_SHL (SH12,PS,X_SC,X,Y,Z,&HX12,&HY12,&HZ12);
        *BX12=FX12+HX12;
        *BY12=FY12+HY12;
        *BZ12=FZ12+HZ12;

        *XKAPPA=*XKAPPA2;//        !  FORWARDED IN BIRK_1N2
        X_SC=*XKAPPA2-1.0e0;//    !  FORWARDED IN BIRK_SHL
        double FX21,FY21,FZ21,HX21,HY21,HZ21;
        BIRK_1N2 (2,1,PS,X,Y,Z,&FX21,&FY21,&FZ21,TS_EXTERNAL_COMMON);//           !  REGION 2, MODE 1
        BIRK_SHL (SH21,PS,X_SC,X,Y,Z,&HX21,&HY21,&HZ21);
        *BX21=FX21+HX21;
        *BY21=FY21+HY21;
        *BZ21=FZ21+HZ21;
        double FX22,FY22,FZ22,HX22,HY22,HZ22;
        BIRK_1N2 (2,2,PS,X,Y,Z,&FX22,&FY22,&FZ22,TS_EXTERNAL_COMMON);//           !  REGION 2, MODE 2
        BIRK_SHL (SH22,PS,X_SC,X,Y,Z,&HX22,&HY22,&HZ22);
        *BX22=FX22+HX22;
        *BY22=FY22+HY22;
        *BZ22=FZ22+HZ22;

        return;
    }//END

    /*
     c ===========================================================================
     C
     SUBROUTINE  SHTBNORM_E (K,L,X,Y,Z,FX,FY,FZ)
     C
     C
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     */
    static void SHTBNORM_E (int const K,int const L,
                            double const X,double const Y,double const Z,
                            double *FX,double *FY,double *FZ,
                            TS07Context const *ctx) {
        //IMPLICIT  REAL * 8  (A - H, O - Z)

        double AK[6];//DIMENSION AK(5)

        double const (*TSE)[6][5]=ctx->TSE_TSE;//COMMON /TSE/ TSE(80,5,4)


        AK[1]=TSE[76][K][L];
        AK[2]=TSE[77][K][L];
        AK[3]=TSE[78][K][L];
        AK[4]=TSE[79][K][L];
        AK[5]=TSE[80][K][L];

        // -------------------------------------------

        double phi=atan2(Y,X);

        // -------------------------------------------

        int L1=0;
        *FX=0.e0;
        *FY=0.e0;
        *FZ=0.e0;

        // ------------------------------------
        double CMP,SMP,RHO,AKN,AKNR,CHZ,SHZ,AKNRI,RHOI,AJM,AJM1,AJMD;
        for (int m=0; m<15; m++) {//DO 2 m1=1,15
            //int m=m1-1;

            CMP=cos(m*phi);
            SMP=sin(m*phi);

            for (int n=1; n<=5; n++) {//DO 2 n=1,5
                // ------------------------------------
                RHO=sqrt(X*X+Y*Y);
                AKN=fabs(AK[n]);
                AKNR=AKN*RHO;
                // -----------------
                CHZ=cosh(Z*AKN);
                SHZ=sinh(Z*AKN);
                // ------------------------------------
                if(AKNR < 1.e-8) {//if(AKNR.lt.1.D-8) then
                    AKNRI=1.e8;
                } else {//else
                    AKNRI=1.e0/AKNR;
                }//end if
                // -----------------
                if(RHO < 1.e-8) {//if(RHO.lt.1.D-8) then
                    RHOI=1.e8;
                } else {//else
                    RHOI=1.e0/RHO;
                }//end if
                // -----------------
                if(m > 2) {//if(m.gt.2) then
                    AJM=bessj(m,AKNR);
                    AJM1=bessj(m-1,AKNR);
                    AJMD=AJM1-m*AJM*AKNRI;
                } else {//else
                    //                  ---------------------
                    if(m == 2) {//if(m.eq.2) then
                        AJM=bessj(2,AKNR);
                        AJM1=bessj1(AKNR);
                        AJMD=AJM1-m*AJM*AKNRI;
                    } else {//else
                        //                  ----------
                        if(m == 1) {//if(m.eq.1) then
                            AJM=bessj1(AKNR);
                            AJM1=bessj0(AKNR);
                            AJMD=AJM1-AJM*AKNRI;
                            //                  -----
                        } else {//else
                            AJM=bessj0(AKNR);
                            AJMD=-bessj1(AKNR);
                            //                  ----------
                        }//end if
                        //                  --------------------
                    }//end if
                    //                  --------------------
                }//endif
                // -----------------
                double DPDX=-Y*RHOI*RHOI;
                double DPDY=X*RHOI*RHOI;
                // ------------------------------------

                double HX1=-m*DPDX*CMP*SHZ*AJM;
                double HX2=-AKN*X*RHOI*SMP*SHZ*AJMD;

                double HX=HX1+HX2;

                double HY1=-m*DPDY*CMP*SHZ*AJM;
                double HY2=-AKN*Y*RHOI*SMP*SHZ*AJMD;

                double HY=HY1+HY2;

                double HZ=-AKN*SMP*CHZ*AJM;


                L1=L1+1;

                *FX=*FX+HX*TSE[L1][K][L];
                *FY=*FY+HY*TSE[L1][K][L];
                *FZ=*FZ+HZ*TSE[L1][K][L];

            }
        }//2       CONTINUE

        return;
    }//END

    /* ===========================================================================
     C
     SUBROUTINE  SHTBNORM_O (K,L,X,Y,Z,FX,FY,FZ)
     C
     C
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     */
    static void SHTBNORM_O (int const K,int const L,
                            double const X,double const Y,double const Z,
                            double *FX,double *FY,double *FZ,
                            TS07Context const *ctx) {
        //IMPLICIT  REAL * 8  (A - H, O - Z)

        double AK[6];//DIMENSION AK(5)

        double const (*TSO)[6][5]=ctx->TSO_TSO;//COMMON /TSO/ TSO(80,5,4)


        AK[1]=TSO[76][K][L];
        AK[2]=TSO[77][K][L];
        AK[3]=TSO[78][K][L];
        AK[4]=TSO[79][K][L];
        AK[5]=TSO[80][K][L];

        // -------------------------------------------

        double phi=atan2(Y,X);

        // -------------------------------------------

        int L1=0;
        *FX=0.e0;
        *FY=0.e0;
        *FZ=0.e0;

        double CMP,SMP,RHO,AKN,AKNR,CHZ,SHZ,AKNRI,RHOI,AJM,AJM1,AJMD;
        for (int m=0; m<15; m++) {//DO 2 m1=1,15
            //int m=m1-1;

            CMP=cos(m*phi);
            SMP=sin(m*phi);

            for (int n=1; n<=5; n++) {//DO 2 n=1,5
                // ------------------------------------
                RHO=sqrt(X*X+Y*Y);
                AKN=fabs(AK[n]);
                AKNR=AKN*RHO;
                // -----------------
                CHZ=cosh(Z*AKN);
                SHZ=sinh(Z*AKN);
                // ------------------------------------
                if(AKNR < 1.e-8) {//if(AKNR.lt.1.D-8) then
                    AKNRI=1.e8;
                } else {//else
                    AKNRI=1.e0/AKNR;
                }//end if
                // -----------------
                if(RHO < 1.e-8) {//if(RHO.lt.1.D-8) then
                    RHOI=1.e8;
                } else {//else
                    RHOI=1.e0/RHO;
                }//end if
                // -----------------
                if(m > 2) {//if(m.gt.2) then
                    AJM=bessj(m,AKNR);
                    AJM1=bessj(m-1,AKNR);
                    AJMD=AJM1-m*AJM*AKNRI;
                } else {//else
                    //                  ---------------------
                    if(m == 2) {//if(m.eq.2) then
                        AJM=bessj(2,AKNR);
                        AJM1=bessj1(AKNR);
                        AJMD=AJM1-m*AJM*AKNRI;
                    } else {//else
                        //                  ----------
                        if(m == 1) {//if(m.eq.1) then
                            AJM=bessj1(AKNR);
                            AJM1=bessj0(AKNR);
                            AJMD=AJM1-AJM*AKNRI;
                            //                  -----
                        } else {//else
                            AJM=bessj0(AKNR);
                            AJMD=-bessj1(AKNR);
                            //                  ----------
                        }//end if
                        //                  --------------------
                    }//end if
                    //                  --------------------
                }//endif
                // -----------------
                double DPDX=-Y*RHOI*RHOI;
                double DPDY=X*RHOI*RHOI;
                // ------------------------------------

                double HX1=m*DPDX*SMP*SHZ*AJM;
                double HX2=-AKN*X*RHOI*CMP*SHZ*AJMD;

                double HX=HX1+HX2;

                double HY1=m*DPDY*SMP*SHZ*AJM;
                double HY2=-AKN*Y*RHOI*CMP*SHZ*AJMD;

                double HY=HY1+HY2;

                double HZ=-AKN*CMP*CHZ*AJM;
                // ------------------------------------

                L1=L1+1;

                *FX=*FX+HX*TSO[L1][K][L];
                *FY=*FY+HY*TSO[L1][K][L];
                *FZ=*FZ+HZ*TSO[L1][K][L];

            }
        }//2       CONTINUE

        return;
    }//END

    // ===================================================================
    //     SUBROUTINE TAILSHT_OE (IEVO,MK,M,X,Y,Z,BX,BY,BZ)

    static void TAILSHT_OE (int const IEVO,int const MK,int const M,
                            double const X,double const Y,double const Z,
                            double *BX,double *BY,double *BZ,
                            struct TS_EXTERNAL_COMMON_ const* TS_EXTERNAL_COMMON) {

        //IMPLICIT REAL*8 (A-H,O-Z)

        double const* D0=&TS_EXTERNAL_COMMON->TAIL_D;//COMMON /TAIL/ D0  ! THE COMMON BLOCKS FORWARDS TAIL SHEET THICKNESS

        //-----------------------------------------------------------------------------------

        double RNOT=20.0;//     !    Rho_0 - scale parameter along the tail axis
        double DLTK=1.0;//      !    step in Km

        // -----------------------------------------------------------------------------------

        double RHO=sqrt(X*X+Y*Y);

        double CSPHI=X/RHO;
        double SNPHI=Y/RHO;

        double phi=atan2(Y,X);
        double CSMPHI=cos(M*phi);
        double SNMPHI=sin(M*phi);

        double DKM=1.e0+(MK-1)*DLTK;
        double RKM=DKM/RNOT;

        //double RKMZ=RKM*Z;
        double RKMR=RKM*RHO;

        double ZD=sqrt(Z*Z+*D0**D0);

        double REX=exp(RKM*ZD);

        // ---- calculating Jm and its derivatives ------
        double AJM,AJM1,AJMD;
        if(M > 2) {//if(m.gt.2) then
            AJM=bessj(M,RKMR);
            AJM1=bessj(M-1,RKMR);
            AJMD=AJM1-M*AJM/RKMR;
        } else {//else
            //                  --------------------
            if(M == 2) {//if(m.eq.2) then
                AJM=bessj(2,RKMR);
                AJM1=bessj1(RKMR);
                AJMD=AJM1-M*AJM/RKMR;
            } else {//else
                //                  --------------------
                AJM=bessj1(RKMR);
                AJM1=bessj0(RKMR);
                AJMD=AJM1-AJM/RKMR;
                //                  --------------------
            }//end if
            //                  --------------------
        }//endif
        // -----------------------------------------
        double BRO,BPHI;
        if(IEVO == 0) {//if(ievo.eq.0) then
            /* -----------------------------------------
             c calculating symmetric modes
             c -----------------------------------------
             */
            BRO=-M*SNMPHI*Z*AJMD/ZD/REX;
            BPHI=-M*M*CSMPHI*Z*AJM/RKMR/ZD/REX;
            *BZ=M*SNMPHI*AJM/REX;

            // -----------------------------------------
        } else {//else
            /* -----------------------------------------
             c calculating asymmetric modes
             c -----------------------------------------
             */
            BRO=M*CSMPHI*Z*AJMD/ZD/REX;
            BPHI=-M*M*SNMPHI*Z*AJM/RKMR/ZD/REX;
            *BZ=-M*CSMPHI*AJM/REX;

            // -----------------------------------------
        }//end if
        /*
         c --- transformation from cylindrical ccordinates to GSM ---
         */
        *BX=BRO*CSPHI-BPHI*SNPHI;
        *BY=BRO*SNPHI+BPHI*CSPHI;
        /*
         C    CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED
         */
        return;
    }//END

    //==========================================================================
    //     SUBROUTINE  SHTBNORM_S(K,X,Y,Z,FX,FY,FZ)
    static void SHTBNORM_S(int const K,
                           double const X,double const Y,double const Z,
                           double *FX,double *FY,double *FZ,
                           TS07Context const *ctx) {

        //IMPLICIT  REAL * 8  (A - H, O - Z)

        double AK[6];//DIMENSION AK(5)

        double const (*TSS)[6]=ctx->TSS_TSS;//COMMON /TSS/ TSS(80,5)


        AK[1]=TSS[76][K];
        AK[2]=TSS[77][K];
        AK[3]=TSS[78][K];
        AK[4]=TSS[79][K];
        AK[5]=TSS[80][K];

        // -------------------------------------------

        double phi=atan2(Y,X);

        // -------------------------------------------

        int L=0;

        *FX=0.e0;
        *FY=0.e0;
        *FZ=0.e0;

        double CMP,SMP,RHO,AKN,AKNR,CHZ,SHZ,AKNRI,RHOI,AJM,AJM1,AJMD;
        for (int m=0; m<15; m++) {//DO 2 m1=1,15
            //int m=m1-1;

            CMP=cos(m*phi);
            SMP=sin(m*phi);

            for (int n=1; n<=5; n++) {//DO 2 n=1,5
                // ------------------------------------
                RHO=sqrt(X*X+Y*Y);
                AKN=fabs(AK[n]);
                AKNR=AKN*RHO;
                // -----------------
                CHZ=cosh(Z*AKN);
                SHZ=sinh(Z*AKN);
                // ------------------------------------
                if(AKNR < 1.e-8) {//if(AKNR.lt.1.D-8) then
                    AKNRI=1.e8;
                } else {//else
                    AKNRI=1.e0/AKNR;
                }//end if
                // -----------------
                if(RHO < 1.e-8) {//if(RHO.lt.1.D-8) then
                    RHOI=1.e8;
                } else {//else
                    RHOI=1.e0/RHO;
                }//end if
                // -----------------
                if(m > 2) {//if(m.gt.2) then
                    AJM=bessj(m,AKNR);
                    AJM1=bessj(m-1,AKNR);
                    AJMD=AJM1-m*AJM*AKNRI;
                } else {//else
                    //                  ---------------------
                    if(m == 2) {//if(m.eq.2) then
                        AJM=bessj(2,AKNR);
                        AJM1=bessj1(AKNR);
                        AJMD=AJM1-m*AJM*AKNRI;
                    } else {//else
                        //                  ----------
                        if(m == 1) {//if(m.eq.1) then
                            AJM=bessj1(AKNR);
                            AJM1=bessj0(AKNR);
                            AJMD=AJM1-AJM*AKNRI;
                            //                  -----
                        } else {//else
                            AJM=bessj0(AKNR);
                            AJMD=-bessj1(AKNR);
                            //                  ----------
                        }//end if
                        //                  --------------------
                    }//end if
                    //                  --------------------
                }//endif
                // -----------------
                double DPDX=-Y*RHOI*RHOI;
                double DPDY=X*RHOI*RHOI;
                // ------------------------------------

                double HX1=m*DPDX*SMP*SHZ*AJM;
                double HX2=-AKN*X*RHOI*CMP*SHZ*AJMD;

                double HX=HX1+HX2;

                double HY1=m*DPDY*SMP*SHZ*AJM;
                double HY2=-AKN*Y*RHOI*CMP*SHZ*AJMD;

                double HY=HY1+HY2;

                double HZ=-AKN*CMP*CHZ*AJM;
                // ------------------------------------

                L=L+1;

                *FX=*FX+HX*TSS[L][K];
                *FY=*FY+HY*TSS[L][K];
                *FZ=*FZ+HZ*TSS[L][K];

            }
        }//2       CONTINUE

        return;
    }//END

    /*
     C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     C
     c
     C
     SUBROUTINE TAILSHT_S (M,X,Y,Z,BX,BY,BZ)
     */
    static void TAILSHT_S (int const M,
                           double const X,double const Y,double const Z,
                           double *BX,double *BY,double *BZ,
                           struct TS_EXTERNAL_COMMON_ const* TS_EXTERNAL_COMMON) {
        //IMPLICIT REAL*8 (A-H,O-Z)

        double const* D=&TS_EXTERNAL_COMMON->TAIL_D;//COMMON /TAIL/ D  ! THE COMMON BLOCKS FORWARDS TAIL SHEET THICKNESS

        //-----------------------------------------------------------------------------------

        double RNOT=20.0;//        !    This can be replaced by introducing them
        double DLTK=1.0;//         !    through the above common block


        double RHO=sqrt(X*X+Y*Y);
        double CSPHI=X/RHO;
        double SNPHI=Y/RHO;

        double DKM=1.e0+(M-1)*DLTK;
        double RKM=DKM/RNOT;

        double RKMZ=RKM*Z;
        double RKMR=RKM*RHO;

        double ZD=sqrt(Z*Z+*D**D);

        double RJ0=bessj0(RKMR);
        double RJ1=bessj1(RKMR);
        double REX=exp(RKM*ZD);

        *BX=RKMZ*RJ1*CSPHI/ZD/REX;
        *BY=RKMZ*RJ1*SNPHI/ZD/REX;
        *BZ=RKM*RJ0/REX;
        /*
         C    CALCULATION OF THE MAGNETOTAIL CURRENT CONTRIBUTION IS FINISHED
         */
        return;
    }//END

    /*
     C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     C
     SUBROUTINE UNWARPED (X,Y,Z,BXS,BYS,BZS,BXO,BYO,BZO,BXE,BYE,BZE)

     C    CALCULATES GSM COMPONENTS OF THE SHIELDED FIELD OF 45 TAIL MODES WITH UNIT
     C    AMPLITUDES,  WITHOUT ANY WARPING OR BENDING.  NONLINEAR PARAMETERS OF THE MODES
     C    ARE FORWARDED HERE VIA A COMMON BLOCK /TAIL/.
     */
    static void UNWARPED (double const X,double const Y,double const Z,
                          double BXS[],double BYS[],double BZS[],
                          double (*BXO)[5],double (*BYO)[5],double (*BZO)[5],
                          double (*BXE)[5],double (*BYE)[5],double (*BZE)[5],
                          struct TS_EXTERNAL_COMMON_ const* TS_EXTERNAL_COMMON, TS07Context const *ctx) {
        //IMPLICIT REAL*8 (A-H,O-Z)

        //double const (*TSS)[6]=ctx.TSS_TSS;//COMMON /TSS/ TSS(80,5)
        //double const (*TSO)[6][5]=ctx.TSO_TSO;//COMMON /TSO/ TSO(80,5,4)
        //double const (*TSE)[6][5]=ctx.TSE_TSE;//COMMON /TSE/ TSE(80,5,4)

        //double const& D0=ctx.TAIL_D;//COMMON /TAIL/ D0

        //DIMENSION BXS(5),BXO(5,4),BXE(5,4)
        //DIMENSION BYS(5),BYO(5,4),BYE(5,4)
        //DIMENSION BZS(5),BZO(5,4),BZE(5,4)

        // --- New tail structure -------------
        double BXSK,BYSK,BZSK,HXSK,HYSK,HZSK;
        for (int K=1; K<=5; K++) {//DO 11 K=1,5

            TAILSHT_S (K,X,Y,Z,&BXSK,&BYSK,&BZSK,TS_EXTERNAL_COMMON);
            SHTBNORM_S (K,X,Y,Z,&HXSK,&HYSK,&HZSK,ctx);

            BXS[K]=BXSK+HXSK;
            BYS[K]=BYSK+HYSK;
            BZS[K]=BZSK+HZSK;

        }//11     CONTINUE

        double BXOKL,BYOKL,BZOKL,HXOKL,HYOKL,HZOKL,BXEKL,BYEKL,BZEKL,HXEKL,HYEKL,HZEKL;
        for (int K=1; K<=5; K++) {//DO 12 K=1,5
            for (int L=1; L<=4; L++) {//DO 13 L=1,4

                TAILSHT_OE (1,K,L,X,Y,Z,&BXOKL,&BYOKL,&BZOKL,TS_EXTERNAL_COMMON);
                SHTBNORM_O (  K,L,X,Y,Z,&HXOKL,&HYOKL,&HZOKL,ctx);

                BXO[K][L]=BXOKL+HXOKL;
                BYO[K][L]=BYOKL+HYOKL;
                BZO[K][L]=BZOKL+HZOKL;

                TAILSHT_OE (0,K,L,X,Y,Z,&BXEKL,&BYEKL,&BZEKL,TS_EXTERNAL_COMMON);
                SHTBNORM_E (  K,L,X,Y,Z,&HXEKL,&HYEKL,&HZEKL,ctx);

                BXE[K][L]=BXEKL+HXEKL;
                BYE[K][L]=BYEKL+HYEKL;
                BZE[K][L]=BZEKL+HZEKL;

            }//13       CONTINUE
        }//12   CONTINUE

        return;
    }//END

    /*
     C------------------------------------------------------------------
     c
     C
     SUBROUTINE WARPED (PS,X,Y,Z,
     BXS,BYS,BZS,BXO,BYO,BZO,BXE,BYE,BZE)
     C
     C   CALCULATES GSM COMPONENTS OF THE WARPED FIELD FOR TWO TAIL UNIT MODES.
     C   THE WARPING DEFORMATION IS IMPOSED ON THE UNWARPED FIELD, COMPUTED
     C   BY THE S/R "UNWARPED".  THE WARPING PARAMETERS WERE TAKEN FROM THE
     C   RESULTS OF GEOTAIL OBSERVATIONS (TSYGANENKO ET AL. [1998]).
     C   NB # 6, P.106, OCT 12, 2000.
     */
    static void WARPED (double const PS,double const X,double const Y,double const Z,
                        double BXS[],double BYS[],double BZS[],
                        double (*BXO)[5],double (*BYO)[5],double (*BZO)[5],
                        double (*BXE)[5],double (*BYE)[5],double (*BZE)[5],
                        struct TS_EXTERNAL_COMMON_ const* TS_EXTERNAL_COMMON, TS07Context const *ctx) {
        //IMPLICIT REAL*8 (A-H,O-Z)

        //DIMENSION BXS(5),BXO(5,4),BXE(5,4)
        //DIMENSION BYS(5),BYO(5,4),BYE(5,4)
        //DIMENSION BZS(5),BZO(5,4),BZE(5,4)

        double BX_ASS[6],BX_ASO[6][5],BX_ASE[6][5];//DIMENSION BX_ASS(5),BX_ASO(5,4),BX_ASE(5,4)
        double BY_ASS[6],BY_ASO[6][5],BY_ASE[6][5];//DIMENSION BY_ASS(5),BY_ASO(5,4),BY_ASE(5,4)
        double BZ_ASS[6],BZ_ASO[6][5],BZ_ASE[6][5];//DIMENSION BZ_ASS(5),BZ_ASO(5,4),BZ_ASE(5,4)

        double const* G=&TS_EXTERNAL_COMMON->G_G, *TW=&TS_EXTERNAL_COMMON->G_TW;//COMMON /G/ G,TW
        double DGDX=0.e0;
        double XL=20.e0;
        double DXLDX=0.e0;

        double SPS=sin(PS);
        double RHO2=Y*Y+Z*Z;
        double RHO=sqrt(RHO2);

        double PHI,CPHI,SPHI;
        if (Y == 0.e0 && Z == 0.e0) {//IF (Y.EQ.0.D0.AND.Z.EQ.0.D0) THEN
            PHI=0.e0;
            CPHI=1.e0;
            SPHI=0.e0;
        } else {//ELSE
            PHI=atan2(Z,Y);
            CPHI=Y/RHO;
            SPHI=Z/RHO;
        }//ENDIF

        double RR4L4=RHO/(RHO2*RHO2+XL*XL*XL*XL);

        double F=PHI+*G*RHO2*RR4L4*CPHI*SPS +*TW*(X/10.e0);
        double DFDPHI=1.e0-*G*RHO2*RR4L4*SPHI*SPS;
        double DFDRHO=*G*RR4L4*RR4L4*(3.e0*XL*XL*XL*XL-RHO2*RHO2)*CPHI*SPS;
        double DFDX=RR4L4*CPHI*SPS*(DGDX*RHO2-*G*RHO*RR4L4*4.e0*XL*XL*XL*DXLDX)+*TW/10.e0;//        !  THE LAST TERM DESCRIBES THE IMF-INDUCED TWISTING (ADDED 04/21/06)

        double CF=cos(F);
        double SF=sin(F);
        double YAS=RHO*CF;
        double ZAS=RHO*SF;

        UNWARPED (X,YAS,ZAS,
                  BX_ASS,BY_ASS,BZ_ASS,
                  BX_ASO,BY_ASO,BZ_ASO,
                  BX_ASE,BY_ASE,BZ_ASE, TS_EXTERNAL_COMMON,ctx);

        double BRHO_AS,BPHI_AS,BRHO_S,BPHI_S;
        for (int K=1; K<=5; K++) {//DO 11 K=1,5
            // ------------------------------------------- Deforming symmetric modules
            BRHO_AS =  BY_ASS[K]*CF+BZ_ASS[K]*SF;
            BPHI_AS = -BY_ASS[K]*SF+BZ_ASS[K]*CF;

            BRHO_S = BRHO_AS*DFDPHI;
            BPHI_S = BPHI_AS-RHO*(BX_ASS[K]*DFDX+BRHO_AS*DFDRHO);

            BXS[K]=BX_ASS[K]*DFDPHI;
            BYS[K]=BRHO_S*CPHI-BPHI_S*SPHI;
            BZS[K]=BRHO_S*SPHI+BPHI_S*CPHI;

        }//11     CONTINUE

        for (int K=1; K<=5; K++) {//DO 12 K=1,5

            for (int L=1; L<=4; L++) {//DO 13 L=1,4
                // -------------------------------------------- Deforming odd modules
                BRHO_AS =  BY_ASO[K][L]*CF+BZ_ASO[K][L]*SF;
                BPHI_AS = -BY_ASO[K][L]*SF+BZ_ASO[K][L]*CF;

                BRHO_S = BRHO_AS*DFDPHI;
                BPHI_S = BPHI_AS-RHO*(BX_ASO[K][L]*DFDX+BRHO_AS*DFDRHO);

                BXO[K][L]=BX_ASO[K][L]*DFDPHI;
                BYO[K][L]=BRHO_S*CPHI-BPHI_S*SPHI;
                BZO[K][L]=BRHO_S*SPHI+BPHI_S*CPHI;
                // ------------------------------------------- Deforming even modules
                BRHO_AS =  BY_ASE[K][L]*CF+BZ_ASE[K][L]*SF;
                BPHI_AS = -BY_ASE[K][L]*SF+BZ_ASE[K][L]*CF;

                BRHO_S = BRHO_AS*DFDPHI;
                BPHI_S = BPHI_AS-RHO*(BX_ASE[K][L]*DFDX+BRHO_AS*DFDRHO);

                BXE[K][L]=BX_ASE[K][L]*DFDPHI;
                BYE[K][L]=BRHO_S*CPHI-BPHI_S*SPHI;
                BZE[K][L]=BRHO_S*SPHI+BPHI_S*CPHI;

            }//13       CONTINUE

        }//12   CONTINUE

        return;
    }//END

    /*
     c############################################################################
     c
     C
     SUBROUTINE DEFORMED (PS,X,Y,Z,
     BXS,BYS,BZS,BXO,BYO,BZO,BXE,BYE,BZE)
     C
     C    CALCULATES GSM COMPONENTS OF 104 UNIT-AMPLITUDE TAIL FIELD MODES,
     C    TAKING INTO ACCOUNT BOTH EFFECTS OF DIPOLE TILT:
     C    WARPING IN Y-Z (DONE BY THE S/R WARPED) AND BENDING IN X-Z (DONE BY THIS SUBROUTINE)
     */
    static void DEFORMED (double const PS,double const X,double const Y,double const Z,
                          double BXS[],double BYS[],double BZS[],
                          double (*BXO)[5],double (*BYO)[5],double (*BZO)[5],
                          double (*BXE)[5],double (*BYE)[5],double (*BZE)[5],
                          struct TS_EXTERNAL_COMMON_ const* TS_EXTERNAL_COMMON, TS07Context const *ctx) {
        //IMPLICIT REAL*8 (A-H,O-Z)

        //DIMENSION BXS(5),BXO(5,4),BXE(5,4)
        //DIMENSION BYS(5),BYO(5,4),BYE(5,4)
        //DIMENSION BZS(5),BZO(5,4),BZE(5,4)

        double BXASS[6],BXASO[6][5],BXASE[6][5];//DIMENSION BXASS(5),BXASO(5,4),BXASE(5,4)
        double BYASS[6],BYASO[6][5],BYASE[6][5];//DIMENSION BYASS(5),BYASO(5,4),BYASE(5,4)
        double BZASS[6],BZASO[6][5],BZASE[6][5];//DIMENSION BZASS(5),BZASO(5,4),BZASE(5,4)

        double const* RH0=&TS_EXTERNAL_COMMON->RH0_RH0;//COMMON /RH0/ RH0

        static double const RH2=-5.2e0;
        static int const IEPS=3;//DATA RH2,IEPS /-5.2D0,3/
        /*
         C  RH0,RH1,RH2, AND IEPS CONTROL THE TILT-RELATED DEFORMATION OF THE TAIL FIELD
         */
        double SPS=sin(PS);
        //double CPS=sqrt(1.e0-SPS*SPS);
        double R2=X*X+Y*Y+Z*Z;
        double R=sqrt(R2);
        double ZR=Z/R;
        double RH=*RH0+RH2*ZR*ZR;
        double DRHDR=-ZR/R*2.e0*RH2*ZR;
        double DRHDZ= 2.e0*RH2*ZR/R;

        double RRH=R/RH;
        double F=1.e0/pow((1.e0+pow(RRH, IEPS)), (1.e0/IEPS));
        double DFDR=-pow(RRH, (IEPS-1))*pow(F, (IEPS+1))/RH;
        double DFDRH=-RRH*DFDR;

        double SPSAS=SPS*F;
        double CPSAS=sqrt(1.e0-SPSAS*SPSAS);

        double XAS=X*CPSAS-Z*SPSAS;
        double ZAS=X*SPSAS+Z*CPSAS;

        double FACPS=SPS/CPSAS*(DFDR+DFDRH*DRHDR)/R;
        double PSASX=FACPS*X;
        double PSASY=FACPS*Y;
        double PSASZ=FACPS*Z+SPS/CPSAS*DFDRH*DRHDZ;

        double DXASDX=CPSAS-ZAS*PSASX;
        double DXASDY=-ZAS*PSASY;
        double DXASDZ=-SPSAS-ZAS*PSASZ;
        double DZASDX=SPSAS+XAS*PSASX;
        double DZASDY=XAS*PSASY;
        double DZASDZ=CPSAS+XAS*PSASZ;
        double FAC1=DXASDZ*DZASDY-DXASDY*DZASDZ;
        double FAC2=DXASDX*DZASDZ-DXASDZ*DZASDX;
        double FAC3=DZASDX*DXASDY-DXASDX*DZASDY;
        /*
         C     DEFORM:
         */
        WARPED(PS,XAS,Y,ZAS,
               BXASS,BYASS,BZASS,BXASO,BYASO,BZASO,BXASE,BYASE,BZASE, TS_EXTERNAL_COMMON,ctx);

        // --- New tail structure -------------

        for (int K=1; K<=5; K++) {//DO 11 K=1,5

            BXS[K]=BXASS[K]*DZASDZ-BZASS[K]*DXASDZ+BYASS[K]*FAC1;
            BYS[K]=BYASS[K]*FAC2;
            BZS[K]=BZASS[K]*DXASDX-BXASS[K]*DZASDX+BYASS[K]*FAC3;

        }//11     CONTINUE

        for (int K=1; K<=5; K++) {//DO 12 K=1,5

            for (int L=1; L<=4; L++) {//DO 13 L=1,4

                BXO[K][L]=BXASO[K][L]*DZASDZ-BZASO[K][L]*DXASDZ
                +BYASO[K][L]*FAC1;
                BYO[K][L]=BYASO[K][L]*FAC2;
                BZO[K][L]=BZASO[K][L]*DXASDX-BXASO[K][L]*DZASDX
                +BYASO[K][L]*FAC3;

                BXE[K][L]=BXASE[K][L]*DZASDZ-BZASE[K][L]*DXASDZ
                +BYASE[K][L]*FAC1;
                BYE[K][L]=BYASE[K][L]*FAC2;
                BZE[K][L]=BZASE[K][L]*DXASDX-BXASE[K][L]*DZASDX
                +BYASE[K][L]*FAC3;

            }//13       CONTINUE
        }//12   CONTINUE
        // ------------------------------------
        return;
    }//END

    /*
     C XXXXXXXXXXXXXXXXXXXXXXXXXXX11/15/05 16:06 XXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXXX
     C
     C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     C
     SUBROUTINE  SHLCAR3X3(X,Y,Z,PS,BX,BY,BZ)
     C
     C   THIS S/R RETURNS THE SHIELDING FIELD FOR THE EARTH'S DIPOLE,
     C   REPRESENTED BY  2x3x3=18 "CARTESIAN" HARMONICS, tilted with respect
     C   to the z=0 plane (see NB#4, p.74-74)
     C
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C  The 36 coefficients enter in pairs in the amplitudes of the "cartesian"
     c    harmonics (A(1)-A(36).
     c  The 14 nonlinear parameters (A(37)-A(50) are the scales Pi,Ri,Qi,and Si
     C   entering the arguments of exponents, sines, and cosines in each of the
     C   18 "Cartesian" harmonics  PLUS TWO TILT ANGLES FOR THE CARTESIAN HARMONICS
     C       (ONE FOR THE PSI=0 MODE AND ANOTHER FOR THE PSI=90 MODE)
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     */
    static void SHLCAR3X3(double const X,double const Y,double const Z,double const PS,
                          double *BX,double *BY,double *BZ) {
        //IMPLICIT  REAL * 8  (A - H, O - Z)

        //DIMENSION A(50)
        static double const A[51] = {NAN/*dummy*/,-901.2327248,895.8011176,817.6208321,-845.5880889,
            -83.73539535,86.58542841,336.8781402,-329.3619944,-311.2947120,
            308.6011161,31.94469304,-31.30824526,125.8739681,-372.3384278,
            -235.4720434,286.7594095,21.86305585,-27.42344605,-150.4874688,
            2.669338538,1.395023949,-.5540427503,-56.85224007,3.681827033,
            -43.48705106,5.103131905,1.073551279,-.6673083508,12.21404266,
            4.177465543,5.799964188,-.3977802319,-1.044652977,.5703560010,
            3.536082962,-3.222069852,9.620648151,6.082014949,27.75216226,
            12.44199571,5.122226936,6.982039615,20.12149582,6.150973118,
            4.663639687,15.73319647,2.303504968,5.840511214,.8385953499e-01,
            .3477844929};

        double const* P1=&A[37];
        double const* P2=&A[38];
        double const* P3=&A[39];
        double const* R1=&A[40];
        double const* R2=&A[41];
        double const* R3=&A[42];
        double const* Q1=&A[43];
        double const* Q2=&A[44];
        double const* Q3=&A[45];
        double const* S1=&A[46];
        double const* S2=&A[47];
        double const* S3=&A[48];

        double const* T1=&A[49];
        double const* T2=&A[50];

        double CPS=cos(PS);
        double SPS=sin(PS);
        double S2PS=2.e0*CPS;//      !   MODIFIED HERE (INSTEAD OF SIN(3*PS) I TRY SIN(2*PS)


        double ST1=sin(PS**T1);
        double CT1=cos(PS**T1);
        double ST2=sin(PS**T2);
        double CT2=cos(PS**T2);

        //     print *,X,Z

        double X1=X*CT1-Z*ST1;

        //         print *,'X1=',X1

        double Z1=X*ST1+Z*CT1;
        double X2=X*CT2-Z*ST2;
        double Z2=X*ST2+Z*CT2;
        /*
         C
         c  MAKE THE TERMS IN THE 1ST SUM ("PERPENDICULAR" SYMMETRY):
         C
         C       I=1
         */
        double SQPR= sqrt(1.e0/(*P1**P1)+1.e0/(*R1**R1));
        double CYP = cos(Y/ *P1);
        double SYP = sin(Y/ *P1);
        double CZR = cos(Z1/ *R1);
        double SZR = sin(Z1/ *R1);
        //       print *,X1
        double EXPR= exp(SQPR*X1);
        double FX1 =-SQPR*EXPR*CYP*SZR;
        double HY1 = EXPR/ *P1*SYP*SZR;
        double FZ1 =-EXPR*CYP/ *R1*CZR;
        double HX1 = FX1*CT1+FZ1*ST1;
        double HZ1 =-FX1*ST1+FZ1*CT1;

        SQPR= sqrt(1.e0/(*P1**P1)+1.e0/(*R2**R2));
        CYP = cos(Y/ *P1);
        SYP = sin(Y/ *P1);
        CZR = cos(Z1/ *R2);
        SZR = sin(Z1/ *R2);
        EXPR= exp(SQPR*X1);
        double FX2 =-SQPR*EXPR*CYP*SZR;
        double HY2 = EXPR/ *P1*SYP*SZR;
        double FZ2 =-EXPR*CYP/ *R2*CZR;
        double HX2 = FX2*CT1+FZ2*ST1;
        double HZ2 =-FX2*ST1+FZ2*CT1;

        SQPR= sqrt(1.e0/(*P1**P1)+1.e0/(*R3**R3));
        CYP = cos(Y/ *P1);
        SYP = sin(Y/ *P1);
        CZR = cos(Z1/ *R3);
        SZR = sin(Z1/ *R3);
        EXPR= exp(SQPR*X1);
        double FX3 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/ *R3*(X1+1.e0/SQPR));
        double HY3 = EXPR/ *P1*SYP*(Z1*CZR+X1/ *R3*SZR/SQPR);
        double FZ3 =-EXPR*CYP*(CZR*(1.e0+X1/(*R3**R3)/SQPR)-Z1/ *R3*SZR);
        double HX3 = FX3*CT1+FZ3*ST1;
        double HZ3 =-FX3*ST1+FZ3*CT1;
        /*
         C       I=2:
         */
        SQPR= sqrt(1.e0/(*P2**P2)+1.e0/(*R1**R1));
        CYP = cos(Y/ *P2);
        SYP = sin(Y/ *P2);
        CZR = cos(Z1/ *R1);
        SZR = sin(Z1/ *R1);
        EXPR= exp(SQPR*X1);
        double FX4 =-SQPR*EXPR*CYP*SZR;
        double HY4 = EXPR/ *P2*SYP*SZR;
        double FZ4 =-EXPR*CYP/ *R1*CZR;
        double HX4 = FX4*CT1+FZ4*ST1;
        double HZ4 =-FX4*ST1+FZ4*CT1;

        SQPR= sqrt(1.e0/(*P2**P2)+1.e0/(*R2**R2));
        CYP = cos(Y/ *P2);
        SYP = sin(Y/ *P2);
        CZR = cos(Z1/ *R2);
        SZR = sin(Z1/ *R2);
        EXPR= exp(SQPR*X1);
        double FX5 =-SQPR*EXPR*CYP*SZR;
        double HY5 = EXPR/ *P2*SYP*SZR;
        double FZ5 =-EXPR*CYP/ *R2*CZR;
        double HX5 = FX5*CT1+FZ5*ST1;
        double HZ5 =-FX5*ST1+FZ5*CT1;

        SQPR= sqrt(1.e0/(*P2**P2)+1.e0/(*R3**R3));
        CYP = cos(Y/ *P2);
        SYP = sin(Y/ *P2);
        CZR = cos(Z1/ *R3);
        SZR = sin(Z1/ *R3);
        EXPR= exp(SQPR*X1);
        double FX6 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/ *R3*(X1+1.e0/SQPR));
        double HY6 = EXPR/ *P2*SYP*(Z1*CZR+X1/ *R3*SZR/SQPR);
        double FZ6 =-EXPR*CYP*(CZR*(1.e0+X1/(*R3**R3)/SQPR)-Z1/ *R3*SZR);
        double HX6 = FX6*CT1+FZ6*ST1;
        double HZ6 =-FX6*ST1+FZ6*CT1;
        /*
         C       I=3:
         */
        SQPR= sqrt(1.e0/(*P3**P3)+1.e0/(*R1**R1));
        CYP = cos(Y/ *P3);
        SYP = sin(Y/ *P3);
        CZR = cos(Z1/ *R1);
        SZR = sin(Z1/ *R1);
        EXPR= exp(SQPR*X1);
        double FX7 =-SQPR*EXPR*CYP*SZR;
        double HY7 = EXPR/ *P3*SYP*SZR;
        double FZ7 =-EXPR*CYP/ *R1*CZR;
        double HX7 = FX7*CT1+FZ7*ST1;
        double HZ7 =-FX7*ST1+FZ7*CT1;

        SQPR= sqrt(1.e0/(*P3**P3)+1.e0/(*R2**R2));
        CYP = cos(Y/ *P3);
        SYP = sin(Y/ *P3);
        CZR = cos(Z1/ *R2);
        SZR = sin(Z1/ *R2);
        EXPR= exp(SQPR*X1);
        double FX8 =-SQPR*EXPR*CYP*SZR;
        double HY8 = EXPR/ *P3*SYP*SZR;
        double FZ8 =-EXPR*CYP/ *R2*CZR;
        double HX8 = FX8*CT1+FZ8*ST1;
        double HZ8 =-FX8*ST1+FZ8*CT1;

        SQPR= sqrt(1.e0/(*P3**P3)+1.e0/(*R3**R3));
        CYP = cos(Y/ *P3);
        SYP = sin(Y/ *P3);
        CZR = cos(Z1/ *R3);
        SZR = sin(Z1/ *R3);
        EXPR= exp(SQPR*X1);
        double FX9 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/ *R3*(X1+1.e0/SQPR));
        double HY9 = EXPR/ *P3*SYP*(Z1*CZR+X1/ *R3*SZR/SQPR);
        double FZ9 =-EXPR*CYP*(CZR*(1.e0+X1/(*R3**R3)/SQPR)-Z1/ *R3*SZR);
        double HX9 = FX9*CT1+FZ9*ST1;
        double HZ9 =-FX9*ST1+FZ9*CT1;


        double A1=A[1]+A[2]*CPS;
        double A2=A[3]+A[4]*CPS;
        double A3=A[5]+A[6]*CPS;
        double A4=A[7]+A[8]*CPS;
        double A5=A[9]+A[10]*CPS;
        double A6=A[11]+A[12]*CPS;
        double A7=A[13]+A[14]*CPS;
        double A8=A[15]+A[16]*CPS;
        double A9=A[17]+A[18]*CPS;
        *BX=A1*HX1+A2*HX2+A3*HX3+A4*HX4+A5*HX5+A6*HX6+A7*HX7+A8*HX8+A9*HX9;
        *BY=A1*HY1+A2*HY2+A3*HY3+A4*HY4+A5*HY5+A6*HY6+A7*HY7+A8*HY8+A9*HY9;
        *BZ=A1*HZ1+A2*HZ2+A3*HZ3+A4*HZ4+A5*HZ5+A6*HZ6+A7*HZ7+A8*HZ8+A9*HZ9;

        /*
         c  MAKE THE TERMS IN THE 2ND SUM ("PARALLEL" SYMMETRY):
         C
         C       I=1
         */
        double SQQS= sqrt(1.e0/(*Q1**Q1)+1.e0/(*S1**S1));
        double CYQ = cos(Y/ *Q1);
        double SYQ = sin(Y/ *Q1);
        double CZS = cos(Z2/ *S1);
        double SZS = sin(Z2/ *S1);
        double EXQS= exp(SQQS*X2);
        FX1 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY1 = EXQS/ *Q1*SYQ*CZS   *SPS;
        FZ1 = EXQS*CYQ/ *S1*SZS   *SPS;
        HX1 = FX1*CT2+FZ1*ST2;
        HZ1 =-FX1*ST2+FZ1*CT2;

        SQQS= sqrt(1.e0/(*Q1**Q1)+1.e0/(*S2**S2));
        CYQ = cos(Y/ *Q1);
        SYQ = sin(Y/ *Q1);
        CZS = cos(Z2/ *S2);
        SZS = sin(Z2/ *S2);
        EXQS= exp(SQQS*X2);
        FX2 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY2 = EXQS/ *Q1*SYQ*CZS   *SPS;
        FZ2 = EXQS*CYQ/ *S2*SZS   *SPS;
        HX2 = FX2*CT2+FZ2*ST2;
        HZ2 =-FX2*ST2+FZ2*CT2;

        SQQS= sqrt(1.e0/(*Q1**Q1)+1.e0/(*S3**S3));
        CYQ = cos(Y/ *Q1);
        SYQ = sin(Y/ *Q1);
        CZS = cos(Z2/ *S3);
        SZS = sin(Z2/ *S3);
        EXQS= exp(SQQS*X2);
        FX3 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY3 = EXQS/ *Q1*SYQ*CZS   *SPS;
        FZ3 = EXQS*CYQ/ *S3*SZS   *SPS;
        HX3 = FX3*CT2+FZ3*ST2;
        HZ3 =-FX3*ST2+FZ3*CT2;
        /*
         C       I=2
         */
        SQQS= sqrt(1.e0/(*Q2**Q2)+1.e0/(*S1**S1));
        CYQ = cos(Y/ *Q2);
        SYQ = sin(Y/ *Q2);
        CZS = cos(Z2/ *S1);
        SZS = sin(Z2/ *S1);
        EXQS= exp(SQQS*X2);
        FX4 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY4 = EXQS/ *Q2*SYQ*CZS   *SPS;
        FZ4 = EXQS*CYQ/ *S1*SZS   *SPS;
        HX4 = FX4*CT2+FZ4*ST2;
        HZ4 =-FX4*ST2+FZ4*CT2;

        SQQS= sqrt(1.e0/(*Q2**Q2)+1.e0/(*S2**S2));
        CYQ = cos(Y/ *Q2);
        SYQ = sin(Y/ *Q2);
        CZS = cos(Z2/ *S2);
        SZS = sin(Z2/ *S2);
        EXQS= exp(SQQS*X2);
        FX5 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY5 = EXQS/ *Q2*SYQ*CZS   *SPS;
        FZ5 = EXQS*CYQ/ *S2*SZS   *SPS;
        HX5 = FX5*CT2+FZ5*ST2;
        HZ5 =-FX5*ST2+FZ5*CT2;

        SQQS= sqrt(1.e0/(*Q2**Q2)+1.e0/(*S3**S3));
        CYQ = cos(Y/ *Q2);
        SYQ = sin(Y/ *Q2);
        CZS = cos(Z2/ *S3);
        SZS = sin(Z2/ *S3);
        EXQS= exp(SQQS*X2);
        FX6 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY6 = EXQS/ *Q2*SYQ*CZS   *SPS;
        FZ6 = EXQS*CYQ/ *S3*SZS   *SPS;
        HX6 = FX6*CT2+FZ6*ST2;
        HZ6 =-FX6*ST2+FZ6*CT2;
        /*
         C       I=3
         */
        SQQS= sqrt(1.e0/(*Q3**Q3)+1.e0/(*S1**S1));
        CYQ = cos(Y/ *Q3);
        SYQ = sin(Y/ *Q3);
        CZS = cos(Z2/ *S1);
        SZS = sin(Z2/ *S1);
        EXQS= exp(SQQS*X2);
        FX7 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY7 = EXQS/ *Q3*SYQ*CZS   *SPS;
        FZ7 = EXQS*CYQ/ *S1*SZS   *SPS;
        HX7 = FX7*CT2+FZ7*ST2;
        HZ7 =-FX7*ST2+FZ7*CT2;

        SQQS= sqrt(1.e0/(*Q3**Q3)+1.e0/(*S2**S2));
        CYQ = cos(Y/ *Q3);
        SYQ = sin(Y/ *Q3);
        CZS = cos(Z2/ *S2);
        SZS = sin(Z2/ *S2);
        EXQS= exp(SQQS*X2);
        FX8 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY8 = EXQS/ *Q3*SYQ*CZS   *SPS;
        FZ8 = EXQS*CYQ/ *S2*SZS   *SPS;
        HX8 = FX8*CT2+FZ8*ST2;
        HZ8 =-FX8*ST2+FZ8*CT2;

        SQQS= sqrt(1.e0/(*Q3**Q3)+1.e0/(*S3**S3));
        CYQ = cos(Y/ *Q3);
        SYQ = sin(Y/ *Q3);
        CZS = cos(Z2/ *S3);
        SZS = sin(Z2/ *S3);
        EXQS= exp(SQQS*X2);
        FX9 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY9 = EXQS/ *Q3*SYQ*CZS   *SPS;
        FZ9 = EXQS*CYQ/ *S3*SZS   *SPS;
        HX9 = FX9*CT2+FZ9*ST2;
        HZ9 =-FX9*ST2+FZ9*CT2;

        A1=A[19]+A[20]*S2PS;
        A2=A[21]+A[22]*S2PS;
        A3=A[23]+A[24]*S2PS;
        A4=A[25]+A[26]*S2PS;
        A5=A[27]+A[28]*S2PS;
        A6=A[29]+A[30]*S2PS;
        A7=A[31]+A[32]*S2PS;
        A8=A[33]+A[34]*S2PS;
        A9=A[35]+A[36]*S2PS;

        *BX=*BX+A1*HX1+A2*HX2+A3*HX3+A4*HX4+A5*HX5+A6*HX6+A7*HX7+A8*HX8
        +A9*HX9;
        *BY=*BY+A1*HY1+A2*HY2+A3*HY3+A4*HY4+A5*HY5+A6*HY6+A7*HY7+A8*HY8
        +A9*HY9;
        *BZ=*BZ+A1*HZ1+A2*HZ2+A3*HZ3+A4*HZ4+A5*HZ5+A6*HZ6+A7*HZ7+A8*HZ8
        +A9*HZ9;

        return;
    }//END

    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     c
     SUBROUTINE EXTERN (IOPGEN,A,NTOT,
     *  PS,PDYN,X,Y,Z,
     *  BXCF,BYCF,BZCF,
     *  BXTS,BYTS,BZTS,BXTO,BYTO,BZTO,BXTE,BYTE,BZTE,
     *  BXR11,BYR11,BZR11,BXR12,BYR12,BZR12,
     *  BXR21a,BYR21a,BZR21a,BXR21s,BYR21s,BZR21s,
     *  BX,BY,BZ)
     C
     C   IOPGEN - GENERAL OPTION FLAG:  IOPGEN=0 - CALCULATE TOTAL FIELD
     C                                  IOPGEN=1 - DIPOLE SHIELDING ONLY
     C                                  IOPGEN=2 - TAIL FIELD ONLY
     C                                  IOPGEN=3 - BIRKELAND FIELD ONLY
     C                                  IOPGEN=4 - RING CURRENT FIELD ONLY
     */
    static void EXTERN (int const IOPGEN,double const A[],int const NTOT,
                        double const PS,double const PDYN,
                        double const X,double const Y,double const Z,
                        double *BXCF,double *BYCF,double *BZCF,
                        double BXTS[],double BYTS[],double BZTS[],
                        double (*BXTO)[5],double (*BYTO)[5],double (*BZTO)[5],double (*BXTE)[5],double (*BYTE)[5],double (*BZTE)[5],
                        double *BXR11,double *BYR11,double *BZR11,double *BXR12,double *BYR12,double *BZR12,
                        double *BXR21a,double *BYR21a,double *BZR21a,double *BXR21s,double *BYR21s,double *BZR21s,
                        double *BX,double *BY,double *BZ,
                        TS07Context const *ctx) {
        //IMPLICIT  REAL * 8  (A - H, O - Z)

        //DIMENSION A(NTOT)

        //DIMENSION BXTS(5),BXTO(5,4),BXTE(5,4)
        //DIMENSION BYTS(5),BYTO(5,4),BYTE(5,4)
        //DIMENSION BZTS(5),BZTO(5,4),BZTE(5,4)

        struct TS_EXTERNAL_COMMON_ TS_EXTERNAL_COMMON;

        double *D=&TS_EXTERNAL_COMMON.TAIL_D; //COMMON /TAIL/ D  ! THE COMMON BLOCK FORWARDS TAIL SHEET THICKNESS
        double *XKAPPA1=&TS_EXTERNAL_COMMON.BIRKPAR_XKAPPA1,*XKAPPA2=&TS_EXTERNAL_COMMON.BIRKPAR_XKAPPA2; //COMMON /BIRKPAR/ XKAPPA1,XKAPPA2   !  SCALING FACTORS FOR BIRKELAND CURRENTS
        double *G=&TS_EXTERNAL_COMMON.G_G,*TW=&TS_EXTERNAL_COMMON.G_TW; //COMMON /G/ G,TW
        double *RH0=&TS_EXTERNAL_COMMON.RH0_RH0; //COMMON /RH0/ RH0

        //double const A0_A=34.586e0,A0_S0=1.1960e0,A0_X0=3.4397e0; // DATA A0_A,A0_S0,A0_X0 /34.586D0,1.1960D0,3.4397D0/   !   SHUE ET AL. PARAMETERS
        //double const DSIG=0.005e0,RH2=-5.2e0; // DATA DSIG /0.005D0/, RH2 /-5.2D0/

        double XAPPA=pow((PDYN/2.e0), 0.155); //   !   0.155 is the value obtained in TS05
        double XAPPA3=XAPPA*XAPPA*XAPPA;

        *D=      A[96];
        *RH0=    A[97];
        *G=      A[98];
        *XKAPPA1=A[99];
        *XKAPPA2=A[100];
        *TW=     A[101];//       !   THIS PARAMETER CONTROLS THE IMF-INDUCED TWISTING (ADDED 04/21/06)

        double XX=X*XAPPA;//  ! pressure scaling has been reinstated here
        double YY=Y*XAPPA;
        double ZZ=Z*XAPPA;

        //     print *,XAPPA,PDYN

        //double SPS=sin(PS);

        //double X0=A0_X0/XAPPA;//   ! pressure scaling has been reinstated, even though these parameters are not used in this code
        //double AM=A0_A/XAPPA;//    ! pressure scaling has been reinstated, even though these parameters are not used in this code
        //double S0=A0_S0;//
        /*
         C   CALCULATE THE IMF CLOCK ANGLE:
         C
         C        IF (BYIMF.EQ.0.D0.AND.BZIMF.EQ.0.D0) THEN
         C            THETA=0.D0
         C         ELSE
         C            THETA=DATAN2(BYIMF,BZIMF)
         C            IF (THETA.LE.0.D0) THETA=THETA+6.283185307D0
         C        ENDIF
         C
         C       CT=COS(THETA)
         C       ST=SIN(THETA)
         C       YS=Y*CT-Z*ST
         C       ZS=Z*CT+Y*ST
         C
         C       STHETAH=SIN(THETA/2.)**2
         C
         C  CALCULATE "IMF" COMPONENTS OUTSIDE THE MAGNETOPAUSE LAYER (HENCE BEGIN WITH "O")
         C  THEY ARE NEEDED ONLY IF THE POINT (X,Y,Z) IS WITHIN THE TRANSITION MAGNETOPAUSE LAYER
         C  OR OUTSIDE THE MAGNETOSPHERE:
         C
         C      FACTIMF=A(24)+A(25)*STHETAH
         C
         C      OIMFX=0.D0
         C      OIMFY=BYIMF*FACTIMF
         C      OIMFZ=BZIMF*FACTIMF
         c
         C =====================================================================
         C  THIS FRAGMENT (BETWEEN THE ===== LINES) DISABLES THE CALCULATION OF THE MAGNETOPAUSE POSITION
         C  IT SHOULD BE USED ONLY FOR THE FITTING (WE ASSUME THAT NO POINTS FROM THE SHEATH ARE PRESENT
         C  IN THE DATASET, WHICH ITSELF IS STILL A QUESTION).
         C
         C  REMOVE IT IN THE FINAL VERSION.
         C
         C      SIGMA=0.D0
         C      GOTO 1111
         C======================================================================
         c
         C      R=SQRT(X**2+Y**2+Z**2)
         C      XSS=X
         C      ZSS=Z
         C
         C  1   XSOLD=XSS      !   BEGIN ITERATIVE SEARCH OF UNWARPED COORDS (TO FIND SIGMA)
         C      ZSOLD=ZSS
         C
         C      RH=RH0+RH2*(ZSS/R)**2
         C      SINPSAS=SPS/(1.D0+(R/RH)**3)**0.33333333D0
         C      COSPSAS=DSQRT(1.D0-SINPSAS**2)
         C      ZSS=X*SINPSAS+Z*COSPSAS
         C      XSS=X*COSPSAS-Z*SINPSAS
         C      DD=DABS(XSS-XSOLD)+DABS(ZSS-ZSOLD)
         C      IF (DD.GT.1.D-6) GOTO 1
         C                                END OF ITERATIVE SEARCH
         C      RHO2=Y**2+ZSS**2
         C      ASQ=AM**2
         C      XMXM=AM+XSS-X0
         C      IF (XMXM.LT.0.) XMXM=0. ! THE BOUNDARY IS A CYLINDER TAILWARD OF X=X0-AM
         C      AXX0=XMXM**2
         C      ARO=ASQ+RHO2
         C      SIGMA=DSQRT((ARO+AXX0+SQRT((ARO+AXX0)**2-4.*ASQ*AXX0))/(2.*ASQ))
         C
         C==================================================================
         C 1111 CONTINUE  !!!!!!!!!!!!  REMOVE IN THE FINAL VERSION
         C==================================================================
         C
         C
         C   NOW, THERE ARE THREE POSSIBLE CASES:
         C    (1) INSIDE THE MAGNETOSPHERE   (SIGMA
         C    (2) IN THE BOUNDARY LAYER
         C    (3) OUTSIDE THE MAGNETOSPHERE AND B.LAYER
         C       FIRST OF ALL, CONSIDER THE CASES (1) AND (2):
         C
         C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         C      IF (SIGMA.LT.S0+DSIG) THEN  !  CASES (1) OR (2); CALCULATE THE MODEL FIELD
         C                                   (WITH THE POTENTIAL "PENETRATED" INTERCONNECTION FIELD):
         C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         */
        if (IOPGEN <= 1) { //IF (IOPGEN.LE.1) THEN

            //     print *,XX,YY,ZZ
            double CFX,CFY,CFZ;
            SHLCAR3X3(XX,YY,ZZ,PS,&CFX,&CFY,&CFZ);//         !  DIPOLE SHIELDING FIELD
            *BXCF=CFX  *XAPPA3;
            *BYCF=CFY  *XAPPA3;
            *BZCF=CFZ  *XAPPA3;
        } else { //ELSE
            *BXCF=0.e0;
            *BYCF=0.e0;
            *BZCF=0.e0;
        } //ENDIF                                              !  DONE

        if (IOPGEN == 0 || IOPGEN == 2) { //IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.2) THEN
            DEFORMED (PS,XX,YY,ZZ,                //  TAIL FIELD (THREE MODES)
                      BXTS,BYTS,BZTS,BXTO,BYTO,BZTO,BXTE,BYTE,BZTE, &TS_EXTERNAL_COMMON,ctx);
        } else {//ELSE
            for (int K=1; K<=5; K++) { //DO 11 K=1,5

                BXTS[K]=0.e0;
                BYTS[K]=0.e0;
                BZTS[K]=0.e0;

            }// 11     CONTINUE

            for (int K=1; K<=5; K++) {//DO 12 K=1,5

                for (int L=1; L<=4; L++) {//DO 13 L=1,4

                    BXTO[K][L]=0.e0;
                    BYTO[K][L]=0.e0;
                    BZTO[K][L]=0.e0;

                    BXTE[K][L]=0.e0;
                    BYTE[K][L]=0.e0;
                    BZTE[K][L]=0.e0;

                }//13       CONTINUE

            } //12   CONTINUE

        }//ENDIF

        double BXR22a,BYR22a,BZR22a,BXR11s,BYR11s,BZR11s,BXR12s,BYR12s,BZR12s,BXR22s,BYR22s,BZR22s;
        if (IOPGEN == 0 || IOPGEN == 3) {//IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.3) THEN
            BIRK_TOT (PS,XX,YY,ZZ,BXR11,BYR11,BZR11,BXR12,BYR12,
                      BZR12,BXR21a,BYR21a,BZR21a,&BXR22a,&BYR22a,&BZR22a,&TS_EXTERNAL_COMMON);//    !   BIRKELAND FIELD (TWO MODES FOR R1 AND TWO MODES FOR R2)

            BIRTOTSY (PS,XX,YY,ZZ,&BXR11s,&BYR11s,&BZR11s,&BXR12s,
                      &BYR12s,&BZR12s,BXR21s,BYR21s,BZR21s,&BXR22s,&BYR22s,&BZR22s,&TS_EXTERNAL_COMMON);//    !   "SYMMETRIC" BIRKELAND FIELD
            //       (TWO MODES FOR R1s AND TWO MODES FOR R2s)
            //       (but we actually use from here only R2s modes)
        } else {//ELSE
            *BXR11=0.e0;
            *BYR11=0.e0;
            *BZR11=0.e0;
            *BXR12=0.e0;
            *BYR12=0.e0;
            *BZR12=0.e0;
            *BXR21a=0.e0;
            *BYR21a=0.e0;
            *BZR21a=0.e0;
            *BXR21s=0.e0;
            *BYR21s=0.e0;
            *BZR21s=0.e0;
        }//ENDIF
        /*
         C-----------------------------------------------------------
         C
         C    NOW, ADD UP ALL THE COMPONENTS:
         */
        double A_R11=A[92];
        double A_R12=A[93];
        double A_R21a=A[94];
        double A_R21s=A[95];

        double TX=0.e0;
        double TY=0.e0;
        double TZ=0.e0;

        // --- New tail structure -------------


        double PDYN_0=2.e0;//   !   AVERAGE PRESSURE USED FOR NORMALIZATION

        double P_FACTOR=sqrt(PDYN/PDYN_0)-1.e0;


        int IND=1;

        for (int K=1; K<=5; K++) {//DO 911 K=1,5
            IND=IND+1;
            TX=TX+(A[IND]+A[IND+45]*P_FACTOR)*BXTS[K];//    !   2 - 6  &  47 - 51
            TY=TY+(A[IND]+A[IND+45]*P_FACTOR)*BYTS[K];
            TZ=TZ+(A[IND]+A[IND+45]*P_FACTOR)*BZTS[K];
        }//911     CONTINUE


        for (int K=1; K<=5; K++) {//DO 912 K=1,5

            for (int L=1; L<=4; L++) {//DO 913 L=1,4

                IND=IND+1;

                TX=TX+(A[IND]+A[IND+45]*P_FACTOR)*BXTO[K][L];//  !   7 -26  &  52 - 71
                TY=TY+(A[IND]+A[IND+45]*P_FACTOR)*BYTO[K][L];
                TZ=TZ+(A[IND]+A[IND+45]*P_FACTOR)*BZTO[K][L];

                TX=TX+(A[IND+20]+A[IND+65]*P_FACTOR)*BXTE[K][L];// !   27 -46  &  72 - 91
                TY=TY+(A[IND+20]+A[IND+65]*P_FACTOR)*BYTE[K][L];
                TZ=TZ+(A[IND+20]+A[IND+65]*P_FACTOR)*BZTE[K][L];

            }//913       CONTINUE

        }//912   CONTINUE

        double BBX=A[1]**BXCF+TX+
        A_R11**BXR11+A_R12**BXR12+A_R21a**BXR21a+A_R21s**BXR21s;

        double BBY=A[1]**BYCF+TY+
        A_R11**BYR11+A_R12**BYR12+A_R21a**BYR21a+A_R21s**BYR21s;

        double BBZ=A[1]**BZCF+TZ+
        A_R11**BZR11+A_R12**BZR12+A_R21a**BZR21a+A_R21s**BZR21s;
        /*
         c   -----------------------------------------------------------
         C
         C   AND WE HAVE THE TOTAL EXTERNAL FIELD.
         C
         C  NOW, LET US CHECK WHETHER WE HAVE THE CASE (1). IF YES - WE ARE DONE:
         C
         C      IF (SIGMA.LT.S0-DSIG) THEN    !  (X,Y,Z) IS INSIDE THE MAGNETOSPHERE
         C
         C       BX=BBX
         C       BY=BBY
         C       BZ=BBZ
         C                     ELSE           !  THIS IS THE MOST COMPLEX CASE: WE ARE INSIDE
         C                                             THE INTERPOLATION REGION
         C       FINT=0.5*(1.-(SIGMA-S0)/DSIG)
         C       FEXT=0.5*(1.+(SIGMA-S0)/DSIG)
         C
         C       CALL DIPOLE (PS,X,Y,Z,QX,QY,QZ)
         C       BX=(BBX+QX)*FINT+OIMFX*FEXT -QX
         C       BY=(BBY+QY)*FINT+OIMFY*FEXT -QY
         C       BZ=(BBZ+QZ)*FINT+OIMFZ*FEXT -QZ
         c
         C        ENDIF  !   THE CASES (1) AND (2) ARE EXHAUSTED; THE ONLY REMAINING
         C                      POSSIBILITY IS NOW THE CASE (3):
         C++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         C        ELSE
         C                CALL DIPOLE (PS,X,Y,Z,QX,QY,QZ)
         C                BX=OIMFX-QX
         C                BY=OIMFY-QY
         C                BZ=OIMFZ-QZ
         C        ENDIF
         */

        *BX=BBX;
        *BY=BBY;
        *BZ=BBZ;

        return;
    }//END

    /***********************************************************************
     c
     SUBROUTINE EXTMODEL (IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
     C
     C  THIS SUBROUTINE SERVES JUST AS AN INTERFACE BETWEEN THE EXTERNAL FIELD
     C  SUBROUTINE EXTERN (WHOSE OUTPUT IS ALL IN DOUBLE PRECISION) AND THE MAIN
     C  PROGRAM (WHERE WE USE ONLY SINGLE PRECISION)
     C*/
    static void TS07D (int const IOPT,double const PARMOD[]/*zero index based 10 element C array*/,double const PS,
                       double const X,double const Y,double const Z,
                       double *BX,double *BY,double *BZ,
                       TS07Context const *ctx) {
        //IMPLICIT REAL*8 (A-H,O-Z)
        //REAL PARMOD(10),PS,X,Y,Z,BX,BY,BZ  ! IN FACT, PARMOD(10) IS NOT USED ANYWHERE IN THIS REALIZATION OF THE MODEL
        //                                   !   BUT KEPT ONLY FOR CONSISTENCY WITH GEOPACK_08 S/W
        int const NTOT=101; //PARAMETER (NTOT=101)

        /* K.MIN: In original fortran code, PARMOD is not used and PDYN is passed via COMMON BLOCK. Here, PDYN is set to the first element of PARMOD.
         */
        double const* PDYN = &PARMOD[0]; //COMMON /INPUT/  PDYN
        double BXTS[6],BXTO[6][5],BXTE[6][5]; //DIMENSION BXTS(5),BXTO(5,4),BXTE(5,4)
        double BYTS[6],BYTO[6][5],BYTE[6][5]; //DIMENSION BYTS(5),BYTO(5,4),BYTE(5,4)
        double BZTS[6],BZTO[6][5],BZTE[6][5]; //DIMENSION BZTS(5),BZTO(5,4),BZTE(5,4)
        double const *A=PARMOD;//double const *A=ctx->PARAM_A; //COMMON /PARAM/ A(NTOT)

        //double const& PSS=PS;
        //double const& XX=X;
        //double const& YY=Y;
        //double const& ZZ=Z;

        double BXCF,BYCF,BZCF,BXR11,BYR11,BZR11,BXR12,BYR12,BZR12,BXR21a,BYR21a,BZR21a,BXR21s,BYR21s,BZR21s;//,BBX,BBY,BBZ;
        EXTERN (0,A,NTOT,PS,*PDYN,X,Y,Z,&BXCF,&BYCF,&BZCF,
                BXTS,BYTS,BZTS,BXTO,BYTO,BZTO,BXTE,BYTE,BZTE,
                &BXR11,&BYR11,&BZR11,&BXR12,&BYR12,&BZR12,&BXR21a,&BYR21a,&BZR21a,&BXR21s,
                &BYR21s,&BZR21s,BX,BY,BZ,ctx);

        //BX=BBX;
        //BY=BBY;
        //BZ=BBZ;

        return;

    }//END


    ////////////////////////////////////////////////////////////////
    // Hard coded shielding field parameters
    ////////////////////////////////////////////////////////////////
#pragma mark - Hard coded shielding field parameters
    namespace __ts07_private {
        // Hard coded shielding field parameters in ROW-MAJOR ORDER
        double const TSS[5][80] = {
            /*tailamebhr1.par*/{13.25539804,255.3217538,-79.89706790,195.0019560,-347.6926018,0.3973046034,174.5316027,-23.35668432,428.5107065,-400.5955764,10.91783806,231.1235658,-64.87034621,451.0944052,-407.5943909,-0.1622864939,20.66583971,-4.106874805,-244.9537855,10.82528900,6.969519442,99.14809699,-35.75424066,-77.11980961,-114.8389102,2.838045646,-112.1078206,-0.4433978304,107.4992199,399.4372386,16.10195878,158.4046678,-74.12900346,-34.55378744,-207.0295882,22.75349480,-233.5493121,-63.84998082,4.779104523,20.77461950,49.20154366,139.3056306,-324.4424523,0.1490629538,21.35258296,-58.42712912,-158.3722924,-82.73035450,-0.3341216639,-23.70307931,-186.0731254,109.0391770,287.9703869,0.1046095319,11.61799927,-42.12650667,-48.03543730,-237.4279218,-0.1993909017E-01,-3.642779050,85.34232763,14.24194314,104.4439846,0.2660826013E-02,0.7942294559,4.770045828,-2.758556969,-23.90987296,-0.2520465807E-03,-0.1196121512,-44.13661853,0.2496181491,-0.4164982984,0.1511020603E-04,0.1041401758E-01,-0.1192315041,0.8089158143E-01,0.1021836245,0.3905866076E-01,0.6230692710E-01},
            /*tailamebhr2.par*/{1413.250347,3388.986757,-3209.156273,1453.189894,-2812.519798,1607.694575,5411.997496,-4048.293303,5763.317904,-6665.784469,2576.897249,8607.738810,-6391.771568,14749.11298,-12115.94957,1582.448621,5464.278466,-3982.413713,10525.51121,-7821.318422,1862.681642,5899.773770,-4629.551410,-11556.80966,-4910.266539,355.3506268,2561.050018,-1192.865402,-125.1672575,-5142.769670,-1178.009932,-5686.161661,3245.503342,1941.225718,12038.37331,482.7613468,3441.049646,-1552.962277,-718.5278467,-6907.875657,-1187.068679,-3333.315354,3125.830980,109.4077536,1220.357614,1116.497337,4142.211324,-2026.672127,2.138247571,530.8706281,241.7507168,-2232.936800,2007.092920,-4.288653843,-374.0412307,1610.945795,597.4109334,-978.1077288,0.9559599404,109.5464428,-679.2994761,106.0523856,1005.194380,-0.1126468710,-15.46511845,-213.4507825,-170.6503933,-805.1737415,0.5126621600E-02,-0.8188571651,226.4973178,81.15966318,386.0208285,0.7942566154E-03,1.008747946,-0.1039004631,0.8415142137E-01,0.9791085633E-01,0.3975123479E-01,0.6242141261E-01},
            /*tailamebhr3.par*/{2561.932600,6639.253742,-5619.538006,2106.965705,-5276.951152,3716.245708,13359.03475,-8923.118125,10533.53613,-14293.45731,5815.101290,22381.07333,-14074.84055,32335.30659,-27692.94504,4448.980831,17216.99977,-10788.01118,31294.44794,-21793.67842,5321.615127,18941.61758,-12814.91867,-28569.99406,-17099.33616,1855.759384,10759.12212,-5136.070414,-2484.659131,-16105.89050,-2614.872817,-18758.20973,7541.820300,5017.661232,37601.32630,1836.091848,13422.80658,-5493.858857,-1513.064830,-19158.24575,-5185.209330,-14665.34093,12908.99484,166.5050017,366.2175899,4360.651257,16175.24228,-7786.512178,21.33920771,4207.829478,785.3331105,-8253.487501,6963.923738,-10.89936917,-2212.570029,5409.724335,2128.850954,-3313.321114,2.007489643,599.6034890,-2769.368086,251.2069477,3175.623710,-0.2078567171,-73.99022951,-125.5725850,-470.7474776,-2459.158838,0.7708521717E-02,-10.04047953,386.1546517,218.8488833,1140.719560,0.1454835647E-02,7.405377840,-0.1037764722,0.8289173533E-01,0.9792480989E-01,0.3781082252E-01,0.6548732570E-01},
            /*tailamebhr4.par*/{1007.772729,5888.099069,-3973.767508,1956.082474,-4544.191741,1497.560875,12793.54815,-7000.144307,10375.63086,-13758.38190,2040.430518,20724.27070,-10177.61292,33633.23739,-27853.52345,879.8366650,11755.54062,-5005.020217,36842.11984,-19740.51455,1052.862248,10185.95498,-5391.129821,-29855.92571,-7588.764721,-50.47223937,10322.53311,-2309.880310,-2485.822262,-21795.94233,-2299.144186,-21662.58363,10855.33649,5318.922966,30622.10372,734.3952846,15304.21767,-5627.390706,-1808.514559,-16316.03841,-2439.048656,-12019.74265,10785.60416,254.5303220,1705.709709,2453.184802,12800.84334,-4914.221808,13.41491883,2640.498918,2787.236046,-6820.734791,4369.534647,-12.71855819,-1581.882978,4506.548208,1932.751530,-2858.809415,2.736241932,464.8436734,-2269.821481,160.0499888,3023.269179,-0.3246745441,-67.41897360,-639.3539647,-411.8486264,-2241.210983,0.1646797120E-01,-3.884378163,761.7182150,203.3498859,1043.233799,0.1833910158E-02,4.828327164,-0.1054748550,0.8332762730E-01,0.9563917472E-01,0.3958834221E-01,0.6476937287E-01},
            /*tailamebhr5.par*/{152.0121860,3065.611859,-1546.755417,1180.965793,-2712.716526,192.5421054,7481.308423,-2990.146935,6452.196295,-9311.657731,67.95145094,10880.48087,-3192.221706,21572.74305,-19604.74342,-564.0789001,-3042.746045,3641.862842,26607.22528,-8773.289193,-845.5116889,-15301.47936,7800.243345,-14824.70914,15603.73369,-1017.444788,-2097.764676,5465.778093,3149.415804,-15491.71049,-1792.790224,-19320.59910,13988.37832,-620.7342789,-4539.431427,-1119.482987,4390.753073,2952.175430,-1327.451528,-7201.850067,-1838.923414,12658.48049,2433.480132,787.0373924,8959.813558,-2356.448334,-4948.893880,5087.281938,-191.6507180,-3241.191775,-3089.997231,-1876.311884,-10734.00506,22.89587098,375.0633295,-8051.379908,2061.080192,6576.094231,-0.9157485222E-01,112.7994803,-1136.646006,-1155.026743,-3787.707002,-0.6234105676,-79.69108605,6758.267975,473.0142648,1907.186619,0.1621568982,26.78301103,-4801.725215,-142.8725456,-713.0206101,-0.2595854277E-01,-6.204706618,-0.1098200713,0.8222202821E-01,0.9327852053E-01,0.4454373405E-01,0.6527850671E-01}
        };
        double const TSO[4][5][80] = {
            {
                /*tailamhr_o_11.par*/{-9129.863693,9641.853031,3747.603145,3037.108159,-6638.170746,-5494.853909,8824.575167,1724.062495,10244.73317,-10467.67027,-4806.552752,7595.605693,1630.259761,18024.50570,-10776.80359,6338.919743,-7339.044323,-2515.715146,-5019.240162,6344.584531,-3549.395103,3706.458035,1468.225176,-2182.240459,-2382.199685,4166.398630,-5403.930543,-1637.177714,2416.777528,6653.826247,-4042.440909,5553.918871,1327.355786,-737.6335425,-3469.342970,4424.410454,-3209.326658,-2290.696676,91.24816922,-174.7841954,-4663.062520,2131.127447,1151.404990,10.74380109,1303.366532,-272.3142702,-3294.939688,-3005.117046,-8.111389242,-1062.616558,-1783.962876,1650.492974,-5006.649181,2.006544201,439.0314168,-1137.968098,-1085.661072,-755.7385646,-0.3450012548,-143.9553017,2423.810861,661.8671906,4684.680422,0.4695092355E-01,40.52429986,-1575.621891,-280.8616573,-3570.856175,-0.5135707926E-02,-9.341890612,632.3829557,85.27455898,1595.531822,0.4409555230E-03,1.698645476,-0.1018427152,0.8925050613E-01,0.1080180596,0.3878781213E-01,0.6827678942E-01},
                /*tailamhr_o_21.par*/{-2448.190539,4577.162745,395.3375337,2445.313203,-4469.467772,-1282.617335,3795.768240,97.10212479,7056.877598,-6484.708415,-1733.990807,4147.161725,239.4483999,11015.09495,-7062.903717,1128.951250,-2464.146093,-163.2242182,-4412.311504,3462.510848,-713.7126717,790.4122576,161.3125259,94.63785933,207.0820406,309.1458490,-1551.964402,3.254226855,1017.344717,4330.276272,-1126.614532,1789.548256,236.6959573,-389.0591302,-3000.442846,-316.3980379,-2483.903914,286.6675517,63.92118498,594.0233590,-3769.165329,828.0529835,792.3829717,-1.613765568,159.6379552,-563.2181443,-2081.344708,-426.9942095,-3.037219513,-334.1890520,2884.057627,1914.181809,-1878.416513,1.048701523,190.0069897,-2448.788404,-999.3005084,-196.5991527,-0.2050759555,-64.96442722,1120.254477,339.2504994,822.9164030,0.2774209983E-01,15.25112053,-291.1824719,-75.41230470,-170.6616652,-0.2694222160E-02,-2.506005170,22.17172045,8.772391885,-201.7001725,0.1759599126E-03,0.2589954877,-0.9990472018E-01,0.8563366860E-01,0.1149467941,0.3805331746E-01,0.6395083260E-01},
                /*tailamhr_o_31.par*/{-1983.036180,3353.126844,404.0772721,1239.412765,-2762.178591,-1771.740909,3691.926018,298.4417785,3464.623744,-4201.145545,-3244.399873,5792.356057,641.8142122,4132.637212,-5499.175722,-2525.956173,4460.035334,479.8670445,-3952.767024,-3200.843936,-1538.600047,2144.437484,348.2367760,-2117.589085,-756.1996288,-318.7476572,-2095.443424,219.1939238,1968.022400,8331.229302,-1210.566337,2514.085434,211.9041245,-480.6560983,-3405.662465,2042.674032,-5145.558794,-238.5911310,17.73149363,-1128.501659,-3097.669999,4153.735415,451.8536868,18.52709019,1461.084753,1576.276317,-2451.180942,-720.0875414,-6.269141648,-722.4183444,-951.9312477,900.7978416,-378.7351041,1.161492520,221.0278985,556.3039260,-179.3093296,-258.3457407,-0.1444096774,-44.59186527,-276.0448308,-9.288649924,726.6149605,0.1176936969E-01,5.008960954,110.8252036,21.96101708,-598.3913039,-0.3820397199E-03,0.2260118395,-32.79896445,-9.262171164,332.9596671,-0.5761191320E-04,-0.2516548038,-0.1002865642,0.8628802713E-01,0.1134301290,0.3833056606E-01,0.6653723847E-01},
                /*tailamhr_o_41.par*/{1826.705558,-1479.990713,-1244.899702,-706.4008455,1482.503327,1232.845822,-1995.831935,-718.7693486,-3032.362759,3293.996874,-841.3280272,-598.4878506,690.0496521,-6489.332802,3162.922678,-4808.359356,4736.569975,3159.632303,870.2924292,-5322.394518,1011.593505,-177.3395405,-785.9187926,-1711.110454,-205.9480940,-643.7476905,2.599803310,456.1045276,569.7606964,1710.536155,-134.7760370,315.8023178,56.09920074,78.99185836,733.8674234,1062.815745,-1448.515623,-457.9606232,-96.69250090,-1491.612426,-988.7777634,324.4832354,1299.347212,25.72567699,605.1679145,2272.136389,1249.143759,-80.36970619,-2.062623937,22.94702844,-59.79745997,-1076.180905,2607.292044,-0.5151698282,-100.3560281,849.4320662,667.6361523,-439.1403591,0.2371229319,55.08652148,-1229.466100,-310.9803574,-1047.846344,-0.5255335718E-01,-19.15514987,720.7698774,105.1232092,820.1709027,0.8048647722E-02,4.782845342,-248.5622378,-25.87986884,-310.9568739,-0.9192116012E-03,-0.8895771861,-0.1054652405,0.8643674931E-01,0.1084689939,0.4188030664E-01,0.6739894547E-01},
                /*tailamhr_o_51.par*/{2039.930169,-3657.152990,-932.4261701,-1285.305972,3563.206455,1693.517871,-5117.703622,-629.7582611,-5334.200434,6919.965569,982.8231926,-3849.271322,-362.6264400,-10568.24625,6780.719768,-1497.262376,4835.021456,508.7309733,2996.249501,-6299.856272,2414.820816,-1808.198491,-1312.694840,-62.03519922,-137.6173128,1368.235641,1035.126599,-916.5284693,-978.5798425,-3502.775259,1991.104620,-56.14098132,-1348.750212,355.8002187,2446.120237,1162.756215,4499.167524,-1190.509226,-34.17805035,963.4747360,3124.088517,-1348.807904,-1646.321051,-4.083140578,-635.7966795,-2745.677996,730.3007416,1746.134370,3.281176733,414.3687093,1073.656627,-667.8283568,408.0016909,-0.9792895713,-230.0755914,363.0374670,398.7030869,277.3100534,0.1906879320,89.66933204,-543.6830042,-150.5452811,-732.7041154,-0.2656629804E-01,-24.40502674,247.0460015,37.09767782,435.9190114,0.2687897894E-02,4.631098862,-55.72347686,-5.361165460,-130.7490548,-0.1846741026E-03,-0.5552524720,-0.1058739688,0.8533246265E-01,0.1125116407,0.3959459877E-01,0.7070985883E-01}
            },{
                /*tailamhr_o_12.par*/{4112.619480,-3191.522090,-2000.985561,-734.2928434,1687.833700,-2169.583355,3001.305403,755.1572971,3908.209658,-3780.589303,12472.46371,-7663.033547,-6376.878401,14073.62974,-2106.952615,15287.55684,-15150.49641,-6577.177415,1818.028834,11029.78470,-530.0748718,340.0648518,269.3543338,7145.670976,-547.3290584,-5022.877387,6697.580247,1923.128252,-4399.712887,-10176.72215,4121.197726,-4599.009732,-1705.798124,901.1001242,3458.964874,-6187.736386,3883.527262,2846.760983,53.39023881,3107.705075,3911.306510,-5083.628309,-2223.678927,-72.79197854,-3188.580080,-4022.130486,1400.477686,282.0067338,18.41285997,1195.376222,1947.419228,776.2562595,-3389.352052,-2.449004952,-208.4816445,-1878.805665,-1260.507729,1575.016195,0.1096921177,-19.15598123,1431.894744,789.5281717,67.46993535,0.2779599703E-01,23.91701567,-624.6257282,-304.3518354,-163.1453235,-0.7790780330E-02,-8.141303468,145.1465191,78.94786376,-51.63402725,0.1133928721E-02,1.777537026,-0.1020916041,0.9053425532E-01,0.1072539447,0.4031689591E-01,0.6707148582E-01},
                /*tailamhr_o_22.par*/{850.0809975,-880.9232758,-307.7253571,-171.2890140,487.3672941,-1922.846158,4073.941404,391.4483187,6883.679591,-6357.457515,5361.073197,-4062.283203,-1974.853860,21263.49411,-5258.754619,8660.478464,-14485.42461,-2220.567015,-2298.481180,15756.77752,-1581.803615,689.3955849,650.3474616,3760.974659,2104.474583,-2144.096377,3320.269363,592.1894329,-1259.145370,-3904.617777,-448.0537745,-57.62718780,46.39230585,192.7157575,1739.601607,-1009.601862,499.2528898,-398.5021208,21.83203775,498.9478962,-1225.359027,-2138.947543,-2115.059087,-16.24468851,-692.4447420,-3172.950837,-184.5580727,-5434.819357,2.678529850,158.0012011,-664.0318771,295.9923171,-7736.865122,-0.1023964450,6.067766538,-1651.607339,-738.6780037,1008.502233,-0.8296768021E-01,-32.58066632,2873.558283,643.9392428,5663.356237,0.2805648632E-01,17.95657029,-1863.610334,-308.5946658,-4997.990717,-0.5196155472E-02,-5.741608922,748.4976117,99.52896237,2363.858460,0.6737227895E-03,1.282353108,-0.1028586142,0.8892722516E-01,0.1120675784,0.3934202720E-01,0.6566983932E-01},
                /*tailamhr_o_32.par*/{1716.388462,-1684.983754,-571.2313512,-45.25119113,587.4683054,-1536.516702,3867.941847,107.6784144,3214.908029,-4099.043790,9137.335392,-7199.661285,-3098.187340,11909.19805,-1881.391726,17462.34944,-23871.46896,-4318.293456,12561.55608,12499.50846,-4706.290470,2656.942759,1829.925668,11900.00216,1263.827007,-8857.775306,11110.31538,2382.655071,-7197.139711,-5069.763504,1581.103360,-5681.329529,11.73426364,2058.567033,8424.277903,-3159.757436,2726.272273,873.2829999,-205.2281571,45.20059060,2860.729502,-4152.754244,-1160.827668,-46.71639990,-3754.845249,-2532.276413,2652.671017,-1201.360077,21.24961153,2322.490281,809.4823185,-800.4137425,-3003.011743,-4.009092902,-711.5586378,-1667.872981,-523.4373089,447.2563516,0.3979628681,51.07400162,1712.366494,629.7124105,1711.714899,-0.1499878222E-02,52.52266623,-927.9710232,-314.0201069,-1371.352108,-0.6821670905E-02,-27.55680667,321.3032345,100.4890838,532.4541281,0.1361956076E-02,7.881648893,-0.9861122058E-01,0.8886723990E-01,0.1066397651,0.4053325320E-01,0.7340035571E-01},
                /*tailamhr_o_42.par*/{-1033.150761,2592.426618,136.2015712,2074.534537,-3529.853924,-560.5291302,2536.485421,20.00868790,5673.223847,-5995.481918,-337.4655817,1754.721247,20.57292817,7905.511676,-5741.530289,1102.268555,-3402.597563,-113.8261619,-2880.787200,5612.361221,-525.2301491,128.2459248,122.2935034,1732.566824,1806.177052,-863.1649157,1180.167879,135.1568399,822.4543083,2136.769747,-1099.121397,1582.160003,221.4072860,-493.1757825,-1395.818851,-1774.834147,-2258.439149,519.3784019,123.8553767,379.3878767,-5012.555099,-262.6222882,816.6330440,-17.83634009,-101.8101103,-845.1251395,-1402.391537,-577.3230491,-7.379652335,-173.6526761,4677.645049,1740.109916,-2066.355834,4.575917512,141.8862246,-3969.079126,-986.4401507,-499.1323639,-1.252854916,-55.69056548,1780.286223,344.6494814,370.7032179,0.2222324271,14.08741626,-444.5485206,-78.29615376,942.3174431,-0.2763928350E-01,-2.451893649,23.19990650,9.865845931,-1184.390198,0.2336934375E-02,0.2768368948,-0.1010849109,0.8402875949E-01,0.1190039198,0.4496931167E-01,0.6316419903E-01},
                /*tailamhr_o_52.par*/{-4174.249290,5891.414468,1233.386211,2690.798139,-5290.091891,-2791.382165,5180.922970,673.1626127,5790.963933,-7063.717085,-5314.138089,7648.455351,1552.394458,4016.478610,-7182.434316,-4081.658384,6618.398636,1025.399384,-7210.550941,-4454.574198,-686.8918125,2374.104703,2.925346772,-5327.385168,-3200.105813,2938.101913,-5424.773649,-899.8348234,5318.875620,11095.93108,-1119.419559,6776.820191,-551.4329015,-1518.168169,-4415.156376,7195.165929,-4044.928926,-3146.968207,91.26409115,-721.5517893,-3704.214566,6116.637135,-501.4272331,90.09838260,1946.569798,-1275.668463,-4880.059751,-2674.778292,-37.26165787,-1071.065965,-1457.023612,1484.561637,-4037.699752,7.650383355,306.6802765,438.0405147,-316.3708883,-378.7039795,-1.055747962,-59.51261238,713.5100161,114.6889283,4109.182555,0.1062234479,8.961144330,-681.6470109,-55.53486777,-3499.809100,-0.7426226205E-02,-1.131034205,333.2626864,21.37245578,1762.511873,0.2129725616E-03,0.1261987603,-0.1011903824,0.8731514653E-01,0.1106323079,0.4426842668E-01,0.6570691743E-01}
            },{
                /*tailamhr_o_13.par*/{-747.3405783,-4123.025295,2350.789854,-1520.747603,3780.326776,-915.1420786,-8046.593981,3425.218525,-7073.004867,9960.266748,-665.6546382,-7726.502029,2731.378612,-14080.98792,12001.07544,162.3554481,1482.339875,-542.1660762,9170.448901,-3284.131035,237.2089101,340.5840893,-645.5251500,592.5781019,1033.235920,53.44901901,4390.063357,-1031.270732,-2547.677914,-5775.802281,-953.1304780,-1634.260632,2005.799311,1082.247132,5579.810538,-1580.672039,4464.988850,2670.401508,-182.4325605,-774.3994929,-1747.657450,-918.6411901,5396.625143,5.218899888,-328.9424825,1583.940095,1947.795616,-2314.955599,9.417297454,676.9760929,961.0389758,-2008.446788,-1213.711830,-3.563270557,-428.9587528,-744.3860538,1054.393446,1675.842357,0.7515489713,157.7204148,675.5273859,-334.6558119,-684.1046656,-0.1077550550,-38.30993570,-750.5230367,61.41947982,30.87356713,0.1077200034E-01,6.121671808,542.4310536,-1.580812072,105.6981982,-0.6597319616E-03,-0.4785960429,-0.1125546762,0.8601082946E-01,0.1030046711,0.4181406127E-01,0.6880751127E-01},
                /*tailamhr_o_23.par*/{3149.982502,-5881.471562,-1005.921983,-1849.714655,5169.812303,4816.546627,-11985.51060,-1302.758075,-9642.185257,13668.44085,4195.675652,-12383.60105,-1032.105313,-20591.78288,16793.67508,1030.180443,-3394.349858,-151.5885288,11764.25935,2779.443012,-1780.093772,-1556.318632,873.6338910,5353.561077,6996.131315,-5194.923109,13337.99281,1272.945805,-4523.358077,-13147.71665,3207.003284,-5436.918522,-1275.552185,1425.692287,11371.21629,1443.665735,7890.484102,-1257.387925,-193.3170072,-1096.501173,8990.011804,-4189.111794,-2719.813678,-3.753942373,-1711.300073,-3599.158236,4734.844247,2729.438618,9.425977652,1822.415727,-809.9365239,-3382.471114,1725.556570,-2.618128183,-929.1039863,1701.080718,1490.273431,-579.7112608,0.4417419756,302.6573553,-755.9755058,-417.3293251,-123.3866386,-0.5183242629E-01,-66.80196739,44.45557788,63.35665508,-157.1634002,0.4169951244E-02,9.278708330,111.9199093,3.172074729,264.4968832,-0.1744176532E-03,-0.3402525020,-0.1031306538,0.8573209679E-01,0.1126062280,0.3810932879E-01,0.7022217216E-01},
                /*tailamhr_o_33.par*/{164.2882058,-738.2555780,66.74224445,-527.5043108,975.2338168,-1417.591412,553.9189640,734.5176323,-844.7579035,734.8013142,-1122.311926,-391.1753435,809.5595475,-333.6182845,1320.519690,858.9727499,-5484.686432,710.8949682,4536.493597,5445.070140,-8012.470189,-2103.995131,4921.035294,12106.85689,11622.69348,-18697.16687,10850.91231,8055.422463,-1203.391517,13197.44701,3905.151273,-3601.889994,-1876.309522,3829.496566,6130.947540,6357.903659,-2157.403954,-3034.111166,-2209.200884,-9257.817289,-1777.686283,5333.569747,-109.1450799,656.9222907,4757.612753,2498.611280,252.2937804,-996.3892468,-80.76264275,-683.8306227,-3358.578607,-2659.921971,2077.327583,-10.84357877,-355.7327004,2296.238617,1818.097700,-1479.110999,7.077872649,241.7199063,-658.0776411,-644.2183761,1438.454616,-1.666731904,-74.96094662,-109.4690019,119.4818219,-1236.454782,0.2469929523,14.38794034,192.1543978,3.077389138,738.9952692,-0.2355900811E-01,-1.530130403,-0.9760666723E-01,0.8703574212E-01,0.1037370794,0.5026492434E-01,0.6947108465E-01},
                /*tailamhr_o_43.par*/{3694.056653,6028.165497,-1257.144559,-3467.064385,-5863.546408,1103.768941,3571.878327,-305.9173743,-3783.848013,-2411.192440,6643.083851,6377.458977,-2273.805963,8966.442198,-9724.533952,6694.421076,25605.60790,-1583.354773,3337.764363,-17271.72011,-4689.714087,18712.34502,2121.447228,6037.305642,-392.5432262,-12441.11041,-14058.79382,3664.392555,-2443.984730,22427.23399,1368.146663,12371.76183,-1224.383405,685.2871630,2305.248583,5567.968212,-3690.853950,-2674.975537,-98.38907688,-2028.102246,1047.796554,276.2228081,-1277.516996,5.005712704,-69.59505857,-1918.170688,246.6437529,-709.5674767,1.107500863,1408.168198,-114.8952505,-123.8746586,-768.3685137,-0.3184791937,-1069.113141,197.6572713,30.25211428,757.1350615,0.4276908342E-01,360.8733635,237.0015625,-4.038553054,193.9058863,-0.3454891983E-02,-40.96290632,-269.6350960,0.5757236927E-01,-461.0196540,0.1225463039E-03,-17.22942988,134.6753861,0.1113690275,272.0367494,0.1022664961E-04,10.63647008,-0.9608831742E-01,0.5824998334E-01,0.1050748382,0.3189037952E-01,0.7905610921E-01},
                /*tailamhr_o_53.par*/{377.7535971,-3255.571500,-94.42445595,-2326.009022,5000.477118,98.76658052,-321.6878644,-27.28269894,527.8218150,22.26341612,944.0954899,-5201.722239,-232.2265331,8252.878279,2966.816450,672.3506920,-12722.80630,-60.29323661,8360.836372,18609.65566,-893.5678295,-2360.672368,298.5172900,9274.505912,9607.477142,-2061.692000,12581.79260,517.9009559,-5354.505657,-9363.234437,-2296.376656,6079.439744,594.4062201,1814.112355,8589.916580,-3723.750525,510.1189800,784.1973528,-190.5854668,-988.3436141,-6495.544636,-5081.876411,238.2190999,-56.01715003,-1355.627283,-906.7968273,2716.893844,-3603.555249,22.89435642,682.9092914,8111.006171,-630.7576940,-6539.622893,-3.805494863,-150.2237000,-9972.146255,-4.477313930,303.1266070,0.2792652529,9.380281047,6376.181875,54.40028634,2266.929448,0.1965972222E-01,4.588073500,-2533.191466,-21.19057483,428.4015864,-0.8790210902E-02,-1.827528071,629.7657971,4.972777519,-2082.866740,0.1405981186E-02,0.3849594418,-0.1071564790,0.7227279547E-01,0.1233168996,0.4132417435E-01,0.6047942965E-01}
            },{
                /*tailamhr_o_14.par*/{7512.495926,2486.819952,-9360.696654,563.0578977,-1152.330349,3019.571297,1009.077069,-3766.029864,261.5053184,-488.7746733,-1119.572061,-2234.670328,2081.923051,-5805.233667,4729.356603,3530.511464,-947.9815159,-3720.955602,-13410.40189,7059.328844,-4219.352404,-1273.285480,5318.719617,12542.51332,-2985.697041,4183.796063,-856.2230740,-4407.933881,-2968.659407,4681.011612,-2854.011048,268.3251337,3349.694456,-870.2632267,-7659.351793,3350.035545,-4342.068781,-2326.322973,684.7425731,4419.514188,-1692.161479,69.52622100,786.0231129,-274.0491419,-2333.763102,-997.0036891,1037.339897,-1313.836673,71.70612280,876.4182738,-755.2551057,-897.5493013,-616.1083041,-13.16804120,-232.5482420,823.2389957,292.0610757,663.6925280,1.632661584,39.69900804,-191.7827093,9.726142597,-121.5781918,-0.9959350969E-01,-2.343011288,-48.42161625,-45.45643678,-61.69962169,-0.9089478422E-02,-0.9254394826,45.04229726,20.67652580,44.78727114,0.3534551786E-02,0.3646321691,-0.9673820449E-01,0.8626706627E-01,0.9553040712E-01,0.4594795222E-01,0.6324950410E-01},
                /*tailamhr_o_24.par*/{-4853.710435,7539.283132,1383.617734,3553.580950,-7159.851728,-5685.349470,11177.14552,1381.779719,11682.98175,-15048.68600,-7098.902385,13267.49586,1851.532731,18287.93162,-18689.92051,-2628.727742,2461.506099,858.8676932,-15030.51244,3481.878368,-2755.736191,1966.559300,996.8667831,-527.7609427,2550.317799,-2890.932259,-1276.066217,1259.543729,5926.555942,15176.65908,-2928.987950,3124.936699,1107.569262,-2115.415167,-7949.323711,-226.4469010,-11330.58854,1225.897333,179.2663617,-983.0341935,-9561.028420,5582.070553,2741.940729,69.73826532,1819.231011,4686.865492,-4465.776346,-2119.970808,-39.64162900,-1187.147810,417.5173947,3103.073936,-623.9144335,11.13795716,503.0055461,-1286.814147,-1369.877011,444.0682990,-2.097774502,-145.9827222,382.4438915,365.2261539,-619.7720630,0.2786518782,29.08642493,143.6009285,-40.86552007,814.4980342,-0.2482208127E-01,-3.573753282,-172.2977532,-11.41097139,-582.9693771,0.9586287944E-03,0.6550450097E-01,-0.1023464717,0.8697081774E-01,0.1128488055,0.4418222039E-01,0.6592493733E-01},
                /*tailamhr_o_34.par*/{-6844.723862,5426.891083,3653.847969,1542.870615,-3629.681704,-4377.872568,3667.910862,2298.376834,1512.973594,-2799.229979,-11719.62599,7689.376500,6519.024962,-3466.935063,-1601.662855,-12369.87070,9006.150846,6626.385742,-7887.535681,-1888.790286,-4213.531125,-419.2622393,2801.956900,-6363.998852,8185.623077,-4397.074014,-2387.781358,2870.388849,10927.02776,18452.10519,-3209.427908,7298.937187,235.7335859,-1193.412242,1769.165152,14783.67551,-4863.046503,-9066.955435,-1236.348690,-10586.01126,-10032.60063,9562.597530,4311.024281,720.7743941,8242.326284,5404.977314,-4079.104412,-3512.916509,-190.3882130,-3005.702253,-4957.385310,-1376.148804,2307.672919,27.26955240,481.4295180,2968.117786,2042.987987,-2501.642397,-1.052073432,55.16455543,-624.1601136,-998.4836286,2915.957575,-0.4848664503,-54.23899172,-409.7450544,280.9372885,-2395.181796,0.1399277923,16.73252235,456.4290243,-39.27870263,1413.327013,-0.2213203960E-01,-3.256976958,-0.1045726782,0.9164398599E-01,0.1090811978,0.4929334097E-01,0.7059253890E-01},
                /*tailamhr_o_44.par*/{-5136.588722,-6081.507807,2292.216881,3027.583367,6564.259567,-1997.420290,-2743.151204,877.7907228,1947.113999,2677.928202,-10007.05585,-5514.527538,4604.570151,-9600.955999,10947.64273,-8747.372834,-11529.16432,3693.563050,-2748.191098,12386.74458,-8414.767446,6126.592751,3905.992467,-1467.047432,6543.940311,-10640.24637,16493.52999,3862.497263,1799.628259,13395.55402,4240.195650,2702.443821,-3987.414537,1.624110277,9914.456544,19348.45673,-5958.487456,-9601.275559,-131.6127411,-9609.436587,-4575.685143,4108.754138,2043.595141,47.42614261,14432.92531,393.0228264,-1360.491410,2942.266368,-8.536220496,-7123.058559,198.4659745,233.7573440,4432.873351,0.8749265306,1402.573732,2707.198962,2.757352311,-160.9129747,-0.2607188969E-01,460.3149218,-2728.549928,-13.90858727,-2543.615026,-0.7795262408E-02,-443.0286005,1300.326149,4.476886108,1638.592950,0.1660674813E-02,171.2708980,-372.0598768,-0.8797493459,-474.5423293,-0.1926290493E-03,-42.63059660,-0.9711151673E-01,0.6011736295E-01,0.1038710134,0.3376088362E-01,0.7983628439E-01},
                /*tailamhr_o_54.par*/{524.8063550,8260.303369,-4090.159943,4949.779315,-9230.441878,315.7427688,4477.535580,-2334.216561,2438.113403,-4646.932500,1216.954142,12993.19009,-8307.149124,-14038.34373,-1138.175705,1334.390073,22236.62892,-11157.29925,-9232.920219,-16424.99567,1263.044608,12852.54379,-9023.747360,-5667.900372,3848.945180,1014.988936,15585.64256,-10298.10372,9266.959126,22923.18785,-1117.195431,9060.142885,1732.019060,-342.1829195,1151.717611,-2244.510549,-19442.57627,17843.68608,-1178.268805,-8153.247931,817.7298737,19430.67817,867.9799778,543.8971947,5090.129433,3353.070594,-6392.138818,3944.897938,-116.4193938,-1429.792818,5356.396528,-492.6277671,-8883.223186,11.87148451,152.8530839,-456.4900539,1448.401267,7751.092758,0.6227352683,40.32568299,-1616.930849,-739.1056651,-3959.705839,-0.4743080464,-23.05829013,-74.25697497,228.3995719,1329.166112,0.9822599491E-01,5.955737887,1051.054577,-48.68011219,-293.7303408,-0.1319864421E-01,-1.037585123,-0.1112028593,0.7979730649E-01,0.9351340512E-01,0.4448288312E-01,0.6015772503E-01}
            }
        };
        double const TSE[4][5][80] = {
            {
                /*tailamhr_e_11.par*/{0.,0.,0.,0.,0.,-1635.691948,-28930.97135,8245.823066,-29925.08403,42275.94074,806.6684224,560.8077583,-2513.682257,-23549.41085,10103.35539,-76.23276554,5803.489808,-470.6579522,21658.84980,-14132.13783,14.06899727,-2243.113362,246.3587424,-4691.478961,4769.586110,-43.61281785,4982.765351,-367.8342988,-3406.973074,-11915.00792,-241.9393467,-8765.720963,2052.948947,1967.319660,7815.662746,944.6143675,2255.108158,-3128.529423,-585.8194200,-4187.529690,700.0942102,-2454.797575,-3353.939057,94.65725868,816.5296898,-724.3606196,505.8595172,-2684.782259,-8.917125328,-78.76697889,-3566.401210,198.2729446,4142.089169,0.1210201393,-1.508721438,-495.4001844,-238.1479271,-3492.356398,0.1048599843,-2.986583329,1261.504771,109.1306061,1736.070148,-0.1903775245E-01,2.390145712,2.393483005,-30.38754334,-512.8640180,0.1978701967E-02,-0.7276602613,-531.3895168,5.489599954,69.67497715,-0.1298054520E-03,0.1357890368,-0.1099674919,0.7873252974E-01,0.9762197017E-01,0.4071780300E-01,0.6343778344E-01},
                /*tailamhr_e_21.par*/{0.,0.,0.,0.,0.,471.0648177,4657.810326,-2682.872430,1934.482342,-4064.487443,292.9248911,1172.657060,-1373.572228,-6405.273827,2655.244934,718.4082435,9941.585815,-5030.852873,-4407.576868,-7740.864200,-125.4509581,549.5241683,268.2114495,-3307.277386,-1579.583299,-85.57463696,-3045.659611,1040.287938,2473.333006,4384.731947,-88.50386591,-764.2691265,657.5020888,-741.7227357,-2953.919108,203.3355925,-2686.487179,-40.33760783,-63.23114168,-954.7897018,422.6821176,2295.395291,-3440.322512,85.54260840,1193.498297,-783.4682277,-2240.540429,164.8357536,-30.64982725,-656.6147322,-1420.551966,1116.889717,566.4276397,6.608928299,216.3198419,-347.2077398,-342.9546493,-397.7159324,-0.9596426309,-46.91722842,913.9943763,52.52799308,92.24493413,0.8734504994E-01,5.953803684,-458.4126006,7.258542501,31.91255699,-0.1761306981E-02,0.4661200516E-01,126.6857465,-7.150982930,-32.25479984,-0.9849341433E-03,-0.2316652686,-0.1094765838,0.8252055713E-01,0.9538712160E-01,0.4419781741E-01,0.6415453686E-01},
                /*tailamhr_e_31.par*/{0.,0.,0.,0.,0.,3618.057554,12788.95907,-6347.761316,13277.10046,-20216.73390,-995.5349461,736.5059856,1393.837460,7892.701526,-5303.617599,624.7468943,-245.3916845,-965.9764341,-12351.01328,5189.716871,-324.8626611,3448.320787,144.6730679,-1680.008732,-7174.864126,-2103.608289,-7366.433303,3752.251737,5137.922758,8855.068755,1621.214948,4578.122803,-2700.761714,-2660.858072,-6504.927202,-318.0217215,-3744.299081,898.8715105,657.4500435,1886.368886,1985.300583,1345.445613,-3123.880668,-82.29435964,-156.0379084,-2065.744601,-623.4235346,3157.439408,-6.028437820,-194.9999960,1337.937018,328.8171793,-1776.100085,5.410842887,121.0046101,-394.3094848,-91.93484016,977.5219612,-1.312307290,-36.21460594,86.47007040,0.2111675670,-475.7275182,0.1892168108,6.142636332,-105.3983702,10.39300922,156.1564632,-0.1592581355E-01,-0.3216640187,114.3749831,-4.776368468,-18.00949142,0.1012542556E-03,-0.1445843096,-0.1080369027,0.8107196419E-01,0.1037153332,0.4674820089E-01,0.6560415087E-01},
                /*tailamhr_e_41.par*/{0.,0.,0.,0.,0.,-13783.44671,9436.537277,11954.50740,8574.185406,-13887.96148,20251.62627,223.4185485,-18821.27226,12801.70542,-7967.052160,19519.12475,-4415.200082,-17733.48937,6941.946147,1180.684120,-5284.295469,7977.130903,4229.409860,5099.938243,-13617.29398,7397.565801,6028.456646,-7229.457645,-9778.586776,-14861.94084,-1272.927462,1211.387402,1276.859015,-1778.313739,-5354.394491,-6817.274613,-4267.716007,6611.480108,3039.039229,14623.41517,6667.434161,-3193.023802,-5406.321705,-1257.910459,-10318.85727,-3841.704236,-2094.667622,4156.382573,253.1418933,2580.543648,2848.227648,5123.629199,-2285.516804,-16.23382865,370.3352639,-1680.252112,-3527.387542,1200.767220,-5.592823234,-505.2690332,-756.2850263,1411.024991,-2214.309718,2.093724489,198.8499888,1822.282014,-372.3057903,2527.865683,-0.3969072342,-48.89640678,-1414.127983,61.99821648,-1733.600215,0.5141101420E-01,8.244291367,-0.1040991681,0.8302302776E-01,0.1048893768,0.4795465931E-01,0.6940868945E-01},
                /*tailamhr_e_51.par*/{0.,0.,0.,0.,0.,-932.5379008,1687.854879,928.0484479,2044.493740,-3131.196405,-4029.967093,-3113.781590,5803.433982,3491.313276,-64.94642303,1104.374977,-974.9986165,-903.5218828,6630.650789,-851.1588278,3852.320261,-2271.458417,-4318.320557,1387.600728,4275.535661,3157.721167,-1806.281992,-3360.695160,617.0093913,2241.594698,43.28081870,-3434.009440,1260.223210,352.1218913,-1955.946766,1231.661652,-2759.871499,-382.3254070,-555.4112023,-2623.210024,467.1455431,-466.8878517,-435.6675177,95.62242530,284.2102454,-1950.315771,2576.014521,695.0952818,28.33918165,926.3882280,-570.6648413,-2072.408938,-2145.624612,-18.35167592,-672.8803527,-749.4055562,722.0969576,226.2136780,4.523179151,219.5046613,1691.984399,-78.19723630,979.5288436,-0.6267012592,-34.66611751,-1237.632490,-41.26640579,-819.7505876,0.3709737656E-01,-1.343441717,562.6521210,25.73706309,380.0171718,0.4956993846E-02,2.395797440,-0.9819136695E-01,0.8164669390E-01,0.9590489566E-01,0.4920865643E-01,0.7089769123E-01}
            },{
                /*tailamhr_e_12.par*/{0.,0.,0.,0.,0.,-153.3676502,-2291.166780,256.7201787,-4920.613115,5111.351218,39.00569341,-5955.430382,337.7489222,-25306.27620,16080.37071,-10108.10163,22417.75066,1781.631639,3099.357366,-21861.57088,5604.657205,-9140.400932,-1276.738046,-1823.300973,5989.453327,-569.1364939,2655.870640,52.92490160,-1284.065123,-5438.169614,6252.048930,-10428.27438,-1176.406730,541.0453412,-1575.603839,-5727.252131,3063.184878,1713.531964,-273.0630735,-958.2825689,-3811.487979,-1922.775363,710.8478819,77.11849915,738.1353454,-2643.044701,-1882.233677,-912.0519091,-19.71186338,-708.9974082,3985.834475,1923.044863,-4062.157976,4.001611671,346.3415977,-3130.859235,-982.6165285,-225.7216245,-0.6255066503,-112.2149908,1422.281741,319.6644624,1425.142115,0.7163960360E-01,25.14557958,-355.2243057,-64.39573767,-351.8599166,-0.5125067970E-02,-3.624150323,19.50975308,5.072422786,-224.3906922,0.2149942478E-04,0.1924280756,-0.9674927775E-01,0.8299548447E-01,0.1086388415,0.4196060660E-01,0.6684130935E-01},
                /*tailamhr_e_22.par*/{0.,0.,0.,0.,0.,-160.5346318,-5405.836237,2176.225259,-15306.46149,14086.26259,2054.980104,10761.22876,-8688.621863,-7075.372738,-3913.115650,1072.471941,10194.25852,-6074.895647,-351.0399978,-12452.34540,337.9263571,4210.746945,-2236.824494,-6784.725503,-5559.931347,-122.3033383,-99.68787985,-86.45423979,3206.265844,8013.281241,-451.8515118,8170.767364,-2208.624252,394.8105037,3024.089465,-4289.504144,-7915.835569,13434.59805,-568.8148415,-4403.387263,-1048.969071,11791.50685,-786.8201020,218.0944511,2445.320683,-353.3373418,-5346.666841,-788.0084114,-44.06165935,-691.8949199,-1618.950332,-145.6534115,-3565.785956,4.491220450,86.05910564,-2404.442503,936.4748611,2680.044311,0.3547843381E-01,7.333647575,4787.518187,-408.7727924,-686.7160060,-0.9420965270E-01,-5.600462014,-3791.618296,89.08115933,-112.2127859,0.1803612648E-01,1.315611938,1980.813912,-4.812157260,174.0079503,-0.2097142534E-02,-0.1919935469,-0.1083953570,0.8528017517E-01,0.9723320957E-01,0.4221285305E-01,0.5830108718E-01},
                /*tailamhr_e_32.par*/{0.,0.,0.,0.,0.,-956.9159973,-12108.15199,4251.469821,-21139.17184,17368.22578,641.8654077,3932.128542,-2296.150160,-11184.86910,-1644.613179,-1081.336296,-7008.301772,3728.485337,6258.849986,6090.003027,-251.7646059,-7981.447780,1697.170133,1935.097152,16136.38716,1269.698288,17924.55797,-6664.665731,-1170.854394,-14743.58864,-1479.249159,14505.18225,-540.9172571,388.3795249,17804.66889,-7953.394730,875.8009245,20317.53458,-50.61437331,-4246.442977,-4786.405305,3981.494912,16332.30961,3.368467890,814.9491747,5578.575168,496.3125703,-10778.58171,0.5021401607,354.2374592,886.3949096,-2048.317963,224.9749791,-0.1869855573,-313.5133707,-1928.355906,1098.248667,1835.670791,0.2882852258E-01,108.5305517,2986.111460,-285.1570650,-330.0170117,-0.2775056713E-02,-21.92505308,-2859.789608,24.63761222,-416.2150298,0.1679870467E-03,2.398503812,1695.127213,10.61347613,339.5039575,-0.3783968109E-05,0.5014646290E-01,-0.1090329330,0.7974339968E-01,0.9797696311E-01,0.2940324022E-01,0.6155823157E-01},
                /*tailamhr_e_42.par*/{0.,0.,0.,0.,0.,-357.2107920,-21779.34849,2344.338145,-11517.58347,27273.30578,252.5834848,-2221.374481,-852.9594572,-11541.58617,7114.454121,-36.22575206,10169.79128,-398.2684635,13184.20308,-15349.12472,-460.6763511,-5961.473660,1809.118875,-3950.208442,4612.492705,-1008.441855,-1630.201605,3866.112272,-3789.011794,-11927.94152,-1324.929732,-16870.32361,6198.241555,1803.877769,-2350.457635,-570.0633956,5537.960149,1601.791156,-713.8468942,-1710.360176,-1065.696994,2594.341735,5441.339412,199.9289232,3049.669055,1569.277545,-466.8358084,367.8775081,-28.52832225,-655.2954163,4456.588863,-1201.321971,-6693.978326,-0.4695404475,-387.3655013,1586.281669,984.1991720,7773.341327,1.192836442,328.9192972,-625.2965810,-394.0873774,-4036.101526,-0.2969851570,-121.0580425,-2926.034607,97.26037666,974.1650803,0.4394252805E-01,27.64123627,3431.281860,-13.94348656,76.22688987,-0.4346068044E-02,-3.918900319,-0.1186105318,0.7949697534E-01,0.1037116210,0.4350625618E-01,0.7145210367E-01},
                /*tailamhr_e_52.par*/{0.,0.,0.,0.,0.,548.7059122,-1221.778391,-498.2371056,-1290.427153,2152.768253,1833.408335,2412.287085,-3087.080045,771.2911136,-1639.795674,-1009.048421,-6591.314865,3358.282077,5467.774993,3828.156667,7222.456913,-4014.369878,-9056.505074,11625.85902,7949.112961,11799.74230,-978.5540825,-16305.14247,399.2419484,10482.67607,3973.063929,-13997.28121,469.0496744,1534.038972,-5907.105565,10769.24528,-20261.64830,-7741.697836,-2392.083464,-13566.53796,2806.396473,-6772.403489,-4432.025772,11.57873667,-2063.497053,-14497.02896,18265.00380,9522.048215,359.7159760,6810.853224,903.7382632,-12961.22161,-13429.79931,-153.1846769,-3997.157623,-5595.408071,4273.832343,2560.812039,34.02675689,1199.423524,9148.791839,-418.6882611,3656.632058,-4.451482025,-180.6403105,-6377.914157,-258.7268662,-3270.904141,0.2294848080,-6.812457862,2845.520419,154.9974731,1513.125122,0.4336022530E-01,11.79623091,-0.9720390077E-01,0.8166535722E-01,0.9359642020E-01,0.4996174383E-01,0.6999451995E-01}
            },{
                /*tailamhr_e_13.par*/{0.,0.,0.,0.,0.,-12510.99836,23741.30362,1471.271517,14532.75294,-23557.17662,-4789.916584,10665.99003,443.4387757,12564.74323,-13900.22009,-4587.563439,5691.880396,666.9007434,-24883.73605,4927.655675,-5747.861409,13803.78991,399.9872563,-10404.86842,-14693.09545,10237.57775,-21377.06693,-1002.464267,12731.90847,20804.56003,-6892.208596,12903.10078,830.3210808,-4807.899188,-12068.07855,6272.637643,-12369.79829,-584.7546086,622.1424912,-591.8245198,-8155.726091,8837.785118,915.8422888,117.5963947,3248.713269,4777.079731,-5064.496243,-1875.322617,-75.95799421,-1894.536729,-3596.249880,1480.096313,-824.4282862,18.63763084,597.2662867,2409.904869,8.835158823,-571.9380738,-2.791790078,-110.9033645,-1207.272145,-216.5451881,1883.133328,0.2445648078,7.025073389,442.8589109,112.9851745,-1775.453457,-0.1244727200E-02,2.855928844,-114.1302259,-35.33544733,1113.666197,-0.3647561303E-02,-1.199840263,-0.9461100060E-01,0.8382703222E-01,0.1093372916,0.4489173259E-01,0.6570948143E-01},
                /*tailamhr_e_23.par*/{0.,0.,0.,0.,0.,387.4572241,7806.929525,-1402.781098,10224.23225,-13158.34457,-601.6365433,-7211.893510,1671.217958,-15648.57677,13167.90009,377.8157383,13321.25180,-1899.731125,-1013.014330,-20436.21140,-2155.420381,-4086.208958,5059.809786,-5442.210698,-7848.208719,-2059.092636,-13511.09786,5898.876091,1410.982746,357.7801407,-857.1215503,-12204.18506,3330.023529,-611.8258485,-7399.158959,2875.561408,2665.404552,-7137.229225,27.02489231,945.1906461,2401.423164,1120.764298,-9237.702631,28.48141724,780.1917799,-2674.698329,-3361.575154,-2685.006266,-13.14073883,-845.8844136,-7475.091079,2093.981054,5503.505028,3.065952618,365.9986405,5.181054424,-785.3582433,-5225.821957,-0.4783455054,-101.1137987,2778.167173,196.6910080,2619.448000,0.5242088187E-01,19.16778579,-1104.886585,-30.43601946,-662.7053491,-0.3736829299E-02,-2.309121571,-54.00018450,1.172106403,11.72438395,0.8456801975E-04,0.7960810763E-01,-0.1098948178,0.7501279604E-01,0.1001943421,0.3774202875E-01,0.6107579857E-01},
                /*tailamhr_e_33.par*/{0.,0.,0.,0.,0.,1114.404110,23849.81449,-6592.016308,59570.34212,-51352.25184,-27.97272393,6030.981398,-655.5342866,45719.50403,-22122.59993,414.3043798,13991.43345,-3319.505243,-34567.43688,-22798.44697,-1076.820762,1814.098929,3586.445885,-3576.021184,-35194.84425,-1802.697685,-31997.32604,10209.73294,5425.508492,42578.81157,-1197.493975,-13147.07361,6182.394052,-1965.284079,-29048.45815,-658.0380304,-1191.075734,2779.224069,305.0984266,7094.265251,-188.7611992,16979.57837,790.3716351,-8.658435152,529.2306055,1225.873734,-12152.63985,2700.287649,-5.543880407,-822.1238026,5422.143355,4173.488719,-9302.055657,1.257346292,260.1128832,2856.449668,-556.1596653,11604.79852,-0.1455880270,-44.00512662,-348.4839046,-141.8479933,-7168.882408,0.9621479085E-02,3.250042257,-5408.375465,97.22597454,2643.700103,-0.1261227320E-03,0.4039691189,6170.977473,-28.19126115,-567.1831438,-0.5445297726E-04,-0.1741798789,-0.1173227513,0.7812928016E-01,0.1007329280,0.3300023726E-01,0.5562665838E-01},
                /*tailamhr_e_43.par*/{0.,0.,0.,0.,0.,1905.715873,10618.69862,-9774.232125,12407.89971,-11269.92332,553.2835823,3735.279433,-3240.802524,3630.548093,-3869.569152,1754.944065,12143.92944,-10413.06100,-11108.41131,-8345.577623,-1161.956926,1919.202869,1642.674775,-6173.139534,-18652.64479,-4063.848935,-21761.04058,20646.65433,3040.300954,7532.016882,-686.8912691,-2989.299966,3319.529545,-2006.839954,-12985.02912,196.0217172,7604.865620,-4076.097952,371.9777937,3857.053692,-1517.183984,3879.195756,3760.312802,16.21805076,276.9667006,840.6119053,-379.3341626,4518.844056,-16.21971941,-355.0879416,8004.182067,-4697.461870,-8612.312372,2.431803057,61.87395896,400.5148583,5941.983460,10262.79853,0.4155617423E-01,11.66296226,-81.18076384,-3240.417914,-5650.909145,-0.7649688510E-01,-8.163209939,-4007.637199,911.2064004,1266.127154,0.1525775799E-01,2.051540656,4519.164304,-56.09611767,308.5174151,-0.1737344839E-02,-0.3035047160,-0.1189154977,0.9992129127E-01,0.1079711756,0.4232104853E-01,0.6049450374E-01},
                /*tailamhr_e_53.par*/{0.,0.,0.,0.,0.,-2055.990625,3973.419038,-3804.537867,-3785.871973,5360.032464,3924.741065,-8553.693203,8843.626789,1796.299125,-6994.231884,-1406.876206,-644.6175509,6287.783493,-13527.25889,1169.479219,-10430.31429,12696.14242,6666.844032,-12326.69059,-15077.93692,-14323.67665,21109.36965,-287.7744305,-10589.09900,-18952.72128,1785.979142,-6533.920956,9256.984700,-931.2362851,3095.239490,-5061.471682,3617.798178,9849.237239,4822.453061,8336.650897,-1422.824407,6168.532565,6315.256959,-506.8608713,-251.6057439,16121.47641,-3088.441738,-11645.43294,-652.5349848,-2355.517484,3405.758554,11875.44498,10116.15215,398.6206081,1566.883439,-2969.790597,-6026.299586,-4135.478463,-114.2671122,-510.6456289,-2212.630399,-18.24079465,755.0655355,19.22578627,94.61572772,2777.061626,1282.370964,74.36348528,-1.586580984,-6.916705331,-1465.270451,-749.8913404,-95.32802084,-0.1224748672,-1.644761412,-0.9469716563E-01,0.9107284134E-01,0.8084435769E-01,0.5588452995E-01,0.6430659167E-01}
            },{
                /*tailamhr_e_14.par*/{0.,0.,0.,0.,0.,3968.169109,-8129.605534,-1234.785613,-5324.362057,9809.099364,503.3910459,-474.9596156,-191.0029229,1020.994781,-341.9188240,2944.203970,-4100.977662,-1002.655911,6571.693376,274.2527157,961.9400286,-2194.237486,-267.8497029,2055.055210,1936.310585,-1216.651434,4193.800190,203.9428018,-4069.466415,-5109.422115,567.4301756,-502.6431810,-401.1933806,2022.075464,4400.871928,1239.900044,3287.277452,-1056.272167,-293.6256440,-122.6689234,3151.639774,-1321.874528,-1593.468872,-32.93885185,-530.0775262,-3018.082970,765.8213684,291.6584405,35.73888961,402.1813290,593.7039368,-549.7190023,-1555.084524,-12.18164510,-182.5095370,-485.3062851,144.8446129,-121.4387088,2.376904969,45.86651966,670.7525145,29.56484800,1360.052389,-0.2536461701,-4.322915405,-448.6731276,-36.87766412,-1111.243899,-0.2390353640E-03,-1.242397877,187.7281341,15.37891652,509.6520008,0.5782472177E-02,0.6512897010,-0.9901923343E-01,0.8276619165E-01,0.1078751443,0.4952637425E-01,0.6713525932E-01},
                /*tailamhr_e_24.par*/{0.,0.,0.,0.,0.,-10078.14225,-29426.86961,25082.98545,-23253.02364,30725.03575,-526.8927921,408.7001197,851.0795171,7854.513212,-4116.467712,-8842.314975,-22013.53360,21267.19114,18675.47393,11060.76519,-37.93926918,1242.252998,-158.6290427,-219.8867251,-3879.951117,-3118.350507,-130.1528642,5956.137418,-5724.730567,-20430.71926,946.2629542,-1753.231782,-1187.375300,1645.831088,8355.346344,1162.096889,9065.033297,-4408.112151,13.92261950,2114.674352,-1622.994128,-7360.439678,6332.193292,-105.5783064,-2585.074669,4223.430502,2054.681311,-5410.973304,29.57546567,939.3832765,-509.7225707,1091.244102,5961.575933,-4.154189194,-151.0047801,1682.353752,-1052.246540,-2995.985544,0.2521868710,-5.190741448,-2943.535526,409.0241063,483.7412290,0.2449791415E-01,9.202939601,2354.414273,-91.34426338,301.6918549,-0.8708658753E-02,-2.667584115,-1248.470310,7.642033421,-274.4283732,0.1294421415E-02,0.4752397742,-0.1036609460,0.8403060843E-01,0.9727927609E-01,0.4084135690E-01,0.6184064129E-01},
                /*tailamhr_e_34.par*/{0.,0.,0.,0.,0.,-1275.648990,-16634.06866,5457.084923,-44566.29953,39953.33906,-382.6256178,-5412.946644,1665.083584,-23637.92660,15485.72390,-1221.847867,-11933.79394,4854.403342,26511.25653,13365.17784,-890.9856459,-5139.316012,3235.687519,-4196.229995,-3541.229228,-2041.693204,-16810.58151,8407.726119,-3140.169790,-22855.40133,-1145.104310,-20395.21396,5995.349028,1122.172198,9353.779118,-389.8368869,14793.43124,60.63135053,-253.7280804,-3105.088874,-1616.715476,-3789.710913,8934.764172,40.11458783,807.1039008,2151.727150,3524.180699,1446.685231,-2.311318447,-7.193271953,7032.966745,-3204.893454,-9345.213243,-0.6165550082,-80.43086687,2844.492232,1753.063746,10152.83298,0.2049177237,34.52739845,-1404.606742,-597.7522241,-4818.801646,-0.3286595756E-01,-8.340759608,-4100.460991,127.0571096,800.7150942,0.3455058546E-02,1.325885625,4913.065521,-12.28173058,369.7121028,-0.2398171364E-03,-0.1317127542,-0.1189225949,0.8021135089E-01,0.1042692594,0.3527838589E-01,0.5507730871E-01},
                /*tailamhr_e_44.par*/{0.,0.,0.,0.,0.,-1035.340294,-19041.12744,4306.112238,-30205.58228,41533.73028,-546.6152490,-9506.002391,2237.215734,-15228.45639,20453.68449,-1054.319469,-8499.257654,3837.043445,33980.55720,-7761.867948,-1170.785353,-8887.184949,4265.679002,-6106.694456,2653.061769,-1788.086226,-11564.17668,7291.938596,-15543.91438,-32071.31554,-779.3392626,-32746.54341,6004.235763,5966.735179,9174.021916,1698.635135,19566.14776,-5918.326066,-1821.610571,-3188.337739,176.6932945,-6638.170696,4280.015086,389.6698080,887.6286850,3148.183908,3599.391359,980.7502289,-21.74170518,127.5655223,6254.196671,-2386.314848,-4057.861070,-16.81658724,-216.4139974,2478.674252,1209.846264,5922.630295,7.033170544,97.71946530,-2711.696967,-408.0433544,-3177.573562,-1.548487607,-26.41705085,-1822.484830,87.89234258,507.0666132,0.2249035314,4.725278645,3180.746106,-9.218895112,330.9399830,-0.2123780234E-01,-0.5198722176,-0.1196535072,0.7900564616E-01,0.1051433708,0.4855125877E-01,0.6153614083E-01},
                /*tailamhr_e_54.par*/{0.,0.,0.,0.,0.,3756.093256,-4635.142096,1979.391861,2266.916067,-2899.035622,-4246.522389,5428.713032,-2956.981230,-4796.011338,5179.651742,8592.615258,-10355.00526,3004.717657,-7470.282026,791.3326517,-9598.482229,11176.94867,-735.5418407,2830.452813,-6622.809340,902.1370448,-676.0796821,1200.411423,-5989.671363,-11251.20486,499.1976667,416.8553631,-3190.578409,89.25759967,-496.7378854,166.1082826,-1057.165205,4709.051329,768.9340365,3674.798002,-2650.339850,3730.556600,11227.34700,-21.35099478,629.8509760,9257.962434,1103.380137,-13387.60881,-103.1094049,-1444.007275,2438.343344,6628.777118,9178.773828,46.68221869,786.1464851,233.4319721,-1787.785763,-3377.783281,-10.92071600,-233.5516877,-4241.829329,-2927.122153,530.1582795,1.556760991,40.64940999,3730.517702,2940.747482,110.6329197,-0.1155552762,-2.695074540,-1864.628869,-1499.559814,-99.01627592,-0.5230514058E-02,-0.7336180231,-0.1011139088,0.9985287476E-01,0.8542176011E-01,0.4916933442E-01,0.6462942234E-01}
            }
        };
    }
    
    TS07Context::TS07Context() {
        // Copy parameters in COLUMN-MAJOR ORDER
        for (int i=0; i<80; i++) {
            for (int j=0; j<5; j++) {
                TSS_TSS[i+1][j+1] = __ts07_private::TSS[j][i];
                for (int k=0; k<4; k++) {
                    TSE_TSE[i+1][j+1][k+1] = __ts07_private::TSE[k][j][i];
                    TSO_TSO[i+1][j+1][k+1] = __ts07_private::TSO[k][j][i];
                }
            }
        }
    }
    
}