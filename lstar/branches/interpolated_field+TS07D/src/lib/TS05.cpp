//
//  TS05.cpp
//  UBJDevelopment
//
//  Translated by Kyungguk Min on 4/6/13.
//  Original code by N. Tsyganenko.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#include "TS05.h"
#include "Geopack.h"
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace UBK {
    using namespace std;

    TS05::TS05 (Geopack const* geopack, double const parmod[]) : TSExternalField(geopack, 1, parmod)
    {
        if (NULL == parmod) {
            throw invalid_argument("PARMOD must be 10-element array.");
        }
    }

    static void T04_s (int const IOPT,double const PARMOD[],double const PS,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ);
    void TS05::getFieldInGSW_atPoint(Point *bOut, const Point ptgsw) const
    {
        T04_s(this->iopt(), this->parmod(), this->geopack()->psi(), ptgsw.x, ptgsw.y, ptgsw.z, &bOut->x, &bOut->y, &bOut->z);
    }

    /*-----------------------------Model---------------------*/
    /* Common blocks
     */
    struct TS_EXTERNAL_COMMON_ {
        // COMMON /TAIL/
        double DXSHIFT1,DXSHIFT2,D,DELTADY;
        // COMMON /BIRKPAR/
        double XKAPPA1,XKAPPA2;
        // COMMON /RCPAR/
        double SC_SY,SC_AS,PHI;
        // COMMON /G/
        double G;
        // COMMON /RH0/
        double RH0;
        // COMMON /DPHI_B_RHO0/
        double DPHI,B,RHO_0,XKAPPA;
        // COMMON /MODENUM/
        int M;
        // COMMON /DTHETA/
        double DTHETA;
    };

    /*
     c===========================================================================
     c
     SUBROUTINE DIPOLE (PS,X,Y,Z,BX,BY,BZ)
     C
     C      A DOUBLE PRECISION ROUTINE
     C
     C  CALCULATES GSM COMPONENTS OF A GEODIPOLE FIELD WITH THE DIPOLE MOMENT
     C  CORRESPONDING TO THE EPOCH OF 2000.
     C
     C----INPUT PARAMETERS:
     C     PS - GEODIPOLE TILT ANGLE IN RADIANS,
     C     X,Y,Z - GSM COORDINATES IN RE (1 RE = 6371.2 km)
     C
     C----OUTPUT PARAMETERS:
     C     BX,BY,BZ - FIELD COMPONENTS IN GSM SYSTEM, IN NANOTESLA.
     C
     */
    static void DIPOLE (double const PS,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {
        double SPS=sin(PS);
        double CPS=cos(PS);
        double P=X*X;
        double U=Z*Z;
        double V=3.e0*Z*X;
        double T=Y*Y;
        double Q=30115.e0/pow(sqrt(P+T+U), 5);
        *BX=Q*((T+U-2.e0*P)*SPS-V*CPS);
        *BY=-3.e0*Y*Q*(X*SPS+Z*CPS);
        *BZ=Q*((P+T-2.e0*U)*CPS-V*SPS);
        return;
    }//END

    /*
     C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     C
     SUBROUTINE TAILDISK(D0,DELTADX,DELTADY,X,Y,Z,BX,BY,BZ)
     c
     c      THIS SUBROUTINE COMPUTES THE COMPONENTS OF THE TAIL CURRENT FIELD,
     C       SIMILAR TO THAT DESCRIBED BY TSYGANENKO AND PEREDO (1994).  THE
     C        DIFFERENCE IS THAT NOW WE USE SPACEWARPING, AS DESCRIBED IN OUR
     C         PAPER ON MODELING BIRKELAND CURRENTS (TSYGANENKO AND STERN, 1996)
     C          INSTEAD OF SHEARING IT IN THE SPIRIT OF T89 TAIL MODEL.
     C
     */
    static void TAILDISK(double const D0,double const DELTADX,double const DELTADY,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {
        //DIMENSION F(5),B(5),C(5)

        static double const F[6] = {0./**/,-71.09346626e0,-1014.308601e0,-1272.939359e0,
            -3224.935936e0,-44546.86232e0};
        static double const B[6] = {0./**/,10.90101242e0,12.68393898e0,13.51791954e0,14.86775017e0,
            15.12306404e0};
        static double const C[6] = {0./**/,.7954069972e0,.6716601849e0,1.174866319e0,2.565249920e0,
            10.01986790e0};

        double RHO=sqrt(X*X+Y*Y);
        double DRHODX=X/RHO;
        double DRHODY=Y/RHO;

        double DEX=exp(X/7.e0);
        double D=D0+DELTADY*pow((Y/20.e0), 2)  +DELTADX*DEX;// !   THE LAST TERM (INTRODUCED 10/11/2000) MAKES THE SHEET
        double DDDY=DELTADY*Y*0.005e0;//                  !   THICKEN SUNWARD, TO AVOID PROBLEMS IN THE SUBSOLAR REGION
        double DDDX=DELTADX/7.e0*DEX;

        double DZETA=sqrt(Z*Z+D*D);//  !  THIS IS THE SAME SIMPLE WAY TO SPREAD
        //C                                        OUT THE SHEET, AS THAT USED IN T89
        double DDZETADX=D*DDDX/DZETA;
        double DDZETADY=D*DDDY/DZETA;
        double DDZETADZ=Z/DZETA;

        double DBX=0.e0;
        double DBY=0.e0;
        double DBZ=0.e0;

        for (int I=1; I<=5; I++) {//DO 1 I=1,5
            double BI=B[I];
            double CI=C[I];

            double S1=sqrt(pow((RHO+BI), 2)+pow((DZETA+CI), 2));
            double S2=sqrt(pow((RHO-BI), 2)+pow((DZETA+CI), 2));

            double DS1DRHO=(RHO+BI)/S1;
            double DS2DRHO=(RHO-BI)/S2;
            double DS1DDZ=(DZETA+CI)/S1;
            double DS2DDZ=(DZETA+CI)/S2;

            double DS1DX=DS1DRHO*DRHODX  +DS1DDZ*DDZETADX;
            double DS1DY=DS1DRHO*DRHODY  +DS1DDZ*DDZETADY;
            double DS1DZ=                 DS1DDZ*DDZETADZ;

            double DS2DX=DS2DRHO*DRHODX  +DS2DDZ*DDZETADX;
            double DS2DY=DS2DRHO*DRHODY  +DS2DDZ*DDZETADY;
            double DS2DZ=                 DS2DDZ*DDZETADZ;

            double S1TS2=S1*S2;
            double S1PS2=S1+S2;
            double S1PS2SQ=S1PS2*S1PS2;

            double FAC1=sqrt(S1PS2SQ-pow((2.e0*BI), 2));
            double AS=FAC1/(S1TS2*S1PS2SQ);
            double DASDS1=(1.e0/(FAC1*S2)-AS/S1PS2*(S2*S2+S1*(3.e0*S1+4.e0*S2)))
            /(S1*S1PS2);
            double DASDS2=(1.e0/(FAC1*S1)-AS/S1PS2*(S1*S1+S2*(3.e0*S2+4.e0*S1)))
            /(S2*S1PS2);

            double DASDX=DASDS1*DS1DX+DASDS2*DS2DX;
            double DASDY=DASDS1*DS1DY+DASDS2*DS2DY;
            double DASDZ=DASDS1*DS1DZ+DASDS2*DS2DZ;

            DBX=DBX-F[I]*X*DASDZ;
            DBY=DBY-F[I]*Y*DASDZ;
            DBZ=DBZ+F[I]*(2.e0*AS+X*DASDX+Y*DASDY); // 1
        }

        *BX=DBX;
        *BY=DBY;
        *BZ=DBZ;

        return;
    }//END

    /*
     C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     C
     C THIS CODE RETURNS THE SHIELDING FIELD REPRESENTED BY  5x5=25 "CARTESIAN"
     C    HARMONICS
     C
     SUBROUTINE  SHLCAR5X5(A,X,Y,Z,DSHIFT,HX,HY,HZ)
     C
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C  The NLIN coefficients are the amplitudes of the "cartesian"
     c    harmonics (A(1)-A(NLIN).
     c  The NNP nonlinear parameters (A(NLIN+1)-A(NTOT) are the scales Pi and Ri
     C   entering the arguments of exponents, sines, and cosines in each of the
     C   NLIN "Cartesian" harmonics
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C
     */
    static void SHLCAR5X5(double const A[],double const X,double const Y,double const Z,double const DSHIFT,double *HX,double *HY,double *HZ) {
        //DIMENSION A(60)

        double DHX=0.e0;
        double DHY=0.e0;
        double DHZ=0.e0;

        int L=0;

        for (int I=1; I<=5; I++) {//DO 2 I=1,5
            double RP=1.e0/A[50+I];
            double CYPI=cos(Y*RP);
            double SYPI=sin(Y*RP);

            for (int K=1; K<=5; K++) {//DO 2 K=1,5
                double RR=1.e0/A[55+K];
                double SZRK=sin(Z*RR);
                double CZRK=cos(Z*RR);
                double SQPR=sqrt(RP*RP+RR*RR);
                double EPR=exp(X*SQPR);

                double DBX=-SQPR*EPR*CYPI*SZRK;
                double DBY= RP*EPR*SYPI*SZRK;
                double DBZ=-RR*EPR*CYPI*CZRK;

                L=L+2;
                double COEF=A[L-1]+A[L]*DSHIFT;

                DHX=DHX+COEF*DBX;
                DHY=DHY+COEF*DBY;
                DHZ=DHZ+COEF*DBZ;

            }
        }//2      CONTINUE

        *HX=DHX;
        *HY=DHY;
        *HZ=DHZ;

        return;
    }//END

    /*
     C=====================================================================================
     DOUBLE PRECISION FUNCTION R_S(A,R,THETA)
     */
    static double R_S(double const A[],double const R,double const THETA) {
        //DIMENSION A(31)

        return /*R_S=*/R+A[2]/R+A[3]*R/sqrt(R*R+A[11]*A[11])+A[4]*R/(R*R+A[12]*A[12])
        +(A[5]+A[6]/R+A[7]*R/sqrt(R*R+A[13]*A[13])+A[8]*R/(R*R+A[14]*A[14]))*
        cos(THETA)
        +(A[9]*R/sqrt(R*R+A[15]*A[15])+A[10]*R/pow((R*R+A[16]*A[16]), 2))
        *cos(2.e0*THETA);

    }//END

    /*
     C-----------------------------------------------------------------------------
     C
     DOUBLE PRECISION FUNCTION THETA_S(A,R,THETA)
     */
    static double THETA_S(double const A[],double const R,double const THETA) {
        //DIMENSION A(31)

        return /*THETA_S=*/THETA+(A[17]+A[18]/R+A[19]/(R*R)
                                  +A[20]*R/sqrt(R*R+A[27]*A[27]))*sin(THETA)
        +(A[21]+A[22]*R/sqrt(R*R+A[28]*A[28])
          +A[23]*R/(R*R+A[29]*A[29]))*sin(2.e0*THETA)
        +(A[24]+A[25]/R+A[26]*R/(R*R+A[30]*A[30]))*sin(3.e0*THETA);

    }//END

    /*
     c!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
     c
     SUBROUTINE FIALCOS(R,THETA,PHI,BTHETA,BPHI,N,THETA0,DT)
     C
     C  CONICAL MODEL OF BIRKELAND CURRENT FIELD; BASED ON THE OLD S/R FIALCO (OF 1990-91)

     C  BTN, AND BPN ARE THE ARRAYS OF BTHETA AND BPHI (BTN(i), BPN(i) CORRESPOND TO i-th MODE).
     C   ONLY FIRST  N  MODE AMPLITUDES ARE COMPUTED (N<=10).
     C    THETA0 IS THE ANGULAR HALF-WIDTH OF THE CONE, DT IS THE ANGULAR H.-W. OF THE CURRENT LAYER

     C   NOTE:  BR=0  (BECAUSE ONLY RADIAL CURRENTS ARE PRESENT IN THIS MODEL)
     C
     */
    static void FIALCOS(double const R,double const THETA,double const PHI,double *BTHETA,double *BPHI,int const N,double const THETA0,double const DT) {
        //DIMENSION  BTN(10),BPN(10),CCCOS(10),SSIN(10)
        double BTN[11],BPN[11],CCCOS[11],SSIN[11];

        double SINTE=sin(THETA);
        double RO=R*SINTE;
        double COSTE=cos(THETA);
        double SINFI=sin(PHI);
        double COSFI=cos(PHI);
        double TG=SINTE/(1.e0+COSTE);//   !        TAN(THETA/2)
        double CTG=SINTE/(1.e0-COSTE);//  !        CTG(THETA/2)

        double TETANP=THETA0+DT;
        double TETANM=THETA0-DT;
        double TGP = 0.0, TGM = 0.0, TGM2 = 0.0, TGP2 = 0.0;
        if (THETA >= TETANM) {//GOTO 1IF(THETA.LT.TETANM) GOTO 1
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

        for (int M=1; M<=N; M++) {//DO 2 M=1,N
            TM=TM*TG;
            CCCOS[M]=COSM1*COSFI-SINM1*SINFI;
            SSIN[M]=SINM1*COSFI+COSM1*SINFI;
            COSM1=CCCOS[M];
            SINM1=SSIN[M];
            double T = 0.0, DTT = 0.0, DTT0 = 0.0, FC = 0.0, FC1 = 0.0;
            if (THETA < TETANM) {//IF(THETA.LT.TETANM) THEN
                T=TM;
                DTT=0.5e0*M*TM*(TG+CTG);
                DTT0=0.e0;
            } else if (THETA < TETANP) {//ELSE IF(THETA.LT.TETANP) THEN
                TGM2M=TGM2M*TGM2;
                FC=1.e0/(TGP-TGM);
                FC1=1.e0/(2*M+1);
                double TGM2M1=TGM2M*TGM;
                double TG21=1.e0+TG*TG;
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

            BTN[M]=M*T*CCCOS[M]/RO;
            BPN[M]=-DTT*SSIN[M]/R; // 2
        }

        *BTHETA=BTN[N] *800e0;
        *BPHI  =BPN[N] *800e0;

        return;
    }//END

    /*
     C-------------------------------------------------------------------------
     C
     SUBROUTINE ONE_CONE(A,X,Y,Z,BX,BY,BZ)
     c
     c  RETURNS FIELD COMPONENTS FOR A DEFORMED CONICAL CURRENT SYSTEM, FITTED TO A BIOSAVART FIELD
     c    BY SIM_14.FOR.  HERE ONLY THE NORTHERN CONE IS TAKEN INTO ACCOUNT.
     c
     */
    static void ONE_CONE(double const A[],double const X,double const Y,double const Z,double *BX,double *BY,double *BZ, struct TS_EXTERNAL_COMMON_ const* TS_EXTERNAL_COMMON) {
        //DIMENSION A(31)

        //COMMON /DTHETA/ DTHETA
        //COMMON /MODENUM/ M
        double const* DTHETA = &TS_EXTERNAL_COMMON->DTHETA;
        int const* M = &TS_EXTERNAL_COMMON->M;

        static double const DR = 1.e-6;
        static double const DT = 1.e-6;//  !   JUST FOR NUMERICAL DIFFERENTIATION

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
        double BTAST = 0.0, BFAST = 0.0;
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

        double BR     =-RSR/R*STSST*BTAST*DRSDT;
        double BTHETA = RSR*STSST*BTAST*DRSDR;
        double BPHI   = RSR*BFAST*(DRSDR*DTSDT-DRSDT*DTSDR);

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
     SUBROUTINE TWOCONES (A,X,Y,Z,BX,BY,BZ)
     C
     C   ADDS FIELDS FROM TWO CONES (NORTHERN AND SOUTHERN), WITH A PROPER SYMMETRY OF THE CURRENT AND FIELD,
     C     CORRESPONDING TO THE REGION 1 BIRKELAND CURRENTS.
     C
     */
    static void TWOCONES (double const A[],double const X,double const Y,double const Z,double *BX,double *BY,double *BZ, struct TS_EXTERNAL_COMMON_ const* TS_EXTERNAL_COMMON) {
        //DIMENSION A(31)

        double BXN = 0.0, BYN = 0.0, BZN = 0.0, BXS = 0.0, BYS = 0.0, BZS = 0.0;
        ONE_CONE (A,X,Y,Z,&BXN,&BYN,&BZN, TS_EXTERNAL_COMMON);
        ONE_CONE (A,X,-Y,-Z,&BXS,&BYS,&BZS, TS_EXTERNAL_COMMON);
        *BX=BXN-BXS;
        *BY=BYN+BYS;
        *BZ=BZN+BZS;

        return;
    }//END

    /*
     c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     c
     c
     SUBROUTINE BIRK_1N2 (NUMB,MODE,PS,X,Y,Z,BX,BY,BZ)
     C
     C  CALCULATES COMPONENTS  OF REGION 1/2 FIELD IN SPHERICAL COORDS.  DERIVED FROM THE S/R DIPDEF2C (WHICH
     C    DOES THE SAME JOB, BUT INPUT/OUTPUT THERE WAS IN SPHERICAL COORDS, WHILE HERE WE USE CARTESIAN ONES)
     C
     C   INPUT:  NUMB=1 (2) FOR REGION 1 (2) CURRENTS
     C           MODE=1 YIELDS SIMPLE SINUSOIDAL MLT VARIATION, WITH MAXIMUM CURRENT AT DAWN/DUSK MERIDIAN
     C     WHILE MODE=2 YIELDS THE SECOND HARMONIC.
     C
     C
     */
    static void BIRK_1N2 (int const NUMB,int const MODE,double const PS,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ, struct TS_EXTERNAL_COMMON_ * TS_EXTERNAL_COMMON) {
        //DIMENSION A11(31),A12(31),A21(31),A22(31)
        //COMMON /MODENUM/ M
        //COMMON /DTHETA/ DTHETA
        int * M = &TS_EXTERNAL_COMMON->M;
        double * DTHETA = &TS_EXTERNAL_COMMON->DTHETA;

        //COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! THESE PARAMETERS CONTROL DAY-NIGHT ASYMMETRY OF F.A.C., AS FOLLOWS:
        double * DPHI = &TS_EXTERNAL_COMMON->DPHI;
        double * B = &TS_EXTERNAL_COMMON->B;
        double * RHO_0 = &TS_EXTERNAL_COMMON->RHO_0;
        double * XKAPPA = &TS_EXTERNAL_COMMON->XKAPPA;
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
        static double const BETA = 0.9e0;
        static double const RH = 10.e0;
        static double const EPS = 3.e0;// ! parameters of the tilt-dependent deformation of the untilted F.A.C. field

        static double const A11[32] = {0./**/,
            .1618068350e0,-.1797957553e0,2.999642482e0,-.9322708978e0,
            -.6811059760e0,.2099057262e0,-8.358815746e0,
            -14.86033550e0,.3838362986e0,
            -16.30945494e0,4.537022847e0,2.685836007e0,
            27.97833029e0,6.330871059e0,
            1.876532361e0,18.95619213e0,.9651528100e0,
            .4217195118e0,-.08957770020e0,
            -1.823555887e0,.7457045438e0,-.5785916524e0,
            -1.010200918e0,.01112389357e0,
            .09572927448e0,-.3599292276e0,8.713700514e0,
            .9763932955e0,3.834602998e0,
            2.492118385e0,.7113544659e0};

        static double const A12[32] = {0./**/,
            .7058026940e0,-.2845938535e0,5.715471266e0,-2.472820880e0,
            -.7738802408e0,.3478293930e0,-11.37653694e0,
            -38.64768867e0,.6932927651e0,
            -212.4017288e0,4.944204937e0,3.071270411e0,
            33.05882281e0,7.387533799e0,
            2.366769108e0,79.22572682e0,.6154290178e0,
            .5592050551e0,-.1796585105e0,
            -1.654932210e0,.7309108776e0,-.4926292779e0,
            -1.130266095e0,-.009613974555e0,
            .1484586169e0,-.2215347198e0,7.883592948e0,
            .02768251655e0,2.950280953e0,
            1.212634762e0,.5567714182e0};

        static double const A21[32] = {0./**/,
            .1278764024e0,-.2320034273e0,1.805623266e0,-32.37241440e0,
            -.9931490648e0,.3175085630e0,-2.492465814e0,
            -16.21600096e0,.2695393416e0,
            -6.752691265e0,3.971794901e0,14.54477563e0,
            41.10158386e0,7.912889730e0,
            1.258297372e0,9.583547721e0,1.014141963e0,
            .5104134759e0,-.1790430468e0,
            -1.756358428e0,.7561986717e0,-.6775248254e0,
            -.04014016420e0,.01446794851e0,
            .1200521731e0,-.2203584559e0,4.508963850e0,
            .8221623576e0,1.779933730e0,
            1.102649543e0,.8867880020e0};

        static double const A22[32] = {0./**/,
            .4036015198e0,-.3302974212e0,2.827730930e0,-45.44405830e0,
            -1.611103927e0,.4927112073e0,-.003258457559e0,
            -49.59014949e0,.3796217108e0,
            -233.7884098e0,4.312666980e0,18.05051709e0,
            28.95320323e0,11.09948019e0,
            .7471649558e0,67.10246193e0,.5667096597e0,
            .6468519751e0,-.1560665317e0,
            -1.460805289e0,.7719653528e0,-.6658988668e0,.2515179349e-5,
            .02426021891e0,.1195003324e0,-.2625739255e0,
            4.377172556e0,.2421190547e0,
            2.503482679e0,1.071587299e0,.7247997430e0};

        *B=0.5e0;
        *RHO_0=7.e0;

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

        double PHI = 0.0;
        if (Xsc == 0.e0 && Zsc == 0.e0) {// IF (Xsc.EQ.0.D0.AND.Zsc.EQ.0.D0) THEN
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
        double DPHISRHO=-2.e0**B*RHO2*RHO/pow((RHO2+RHO*RHO), 2) *sin(PHI)
        +BETA*PS*pow(R1RH, (EPS-1.e0))*RHO/(RH*Rsc*
                                            pow((1.e0+pow(R1RH, EPS)), (1.e0/EPS+1.e0)));
        double DPHISDY= BETA*PS*pow(R1RH, (EPS-1.e0))*Ysc/(RH*Rsc*
                                                           pow((1.e0+pow(R1RH, EPS)), (1.e0/EPS+1.e0)));

        double SPHICS=sin(PHIS);
        double CPHICS=cos(PHIS);

        double XS= RHO*CPHICS;
        double ZS=-RHO*SPHICS;

        double BXS = 0.0, BYAS = 0.0, BZS = 0.0;
        if (NUMB == 1) {//IF (NUMB.EQ.1) THEN
            if (MODE == 1) TWOCONES (A11,XS,Ysc,ZS,&BXS,&BYAS,&BZS, TS_EXTERNAL_COMMON); //IF (MODE.EQ.1) CALL TWOCONES (A11,XS,Ysc,ZS,BXS,BYAS,BZS)
            if (MODE == 2) TWOCONES (A12,XS,Ysc,ZS,&BXS,&BYAS,&BZS, TS_EXTERNAL_COMMON);// IF (MODE.EQ.2) CALL TWOCONES (A12,XS,Ysc,ZS,BXS,BYAS,BZS)
        } else {//ELSE
            if (MODE == 1) TWOCONES (A21,XS,Ysc,ZS,&BXS,&BYAS,&BZS, TS_EXTERNAL_COMMON);// IF (MODE.EQ.1) CALL TWOCONES (A21,XS,Ysc,ZS,BXS,BYAS,BZS)
            if (MODE == 2) TWOCONES (A22,XS,Ysc,ZS,&BXS,&BYAS,&BZS, TS_EXTERNAL_COMMON);// IF (MODE.EQ.2) CALL TWOCONES (A22,XS,Ysc,ZS,BXS,BYAS,BZS)
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
     C
     SUBROUTINE BIRK_SHL (A,PS,X_SC,X,Y,Z,BX,BY,BZ)
     C
     */
    static void BIRK_SHL (double const A[],double const PS,double const X_SC,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {
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

        for (int M=1; M<=2; M++) {//DO 1 M=1,2     !    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
            //C                          AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
            for (int I=1; I<=3; I++) {//DO 2 I=1,3
                double P=A[72+I];
                double Q=A[78+I];
                double CYPI=cos(Y/P);
                double CYQI=cos(Y/Q);
                double SYPI=sin(Y/P);
                double SYQI=sin(Y/Q);

                for (int K=1; K<=3; K++) {//DO 3 K=1,3
                    double R=A[75+K];
                    double S=A[81+K];
                    double SZRK=sin(Z1/R);
                    double CZSK=cos(Z2/S);
                    double CZRK=cos(Z1/R);
                    double SZSK=sin(Z2/S);
                    double SQPR=sqrt(1.e0/(P*P)+1.e0/(R*R));
                    double SQQS=sqrt(1.e0/(Q*Q)+1.e0/(S*S));
                    double EPR=exp(X1*SQPR);
                    double EQS=exp(X2*SQQS);

                    double HX, HY, HZ, FX, FY, FZ;
                    for (int N=1; N<=2; N++) {//DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
                        //C                                AND N=2 IS FOR THE SECOND ONE

                        for (int NN=1; NN<=2; NN++) {//DO 5 NN=1,2 !   NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
                            //C                                         TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE

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

                            double HXR, HZR;
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

    /*
     c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     C
     DOUBLE PRECISION FUNCTION AP(R,SINT,COST)
     C
     C      Calculates azimuthal component of the vector potential of the symmetric
     c  part of the model ring current.
     C
     */
    static double AP(double const R,double const SINT,double const COST) {
        bool PROX;//   !  INDICATES WHETHER WE ARE TOO CLOSE TO THE AXIS OF SYMMETRY, WHERE THE INVERSION
        //C                                                             OF DIPOLAR COORDINATES BECOMES INACCURATE
        static double const A1 = -456.5289941e0;
        static double const A2 = 375.9055332e0;
        static double const RRC1 = 4.274684950e0;
        static double const DD1 = 2.439528329e0;
        static double const RRC2 = 3.367557287e0;
        static double const DD2 = 3.146382545e0;
        static double const P1 = -0.2291904607e0;
        static double const R1 = 3.746064740e0;
        static double const DR1 = 1.508802177e0;
        static double const DLA1 = 0.5873525737e0;
        static double const P2 = 0.1556236119e0;
        static double const R2 = 4.993638842e0;
        static double const DR2 = 3.324180497e0;
        static double const DLA2 = 0.4368407663e0;
        static double const P3 = 0.1855957207e0;
        static double const R3 = 2.969226745e0;
        static double const DR3 = 2.243367377e0;

        PROX=false;
        double SINT1=SINT;
        double COST1=COST;
        if (SINT1 < 1e-2) {//IF (SINT1.LT.1.D-2) THEN  !  TOO CLOSE TO Z-AXIS;  USE LINEAR INTERPOLATION BETWEEN SINT=0 & SINT=0.01
            SINT1=1.e-2;
            COST1=.99994999875e0;
            PROX=true;
        }//ENDIF

        double ALPHA=SINT1*SINT1/R;//         !  R,THETA -> ALPHA,GAMMA
        double GAMMA=COST1/(R*R);

        double ARG1=-pow(((R-R1)/DR1), 2)-pow((COST1/DLA1), 2);
        double ARG2=-pow(((R-R2)/DR2), 2)-pow((COST1/DLA2), 2);
        double ARG3=-pow(((R-R3)/DR3), 2);

        double DEXP1, DEXP2, DEXP3;
        if (ARG1 < -500.) {//IF (ARG1.LT.-500.D0) THEN        !   TO PREVENT "FLOATING UNDERFLOW" CRASHES
            DEXP1=0.e0;
        } else {//ELSE
            DEXP1=exp(ARG1);
        }//ENDIF

        if (ARG2 < -500.) {//IF (ARG2.LT.-500.D0) THEN
            DEXP2=0.e0;
        } else {//ELSE
            DEXP2=exp(ARG2);
        }//ENDIF

        if (ARG3 < -500.) {//IF (ARG3.LT.-500.D0) THEN
            DEXP3=0.e0;
        } else {//ELSE
            DEXP3=exp(ARG3);
        }//ENDIF


        double ALPHA_S=ALPHA*(1.e0+P1*DEXP1+P2*DEXP2+P3*DEXP3);//     !  ALPHA -> ALPHA_S  (DEFORMED)

        double GAMMA_S=GAMMA;
        double GAMMAS2=GAMMA_S*GAMMA_S;


        double ALSQH=ALPHA_S*ALPHA_S/2.e0;//            !  ALPHA_S,GAMMA_S -> RS,SINTS,COSTS
        double F=64.e0/27.e0*GAMMAS2+ALSQH*ALSQH;
        double Q=pow((sqrt(F)+ALSQH), (1.e0/3.e0));
        double C=Q-4.e0*pow(GAMMAS2, (1.e0/3.e0))/(3.e0*Q);
        if (C < 0.) C = 0.; //IF (C.LT.0.D0) C=0.D0
        double G=sqrt(C*C+4.e0*pow(GAMMAS2, (1.e0/3.e0)));
        double RS=4.e0/((sqrt(2.e0*G-C)+sqrt(C))*(G+C));
        double COSTS=GAMMA_S*RS*RS;
        double SINTS=sqrt(1.e0-COSTS*COSTS);
        double RHOS=RS*SINTS;
        //double RHOS2=RHOS*RHOS;
        double ZS=RS*COSTS;
        /*
         c  1st loop:
         */
        double P=pow((RRC1+RHOS), 2)+ZS*ZS+DD1*DD1;
        double XK2=4.e0*RRC1*RHOS/P;
        double XK=sqrt(XK2);
        double XKRHO12=XK*sqrt(RHOS);

        double XK2S=1.e0-XK2;
        double DL=log(1.e0/XK2S);
        double ELK=1.38629436112e0+XK2S*(0.09666344259e0+XK2S*(0.03590092383e0+
                                                               XK2S*(0.03742563713e0+XK2S*0.01451196212e0))) +DL*
        (0.5e0+XK2S*(0.12498593597e0+XK2S*(0.06880248576e0+
                                           XK2S*(0.03328355346e0+XK2S*0.00441787012e0))));
        double ELE=1.e0+XK2S*(0.44325141463e0+XK2S*(0.0626060122e0+XK2S*
                                                    (0.04757383546e0+XK2S*0.01736506451e0))) +DL*
        XK2S*(0.2499836831e0+XK2S*(0.09200180037e0+XK2S*
                                   (0.04069697526e0+XK2S*0.00526449639e0)));

        double APHI1=((1.e0-XK2*0.5e0)*ELK-ELE)/XKRHO12;
        /*
         c  2nd loop:
         */
        P=pow((RRC2+RHOS), 2)+ZS*ZS+DD2*DD2;
        XK2=4.e0*RRC2*RHOS/P;
        XK=sqrt(XK2);
        XKRHO12=XK*sqrt(RHOS);

        XK2S=1.e0-XK2;
        DL=log(1.e0/XK2S);
        ELK=1.38629436112e0+XK2S*(0.09666344259e0+XK2S*(0.03590092383e0+
                                                        XK2S*(0.03742563713e0+XK2S*0.01451196212e0))) +DL*
        (0.5e0+XK2S*(0.12498593597e0+XK2S*(0.06880248576e0+
                                           XK2S*(0.03328355346e0+XK2S*0.00441787012e0))));
        ELE=1.e0+XK2S*(0.44325141463e0+XK2S*(0.0626060122e0+XK2S*
                                             (0.04757383546e0+XK2S*0.01736506451e0))) +DL*
        XK2S*(0.2499836831e0+XK2S*(0.09200180037e0+XK2S*
                                   (0.04069697526e0+XK2S*0.00526449639e0)));

        double APHI2=((1.e0-XK2*0.5e0)*ELK-ELE)/XKRHO12;

        double AP_RET=A1*APHI1+A2*APHI2;
        if (PROX) AP_RET=AP_RET*SINT/SINT1;//IF (PROX) AP=AP*SINT/SINT1   !   LINEAR INTERPOLATION, IF TOO CLOSE TO THE Z-AXIS

        return AP_RET;
    }//END

    /*
     C&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     C
     SUBROUTINE RC_SYMM (X,Y,Z,BX,BY,BZ)
     */
    static void RC_SYMM (double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {
        static double const DS = 1.e-2;
        static double const DC = 0.99994999875e0;
        static double const D = 1.e-4;
        static double const DRD = 5.e3;//  ! DS=SIN(THETA) AT THE BOUNDARY OF THE LINEARITY
        //c                                                                          REGION; DC=SQRT(1-DS**2);  DRD=1/(2*D)
        double RHO2=X*X+Y*Y;
        double R2=RHO2+Z*Z;
        double R=sqrt(R2);
        double RP=R+D;
        double RM=R-D;
        double SINT=sqrt(RHO2)/R;
        double COST=Z/R;

        double FXY;
        if (SINT < DS) {//IF (SINT.LT.DS) THEN  !  TOO CLOSE TO THE Z-AXIS; USING A LINEAR APPROXIMATION A_PHI~SINT,
            //C                                    TO AVOID THE SINGULARITY PROBLEM
            double A=AP(R,DS,DC)/DS;
            double DARDR=(RP*AP(RP,DS,DC)-RM*AP(RM,DS,DC))*DRD;
            FXY=Z*(2.e0*A-DARDR)/(R*R2);
            *BX=FXY*X;
            *BY=FXY*Y;
            *BZ=(2.e0*A*COST*COST+DARDR*SINT*SINT)/R;

        } else {//ELSE

            double THETA=atan2(SINT,COST);
            double TP=THETA+D;
            double TM=THETA-D;
            double SINTP=sin(TP);
            double SINTM=sin(TM);
            double COSTP=cos(TP);
            double COSTM=cos(TM);
            double BR=(SINTP*AP(R,SINTP,COSTP)-SINTM*AP(R,SINTM,COSTM))
            /(R*SINT)*DRD;
            double BT=(RM*AP(RM,SINT,COST)-RP*AP(RP,SINT,COST))/R*DRD;
            FXY=(BR+BT*COST/SINT)/R;
            *BX=FXY*X;
            *BY=FXY*Y;
            *BZ=BR*COST-BT*SINT;

        }//ENDIF

        return;
    }//END

    /*
     c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     C
     C
     DOUBLE PRECISION FUNCTION APPRC(R,SINT,COST)
     C
     C      Calculates azimuthal component of the vector potential of the symmetric
     c  part of the model PARTIAL ring current.
     C
     */
    static double APPRC(double const R,double const SINT,double const COST) {

        bool PROX;
        static double const A1 = -80.11202281e0;
        static double const A2 = 12.58246758e0;
        static double const RRC1 = 6.560486035e0;
        static double const DD1 = 1.930711037e0;
        static double const RRC2 = 3.827208119e0;
        static double const DD2 = .7789990504e0;
        static double const P1 = .3058309043e0;
        static double const ALPHA1 = .1817139853e0;
        static double const DAL1 = .1257532909e0;
        static double const BETA1 = 3.422509402e0;
        static double const DG1 = .04742939676e0;
        static double const P2 = -4.800458958e0;
        static double const ALPHA2 = -.02845643596e0;
        static double const DAL2 = .2188114228e0;
        static double const BETA2 = 2.545944574e0;
        static double const DG2 = .00813272793e0;
        static double const BETA3 = .35868244e0;
        static double const P3 = 103.1601001e0;
        static double const ALPHA3 = -.00764731187e0;
        static double const DAL3 = .1046487459e0;
        static double const BETA4 = 2.958863546e0;
        static double const DG3 = .01172314188e0;
        static double const BETA5 = .4382872938e0;
        static double const Q0 = .01134908150e0;
        static double const Q1 = 14.51339943e0;
        static double const ALPHA4 = .2647095287e0;
        static double const DAL4 = .07091230197e0;
        static double const DG4 = .01512963586e0;
        static double const Q2 = 6.861329631e0;
        static double const ALPHA5 = .1677400816e0;
        static double const DAL5 = .04433648846e0;
        static double const DG5 = .05553741389e0;
        static double const BETA6 = .7665599464e0;
        static double const BETA7 = .7277854652e0;

        PROX=false;
        double SINT1=SINT;
        double COST1=COST;
        if (SINT1 <  1e-2) {//IF (SINT1.LT.1.D-2) THEN  !  TOO CLOSE TO Z-AXIS;  USE LINEAR INTERPOLATION BETWEEN SINT=0 & SINT=0.01
            SINT1=1.e-2;
            COST1=.99994999875e0;
            PROX=true;
        }//ENDIF

        double ALPHA=SINT1*SINT1/R;//         !  R,THETA -> ALPHA,GAMMA
        double GAMMA=COST1/(R*R);

        double ARG1=-pow((GAMMA/DG1), 2);
        double ARG2=-pow(((ALPHA-ALPHA4)/DAL4), 2)-pow((GAMMA/DG4), 2);

        double DEXP1, DEXP2;
        if (ARG1 < -500.) {//IF (ARG1.LT.-500.D0) THEN        !   TO PREVENT "FLOATING UNDERFLOW" CRASHES
            DEXP1=0.e0;
        }else {//ELSE
            DEXP1=exp(ARG1);
        }//ENDIF

        if (ARG2 < -500.) {//IF (ARG2.LT.-500.D0) THEN        !   TO PREVENT "FLOATING UNDERFLOW" CRASHES
            DEXP2=0.e0;
        } else {//ELSE
            DEXP2=exp(ARG2);
        }//ENDIF

        double ALPHA_S=ALPHA*(1.e0+P1/pow((1.e0+pow(((ALPHA-ALPHA1)/DAL1), 2)), BETA1)
                              *DEXP1+P2*(ALPHA-ALPHA2)/pow((1.e0+pow(((ALPHA-ALPHA2)/DAL2), 2)), BETA2)
                              /pow((1.e0+pow((GAMMA/DG2), 2)), BETA3)
                              +P3*pow((ALPHA-ALPHA3), 2)/pow((1.e0+pow(((ALPHA-ALPHA3)/DAL3), 2)), BETA4)
                              /pow((1.e0+pow((GAMMA/DG3), 2)), BETA5));//     !  ALPHA -> ALPHA_S  (DEFORMED)

        double GAMMA_S=GAMMA*(1.e0+Q0+Q1*(ALPHA-ALPHA4)*DEXP2 //             !  GAMMA -> GAMMA_  (DEFORMED)
                              +Q2*(ALPHA-ALPHA5)/pow((1.e0+pow(((ALPHA-ALPHA5)/DAL5), 2)), BETA6)
                              /pow((1.e0+pow((GAMMA/DG5), 2)), BETA7));

        double GAMMAS2=GAMMA_S*GAMMA_S;

        double ALSQH=ALPHA_S*ALPHA_S/2.e0;//                            !  ALPHA_S,GAMMA_S -> RS,SINTS,COSTS
        double F=64.e0/27.e0*GAMMAS2+ALSQH*ALSQH;
        double Q=pow((sqrt(F)+ALSQH), (1.e0/3.e0));
        double C=Q-4.e0*pow(GAMMAS2, (1.e0/3.e0))/(3.e0*Q);
        if (C < 0.) C=0.e0;//IF (C.LT.0.D0) C=0.D0
        double G=sqrt(C*C+4.e0*pow(GAMMAS2, (1.e0/3.e0)));
        double RS=4.e0/((sqrt(2.e0*G-C)+sqrt(C))*(G+C));
        double COSTS=GAMMA_S*RS*RS;
        double SINTS=sqrt(1.e0-COSTS*COSTS);
        double RHOS=RS*SINTS;
        //double RHOS2=RHOS*RHOS;
        double ZS=RS*COSTS;
        /*
         c  1st loop:
         */
        double P=pow((RRC1+RHOS), 2)+ZS*ZS+DD1*DD1;
        double XK2=4.e0*RRC1*RHOS/P;
        double XK=sqrt(XK2);
        double XKRHO12=XK*sqrt(RHOS);

        double XK2S=1.e0-XK2;
        double DL=log(1.e0/XK2S);
        double ELK=1.38629436112e0+XK2S*(0.09666344259e0+XK2S*(0.03590092383e0+
                                                               XK2S*(0.03742563713e0+XK2S*0.01451196212e0))) +DL*
        (0.5e0+XK2S*(0.12498593597e0+XK2S*(0.06880248576e0+
                                           XK2S*(0.03328355346e0+XK2S*0.00441787012e0))));
        double ELE=1.e0+XK2S*(0.44325141463e0+XK2S*(0.0626060122e0+XK2S*
                                                    (0.04757383546e0+XK2S*0.01736506451e0))) +DL*
        XK2S*(0.2499836831e0+XK2S*(0.09200180037e0+XK2S*
                                   (0.04069697526e0+XK2S*0.00526449639e0)));

        double APHI1=((1.e0-XK2*0.5e0)*ELK-ELE)/XKRHO12;
        /*
         c  2nd loop:
         */
        P=pow((RRC2+RHOS), 2)+ZS*ZS+DD2*DD2;
        XK2=4.e0*RRC2*RHOS/P;
        XK=sqrt(XK2);
        XKRHO12=XK*sqrt(RHOS);

        XK2S=1.e0-XK2;
        DL=log(1.e0/XK2S);
        ELK=1.38629436112e0+XK2S*(0.09666344259e0+XK2S*(0.03590092383e0+
                                                        XK2S*(0.03742563713e0+XK2S*0.01451196212e0))) +DL*
        (0.5e0+XK2S*(0.12498593597e0+XK2S*(0.06880248576e0+
                                           XK2S*(0.03328355346e0+XK2S*0.00441787012e0))));
        ELE=1.e0+XK2S*(0.44325141463e0+XK2S*(0.0626060122e0+XK2S*
                                             (0.04757383546e0+XK2S*0.01736506451e0))) +DL*
        XK2S*(0.2499836831e0+XK2S*(0.09200180037e0+XK2S*
                                   (0.04069697526e0+XK2S*0.00526449639e0)));

        double APHI2=((1.e0-XK2*0.5e0)*ELK-ELE)/XKRHO12;

        double APPRC_RET=A1*APHI1+A2*APHI2;
        if (PROX) APPRC_RET=APPRC_RET*SINT/SINT1;//IF (PROX) APPRC=APPRC*SINT/SINT1   !   LINEAR INTERPOLATION, IF TOO CLOSE TO THE Z-AXIS

        return APPRC_RET;
    }//END

    /*
     c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     C
     SUBROUTINE PRC_SYMM (X,Y,Z,BX,BY,BZ)
     */
    static void PRC_SYMM (double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {
        static double const DS = 1.e-2;
        static double const DC = 0.99994999875e0;
        static double const D = 1.e-4;
        static double const DRD = 5.e3;//  ! DS=SIN(THETA) AT THE BOUNDARY OF THE LINEARITY
        //                          REGION; DC=SQRT(1-DS**2);  DRD=1/(2*D)
        double RHO2=X*X+Y*Y;
        double R2=RHO2+Z*Z;
        double R=sqrt(R2);
        double RP=R+D;
        double RM=R-D;
        double SINT=sqrt(RHO2)/R;
        double COST=Z/R;

        if (SINT < DS) {//IF (SINT.LT.DS) THEN  !  TOO CLOSE TO THE Z-AXIS; USING A LINEAR APPROXIMATION A_PHI~SINT,
            //C                                    TO AVOID THE SINGULARITY PROBLEM
            double A=APPRC(R,DS,DC)/DS;
            double DARDR=(RP*APPRC(RP,DS,DC)-RM*APPRC(RM,DS,DC))*DRD;
            double FXY=Z*(2.e0*A-DARDR)/(R*R2);
            *BX=FXY*X;
            *BY=FXY*Y;
            *BZ=(2.e0*A*COST*COST+DARDR*SINT*SINT)/R;

        }else {//ELSE

            double THETA=atan2(SINT,COST);
            double TP=THETA+D;
            double TM=THETA-D;
            double SINTP=sin(TP);
            double SINTM=sin(TM);
            double COSTP=cos(TP);
            double COSTM=cos(TM);
            double BR=(SINTP*APPRC(R,SINTP,COSTP)-SINTM*APPRC(R,SINTM,COSTM))
            /(R*SINT)*DRD;
            double BT=(RM*APPRC(RM,SINT,COST)-RP*APPRC(RP,SINT,COST))/R*DRD;
            double FXY=(BR+BT*COST/SINT)/R;
            *BX=FXY*X;
            *BY=FXY*Y;
            *BZ=BR*COST-BT*SINT;

        }//ENDIF

        return;
    }//END

    /*
     c@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

     SUBROUTINE FFS(A,A0,DA,F,FA,FS)
     */
    static void FFS(double const A,double const A0,double const DA,double *F,double *FA,double *FS) {
        double SQ1=sqrt(pow((A+A0), 2)+DA*DA);
        double SQ2=sqrt(pow((A-A0), 2)+DA*DA);
        *FA=2.e0/(SQ1+SQ2);
        *F=*FA*A;
        *FS=0.5e0*(SQ1+SQ2)/(SQ1*SQ2)*(1.e0-*F**F);
        return;
    }//END

    /*
     c&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&
     C
     DOUBLE PRECISION FUNCTION BR_PRC_Q (R,SINT,COST)
     C
     Calculates the radial component of the "quadrupole" part of the model partial ring current.
     C
     */
    static double BR_PRC_Q (double const R,double const SINT,double const COST) {
        static double const A1 = -21.2666329e0;
        static double const A2 = 32.24527521e0;
        static double const A3 = -6.062894078e0;
        static double const A4 = 7.515660734e0;
        static double const A5 = 233.7341288e0;
        static double const A6 = -227.1195714e0;
        static double const A7 = 8.483233889e0;
        static double const A8 = 16.80642754e0;
        static double const A9 = -24.63534184e0;
        static double const A10 = 9.067120578e0;
        static double const A11 = -1.052686913e0;
        static double const A12 = -12.08384538e0;
        static double const A13 = 18.61969572e0;
        static double const A14 = -12.71686069e0;
        static double const A15 = 47017.35679e0;
        static double const A16 = -50646.71204e0;
        static double const A17 = 7746.058231e0;
        static double const A18 = 1.531069371e0;
        static double const XK1 = 2.318824273e0;
        static double const AL1 = .1417519429e0;
        static double const DAL1 = .6388013110e-02;
        static double const B1 = 5.303934488e0;
        static double const BE1 = 4.213397467e0;
        static double const XK2 = .7955534018e0;
        static double const AL2 = .1401142771e0;
        static double const DAL2 = .2306094179e-01;
        static double const B2 = 3.462235072e0;
        static double const BE2 = 2.568743010e0;
        static double const XK3 = 3.477425908e0;
        static double const XK4 = 1.922155110e0;
        static double const AL3 = .1485233485e0;
        static double const DAL3 = .2319676273e-01;
        static double const B3 = 7.830223587e0;
        static double const BE3 = 8.492933868e0;
        static double const AL4 = .1295221828e0;
        static double const DAL4 = .01753008801e0;
        static double const DG1 = .01125504083e0;
        static double const AL5 = .1811846095e0;
        static double const DAL5 = .04841237481e0;
        static double const DG2 = .01981805097e0;
        static double const C1 = 6.557801891e0;
        static double const C2 = 6.348576071e0;
        static double const C3 = 5.744436687e0;
        static double const AL6 = .2265212965e0;
        static double const DAL6 = .1301957209e0;
        static double const DRM = .5654023158e0;

        double SINT2=SINT*SINT;
        double COST2=COST*COST;
        double SC=SINT*COST;
        double ALPHA=SINT2/R;
        double GAMMA=COST/(R*R);

        double F, FA, FS;
        FFS(ALPHA,AL1,DAL1,&F,&FA,&FS);
        double D1=SC*pow(F, XK1)/(pow((R/B1), BE1)+1.e0);
        double D2=D1*COST2;

        FFS(ALPHA,AL2,DAL2,&F,&FA,&FS);
        double D3=SC*pow(FS, XK2)/(pow((R/B2), BE2)+1.e0);
        double D4=D3*COST2;

        FFS(ALPHA,AL3,DAL3,&F,&FA,&FS);
        double D5=SC*(pow(ALPHA, XK3))*(pow(FS, XK4))/(pow((R/B3), BE3)+1.e0);
        double D6=D5*COST2;

        double ARGA=pow(((ALPHA-AL4)/DAL4), 2)+1.e0;
        double ARGG=1.e0+pow((GAMMA/DG1), 2);

        double D7=SC/ARGA/ARGG;
        double D8=D7/ARGA;
        double D9=D8/ARGA;
        double D10=D9/ARGA;

        ARGA=pow(((ALPHA-AL5)/DAL5), 2)+1.e0;
        ARGG=1.e0+pow((GAMMA/DG2), 2);

        double D11=SC/ARGA/ARGG;
        double D12=D11/ARGA;
        double D13=D12/ARGA;
        double D14=D13/ARGA;


        double D15=SC/(pow(R, 4)+pow(C1, 4));
        double D16=SC/(pow(R, 4)+pow(C2, 4))*COST2;
        double D17=SC/(pow(R, 4)+pow(C3, 4))*COST2*COST2;

        FFS(ALPHA,AL6,DAL6,&F,&FA,&FS);
        double D18=SC*FS/(1.e0+pow(((R-1.2e0)/DRM), 2));

        double BR_PRC_Q_RET=A1*D1+A2*D2+A3*D3+A4*D4+A5*D5+A6*D6+A7*D7+A8*D8+A9*D9+
        A10*D10+A11*D11+A12*D12+A13*D13+A14*D14+A15*D15+A16*D16+A17*D17+
        A18*D18;

        return BR_PRC_Q_RET;
    }//END

    /*
     C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     C
     DOUBLE PRECISION FUNCTION BT_PRC_Q (R,SINT,COST)
     C
     Calculates the Theta component of the "quadrupole" part of the model partial ring current.
     C
     */
    static double BT_PRC_Q (double const R,double const SINT,double const COST) {
        static double const A1 = 12.74640393e0;
        static double const A2 = -7.516393516e0;
        static double const A3 = -5.476233865e0;
        static double const A4 = 3.212704645e0;
        static double const A5 = -59.10926169e0;
        static double const A6 = 46.62198189e0;
        static double const A7 = -.01644280062e0;
        static double const A8 = .1234229112e0;
        static double const A9 = -.08579198697e0;
        static double const A10 = .01321366966e0;
        static double const A11 = .8970494003e0;
        static double const A12 = 9.136186247e0;
        static double const A13 = -38.19301215e0;
        static double const A14 = 21.73775846e0;
        static double const A15 = -410.0783424e0;
        static double const A16 = -69.90832690e0;
        static double const A17 = -848.8543440e0;
        static double const XK1 = 1.243288286e0;
        static double const AL1 = .2071721360e0;
        static double const DAL1 = .05030555417e0;
        static double const B1 = 7.471332374e0;
        static double const BE1 = 3.180533613e0;
        static double const XK2 = 1.376743507e0;
        static double const AL2 = .1568504222e0;
        static double const DAL2 = .02092910682e0;
        static double const BE2 = 1.985148197e0;
        static double const XK3 = .3157139940e0;
        static double const XK4 = 1.056309517e0;
        static double const AL3 = .1701395257e0;
        static double const DAL3 = .1019870070e0;
        static double const B3 = 6.293740981e0;
        static double const BE3 = 5.671824276e0;
        static double const AL4 = .1280772299e0;
        static double const DAL4 = .02189060799e0;
        static double const DG1 = .01040696080e0;
        static double const AL5 = .1648265607e0;
        static double const DAL5 = .04701592613e0;
        static double const DG2 = .01526400086e0;
        static double const C1 = 12.88384229e0;
        static double const C2 = 3.361775101e0;
        static double const C3 = 23.44173897e0;

        double SINT2=SINT*SINT;
        double COST2=COST*COST;
        //double SC=SINT*COST;
        double ALPHA=SINT2/R;
        double GAMMA=COST/(R*R);

        double F, FA, FS;
        FFS(ALPHA,AL1,DAL1,&F,&FA,&FS);
        double D1=pow(F, XK1)/(pow((R/B1), BE1)+1.e0);
        double D2=D1*COST2;

        FFS(ALPHA,AL2,DAL2,&F,&FA,&FS);
        double D3=pow(FA, XK2)/pow(R, BE2);
        double D4=D3*COST2;

        FFS(ALPHA,AL3,DAL3,&F,&FA,&FS);
        double D5=pow(FS, XK3)*pow(ALPHA, XK4)/(pow((R/B3), BE3)+1.e0);
        double D6=D5*COST2;

        FFS(GAMMA,0.e0,DG1,&F,&FA,&FS);
        double FCC=(1.e0+pow(((ALPHA-AL4)/DAL4), 2));
        double D7 =1.e0/FCC*FS;
        double D8 =D7/FCC;
        double D9 =D8/FCC;
        double D10=D9/FCC;

        double ARG=1.e0+pow(((ALPHA-AL5)/DAL5), 2);
        double D11=1.e0/ARG/(1.e0+pow((GAMMA/DG2), 2));
        double D12=D11/ARG;
        double D13=D12/ARG;
        double D14=D13/ARG;

        double D15=1.e0/(pow(R, 4)+C1*C1);
        double D16=COST2/(pow(R, 4)+C2*C2);
        double D17=COST2*COST2/(pow(R, 4)+C3*C3);

        double BT_PRC_Q_RET=A1*D1+A2*D2+A3*D3+A4*D4+A5*D5+A6*D6+A7*D7+A8*D8+A9*D9+
        A10*D10+A11*D11+A12*D12+A13*D13+A14*D14+A15*D15+A16*D16+A17*D17;

        return BT_PRC_Q_RET;
    }//END

    /*
     C@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     C
     C
     SUBROUTINE PRC_QUAD (X,Y,Z,BX,BY,BZ)
     C
     */
    static void PRC_QUAD (double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {
        static double const D = 1.e-4;
        static double const DD = 2.e-4;
        static double const DS = 1.e-2;
        static double const DC = 0.99994999875e0;

        double RHO2=X*X+Y*Y;
        double R=sqrt(RHO2+Z*Z);
        double RHO=sqrt(RHO2);
        double SINT=RHO/R;
        double COST=Z/R;
        double RP=R+D;
        double RM=R-D;

        if (SINT > DS) {//IF (SINT.GT.DS) THEN
            double CPHI=X/RHO;
            double SPHI=Y/RHO;
            double BR=BR_PRC_Q(R,SINT,COST);
            double BT=BT_PRC_Q(R,SINT,COST);
            double DBRR=(BR_PRC_Q(RP,SINT,COST)-BR_PRC_Q(RM,SINT,COST))/DD;
            double THETA=atan2(SINT,COST);
            double TP=THETA+D;
            double TM=THETA-D;
            double SINTP=sin(TP);
            double COSTP=cos(TP);
            double SINTM=sin(TM);
            double COSTM=cos(TM);
            double DBTT=(BT_PRC_Q(R,SINTP,COSTP)-BT_PRC_Q(R,SINTM,COSTM))/DD;
            *BX=SINT*(BR+(BR+R*DBRR+DBTT)*SPHI*SPHI)+COST*BT;
            *BY=-SINT*SPHI*CPHI*(BR+R*DBRR+DBTT);
            *BZ=(BR*COST-BT*SINT)*CPHI;
        } else {//ELSE
            double ST=DS;
            double CT=DC;
            if (Z < 0.) CT = -DC;//IF (Z.LT.0.D0) CT=-DC
            double THETA=atan2(ST,CT);
            double TP=THETA+D;
            double TM=THETA-D;
            double SINTP=sin(TP);
            double COSTP=cos(TP);
            double SINTM=sin(TM);
            double COSTM=cos(TM);
            double BR=BR_PRC_Q(R,ST,CT);
            double BT=BT_PRC_Q(R,ST,CT);
            double DBRR=(BR_PRC_Q(RP,ST,CT)-BR_PRC_Q(RM,ST,CT))/DD;
            double DBTT=(BT_PRC_Q(R,SINTP,COSTP)-BT_PRC_Q(R,SINTM,COSTM))/DD;
            double FCXY=R*DBRR+DBTT;
            *BX=(BR*(X*X+2.e0*Y*Y)+FCXY*Y*Y)/pow((R*ST), 2)+BT*COST;
            *BY=-(BR+FCXY)*X*Y/pow((R*ST), 2);
            *BZ=(BR*COST/ST-BT)*X/R;
        }//ENDIF

        return;
    }//END

    /*
     C---------------------------------------------------------------------------------------
     C
     SUBROUTINE SRC_PRC (IOPR,SC_SY,SC_PR,PHI,PS,X,Y,Z,BXSRC,BYSRC,
     *    BZSRC,BXPRC,BYPRC,BZPRC)
     C
     C   RETURNS FIELD COMPONENTS FROM A MODEL RING CURRENT, INCLUDING ITS SYMMETRIC PART
     C     AND A PARTIAL RING CURRENT, CLOSED VIA BIRKELAND CURRENTS. BASED ON RESULTS, DESCRIBED
     C     IN A PAPER "MODELING THE INNER MAGNETOSPHERE: ASYMMETRIC RING CURRENT AND REGION 2
     C     BIRKELAND CURRENTS REVISITED" (JGR, DEC.2000).
     C
     C     IOPR -  A RING CURRENT CALCULATION FLAG (FOR LEAST-SQUARES FITTING ONLY):
     C             IOPR=0 - BOTH SRC AND PRC FIELDS ARE CALCULATED
     C             IOPR=1 - SRC ONLY
     C             IOPR=2 - PRC ONLY
     C
     C     SC_SY &  SC_PR ARE SCALE FACTORS FOR THE ABOVE COMPONENTS;  TAKING SC<1 OR SC>1 MAKES THE CURRENTS
     C                      SHRINK OR EXPAND, RESPECTIVELY.
     C
     C   PHI IS THE ROTATION ANGLE (RADIANS) OF THE PARTIAL RING CURRENT (MEASURED FROM MIDNIGHT TOWARD DUSK)
     C
     */
    static void SRC_PRC (int const IOPR,double const SC_SY,double const SC_PR,double const PHI,double const PS,double const X,double const Y,double const Z,double *BXSRC,double *BYSRC,double *BZSRC,double *BXPRC,double *BYPRC,double *BZPRC) {
        /*
         c   1.  TRANSFORM TO TILTED COORDINATES (i.e., SM coordinates):
         */
        double CPS=cos(PS);
        double SPS=sin(PS);

        double XT=X*CPS-Z*SPS;
        double ZT=Z*CPS+X*SPS;
        /*
         C   2.  SCALE THE COORDINATES FOR THE SYMMETRIC AND PARTIAL RC COMPONENTS:
         */
        double XTS=XT/SC_SY;//    !  SYMMETRIC
        double YTS=Y /SC_SY;
        double ZTS=ZT/SC_SY;

        double XTA=XT/SC_PR;//   !  PARTIAL
        double YTA=Y /SC_PR;
        double ZTA=ZT/SC_PR;
        /*
         C   3.  CALCULATE COMPONENTS OF THE TOTAL FIELD IN THE TILTED (SOLAR-MAGNETIC) COORDINATE SYSTEM:
         C
         C
         C    3a. SYMMETRIC FIELD:
         */
        double BXS = 0.0, BYS = 0.0, BZS = 0.0;
        double BXA_S = 0.0, BYA_S = 0.0, BZA_S = 0.0;
        if (IOPR <= 1) RC_SYMM(XTS,YTS,ZTS,&BXS,&BYS,&BZS);//IF (IOPR.LE.1) CALL RC_SYMM(XTS,YTS,ZTS,BXS,BYS,BZS)
        if (IOPR == 0 || IOPR == 2) PRC_SYMM(XTA,YTA,ZTA,&BXA_S,&BYA_S,&BZA_S);//IF (IOPR.EQ.0.OR.IOPR.EQ.2)
        /*
         C    3b. ROTATE THE SCALED SM COORDINATES BY PHI AROUND ZSM AXIS AND CALCULATE QUADRUPOLE PRC FIELD
         C         IN THOSE COORDS:
         */

        double CP=cos(PHI);
        double SP=sin(PHI);
        double XR=XTA*CP-YTA*SP;
        double YR=XTA*SP+YTA*CP;

        double BXA_QR = 0.0, BYA_QR = 0.0, BZA_Q = 0.0;
        if (IOPR==0 || IOPR == 2) PRC_QUAD(XR,YR,ZTA,&BXA_QR,&BYA_QR,&BZA_Q);//IF (IOPR.EQ.0.OR.IOPR.EQ.2)
        /*
         C    3c. TRANSFORM THE QUADRUPOLE FIELD COMPONENTS BACK TO THE SM COORDS:
         */
        double BXA_Q= BXA_QR*CP+BYA_QR*SP;
        double BYA_Q=-BXA_QR*SP+BYA_QR*CP;
        /*
         C    3d. FIND THE TOTAL FIELD OF PRC (SYMM.+QUADR.) IN THE SM COORDS:
         */
        double BXP=BXA_S+BXA_Q;
        double BYP=BYA_S+BYA_Q;
        double BZP=BZA_S+BZA_Q;
        /*
         C   4.  TRANSFORM THE FIELDS OF BOTH PARTS OF THE RING CURRENT BACK TO THE GSM SYSTEM:
         */
        *BXSRC=BXS*CPS+BZS*SPS;//   !    SYMMETRIC RC
        *BYSRC=BYS;
        *BZSRC=BZS*CPS-BXS*SPS;

        *BXPRC=BXP*CPS+BZP*SPS;//   !    PARTIAL RC
        *BYPRC=BYP;
        *BZPRC=BZP*CPS-BXP*SPS;

        return;
    }//END

    /*
     C||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||
     C
     C-------------------------------------------------------------------------
     C
     C
     SUBROUTINE RC_SHIELD (A,PS,X_SC,X,Y,Z,BX,BY,BZ)
     C
     */
    static void RC_SHIELD (double const A[],double const PS,double const X_SC,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {
        //DIMENSION A(86)

        double FAC_SC=pow((X_SC+1.e0), 3);

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

        for (int M=1; M<=2; M++) {//DO 1 M=1,2     !    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
            //C                          AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
            for (int I=1; I<=3; I++) {//DO 2 I=1,3
                double P=A[72+I];
                double Q=A[78+I];
                double CYPI=cos(Y/P);
                double CYQI=cos(Y/Q);
                double SYPI=sin(Y/P);
                double SYQI=sin(Y/Q);

                for (int K=1; K<=3; K++) {//DO 3 K=1,3
                    double R=A[75+K];
                    double S=A[81+K];
                    double SZRK=sin(Z1/R);
                    double CZSK=cos(Z2/S);
                    double CZRK=cos(Z1/R);
                    double SZSK=sin(Z2/S);
                    double SQPR=sqrt(1.e0/(P*P)+1.e0/(R*R));
                    double SQQS=sqrt(1.e0/(Q*Q)+1.e0/(S*S));
                    double EPR=exp(X1*SQPR);
                    double EQS=exp(X2*SQQS);

                    double FX, FY, FZ;
                    double HX, HY, HZ;
                    for (int N=1; N<=2; N++) {//DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
                        //C                                AND N=2 IS FOR THE SECOND ONE

                        for (int NN=1; NN<=2; NN++) {//DO 5 NN=1,2 !   NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
                            //C                                         TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE

                            if (M == 1) {//IF (M.EQ.1) THEN
                                FX=-SQPR*EPR*CYPI*SZRK  *FAC_SC;
                                FY=EPR*SYPI*SZRK/P   *FAC_SC;
                                FZ=-EPR*CYPI*CZRK/R  *FAC_SC;
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
                                FX=-SPS*SQQS*EQS*CYQI*CZSK  *FAC_SC;
                                FY=SPS/Q*EQS*SYQI*CZSK   *FAC_SC;
                                FZ=SPS/S*EQS*CYQI*SZSK   *FAC_SC;
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

                            double HXR, HZR;
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

    /*
     C************************************************************************************
     C
     SUBROUTINE FULL_RC (IOPR,PS,X,Y,Z,BXSRC,BYSRC,BZSRC,BXPRC,BYPRC,
     *  BZPRC)
     C
     C   CALCULATES GSM FIELD COMPONENTS OF THE SYMMETRIC (SRC) AND PARTIAL (PRC) COMPONENTS OF THE RING CURRENT
     C   SRC  PROVIDES A DEPRESSION OF -28 nT AT EARTH
     C   PRC  CORRESPONDS TO THE PRESSURE DIFFERENCE OF 2 nPa BETWEEN MIDNIGHT AND NOON RING CURRENT
     C             PARTICLE PRESSURE AND YIELDS A DEPRESSION OF -17 nT AT X=-6Re
     C
     C   SC_SY AND SC_PR ARE SCALING FACTORS FOR THE SYMMETRIC AND PARTIAL COMPONENTS:
     C          VALUES LARGER THAN 1 RESULT IN SPATIALLY LARGER CURRENTS
     C
     C   PHI IS THE ROTATION ANGLE IN RADIANS OF THE PARTIAL RING CURRENT (MEASURED FROM MIDNIGHT TOWARD DUSK)
     C
     C     IOPR -  A RING CURRENT CALCULATION FLAG (FOR LEAST-SQUARES FITTING ONLY):
     C             IOPR=0 - BOTH SRC AND PRC FIELDS ARE CALCULATED
     C             IOPR=1 - SRC ONLY
     C             IOPR=2 - PRC ONLY
     C
     */
    static void FULL_RC (int const IOPR,double const PS,double const X,double const Y,double const Z,double *BXSRC,double *BYSRC,double *BZSRC,double *BXPRC,double *BYPRC,double *BZPRC, struct TS_EXTERNAL_COMMON_ const* TS_EXTERNAL_COMMON) {
        //DIMENSION C_SY(86),C_PR(86)
        //COMMON /RCPAR/ SC_SY,SC_PR,PHI
        double const* SC_SY = &TS_EXTERNAL_COMMON->SC_SY;
        double const* SC_PR = &TS_EXTERNAL_COMMON->SC_AS;
        double const* PHI = &TS_EXTERNAL_COMMON->PHI;

        static double const C_SY[87] = {0./**/,
            -957.25349e0,-817.5450246e0,583.2991249e0,758.856827e0,
            13.17029064e0,68.94173502e0,-15.29764089e0,
            -53.43151590e0,27.34311724e0,
            149.5252826e0,-11.00696044e0,-179.7031814e0,
            953.0914774e0,817.2340042e0,
            -581.0791366e0,-757.5387665e0,-13.10602697e0,
            -68.58155678e0,15.22447386e0,
            53.15535633e0,-27.07982637e0,-149.1413391e0,
            10.91433279e0,179.3251739e0,
            -6.028703251e0,1.303196101e0,-1.345909343e0,
            -1.138296330e0,-0.06642634348e0,
            -0.3795246458e0,.07487833559e0,.2891156371e0,
            -.5506314391e0,-.4443105812e0,
            0.2273682152e0,0.01086886655e0,-9.130025352e0,
            1.118684840e0,1.110838825e0,
            .1219761512e0,-.06263009645e0,-.1896093743e0,
            .03434321042e0,.01523060688e0,
            -.4913171541e0,-.2264814165e0,-.04791374574e0,
            .1981955976e0,-68.32678140e0,
            -48.72036263e0,14.03247808e0,16.56233733e0,
            2.369921099e0,6.200577111e0,
            -1.415841250e0,-0.8184867835e0,-3.401307527e0,
            -8.490692287e0,3.217860767e0,
            -9.037752107e0,66.09298105e0,48.23198578e0,
            -13.67277141e0,-16.27028909e0,
            -2.309299411e0,-6.016572391e0,1.381468849e0,
            0.7935312553e0,3.436934845e0,
            8.260038635e0,-3.136213782e0,8.833214943e0,
            8.041075485e0,8.024818618e0,
            35.54861873e0,12.55415215e0,1.738167799e0,
            3.721685353e0,23.06768025e0,
            6.871230562e0,6.806229878e0,21.35990364e0,
            1.687412298e0,3.500885177e0,
            0.3498952546e0,0.6595919814e0};

        static double const C_PR[87] = {0./**/,
            -64820.58481e0,-63965.62048e0,66267.93413e0,135049.7504e0,
            -36.56316878e0,124.6614669e0,56.75637955e0,
            -87.56841077e0,5848.631425e0,
            4981.097722e0,-6233.712207e0,-10986.40188e0,
            68716.52057e0,65682.69473e0,
            -69673.32198e0,-138829.3568e0,43.45817708e0,
            -117.9565488e0,-62.14836263e0,
            79.83651604e0,-6211.451069e0,-5151.633113e0,
            6544.481271e0,11353.03491e0,
            23.72352603e0,-256.4846331e0,25.77629189e0,
            145.2377187e0,-4.472639098e0,
            -3.554312754e0,2.936973114e0,2.682302576e0,
            2.728979958e0,26.43396781e0,
            -9.312348296e0,-29.65427726e0,-247.5855336e0,
            -206.9111326e0,74.25277664e0,
            106.4069993e0,15.45391072e0,16.35943569e0,
            -5.965177750e0,-6.079451700e0,
            115.6748385e0,-35.27377307e0,-32.28763497e0,
            -32.53122151e0,93.74409310e0,
            84.25677504e0,-29.23010465e0,-43.79485175e0,
            -6.434679514e0,-6.620247951e0,
            2.443524317e0,2.266538956e0,-43.82903825e0,
            6.904117876e0,12.24289401e0,
            17.62014361e0,152.3078796e0,124.5505289e0,
            -44.58690290e0,-63.02382410e0,
            -8.999368955e0,-9.693774119e0,3.510930306e0,
            3.770949738e0,-77.96705716e0,
            22.07730961e0,20.46491655e0,18.67728847e0,
            9.451290614e0,9.313661792e0,
            644.7620970e0,418.2515954e0,7.183754387e0,
            35.62128817e0,19.43180682e0,
            39.57218411e0,15.69384715e0,7.123215241e0,
            2.300635346e0,21.90881131e0,
            -.01775839370e0,.3996346710e0};

        double HXSRC, HYSRC, HZSRC, HXPRC, HYPRC, HZPRC;
        SRC_PRC (IOPR,*SC_SY,*SC_PR,*PHI,PS,X,Y,Z,&HXSRC,&HYSRC,&HZSRC,
                 &HXPRC,&HYPRC,&HZPRC);

        double X_SC=*SC_SY-1.e0;
        double FSX, FSY, FSZ;
        if (IOPR == 0 || IOPR == 1) {//IF (IOPR.EQ.0.OR.IOPR.EQ.1) THEN
            RC_SHIELD (C_SY,PS,X_SC,X,Y,Z,&FSX,&FSY,&FSZ);
        } else {//ELSE
            FSX=0.e0;
            FSY=0.e0;
            FSZ=0.e0;
        }//ENDIF

        X_SC=*SC_PR-1.e0;
        double FPX = 0., FPY = 0., FPZ = 0.;
        if (IOPR == 0 || IOPR == 2) {//IF (IOPR.EQ.0.OR.IOPR.EQ.2) THEN
            RC_SHIELD (C_PR,PS,X_SC,X,Y,Z,&FPX,&FPY,&FPZ);
        } else {//ELSE
            FPX=0.e0;
            FPY=0.e0;
            FPZ=0.e0;
        }//ENDIF

        *BXSRC=HXSRC+FSX;
        *BYSRC=HYSRC+FSY;
        *BZSRC=HZSRC+FSZ;

        *BXPRC=HXPRC+FPX;
        *BYPRC=HYPRC+FPY;
        *BZPRC=HZPRC+FPZ;

        return;
    }//END

    /*
     c %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     C
     SUBROUTINE BIRK_TOT (IOPB,PS,X,Y,Z,BX11,BY11,BZ11,BX12,BY12,BZ12,
     *                          BX21,BY21,BZ21,BX22,BY22,BZ22)
     C
     C      IOPB -  BIRKELAND FIELD MODE FLAG:
     C         IOPB=0 - ALL COMPONENTS
     C         IOPB=1 - REGION 1, MODES 1 & 2
     C         IOPB=2 - REGION 2, MODES 1 & 2
     C
     */
    static void BIRK_TOT (int const IOPB,double const PS,double const X,double const Y,double const Z,
                          double *BX11,double *BY11,double *BZ11, double *BX12,double *BY12,double *BZ12,
                          double *BX21,double *BY21,double *BZ21,double *BX22,double *BY22,double *BZ22, struct TS_EXTERNAL_COMMON_ * TS_EXTERNAL_COMMON) {
        //DIMENSION SH11(86),SH12(86),SH21(86),SH22(86)
        //COMMON /BIRKPAR/ XKAPPA1,XKAPPA2   !  INPUT PARAMETERS, SPECIFIED FROM A MAIN PROGRAM
        //COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! PARAMETERS, CONTROL DAY-NIGHT ASYMMETRY OF F.A.C.
        double const* XKAPPA1 = &TS_EXTERNAL_COMMON->XKAPPA1;
        double const* XKAPPA2 = &TS_EXTERNAL_COMMON->XKAPPA2;
        //double &DPHI = TS_EXTERNAL_COMMON.DPHI;
        //double &B = TS_EXTERNAL_COMMON.B;
        //double &RHO_0 = TS_EXTERNAL_COMMON.RHO_0;
        double * XKAPPA = &TS_EXTERNAL_COMMON->XKAPPA;

        static double const SH11[87] = {0./**/,46488.84663e0,-15541.95244e0,
            -23210.09824e0,-32625.03856e0,
            -109894.4551e0,-71415.32808e0,58168.94612e0,
            55564.87578e0,-22890.60626e0,
            -6056.763968e0,5091.3681e0,239.7001538e0,
            -13899.49253e0,4648.016991e0,
            6971.310672e0,9699.351891e0,32633.34599e0,
            21028.48811e0,-17395.9619e0,
            -16461.11037e0,7447.621471e0,2528.844345e0,
            -1934.094784e0,-588.3108359e0,
            -32588.88216e0,10894.11453e0,16238.25044e0,
            22925.60557e0,77251.11274e0,
            50375.97787e0,-40763.78048e0,-39088.6066e0,
            15546.53559e0,3559.617561e0,
            -3187.730438e0,309.1487975e0,88.22153914e0,
            -243.0721938e0,-63.63543051e0,
            191.1109142e0,69.94451996e0,-187.9539415e0,
            -49.89923833e0,104.0902848e0,
            -120.2459738e0,253.5572433e0,89.25456949e0,
            -205.6516252e0,-44.93654156e0,
            124.7026309e0,32.53005523e0,-98.85321751e0,
            -36.51904756e0,98.8824169e0,
            24.88493459e0,-55.04058524e0,61.14493565e0,
            -128.4224895e0,-45.3502346e0,
            105.0548704e0,-43.66748755e0,119.3284161e0,
            31.38442798e0,-92.87946767e0,
            -33.52716686e0,89.98992001e0,25.87341323e0,
            -48.86305045e0,59.69362881e0,
            -126.5353789e0,-44.39474251e0,101.5196856e0,
            59.41537992e0,41.18892281e0,
            80.861012e0,3.066809418e0,7.893523804e0,
            30.56212082e0,10.36861082e0,
            8.222335945e0,19.97575641e0,2.050148531e0,
            4.992657093e0,2.300564232e0,
            .2256245602e0,-.05841594319e0};

        static double const SH12[87] = {0./**/,
            210260.4816e0,-1443587.401e0,-1468919.281e0,281939.2993e0,
            -1131124.839e0,729331.7943e0,2573541.307e0,
            304616.7457e0,468887.5847e0,
            181554.7517e0,-1300722.65e0,-257012.8601e0,
            645888.8041e0,-2048126.412e0,
            -2529093.041e0,571093.7972e0,-2115508.353e0,
            1122035.951e0,4489168.802e0,
            75234.22743e0,823905.6909e0,147926.6121e0,
            -2276322.876e0,-155528.5992e0,
            -858076.2979e0,3474422.388e0,3986279.931e0,
            -834613.9747e0,3250625.781e0,
            -1818680.377e0,-7040468.986e0,-414359.6073e0,
            -1295117.666e0,-346320.6487e0,
            3565527.409e0,430091.9496e0,-.1565573462e0,
            7.377619826e0,.4115646037e0,
            -6.14607888e0,3.808028815e0,-.5232034932e0,
            1.454841807e0,-12.32274869e0,
            -4.466974237e0,-2.941184626e0,-.6172620658e0,
            12.6461349e0,1.494922012e0,
            -21.35489898e0,-1.65225696e0,16.81799898e0,
            -1.404079922e0,-24.09369677e0,
            -10.99900839e0,45.9423782e0,2.248579894e0,
            31.91234041e0,7.575026816e0,
            -45.80833339e0,-1.507664976e0,14.60016998e0,
            1.348516288e0,-11.05980247e0,
            -5.402866968e0,31.69094514e0,12.28261196e0,
            -37.55354174e0,4.155626879e0,
            -33.70159657e0,-8.437907434e0,36.22672602e0,
            145.0262164e0,70.73187036e0,
            85.51110098e0,21.47490989e0,24.34554406e0,
            31.34405345e0,4.655207476e0,
            5.747889264e0,7.802304187e0,1.844169801e0,
            4.86725455e0,2.941393119e0,
            .1379899178e0,.06607020029e0};

        static double const SH21[87] = {0./**/,
            162294.6224e0,503885.1125e0,-27057.67122e0,-531450.1339e0,
            84747.05678e0,-237142.1712e0,84133.61490e0,
            259530.0402e0,69196.05160e0,
            -189093.5264e0,-19278.55134e0,195724.5034e0,
            -263082.6367e0,-818899.6923e0,
            43061.10073e0,863506.6932e0,-139707.9428e0,
            389984.8850e0,-135167.5555e0,
            -426286.9206e0,-109504.0387e0,295258.3531e0,
            30415.07087e0,-305502.9405e0,
            100785.34e0,315010.9567e0,-15999.50673e0,
            -332052.2548e0,54964.34639e0,
            -152808.375e0,51024.67566e0,166720.0603e0,
            40389.67945e0,-106257.7272e0,
            -11126.14442e0,109876.2047e0,2.978695024e0,
            558.6019011e0,2.685592939e0,
            -338.000473e0,-81.9972409e0,-444.1102659e0,
            89.44617716e0,212.0849592e0,
            -32.58562625e0,-982.7336105e0,-35.10860935e0,
            567.8931751e0,-1.917212423e0,
            -260.2023543e0,-1.023821735e0,157.5533477e0,
            23.00200055e0,232.0603673e0,
            -36.79100036e0,-111.9110936e0,18.05429984e0,
            447.0481e0,15.10187415e0,
            -258.7297813e0,-1.032340149e0,-298.6402478e0,
            -1.676201415e0,180.5856487e0,
            64.52313024e0,209.0160857e0,-53.8557401e0,
            -98.5216429e0,14.35891214e0,
            536.7666279e0,20.09318806e0,-309.7349530e0,
            58.54144539e0,67.45226850e0,
            97.92374406e0,4.752449760e0,10.46824379e0,
            32.91856110e0,12.05124381e0,
            9.962933904e0,15.91258637e0,1.804233877e0,
            6.578149088e0,2.515223491e0,
            .1930034238e0,-.02261109942e0};

        static double const SH22[87] = {0./**/,
            -131287.8986e0,-631927.6885e0,-318797.4173e0,616785.8782e0,
            -50027.36189e0,863099.9833e0,47680.20240e0,
            -1053367.944e0,-501120.3811e0,
            -174400.9476e0,222328.6873e0,333551.7374e0,
            -389338.7841e0,-1995527.467e0,
            -982971.3024e0,1960434.268e0,297239.7137e0,
            2676525.168e0,-147113.4775e0,
            -3358059.979e0,-2106979.191e0,-462827.1322e0,
            1017607.960e0,1039018.475e0,
            520266.9296e0,2627427.473e0,1301981.763e0,
            -2577171.706e0,-238071.9956e0,
            -3539781.111e0,94628.16420e0,4411304.724e0,
            2598205.733e0,637504.9351e0,
            -1234794.298e0,-1372562.403e0,-2.646186796e0,
            -31.10055575e0,2.295799273e0,
            19.20203279e0,30.01931202e0,-302.1028550e0,
            -14.78310655e0,162.1561899e0,
            .4943938056e0,176.8089129e0,-.2444921680e0,
            -100.6148929e0,9.172262228e0,
            137.4303440e0,-8.451613443e0,-84.20684224e0,
            -167.3354083e0,1321.830393e0,
            76.89928813e0,-705.7586223e0,18.28186732e0,
            -770.1665162e0,-9.084224422e0,
            436.3368157e0,-6.374255638e0,-107.2730177e0,
            6.080451222e0,65.53843753e0,
            143.2872994e0,-1028.009017e0,-64.22739330e0,
            547.8536586e0,-20.58928632e0,
            597.3893669e0,10.17964133e0,-337.7800252e0,
            159.3532209e0,76.34445954e0,
            84.74398828e0,12.76722651e0,27.63870691e0,
            32.69873634e0,5.145153451e0,
            6.310949163e0,6.996159733e0,1.971629939e0,
            4.436299219e0,2.904964304e0,
            .1486276863e0,.06859991529e0};

        *XKAPPA=*XKAPPA1;//        !  FORWARDED IN BIRK_1N2
        double X_SC=*XKAPPA1-1.1e0;//    !  FORWARDED IN BIRK_SHL

        if (IOPB == 0 || IOPB == 1) {//IF (IOPB.EQ.0.OR.IOPB.EQ.1) THEN

            double FX11, FY11, FZ11;
            BIRK_1N2 (1,1,PS,X,Y,Z,&FX11,&FY11,&FZ11, TS_EXTERNAL_COMMON);//           !  REGION 1, MODE 1
            double HX11, HY11, HZ11;
            BIRK_SHL (SH11,PS,X_SC,X,Y,Z,&HX11,&HY11,&HZ11);
            *BX11=FX11+HX11;
            *BY11=FY11+HY11;
            *BZ11=FZ11+HZ11;

            double FX12, FY12, FZ12;
            BIRK_1N2 (1,2,PS,X,Y,Z,&FX12,&FY12,&FZ12, TS_EXTERNAL_COMMON);//           !  REGION 1, MODE 2
            double HX12, HY12, HZ12;
            BIRK_SHL (SH12,PS,X_SC,X,Y,Z,&HX12,&HY12,&HZ12);
            *BX12=FX12+HX12;
            *BY12=FY12+HY12;
            *BZ12=FZ12+HZ12;

        }//ENDIF

        *XKAPPA=*XKAPPA2;//        !  FORWARDED IN BIRK_1N2
        X_SC=*XKAPPA2-1.0e0;//    !  FORWARDED IN BIRK_SHL

        if (IOPB == 0 || IOPB == 2) {//IF (IOPB.EQ.0.OR.IOPB.EQ.2) THEN

            double FX21, FY21, FZ21;
            BIRK_1N2 (2,1,PS,X,Y,Z,&FX21,&FY21,&FZ21, TS_EXTERNAL_COMMON);//           !  REGION 2, MODE 1
            double HX21, HY21, HZ21;
            BIRK_SHL (SH21,PS,X_SC,X,Y,Z,&HX21,&HY21,&HZ21);
            *BX21=FX21+HX21;
            *BY21=FY21+HY21;
            *BZ21=FZ21+HZ21;

            double FX22, FY22, FZ22;
            BIRK_1N2 (2,2,PS,X,Y,Z,&FX22,&FY22,&FZ22, TS_EXTERNAL_COMMON);//           !  REGION 2, MODE 2
            double HX22, HY22, HZ22;
            BIRK_SHL (SH22,PS,X_SC,X,Y,Z,&HX22,&HY22,&HZ22);
            *BX22=FX22+HX22;
            *BY22=FY22+HY22;
            *BZ22=FZ22+HZ22;

        }//ENDIF

        return;
    }//END

    /*
     C%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     C
     SUBROUTINE UNWARPED (IOPT,X,Y,Z,BX1,BY1,BZ1,BX2,BY2,BZ2)

     C   IOPT - TAIL FIELD MODE FLAG:   IOPT=0 - THE TWO TAIL MODES ARE ADDED UP
     C                                  IOPT=1 - MODE 1 ONLY
     C                                  IOPT=2 - MODE 2 ONLY
     C
     C    CALCULATES GSM COMPONENTS OF THE SHIELDED FIELD OF TWO TAIL MODES WITH UNIT
     C    AMPLITUDES,  WITHOUT ANY WARPING OR BENDING.  NONLINEAR PARAMETERS OF THE MODES
     C    ARE FORWARDED HERE VIA A COMMON BLOCK /TAIL/.
     C
     */
    static void UNWARPED (int const IOPT,double const X,double const Y,double const Z,double *BX1,double *BY1,double *BZ1,double *BX2,double *BY2,double *BZ2, struct TS_EXTERNAL_COMMON_ const* TS_EXTERNAL_COMMON) {
        //DIMENSION A1(60),A2(60)  !   TAIL SHIELDING FIELD PARAMETERS FOR THE MODES #1 & #2
        //COMMON /TAIL/ DXSHIFT1,DXSHIFT2,D0,DELTADY  ! ATTENTION:  HERE D0 & DELTADY ARE INCLUDED IN /TAIL/
        double const* DXSHIFT1 = &TS_EXTERNAL_COMMON->DXSHIFT1;
        double const* DXSHIFT2 = &TS_EXTERNAL_COMMON->DXSHIFT2;
        double const* D0 = &TS_EXTERNAL_COMMON->D;
        double const* DELTADY = &TS_EXTERNAL_COMMON->DELTADY;

        static double const DELTADX1 = 1.e0;
        static double const ALPHA1 = 1.1e0;
        static double const XSHIFT1 = 6.e0;

        static double const DELTADX2 = 0.e0;
        static double const ALPHA2 = .25e0;
        static double const XSHIFT2 = 4.e0;

        static double const A1[61] = {0./**/,-25.45869857e0,57.3589908e0,317.5501869e0,-2.626756717e0,
            -93.38053698e0,-199.6467926e0,-858.8129729e0,
            34.09192395e0,845.4214929e0,
            -29.07463068e0,47.10678547e0,
            -128.9797943e0,-781.7512093e0,6.165038619e0,
            167.8905046e0,492.068041e0,1654.724031e0,
            -46.7733792e0,-1635.922669e0,
            40.86186772e0,-.1349775602e0,-.9661991179e-1,-.1662302354e0,
            .002810467517e0,.2487355077e0,.1025565237e0,
            -14.41750229e0,-.8185333989e0,
            11.07693629e0,.7569503173e0,-9.655264745e0,
            112.2446542e0,777.5948964e0,
            -5.745008536e0,-83.03921993e0,
            -490.2278695e0,-1155.004209e0,39.0802332e0,
            1172.780574e0,-39.44349797e0,-14.07211198e0,
            -40.41201127e0,-313.2277343e0,
            2.203920979e0,8.232835341e0,197.7065115e0,
            391.2733948e0,-18.57424451e0,
            -437.2779053e0,23.04976898e0,
            11.75673963e0,13.60497313e0,4.69192706e0,
            18.20923547e0,27.59044809e0,6.677425469e0,
            1.398283308e0,2.839005878e0,
            31.24817706e0,24.53577264e0};

        static double const A2[61] = {0./**/,
            -287187.1962e0,4970.499233e0,410490.1952e0,-1347.839052e0,
            -386370.324e0,3317.98375e0,-143462.3895e0,
            5706.513767e0,171176.2904e0,
            250.888275e0,-506570.8891e0,5733.592632e0,
            397975.5842e0,9771.762168e0,
            -941834.2436e0,7990.97526e0,54313.10318e0,
            447.538806e0,528046.3449e0,
            12751.04453e0,-21920.98301e0,-21.05075617e0,
            31971.07875e0,3012.641612e0,
            -301822.9103e0,-3601.107387e0,1797.577552e0,
            -6.315855803e0,142578.8406e0,
            13161.9364e0,804184.841e0,-14168.99698e0,
            -851926.636e0,-1890.885671e0,
            972475.6869e0,-8571.862853e0,26432.49197e0,
            -2554.752298e0,-482308.3431e0,
            -4391.473324e0,105155.916e0,-1134.62205e0,
            -74353.53091e0,-5382.670711e0,
            695055.0788e0,-916.3365144e0,-12111.06667e0,
            67.20923358e0,-367200.9285e0,
            -21414.14421e0,14.75567902e0,20.7563819e0,
            59.78601609e0,16.86431444e0,
            32.58482365e0,23.69472951e0,17.24977936e0,
            13.64902647e0,68.40989058e0,
            11.67828167e0};

        static double const XM1 = -12.e0;
        static double const XM2 = -12.e0;

        if (IOPT != 2) {//IF (IOPT.EQ.2) GOTO 1

            double XSC1=(X-XSHIFT1-*DXSHIFT1)*ALPHA1-XM1*(ALPHA1-1.e0);
            double YSC1=Y*ALPHA1;
            double ZSC1=Z*ALPHA1;
            double D0SC1=*D0*ALPHA1;//   ! HERE WE USE A SINGLE VALUE D0 OF THE THICKNESS FOR BOTH MODES

            double FX1, FY1, FZ1, HX1, HY1, HZ1;
            TAILDISK(D0SC1,DELTADX1,*DELTADY,XSC1,YSC1,ZSC1,&FX1,&FY1,&FZ1);
            SHLCAR5X5(A1,X,Y,Z,*DXSHIFT1,&HX1,&HY1,&HZ1);

            *BX1=FX1+HX1;
            *BY1=FY1+HY1;
            *BZ1=FZ1+HZ1;

            if (IOPT == 1) {//IF (IOPT.EQ.1) THEN
                *BX2=0.e0;
                *BY2=0.e0;
                *BZ2=0.e0;
                return;
            }//ENDIF

        }
        double XSC2=(X-XSHIFT2-*DXSHIFT2)*ALPHA2-XM2*(ALPHA2-1.e0); // 1
        double YSC2=Y*ALPHA2;
        double ZSC2=Z*ALPHA2;
        double D0SC2=*D0*ALPHA2;//   ! HERE WE USE A SINGLE VALUE D0 OF THE THICKNESS FOR BOTH MODES

        double FX2, FY2, FZ2, HX2, HY2, HZ2;
        TAILDISK(D0SC2,DELTADX2,*DELTADY,XSC2,YSC2,ZSC2,&FX2,&FY2,&FZ2);
        SHLCAR5X5(A2,X,Y,Z,*DXSHIFT2,&HX2,&HY2,&HZ2);

        *BX2=FX2+HX2;
        *BY2=FY2+HY2;
        *BZ2=FZ2+HZ2;

        if (IOPT == 2) {//IF (IOPT.EQ.2) THEN
            *BX1=0.e0;
            *BY1=0.e0;
            *BZ1=0.e0;
            return;
        }//ENDIF

        return;
    }//END

    /*
     C------------------------------------------------------------------
     c
     C
     SUBROUTINE WARPED (IOPT,PS,X,Y,Z,BX1,BY1,BZ1,BX2,BY2,BZ2)
     C
     C   CALCULATES GSM COMPONENTS OF THE WARPED FIELD FOR TWO TAIL UNIT MODES.
     C   THE WARPING DEFORMATION IS IMPOSED ON THE UNWARPED FIELD, COMPUTED
     C   BY THE S/R "UNWARPED".  THE WARPING PARAMETERS WERE TAKEN FROM THE
     C   RESULTS OF GEOTAIL OBSERVATIONS (TSYGANENKO ET AL. [1998]).
     C   NB # 6, P.106, OCT 12, 2000.
     C
     C   IOPT - TAIL FIELD MODE FLAG:   IOPT=0 - THE TWO TAIL MODES ARE ADDED UP
     C                                  IOPT=1 - MODE 1 ONLY
     C                                  IOPT=2 - MODE 2 ONLY
     C
     */
    static void WARPED (int const IOPT,double const PS,double const X,double const Y,double const Z,double *BX1,double *BY1,double *BZ1,double *BX2,double *BY2,double *BZ2, struct TS_EXTERNAL_COMMON_ const* TS_EXTERNAL_COMMON) {
        //COMMON /G/ G
        double const* G = &TS_EXTERNAL_COMMON->G;
        double DGDX=0.e0;
        double XL=20.e0;
        double DXLDX=0.e0;

        double SPS=sin(PS);
        double RHO2=Y*Y+Z*Z;
        double RHO=sqrt(RHO2);

        double PHI, CPHI, SPHI;
        if (Y == 0. && Z == 0.) {//IF (Y.EQ.0.D0.AND.Z.EQ.0.D0) THEN
            PHI=0.e0;
            CPHI=1.e0;
            SPHI=0.e0;
        } else {//ELSE
            PHI=atan2(Z,Y);
            CPHI=Y/RHO;
            SPHI=Z/RHO;
        }//ENDIF

        double RR4L4=RHO/(RHO2*RHO2+pow(XL, 4));

        double F=PHI+*G*RHO2*RR4L4*CPHI*SPS;
        double DFDPHI=1.e0-*G*RHO2*RR4L4*SPHI*SPS;
        double DFDRHO=*G*RR4L4*RR4L4*(3.e0*pow(XL, 4)-RHO2*RHO2)*CPHI*SPS;
        double DFDX=RR4L4*CPHI*SPS*(DGDX*RHO2-*G*RHO*RR4L4*4.e0*pow(XL, 3)*DXLDX);

        double CF=cos(F);
        double SF=sin(F);
        double YAS=RHO*CF;
        double ZAS=RHO*SF;

        double BX_AS1, BX_AS2, BY_AS1, BY_AS2, BZ_AS1, BZ_AS2;
        UNWARPED (IOPT,X,YAS,ZAS,&BX_AS1,&BY_AS1,&BZ_AS1,
                  &BX_AS2,&BY_AS2,&BZ_AS2, TS_EXTERNAL_COMMON);

        double BRHO_AS =  BY_AS1*CF+BZ_AS1*SF;//      !   DEFORM THE 1ST MODE
        double BPHI_AS = -BY_AS1*SF+BZ_AS1*CF;

        double BRHO_S = BRHO_AS*DFDPHI;
        double BPHI_S = BPHI_AS-RHO*(BX_AS1*DFDX+BRHO_AS*DFDRHO);
        *BX1    = BX_AS1*DFDPHI;

        *BY1    = BRHO_S*CPHI-BPHI_S*SPHI;
        *BZ1    = BRHO_S*SPHI+BPHI_S*CPHI;//   !   DONE

        BRHO_AS =  BY_AS2*CF+BZ_AS2*SF;//      !   DEFORM THE 2ND MODE
        BPHI_AS = -BY_AS2*SF+BZ_AS2*CF;

        BRHO_S = BRHO_AS*DFDPHI;
        BPHI_S = BPHI_AS-RHO*(BX_AS2*DFDX+BRHO_AS*DFDRHO);
        *BX2    = BX_AS2*DFDPHI;

        *BY2    = BRHO_S*CPHI-BPHI_S*SPHI;
        *BZ2    = BRHO_S*SPHI+BPHI_S*CPHI;//   !   DONE

        return;
    }//END

    /*
     c############################################################################
     c
     C
     SUBROUTINE DEFORMED (IOPT,PS,X,Y,Z,BX1,BY1,BZ1,BX2,BY2,BZ2)
     C
     C   IOPT - TAIL FIELD MODE FLAG:   IOPT=0 - THE TWO TAIL MODES ARE ADDED UP
     C                                  IOPT=1 - MODE 1 ONLY
     C                                  IOPT=2 - MODE 2 ONLY
     C
     C   CALCULATES GSM COMPONENTS OF TWO UNIT-AMPLITUDE TAIL FIELD MODES,
     C    TAKING INTO ACCOUNT BOTH EFFECTS OF DIPOLE TILT:
     C    WARPING IN Y-Z (DONE BY THE S/R WARPED) AND BENDING IN X-Z (DONE BY THIS SUBROUTINE)
     C
     */
    static void DEFORMED (int const IOPT,double const PS,double const X,double const Y,double const Z,double *BX1,double *BY1,double *BZ1,double *BX2,double *BY2,double *BZ2, struct TS_EXTERNAL_COMMON_ const* TS_EXTERNAL_COMMON) {
        //COMMON /RH0/ RH0
        double const* RH0 = &TS_EXTERNAL_COMMON->RH0;

        static double const RH2 = -5.2e0;
        static int const IEPS = 3;
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
        double DFDR=-pow(RRH, (IEPS-1)) * pow(F, (IEPS+1))/RH;
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
        double BXAS1, BXAS2, BYAS1, BYAS2, BZAS1, BZAS2;
        WARPED(IOPT,PS,XAS,Y,ZAS,&BXAS1,&BYAS1,&BZAS1,&BXAS2,&BYAS2,&BZAS2, TS_EXTERNAL_COMMON);

        *BX1=BXAS1*DZASDZ-BZAS1*DXASDZ +BYAS1*FAC1;
        *BY1=BYAS1*FAC2;
        *BZ1=BZAS1*DXASDX-BXAS1*DZASDX +BYAS1*FAC3;

        *BX2=BXAS2*DZASDZ-BZAS2*DXASDZ +BYAS2*FAC1;
        *BY2=BYAS2*FAC2;
        *BZ2=BZAS2*DXASDX-BXAS2*DZASDX +BYAS2*FAC3;

        return;
    }//END

    /*
     C$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     C
     SUBROUTINE  SHLCAR3X3(X,Y,Z,PS,BX,BY,BZ)
     C
     C   THIS S/R RETURNS THE SHIELDING FIELD FOR THE EARTH'S DIPOLE,
     C   REPRESENTED BY  2x3x3=18 "CARTESIAN" HARMONICS, tilted with respect
     C   to the z=0 plane
     C
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C  The 36 coefficients enter in pairs in the amplitudes of the "cartesian"
     c    harmonics (A(1)-A(36).
     c  The 14 nonlinear parameters (A(37)-A(50) are the scales Pi,Ri,Qi,and Si
     C   entering the arguments of exponents, sines, and cosines in each of the
     C   18 "Cartesian" harmonics  PLUS TWO TILT ANGLES FOR THE CARTESIAN HARMONICS
     C       (ONE FOR THE PSI=0 MODE AND ANOTHER FOR THE PSI=90 MODE)
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C
     */
    static void SHLCAR3X3(double const X,double const Y,double const Z,double const PS,double *BX,double *BY,double *BZ) {
        //DIMENSION A(50)
        static double const A[51] = {0./**/,
            -901.2327248e0,895.8011176e0,817.6208321e0,-845.5880889e0,
            -83.73539535e0,86.58542841e0,336.8781402e0,-329.3619944e0,
            -311.294712e0,
            308.6011161e0,31.94469304e0,
            -31.30824526e0,125.8739681e0,-372.3384278e0,
            -235.4720434e0,286.7594095e0,
            21.86305585e0,-27.42344605e0,-150.4874688e0,
            2.669338538e0,1.395023949e0,-.5540427503e0,
            -56.85224007e0,3.681827033e0,
            -43.48705106e0,5.103131905e0,
            1.073551279e0,-.6673083508e0,12.21404266e0,
            4.177465543e0,5.799964188e0,
            -.3977802319e0,-1.044652977e0,.570356001e0,
            3.536082962e0,-3.222069852e0,
            9.620648151e0,6.082014949e0,27.75216226e0,
            12.44199571e0,5.122226936e0,6.982039615e0,
            20.12149582e0,6.150973118e0,
            4.663639687e0,15.73319647e0,
            2.303504968e0,5.840511214e0,.8385953499e-1,
            .3477844929e0};

        double P1=A[37];
        double P2=A[38];
        double P3=A[39];
        double R1=A[40];
        double R2=A[41];
        double R3=A[42];
        double Q1=A[43];
        double Q2=A[44];
        double Q3=A[45];
        double S1=A[46];
        double S2=A[47];
        double S3=A[48];

        double T1  =A[49];
        double T2  =A[50];

        double CPS=cos(PS);
        double SPS=sin(PS);
        double S2PS=2.e0*CPS;

        double ST1=sin(PS*T1);
        double CT1=cos(PS*T1);
        double ST2=sin(PS*T2);
        double CT2=cos(PS*T2);

        double X1=X*CT1-Z*ST1;
        double Z1=X*ST1+Z*CT1;
        double X2=X*CT2-Z*ST2;
        double Z2=X*ST2+Z*CT2;
        /*
         c  MAKE THE TERMS IN THE 1ST SUM ("PERPENDICULAR" SYMMETRY):
         C
         C       I=1
         */
        double SQPR= sqrt(1.e0/(P1*P1)+1.e0/(R1*R1));
        double CYP = cos(Y/P1);
        double SYP = sin(Y/P1);
        double CZR = cos(Z1/R1);
        double SZR = sin(Z1/R1);
        double EXPR= exp(SQPR*X1);
        double FX1 =-SQPR*EXPR*CYP*SZR;
        double HY1 = EXPR/P1*SYP*SZR;
        double FZ1 =-EXPR*CYP/R1*CZR;
        double HX1 = FX1*CT1+FZ1*ST1;
        double HZ1 =-FX1*ST1+FZ1*CT1;

        SQPR= sqrt(1.e0/(P1*P1)+1.e0/(R2*R2));
        CYP = cos(Y/P1);
        SYP = sin(Y/P1);
        CZR = cos(Z1/R2);
        SZR = sin(Z1/R2);
        EXPR= exp(SQPR*X1);
        double FX2 =-SQPR*EXPR*CYP*SZR;
        double HY2 = EXPR/P1*SYP*SZR;
        double FZ2 =-EXPR*CYP/R2*CZR;
        double HX2 = FX2*CT1+FZ2*ST1;
        double HZ2 =-FX2*ST1+FZ2*CT1;

        SQPR= sqrt(1.e0/(P1*P1)+1.e0/(R3*R3));
        CYP = cos(Y/P1);
        SYP = sin(Y/P1);
        CZR = cos(Z1/R3);
        SZR = sin(Z1/R3);
        EXPR= exp(SQPR*X1);
        double FX3 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.e0/SQPR));
        double HY3 = EXPR/P1*SYP*(Z1*CZR+X1/R3*SZR/SQPR);
        double FZ3 =-EXPR*CYP*(CZR*(1.e0+X1/(R3*R3)/SQPR)-Z1/R3*SZR);
        double HX3 = FX3*CT1+FZ3*ST1;
        double HZ3 =-FX3*ST1+FZ3*CT1;
        /*
         C       I=2:
         */
        SQPR= sqrt(1.e0/(P2*P2)+1.e0/(R1*R1));
        CYP = cos(Y/P2);
        SYP = sin(Y/P2);
        CZR = cos(Z1/R1);
        SZR = sin(Z1/R1);
        EXPR= exp(SQPR*X1);
        double FX4 =-SQPR*EXPR*CYP*SZR;
        double HY4 = EXPR/P2*SYP*SZR;
        double FZ4 =-EXPR*CYP/R1*CZR;
        double HX4 = FX4*CT1+FZ4*ST1;
        double HZ4 =-FX4*ST1+FZ4*CT1;

        SQPR= sqrt(1.e0/(P2*P2)+1.e0/(R2*R2));
        CYP = cos(Y/P2);
        SYP = sin(Y/P2);
        CZR = cos(Z1/R2);
        SZR = sin(Z1/R2);
        EXPR= exp(SQPR*X1);
        double FX5 =-SQPR*EXPR*CYP*SZR;
        double HY5 = EXPR/P2*SYP*SZR;
        double FZ5 =-EXPR*CYP/R2*CZR;
        double HX5 = FX5*CT1+FZ5*ST1;
        double HZ5 =-FX5*ST1+FZ5*CT1;

        SQPR= sqrt(1.e0/(P2*P2)+1.e0/(R3*R3));
        CYP = cos(Y/P2);
        SYP = sin(Y/P2);
        CZR = cos(Z1/R3);
        SZR = sin(Z1/R3);
        EXPR= exp(SQPR*X1);
        double FX6 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.e0/SQPR));
        double HY6 = EXPR/P2*SYP*(Z1*CZR+X1/R3*SZR/SQPR);
        double FZ6 =-EXPR*CYP*(CZR*(1.e0+X1/(R3*R3)/SQPR)-Z1/R3*SZR);
        double HX6 = FX6*CT1+FZ6*ST1;
        double HZ6 =-FX6*ST1+FZ6*CT1;
        /*
         C       I=3:
         */
        SQPR= sqrt(1.e0/(P3*P3)+1.e0/(R1*R1));
        CYP = cos(Y/P3);
        SYP = sin(Y/P3);
        CZR = cos(Z1/R1);
        SZR = sin(Z1/R1);
        EXPR= exp(SQPR*X1);
        double FX7 =-SQPR*EXPR*CYP*SZR;
        double HY7 = EXPR/P3*SYP*SZR;
        double FZ7 =-EXPR*CYP/R1*CZR;
        double HX7 = FX7*CT1+FZ7*ST1;
        double HZ7 =-FX7*ST1+FZ7*CT1;

        SQPR= sqrt(1.e0/(P3*P3)+1.e0/(R2*R2));
        CYP = cos(Y/P3);
        SYP = sin(Y/P3);
        CZR = cos(Z1/R2);
        SZR = sin(Z1/R2);
        EXPR= exp(SQPR*X1);
        double FX8 =-SQPR*EXPR*CYP*SZR;
        double HY8 = EXPR/P3*SYP*SZR;
        double FZ8 =-EXPR*CYP/R2*CZR;
        double HX8 = FX8*CT1+FZ8*ST1;
        double HZ8 =-FX8*ST1+FZ8*CT1;

        SQPR= sqrt(1.e0/(P3*P3)+1.e0/(R3*R3));
        CYP = cos(Y/P3);
        SYP = sin(Y/P3);
        CZR = cos(Z1/R3);
        SZR = sin(Z1/R3);
        EXPR= exp(SQPR*X1);
        double FX9 =-EXPR*CYP*(SQPR*Z1*CZR+SZR/R3*(X1+1.e0/SQPR));
        double HY9 = EXPR/P3*SYP*(Z1*CZR+X1/R3*SZR/SQPR);
        double FZ9 =-EXPR*CYP*(CZR*(1.e0+X1/(R3*R3)/SQPR)-Z1/R3*SZR);
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
        double SQQS= sqrt(1.e0/(Q1*Q1)+1.e0/(S1*S1));
        double CYQ = cos(Y/Q1);
        double SYQ = sin(Y/Q1);
        double CZS = cos(Z2/S1);
        double SZS = sin(Z2/S1);
        double EXQS= exp(SQQS*X2);
        FX1 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY1 = EXQS/Q1*SYQ*CZS   *SPS;
        FZ1 = EXQS*CYQ/S1*SZS   *SPS;
        HX1 = FX1*CT2+FZ1*ST2;
        HZ1 =-FX1*ST2+FZ1*CT2;

        SQQS= sqrt(1.e0/(Q1*Q1)+1.e0/(S2*S2));
        CYQ = cos(Y/Q1);
        SYQ = sin(Y/Q1);
        CZS = cos(Z2/S2);
        SZS = sin(Z2/S2);
        EXQS= exp(SQQS*X2);
        FX2 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY2 = EXQS/Q1*SYQ*CZS   *SPS;
        FZ2 = EXQS*CYQ/S2*SZS   *SPS;
        HX2 = FX2*CT2+FZ2*ST2;
        HZ2 =-FX2*ST2+FZ2*CT2;

        SQQS= sqrt(1.e0/(Q1*Q1)+1.e0/(S3*S3));
        CYQ = cos(Y/Q1);
        SYQ = sin(Y/Q1);
        CZS = cos(Z2/S3);
        SZS = sin(Z2/S3);
        EXQS= exp(SQQS*X2);
        FX3 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY3 = EXQS/Q1*SYQ*CZS   *SPS;
        FZ3 = EXQS*CYQ/S3*SZS   *SPS;
        HX3 = FX3*CT2+FZ3*ST2;
        HZ3 =-FX3*ST2+FZ3*CT2;
        /*
         C       I=2
         */
        SQQS= sqrt(1.e0/(Q2*Q2)+1.e0/(S1*S1));
        CYQ = cos(Y/Q2);
        SYQ = sin(Y/Q2);
        CZS = cos(Z2/S1);
        SZS = sin(Z2/S1);
        EXQS= exp(SQQS*X2);
        FX4 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY4 = EXQS/Q2*SYQ*CZS   *SPS;
        FZ4 = EXQS*CYQ/S1*SZS   *SPS;
        HX4 = FX4*CT2+FZ4*ST2;
        HZ4 =-FX4*ST2+FZ4*CT2;

        SQQS= sqrt(1.e0/(Q2*Q2)+1.e0/(S2*S2));
        CYQ = cos(Y/Q2);
        SYQ = sin(Y/Q2);
        CZS = cos(Z2/S2);
        SZS = sin(Z2/S2);
        EXQS= exp(SQQS*X2);
        FX5 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY5 = EXQS/Q2*SYQ*CZS   *SPS;
        FZ5 = EXQS*CYQ/S2*SZS   *SPS;
        HX5 = FX5*CT2+FZ5*ST2;
        HZ5 =-FX5*ST2+FZ5*CT2;

        SQQS= sqrt(1.e0/(Q2*Q2)+1.e0/(S3*S3));
        CYQ = cos(Y/Q2);
        SYQ = sin(Y/Q2);
        CZS = cos(Z2/S3);
        SZS = sin(Z2/S3);
        EXQS= exp(SQQS*X2);
        FX6 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY6 = EXQS/Q2*SYQ*CZS   *SPS;
        FZ6 = EXQS*CYQ/S3*SZS   *SPS;
        HX6 = FX6*CT2+FZ6*ST2;
        HZ6 =-FX6*ST2+FZ6*CT2;
        /*
         C       I=3
         */
        SQQS= sqrt(1.e0/(Q3*Q3)+1.e0/(S1*S1));
        CYQ = cos(Y/Q3);
        SYQ = sin(Y/Q3);
        CZS = cos(Z2/S1);
        SZS = sin(Z2/S1);
        EXQS= exp(SQQS*X2);
        FX7 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY7 = EXQS/Q3*SYQ*CZS   *SPS;
        FZ7 = EXQS*CYQ/S1*SZS   *SPS;
        HX7 = FX7*CT2+FZ7*ST2;
        HZ7 =-FX7*ST2+FZ7*CT2;

        SQQS= sqrt(1.e0/(Q3*Q3)+1.e0/(S2*S2));
        CYQ = cos(Y/Q3);
        SYQ = sin(Y/Q3);
        CZS = cos(Z2/S2);
        SZS = sin(Z2/S2);
        EXQS= exp(SQQS*X2);
        FX8 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY8 = EXQS/Q3*SYQ*CZS   *SPS;
        FZ8 = EXQS*CYQ/S2*SZS   *SPS;
        HX8 = FX8*CT2+FZ8*ST2;
        HZ8 =-FX8*ST2+FZ8*CT2;

        SQQS= sqrt(1.e0/(Q3*Q3)+1.e0/(S3*S3));
        CYQ = cos(Y/Q3);
        SYQ = sin(Y/Q3);
        CZS = cos(Z2/S3);
        SZS = sin(Z2/S3);
        EXQS= exp(SQQS*X2);
        FX9 =-SQQS*EXQS*CYQ*CZS *SPS;
        HY9 = EXQS/Q3*SYQ*CZS   *SPS;
        FZ9 = EXQS*CYQ/S3*SZS   *SPS;
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
     SUBROUTINE EXTERN (IOPGEN,IOPT,IOPB,IOPR,A,NTOT,
     *  PDYN,DST,BXIMF,BYIMF,BZIMF,W1,W2,W3,W4,W5,W6,PS,X,Y,Z,
     *  BXCF,BYCF,BZCF,BXT1,BYT1,BZT1,BXT2,BYT2,BZT2,
     *  BXSRC,BYSRC,BZSRC,BXPRC,BYPRC,BZPRC, BXR11,BYR11,BZR11,
     *  BXR12,BYR12,BZR12,BXR21,BYR21,BZR21,BXR22,BYR22,BZR22,HXIMF,
     *  HYIMF,HZIMF,BX,BY,BZ)
     C
     C   IOPGEN - GENERAL OPTION FLAG:  IOPGEN=0 - CALCULATE TOTAL FIELD
     C                                  IOPGEN=1 - DIPOLE SHIELDING ONLY
     C                                  IOPGEN=2 - TAIL FIELD ONLY
     C                                  IOPGEN=3 - BIRKELAND FIELD ONLY
     C                                  IOPGEN=4 - RING CURRENT FIELD ONLY
     C                                  IOPGEN=5 - INTERCONNECTION FIELD ONLY
     C
     C   IOPT -  TAIL FIELD FLAG:       IOPT=0  -  BOTH MODES
     C                                  IOPT=1  -  MODE 1 ONLY
     C                                  IOPT=2  -  MODE 2 ONLY
     C
     C   IOPB -  BIRKELAND FIELD FLAG:  IOPB=0  -  ALL 4 TERMS
     C                                  IOPB=1  -  REGION 1, MODES 1 AND 2
     C                                  IOPB=2  -  REGION 2, MODES 1 AND 2
     C
     C   IOPR -  RING CURRENT FLAG:     IOPR=0  -  BOTH SRC AND PRC
     C                                  IOPR=1  -  SRC ONLY
     C                                  IOPR=2  -  PRC ONLY
     C
     */
    static void EXTERN (int const IOPGEN,int const IOPT,int const IOPB,int const IOPR,
                        double const A[],int const NTOT,double const PDYN,double const DST,
                        double const BXIMF,double const BYIMF,double const BZIMF,
                        double const W1,double const W2,double const W3,double const W4,double const W5,double const W6,
                        double const PS,double const X,double const Y,double const Z,
                        double *BXCF,double *BYCF,double *BZCF,
                        double *BXT1,double *BYT1,double *BZT1,
                        double *BXT2,double *BYT2,double *BZT2,
                        double *BXSRC,double *BYSRC,double *BZSRC,
                        double *BXPRC,double *BYPRC,double *BZPRC,
                        double *BXR11,double *BYR11,double *BZR11,
                        double *BXR12,double *BYR12,double *BZR12,
                        double *BXR21,double *BYR21,double *BZR21,
                        double *BXR22,double *BYR22,double *BZR22,
                        double *HXIMF,double *HYIMF,double *HZIMF,
                        double *BX,double *BY,double *BZ) {

        //DIMENSION A(NTOT)
        /*
         COMMON /TAIL/ DXSHIFT1,DXSHIFT2,D,DELTADY  ! THE COMMON BLOCKS FORWARD NONLINEAR PARAMETERS
         COMMON /BIRKPAR/ XKAPPA1,XKAPPA2
         COMMON /RCPAR/ SC_SY,SC_AS,PHI
         COMMON /G/ G
         COMMON /RH0/ RH0
         */
        struct TS_EXTERNAL_COMMON_ TS_EXTERNAL_COMMON;
        double *DXSHIFT1 = &TS_EXTERNAL_COMMON.DXSHIFT1;
        double *DXSHIFT2 = &TS_EXTERNAL_COMMON.DXSHIFT2;
        double *D = &TS_EXTERNAL_COMMON.D;
        double *DELTADY = &TS_EXTERNAL_COMMON.DELTADY;
        double *XKAPPA1 = &TS_EXTERNAL_COMMON.XKAPPA1;
        double *XKAPPA2 = &TS_EXTERNAL_COMMON.XKAPPA2;
        double *SC_SY = &TS_EXTERNAL_COMMON.SC_SY;
        double *SC_AS = &TS_EXTERNAL_COMMON.SC_AS;
        double *PHI = &TS_EXTERNAL_COMMON.PHI;
        double *G = &TS_EXTERNAL_COMMON.G;
        double *RH0 = &TS_EXTERNAL_COMMON.RH0;

        static double const A0_A = 34.586e0;
        static double const A0_S0 = 1.1960e0;
        static double const A0_X0 = 3.4397e0;//   !   SHUE ET AL. PARAMETERS

        static double const DSIG = 0.005e0;
        //double const RH0 = 7.5e0; Already in COMMON block
        static double const RH2 = -5.2e0;

        double XAPPA=pow((PDYN/2e0), A[23]);//   !  OVERALL SCALING PARAMETER
        *RH0=7.5e0;//                  !  TAIL HINGING DISTANCE

        *G=  35e0;//                 !  TAIL WARPING PARAMETER

        double XAPPA3=XAPPA*XAPPA*XAPPA;

        double XX=X*XAPPA;
        double YY=Y*XAPPA;
        double ZZ=Z*XAPPA;

        double SPS=sin(PS);

        double X0=A0_X0/XAPPA;
        double AM=A0_A/XAPPA;
        double S0=A0_S0;
        /*
         C  CALCULATE "IMF" COMPONENTS OUTSIDE THE MAGNETOPAUSE LAYER (HENCE BEGIN WITH "O")
         C  THEY ARE NEEDED ONLY IF THE POINT (X,Y,Z) IS WITHIN THE TRANSITION MAGNETOPAUSE LAYER
         C  OR OUTSIDE THE MAGNETOSPHERE:
         */
        double FACTIMF=A[20];

        double OIMFX=0e0;
        double OIMFY=BYIMF*FACTIMF;
        double OIMFZ=BZIMF*FACTIMF;

        double R=sqrt(X*X+Y*Y+Z*Z);
        double XSS=X;
        double ZSS=Z;

        double DD;
        do {
            double XSOLD=XSS; // 1      !   BEGIN ITERATIVE SEARCH OF UNWARPED COORDS (TO FIND SIGMA)
            double ZSOLD=ZSS;

            double RH=*RH0+RH2 * (ZSS/R)*(ZSS/R);
            double SINPSAS=SPS/pow((1e0+(R/RH)*(R/RH)*(R/RH)), 0.33333333e0);
            double COSPSAS=sqrt(1e0-SINPSAS*SINPSAS);
            ZSS=X*SINPSAS+Z*COSPSAS;
            XSS=X*COSPSAS-Z*SINPSAS;
            DD=fabs(XSS-XSOLD)+fabs(ZSS-ZSOLD);
        } while (DD > 1e-6);//IF (DD.GT.1D-6) GOTO 1
        //C                                END OF ITERATIVE SEARCH
        double RHO2=Y*Y+ZSS*ZSS;
        double ASQ=AM*AM;
        double XMXM=AM+XSS-X0;
        if (XMXM < 0.) XMXM = 0.; //IF (XMXM.LT.0D0) XMXM=0D0 ! THE BOUNDARY IS A CYLINDER TAILWARD OF X=X0-AM
        double AXX0=XMXM*XMXM;
        double ARO=ASQ+RHO2;
        double SIGMA= sqrt((ARO+AXX0+sqrt((ARO+AXX0)*(ARO+AXX0)-4e0*ASQ*AXX0))/(2e0*ASQ));
        /*
         C   NOW, THERE ARE THREE POSSIBLE CASES:
         C    (1) INSIDE THE MAGNETOSPHERE   (SIGMA
         C    (2) IN THE BOUNDARY LAYER
         C    (3) OUTSIDE THE MAGNETOSPHERE AND B.LAYER
         C       FIRST OF ALL, CONSIDER THE CASES (1) AND (2):
         C
         C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
         */
        if (SIGMA < S0+DSIG) {//IF (SIGMA.LT.S0+DSIG) THEN  !  CASES (1) OR (2); CALCULATE THE MODEL FIELD
            /*                                  (WITH THE POTENTIAL "PENETRATED" INTERCONNECTION FIELD):
             C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             */
            if (IOPGEN <= 1) {//IF (IOPGEN.LE.1) THEN
                double CFX, CFY, CFZ;
                SHLCAR3X3(XX,YY,ZZ,PS,&CFX,&CFY,&CFZ);//         !  DIPOLE SHIELDING FIELD
                *BXCF=CFX*XAPPA3;
                *BYCF=CFY*XAPPA3;
                *BZCF=CFZ*XAPPA3;
            } else {//ELSE
                *BXCF=0e0;
                *BYCF=0e0;
                *BZCF=0e0;
            }//ENDIF                                              !  DONE
            
            if (IOPGEN == 0 || IOPGEN == 2) {//IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.2) THEN
                double DSTT=-20e0;
                if (DST < DSTT) DSTT = DST; //IF (DST.LT.DSTT) DSTT=DST
                double ZNAM=pow(fabs(DSTT), 0.37e0);
                *DXSHIFT1=A[24]-A[25]/ZNAM;
                *DXSHIFT2=A[26]-A[27]/ZNAM;
                *D=A[36]*exp(-W1/A[37])  +A[69];
                *DELTADY=4.7e0;
                
                DEFORMED (IOPT,PS,XX,YY,ZZ,//              !  TAIL FIELD (THREE MODES)
                          BXT1,BYT1,BZT1,BXT2,BYT2,BZT2, &TS_EXTERNAL_COMMON);
            } else {//ELSE
                *BXT1=0e0;
                *BYT1=0e0;
                *BZT1=0e0;
                *BXT2=0e0;
                *BYT2=0e0;
                *BZT2=0e0;
            }//ENDIF
            
            if (IOPGEN == 0 || IOPGEN == 3) {//IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.3) THEN
                double ZNAM=fabs(DST);
                if (DST >= -20.) ZNAM = 20.;//IF (DST.GE.-20D0) ZNAM=20D0
                *XKAPPA1=A[32]*pow((ZNAM/20e0), A[33]);
                *XKAPPA2=A[34]*pow((ZNAM/20e0), A[35]);
                
                BIRK_TOT (IOPB,PS,XX,YY,ZZ,BXR11,BYR11,BZR11,BXR12,BYR12,
                          BZR12,BXR21,BYR21,BZR21,BXR22,BYR22,BZR22, &TS_EXTERNAL_COMMON);//    !   BIRKELAND FIELD (TWO MODES FOR R1 AND TWO MODES FOR R2)
            } else { //ELSE
                *BXR11=0e0;
                *BYR11=0e0;
                *BZR11=0e0;
                *BXR21=0e0;
                *BYR21=0e0;
                *BZR21=0e0;
            }//ENDIF
            
            if (IOPGEN == 0 || IOPGEN == 4) {//IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.4) THEN
                *PHI=A[38];
                
                double ZNAM=fabs(DST);
                if (DST >= -20.) ZNAM = 20.; //IF (DST.GE.-20D0) ZNAM=20D0
                *SC_SY=A[28]*pow((20e0/ZNAM), A[29]) *XAPPA;
                *SC_AS=A[30]*pow((20e0/ZNAM), A[31]) *XAPPA;//    !  MULTIPLICATION  BY XAPPA IS MADE IN ORDER TO MAKE THE SRC AND PRC
                //!     SCALING COMPLETELY INDEPENDENT OF THE GENERAL SCALING DUE TO THE
                //C                                                         MAGNETOPAUSE COMPRESSION/EXPANSION                             !
                
                FULL_RC(IOPR,PS,XX,YY,ZZ,BXSRC,BYSRC,BZSRC,BXPRC,BYPRC, BZPRC, &TS_EXTERNAL_COMMON);//  !  SHIELDED RING CURRENT (SRC AND PRC)
            } else {//ELSE
                *BXSRC=0.e0;
                *BYSRC=0.e0;
                *BZSRC=0.e0;
                *BXPRC=0.e0;
                *BYPRC=0.e0;
                *BZPRC=0.e0;
            }//ENDIF
            
            if (IOPGEN == 0 || IOPGEN == 5) {//IF (IOPGEN.EQ.0.OR.IOPGEN.EQ.5) THEN
                *HXIMF=0.e0;
                *HYIMF=BYIMF;
                *HZIMF=BZIMF;//   ! THESE ARE COMPONENTS OF THE PENETRATED FIELD PER UNIT OF THE PENETRATION COEFFICIENT.
                //C                             IN OTHER WORDS, THESE ARE DERIVATIVES OF THE PENETRATION FIELD COMPONENTS WITH RESPECT
                //C                             TO THE PENETRATION COEFFICIENT.   WE ASSUME THAT ONLY TRANSVERSE COMPONENT OF THE
                //C                             FIELD PENETRATES INSIDE.
            } else {//ELSE
                *HXIMF=0.e0;
                *HYIMF=0.e0;
                *HZIMF=0.e0;
            }//ENDIF
            /*
             C-----------------------------------------------------------
             C
             C    NOW, ADD UP ALL THE COMPONENTS:
             */
            double DLP1=pow((PDYN/2.e0), A[21]);
            double DLP2=pow((PDYN/2.e0), A[22]);
            
            double TAMP1=A[2]+A[3]*DLP1+A[4]*A[39]*W1/sqrt(W1*W1+A[39]*A[39])+A[5]*DST;
            double TAMP2=A[6]+A[7]*DLP2+A[8]*A[40]*W2/sqrt(W2*W2+A[40]*A[40])+A[9]*DST;
            double A_SRC=A[10]+A[11]*A[41]*W3/sqrt(W3*W3+A[41]*A[41])
            +A[12]*DST;
            double A_PRC=A[13]+A[14]*A[42]*W4/sqrt(W4*W4+A[42]*A[42])
            +A[15]*DST;
            double A_R11=A[16]+A[17]*A[43]*W5/sqrt(W5*W5+A[43]*A[43]);
            double A_R21=A[18]+A[19]*A[44]*W6/sqrt(W6*W6+A[44]*A[44]);
            
            double BBX=A[1]**BXCF+TAMP1**BXT1+TAMP2**BXT2+A_SRC**BXSRC+A_PRC**BXPRC
            +A_R11**BXR11+A_R21**BXR21+A[20]**HXIMF;
            
            double BBY=A[1]**BYCF+TAMP1**BYT1+TAMP2**BYT2+A_SRC**BYSRC+A_PRC**BYPRC
            +A_R11**BYR11+A_R21**BYR21+A[20]**HYIMF;
            
            double BBZ=A[1]**BZCF+TAMP1**BZT1+TAMP2**BZT2+A_SRC**BZSRC+A_PRC**BZPRC
            +A_R11**BZR11+A_R21**BZR21+A[20]**HZIMF;
            /*
             C   AND WE HAVE THE TOTAL EXTERNAL FIELD.
             C
             C
             C  NOW, LET US CHECK WHETHER WE HAVE THE CASE (1). IF YES - ALL DONE:
             */
            if (SIGMA < S0-DSIG) {//IF (SIGMA.LT.S0-DSIG) THEN    !  (X,Y,Z) IS INSIDE THE MAGNETOSPHERE
                *BX=BBX;
                *BY=BBY;
                *BZ=BBZ;
            } else {//ELSE           !  THIS IS THE MOST COMPLEX CASE: WE ARE INSIDE
                //C                                             THE INTERPOLATION REGION
                double FINT=0.5e0*(1e0-(SIGMA-S0)/DSIG);
                double FEXT=0.5e0*(1e0+(SIGMA-S0)/DSIG);
                
                double QX, QY, QZ;
                DIPOLE (PS,X,Y,Z,&QX,&QY,&QZ);
                *BX=(BBX+QX)*FINT+OIMFX*FEXT -QX;
                *BY=(BBY+QY)*FINT+OIMFY*FEXT -QY;
                *BZ=(BBZ+QZ)*FINT+OIMFZ*FEXT -QZ;
                
            }//ENDIF  !   THE CASES (1) AND (2) ARE EXHAUSTED; THE ONLY REMAINING
            //C                      POSSIBILITY IS NOW THE CASE (3):
            //+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
        } else {//ELSE
            double QX, QY, QZ;
            DIPOLE (PS,X,Y,Z,&QX,&QY,&QZ);
            *BX=OIMFX-QX;
            *BY=OIMFY-QY;
            *BZ=OIMFZ-QZ;
        }//ENDIF
        
    }//END
    
    /*
     c====================================================================================
     c
     c
     SUBROUTINE T04_s (IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
     c
     c  ASSEMBLED:  MARCH 25, 2004;
     C  UPDATED:  AUGUST 2 & 31, DECEMBER 27, 2004.
     c  LATEST MODIFICATIONS/BUGS REMOVED:
     c
     C  (1) MARCH 14, 2005:  79 -> 69  (LINE 94; might cause compilation problems with some Fortran compilers)
     c
     C  (2) JUNE 24, 2006:  REPLACED COEFFICIENTS IN
     c       (i)   DATA statement in FUNCTION AP,
     C       (ii)  DATA C_SY statement in SUBROUTINE FULL_RC, and
     c       (iii) DATA A statement in SUBROUTINE T04_s.
     C  This correction was needed due to a bug found in the symmetric ring current module.
     c   Its impact can be significant (up to ~20 nT) only in the innermost magnetosphere (R<=2)
     c      and only for strongly disturbed conditions; otherwise, the change in the model field
     c      does not exceed a few percent.
     c
     c--------------------------------------------------------------------
     C   A DATA-BASED MODEL OF THE EXTERNAL (I.E., WITHOUT EARTH'S CONTRIBUTION) PART OF THE
     C   MAGNETOSPHERIC MAGNETIC FIELD, CALIBRATED BY
     C    (1) SOLAR WIND PRESSURE PDYN (NANOPASCALS),
     C    (2) DST (NANOTESLA),
     C    (3) BYIMF,
     C    (4) BZIMF (NANOTESLA)
     C    (5-10)   INDICES W1 - W6, CALCULATED AS TIME INTEGRALS FROM THE BEGINNING OF A STORM
     c               SEE THE REFERENCE (3) BELOW, FOR A DETAILED DEFINITION OF THOSE VARIABLES
     C
     c   THE ABOVE 10 INPUT PARAMETERS SHOULD BE PLACED IN THE ELEMENTS
     c   OF THE ARRAY PARMOD(10).
     C
     C   THE REST OF THE INPUT VARIABLES ARE: THE GEODIPOLE TILT ANGLE PS (RADIANS),
     C        X,Y,Z -  GSM POSITION (RE)
     C
     c   IOPT IS A DUMMY INPUT PARAMETER, INCLUDED TO MAKE THIS SUBROUTINE
     C   COMPATIBLE WITH THE TRACING SOFTWARE PACKAGE (GEOPACK). IN THIS MODEL,
     C   THE PARAMETER IOPT DOES NOT AFFECT THE OUTPUT FIELD.
     c
     C*******************************************************************************************
     c** ATTENTION:  THE MODEL IS BASED ON DATA TAKEN SUNWARD FROM X=-15Re, AND HENCE BECOMES   *
     C**              INVALID AT LARGER TAILWARD DISTANCES !!!                                  *
     C*******************************************************************************************
     C
     c   OUTPUT:  GSM COMPONENTS OF THE EXTERNAL MAGNETIC FIELD (BX,BY,BZ, nanotesla)
     C            COMPUTED AS A SUM OF CONTRIBUTIONS FROM PRINCIPAL FIELD SOURCES
     C
     c  (C) Copr. 2004, Nikolai A. Tsyganenko, USRA/Code 612.3, NASA GSFC
     c      Greenbelt, MD 20771, USA
     c
     C                            REFERENCES:
     C
     C  (1)   N. A. Tsyganenko, A new data-based model of the near magnetosphere magnetic field:
     c       1. Mathematical structure.
     c       2. Parameterization and fitting to observations.  JGR v. 107(A8), 1176/1179, doi:10.1029/2001JA000219/220, 2002.
     c
     c  (2)  N. A. Tsyganenko, H. J. Singer, J. C. Kasper, Storm-time distortion of the
     c           inner magnetosphere: How severe can it get ?  JGR v. 108(A5), 1209, doi:10.1029/2002JA009808, 2003.
     
     c   (3)  N. A. Tsyganenko and M. I. Sitnov, Modeling the dynamics of the inner magnetosphere during
     c         strong geomagnetic storms, J. Geophys. Res., v. 110 (A3), A03208, doi: 10.1029/2004JA010798, 2005.
     c----------------------------------------------------------------------
     c
     */
    static void T04_s (int const IOPT,double const PARMOD[],double const PS,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {
        /*A(69),*/
        double PDYN = 0.;
        double DST_AST = 0.;
        double BXIMF = 0.;
        double BYIMF = 0.;
        double BZIMF = 0.;
        double W1 = 0.;
        double W2 = 0.;
        double W3 = 0.;
        double W4 = 0.;
        double W5 = 0.;
        double W6 = 0.;
        double PSS = 0.;
        double XX = 0.;
        double YY = 0.;
        double ZZ = 0.;
        double BXCF = 0.;
        double BYCF = 0.;
        double BZCF = 0.;
        double BXT1 = 0.;
        double BYT1 = 0.;
        double BZT1 = 0.;
        double BXT2 = 0.;
        double BYT2 = 0.;
        double BZT2 = 0.;
        double BXSRC = 0.;
        double BYSRC = 0.;
        double BZSRC = 0.;
        double BXPRC = 0.;
        double BYPRC = 0.;
        double BZPRC = 0.;
        double BXR11 = 0.;
        double BYR11 = 0.;
        double BZR11 = 0.;
        double BXR12 = 0.;
        double BYR12 = 0.;
        double BZR12 = 0.;
        double BXR21 = 0.;
        double BYR21 = 0.;
        double BZR21 = 0.;
        double BXR22 = 0.;
        double BYR22 = 0.;
        double BZR22 = 0.;
        double HXIMF = 0.;
        double HYIMF = 0.;
        double HZIMF = 0.;
        double BBX = 0.;
        double BBY = 0.;
        double BBZ = 0.;
        
        static double const A[70] = {0./**/,1e0,5.44118e0,0.891995e0,9.09684e0,0e0,-7.18972e0,12.27e0,
            -4.89408e0,0e0,0.870536e0,1.36081e0,0e0,0.68865e0,0.60233e0,
            0e0,0.316346e0,1.22728e0,-0.363620e-1,-0.405821e0,0.452536e0,
            0.755831e0,0.215662e0,0.152759e0,5.96235e0,
            23.2036e0,11.2994e0,69.9596e0,
            0.989596e0,-0.132131e-1,0.985681e0,
            0.344212e-1,1.02389e0,0.207867e0,
            1.5122e0,0.682715e-1,1.84714e0,
            1.76977e0,1.3769e0,0.69635e0,0.34328e0,
            3.28846e0,111.293e0,5.82287e0,4.39664e0,
            0.383403e0,0.648176e0,0.318752e-1,
            0.581168e0,1.1507e0,0.843004e0,
            0.394732e0,0.846509e0,0.916555e0,0.55092e0,
            0.180725e0,0.898772e0,0.387365e0,
            2.26596e0,1.29123e0,0.436819e0,1.28211e0,
            1.33199e0,.405553e0,1.6229e0,.699074e0,
            1.26131e0,2.42297e0,.537116e0,.619441e0};
        
        static int const IOPGEN = 0;
        static int const IOPTT = 0;
        static int const IOPB = 0;
        static int const IOPR = 0;
        
        PDYN=PARMOD[0]; // Start from 0!
        DST_AST=PARMOD[1]*0.8e0-13e0*sqrt(PDYN);
        BYIMF=PARMOD[2];
        BZIMF=PARMOD[3];
        
        W1=PARMOD[4];
        W2=PARMOD[5];
        W3=PARMOD[6];
        W4=PARMOD[7];
        W5=PARMOD[8];
        W6=PARMOD[9];
        
        PSS=PS;
        XX=X;
        YY=Y;
        ZZ=Z;
        
        EXTERN (IOPGEN,IOPTT,IOPB,IOPR,A,69+1,PDYN,DST_AST,BXIMF,BYIMF,
                BZIMF,W1,W2,W3,W4,W5,W6,PSS,XX,YY,ZZ,&BXCF,&BYCF,&BZCF,&BXT1,&BYT1,
                &BZT1,&BXT2,&BYT2,&BZT2,&BXSRC,&BYSRC,&BZSRC,&BXPRC,&BYPRC,&BZPRC, &BXR11,
                &BYR11,&BZR11,&BXR12,&BYR12,&BZR12,&BXR21,&BYR21,&BZR21,&BXR22,&BYR22,
                &BZR22,&HXIMF,&HYIMF,&HZIMF,&BBX,&BBY,&BBZ);
        
        *BX=BBX;
        *BY=BBY;
        *BZ=BBZ;
        
        return;
    }//END

}