//
//  T96.cpp
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

#include "T96.h"
#include "Geopack.h"
#include <cassert>
#include <cmath>
#include <cstring>
#include <stdexcept>

namespace UBK {
    using namespace std;

    T96::T96 (Geopack const* geopack, double const parmod[]) : TSExternalField(geopack, 1, parmod)
    {
        if (NULL == parmod) {
            throw invalid_argument("PARMOD must be 10-element array.");
        }
    }

    static void T96_01 (int const IOPT,double const PARMOD[],double const PS, double const X, double const Y, double const Z, double *BX,double *BY,double *BZ);
    void T96::getFieldInGSW_atPoint(Point *bOut, const Point ptgsw) const
    {
        T96_01(this->iopt(), this->parmod(), this->geopack()->psi(), ptgsw.x, ptgsw.y, ptgsw.z, &bOut->x, &bOut->y, &bOut->z);
    }

    /*-----------------------------Model---------------------*/
    /* Common blocks
     */
    struct WARP_COMMON {
        double CPSS,SPSS,DPSRR,RPS,WARP,D,XS,ZS,DXSX,DXSY,DXSZ,DZSX,DZSY,DZSZ,DZETAS,DDZETADX,DDZETADY,DDZETADZ,ZSWW;
    };

    struct COORD11_COMMON {
        double XX1[13],YY1[13];
    };

    struct RHDR_COMMON {
        double RH,DR;
    };

    struct LOOPDIP1_COMMON {
        double TILT,XCENTRE[3],RADIUS[3], DIPX,DIPY;
    };

    struct COORD21_COMMON {
        double XX2[15],YY2[15],ZZ2[15];
    };

    struct DX1_COMMON {
        double DX,SCALEIN,SCALEOUT;
    };

    /*-------------------------------------------------------------------
     c
     DOUBLE PRECISION FUNCTION BES0(X)
     */
    static double BES0(double const X)
    {
        double BES0 = 0.;
        //IF (DABS(X).LT.3.D0) THEN
        if (fabs(X)<3.e0) {
            double X32=pow((X/3.e0),2);
            BES0=1.e0-X32*(2.2499997e0-X32*(1.2656208e0-X32* (0.3163866e0-X32*(0.0444479e0-X32*(0.0039444e0 -X32*0.00021e0)))));
        } else {//ELSE
            double XD3=3.e0/X;
            double F0=0.79788456e0-XD3*(0.00000077e0+XD3*(0.00552740e0+XD3* (0.00009512e0-XD3*(0.00137237e0-XD3*(0.00072805e0 -XD3*0.00014476e0)))));
            double T0=X-0.78539816e0-XD3*(0.04166397e0+XD3*(0.00003954e0-XD3* (0.00262573e0-XD3*(0.00054125e0+XD3*(0.00029333e0 -XD3*0.00013558e0)))));
            BES0=F0/sqrt(X)*cos(T0);
        }//ENDIF
        return BES0;
    }//END

    /*--------------------------------------------------------------------------
     c
     DOUBLE PRECISION FUNCTION BES1(X)
     */
    static double BES1(double const X)
    {
        double BES1;
        //IF (DABS(X).LT.3.D0) THEN
        if (fabs(X)<3.e0) {
            double X32=pow((X/3.e0),2);
            double BES1XM1=0.5e0-X32*(0.56249985e0-X32*(0.21093573e0-X32* (0.03954289e0-X32*(0.00443319e0-X32*(0.00031761e0 -X32*0.00001109e0)))));
            BES1=BES1XM1*X;
        } else { //ELSE
            double XD3=3.e0/X;
            double F1=0.79788456e0+XD3*(0.00000156e0+XD3*(0.01659667e0+XD3* (0.00017105e0-XD3*(0.00249511e0-XD3*(0.00113653e0 -XD3*0.00020033e0)))));
            double T1=X-2.35619449e0+XD3*(0.12499612e0+XD3*(0.0000565e0-XD3* (0.00637879e0-XD3*(0.00074348e0+XD3*(0.00079824e0 -XD3*0.00029166e0)))));
            BES1=F1/sqrt(X)*cos(T1);
        }//ENDIF
        return BES1;
    }//END

    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     C
     DOUBLE PRECISION FUNCTION BES(X,K)
     */
    static double BES(double const X, int K)
    {
        //IF (K.EQ.0) THEN
        if (K==0) {
            return BES0(X);
        }//ENDIF

        //IF (K.EQ.1) THEN
        if (K==1) {
            return BES1(X);
        }//ENDIF

        //IF (X.EQ.0.D0) THEN
        if (X==0.e0) {
            return 0.e0;
        }//ENDIF

        double G=2.e0/X;
        //IF (X.LE.DFLOAT(K)) GOTO 10
        if (X>K) {
            int N=1;
            double XJN=BES1(X);
            double XJNM1=BES0(X);

            while (1) {
                double XJNP1=G*N*XJN-XJNM1; //       1
                N=N+1;
                //IF (N.LT.K) GOTO 2
                if (N>=K) {
                    return XJNP1;
                }//RETURN

                XJNM1=XJN; //         2
                XJN=XJNP1;
            }//GOTO 1
        }

        int N=24; //      10
        double XJN=1.e0;
        double XJNP1=0.e0;
        double SUM=0.e0;
        double BES = 0.;

        //3   IF (MOD(N,2).EQ.0) SUM=SUM+XJN
        do {
            if (N%2 == 0) SUM=SUM+XJN;
            double XJNM1=G*N*XJN-XJNP1;
            N=N-1;

            XJNP1=XJN;
            XJN=XJNM1;
            //IF (N.EQ.K) BES=XJN
            if (N==K) BES=XJN;

            //IF (DABS(XJN).GT.1.D5) THEN
            if (fabs(XJN)>1.e5) {
                XJNP1=XJNP1*1.e-5;
                XJN=XJN*1.e-5;
                SUM=SUM*1.e-5;
                //IF (N.LE.K) BES=BES*1.D-5
                if (N<=K) BES=BES*1.e-5;
            }//ENDIF

            //IF (N.EQ.0) GOTO 4
            //GOTO 3
        } while (N);

        SUM=XJN+2.e0*SUM; //      4
        BES=BES/SUM;
        return BES;
    }//END

    /*------------------------------------------------------------
     C
     SUBROUTINE INTERCON(X,Y,Z,BX,BY,BZ)
     C
     C      Calculates the potential interconnection field inside the magnetosphere,
     c  corresponding to  DELTA_X = 20Re and DELTA_Y = 10Re (NB#3, p.90, 6/6/1996).
     C  The position (X,Y,Z) and field components BX,BY,BZ are given in the rotated
     c   coordinate system, in which the Z-axis is always directed along the BzIMF
     c   (i.e. rotated by the IMF clock angle Theta)
     C   It is also assumed that the IMF Bt=1, so that the components should be
     c     (i) multiplied by the actual Bt, and
     c     (ii) transformed to standard GSM coords by rotating back around X axis
     c              by the angle -Theta.
     c
     C      Description of parameters:
     C
     C     X,Y,Z -   GSM POSITION
     C      BX,BY,BZ - INTERCONNECTION FIELD COMPONENTS INSIDE THE MAGNETOSPHERE
     C        OF A STANDARD SIZE (TO TAKE INTO ACCOUNT EFFECTS OF PRESSURE CHANGES,
     C         APPLY THE SCALING TRANSFORMATION)
     C
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C
     C     The 9 linear parameters are amplitudes of the "cartesian" harmonics
     c     The 6 nonlinear parameters are the scales Pi and Ri entering
     c    the arguments of exponents, sines, and cosines in the 9 "Cartesian"
     c       harmonics (3+3)
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     */
    static void INTERCON(double const X,double const Y,double const Z,double *BX,double *BY,double *BZ)
    {
        //DIMENSION A(15),RP(3),RR(3),P(3),R(3)
        //SAVE
        double RP[4], RR[4], P[4], R[4];

        //DATA A/-8.411078731,5932254.951,-9073284.93,-11.68794634,
        //* 6027598.824,-9218378.368,-6.508798398,-11824.42793,18015.66212,
        //* 7.99754043,13.9669886,90.24475036,16.75728834,1015.645781,
        //* 1553.493216/
        static double const A[16] = {0./*dummy*/, -8.411078731,5932254.951,-9073284.93,-11.68794634,
            6027598.824,-9218378.368,-6.508798398,-11824.42793,18015.66212,
            7.99754043,13.9669886,90.24475036,16.75728834,1015.645781,
            1553.493216};

        //DATA M/0/

        //IF (M.NE.0) GOTO 111
        //M=1

        P[1]=A[10];
        P[2]=A[11];
        P[3]=A[12];
        R[1]=A[13];
        R[2]=A[14];
        R[3]=A[15];

        //DO 11 I=1,3
        for (int I=1; I<=3; I++) {
            RP[I]=1.e0/P[I];
            RR[I]=1.e0/R[I]; //         11
        }

        //111   CONTINUE

        int L=0;

        *BX=0.;
        *BY=0.;
        *BZ=0.;
        double HX, HY, HZ;

        //        "PERPENDICULAR" KIND OF SYMMETRY ONLY

        //DO 2 I=1,3
        for(int I=1; I<=3; I++) {
            double CYPI=cos(Y*RP[I]);
            double SYPI=sin(Y*RP[I]);

            //DO 2 K=1,3
            for(int K=1; K<=3; K++) {
                double SZRK=sin(Z*RR[K]);
                double CZRK=cos(Z*RR[K]);
                double SQPR=sqrt(RP[I]*RP[I]+RR[K]*RR[K]);
                double EPR=exp(X*SQPR);

                HX=-SQPR*EPR*CYPI*SZRK;
                HY=RP[I]*EPR*SYPI*SZRK;
                HZ=-RR[K]*EPR*CYPI*CZRK;
                L=L+1;

                *BX=*BX+A[L]*HX;
                *BY=*BY+A[L]*HY;
                *BZ=*BZ+A[L]*HZ;
            }
        }//2   CONTINUE

    }//END

    /*********************************************************************
     C
     SUBROUTINE RINGCURR96(X,Y,Z,BX,BY,BZ)
     c
     c       THIS SUBROUTINE COMPUTES THE COMPONENTS OF THE RING CURRENT FIELD,
     C        SIMILAR TO THAT DESCRIBED BY TSYGANENKO AND PEREDO (1994).  THE
     C          DIFFERENCE IS THAT NOW WE USE SPACEWARPING, AS DESCRIBED IN THE
     C           PAPER ON MODELING BIRKELAND CURRENTS (TSYGANENKO AND STERN, 1996),
     C            INSTEAD OF SHEARING IT IN THE SPIRIT OF THE T89 TAIL MODEL.
     C
     C          IN  ADDITION, INSTEAD OF 7 TERMS FOR THE RING CURRENT MODEL, WE USE
     C             NOW ONLY 2 TERMS;  THIS SIMPLIFICATION ALSO GIVES RISE TO AN
     C                EASTWARD RING CURRENT LOCATED EARTHWARD FROM THE MAIN ONE,
     C                  IN LINE WITH WHAT IS ACTUALLY OBSERVED
     C
     C             FOR DETAILS, SEE NB #3, PAGES 70-73
     */
    static void RINGCURR96(double const X,double const Y,double const Z,double *BX,double *BY,double *BZ, struct WARP_COMMON const* WARP_C)
    {
        //DIMENSION F(2),BETA(2)
        //COMMON /WARP/ CPSS,SPSS,DPSRR, XNEXT(3),XS,ZSWARPED,DXSX,DXSY,
        //*   DXSZ,DZSX,DZSYWARPED,DZSZ,OTHER(4),ZS  !  ZS HERE IS WITHOUT Y-Z WARP
        double const* CPSS = &WARP_C->CPSS;
        double const* SPSS = &WARP_C->SPSS;
        double const* DPSRR = &WARP_C->DPSRR;
        double const* XS = &WARP_C->XS;
        //double &ZSWARPED = WARP_C.ZS;//ZSWARPED;
        double const* DXSX = &WARP_C->DXSX;
        double const* DXSY = &WARP_C->DXSY;
        double const* DXSZ = &WARP_C->DXSZ;
        double const* DZSX = &WARP_C->DZSX;
        //double &DZSYWARPED = WARP_C.DZSY;//DZSYWARPED;
        double const* DZSZ = &WARP_C->DZSZ;
        double const* ZS = &WARP_C->ZSWW;

        //DATA D0,DELTADX,XD,XLDX /2.,0.,0.,4./  !  ACHTUNG !!  THE RC IS NOW
        //                                           COMPLETELY SYMMETRIC (DELTADX=0)
        static double const D0 = 2.,DELTADX = 0.,XD = 0.,XLDX = 4.;

        //DATA F,BETA /569.895366D0,-1603.386993D0,2.722188D0,3.766875D0/
        static double const F[3] = {0./*dummy*/, 569.895366e0,-1603.386993e0},BETA[3] = {0./*dummy*/, 2.722188e0,3.766875e0};
        /*
         C  THE ORIGINAL VALUES OF F(I) WERE MULTIPLIED BY BETA(I) (TO REDUCE THE
         C     NUMBER OF MULTIPLICATIONS BELOW)  AND BY THE FACTOR -0.43, NORMALIZING
         C      THE DISTURBANCE AT ORIGIN  TO  B=-1nT
         */
        double DZSY=*XS*Y**DPSRR;  // NO WARPING IN THE Y-Z PLANE (ALONG X ONLY), AND
        //                         THIS IS WHY WE DO NOT USE  DZSY FROM THE COMMON-BLOCK
        double XXD=X-XD;
        double FDX=0.5e0*(1.e0+XXD/sqrt(XXD*XXD+XLDX*XLDX));
        double DDDX=DELTADX*0.5e0*XLDX*XLDX/pow(sqrt(XXD*XXD+XLDX*XLDX), 3);
        double D=D0+DELTADX*FDX;

        double DZETAS=sqrt(*ZS**ZS+D*D);  //  THIS IS THE SAME SIMPLE WAY TO SPREAD
        //                                        OUT THE SHEET, AS THAT USED IN T89
        double RHOS=sqrt(*XS**XS+Y*Y);
        double DDZETADX=(*ZS**DZSX+D*DDDX)/DZETAS;
        double DDZETADY=*ZS*DZSY/DZETAS;
        double DDZETADZ=*ZS**DZSZ/DZETAS;

        double DRHOSDX, DRHOSDY, DRHOSDZ;
        //IF (RHOS.LT.1.D-5) THEN
        if (RHOS<1.e-5) {
            DRHOSDX=0.e0;
            DRHOSDY=(Y>=0. ? 1. : -1.);//DRHOSDY=DSIGN(1.D0,Y)
            DRHOSDZ=0.e0;
        } else {//ELSE
            DRHOSDX=*XS**DXSX/RHOS;
            DRHOSDY=(*XS**DXSY+Y)/RHOS;
            DRHOSDZ=*XS**DXSZ/RHOS;
        }//ENDIF

        *BX=0.e0;
        *BY=0.e0;
        *BZ=0.e0;

        //DO 1 I=1,2
        for (int I=1; I<=2; I++) {

            double BI=BETA[I];

            double S1=sqrt(pow((DZETAS+BI), 2)+pow((RHOS+BI), 2));
            double S2=sqrt(pow((DZETAS+BI), 2)+pow((RHOS-BI), 2));
            double DS1DDZ=(DZETAS+BI)/S1;
            double DS2DDZ=(DZETAS+BI)/S2;
            double DS1DRHOS=(RHOS+BI)/S1;
            double DS2DRHOS=(RHOS-BI)/S2;

            double DS1DX=DS1DDZ*DDZETADX+DS1DRHOS*DRHOSDX;
            double DS1DY=DS1DDZ*DDZETADY+DS1DRHOS*DRHOSDY;
            double DS1DZ=DS1DDZ*DDZETADZ+DS1DRHOS*DRHOSDZ;

            double DS2DX=DS2DDZ*DDZETADX+DS2DRHOS*DRHOSDX;
            double DS2DY=DS2DDZ*DDZETADY+DS2DRHOS*DRHOSDY;
            double DS2DZ=DS2DDZ*DDZETADZ+DS2DRHOS*DRHOSDZ;

            double S1TS2=S1*S2;
            double S1PS2=S1+S2;
            double S1PS2SQ=S1PS2*S1PS2;
            double FAC1=sqrt(S1PS2SQ-pow((2.e0*BI), 2));
            double AS=FAC1/(S1TS2*S1PS2SQ);
            double TERM1=1.e0/(S1TS2*S1PS2*FAC1);
            double FAC2=AS/S1PS2SQ;
            double DASDS1=TERM1-FAC2/S1*(S2*S2+S1*(3.e0*S1+4.e0*S2));
            double DASDS2=TERM1-FAC2/S2*(S1*S1+S2*(3.e0*S2+4.e0*S1));

            double DASDX=DASDS1*DS1DX+DASDS2*DS2DX;
            double DASDY=DASDS1*DS1DY+DASDS2*DS2DY;
            double DASDZ=DASDS1*DS1DZ+DASDS2*DS2DZ;

            *BX=*BX+F[I]*((2.e0*AS+Y*DASDY)**SPSS-*XS*DASDZ +AS**DPSRR*(Y*Y**CPSS+Z**ZS));
            *BY=*BY-F[I]*Y*(AS**DPSRR**XS+DASDZ**CPSS+DASDX**SPSS);
            *BZ=*BZ+F[I]*((2.e0*AS+Y*DASDY)**CPSS+*XS*DASDX -AS**DPSRR*(X**ZS+Y*Y**SPSS)); //        1
        }

    }//END

    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     C
     C THIS CODE RETURNS THE SHIELDING FIELD REPRESENTED BY  2x3x3=18 "CARTESIAN"
     C    HARMONICS
     C
     SUBROUTINE  SHLCAR3X3(A,X,Y,Z,SPS,HX,HY,HZ)
     C
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C  The 36 coefficients enter in pairs in the amplitudes of the "cartesian"
     c    harmonics (A(1)-A(36).
     c  The 12 nonlinear parameters (A(37)-A(48) are the scales Pi,Ri,Qi,and Si
     C   entering the arguments of exponents, sines, and cosines in each of the
     C   18 "Cartesian" harmonics
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     */
    static void SHLCAR3X3(double const A[],double const X,double const Y,double const Z,double const SPS,double *HX,double *HY,double *HZ)
    {
        //DIMENSION A(48)

        double CPS=sqrt(1.e0-SPS*SPS);
        double S3PS=4.e0*CPS*CPS-1.e0;   //  THIS IS SIN(3*PS)/SIN(PS)

        *HX=0.e0;
        *HY=0.e0;
        *HZ=0.e0;
        int L=0;

        //DO 1 M=1,2     !    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
        //                           AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
        for (int M=1; M<=2; M++) {
            //DO 2 I=1,3
            for (int I=1; I<=3; I++) {
                double P=A[36+I];
                double Q=A[42+I];
                double CYPI=cos(Y/P);
                double CYQI=cos(Y/Q);
                double SYPI=sin(Y/P);
                double SYQI=sin(Y/Q);

                //DO 3 K=1,3
                for (int K=1; K<=3; K++) {
                    double R=A[39+K];
                    double S=A[45+K];
                    double SZRK=sin(Z/R);
                    double CZSK=cos(Z/S);
                    double CZRK=cos(Z/R);
                    double SZSK=sin(Z/S);
                    double SQPR=sqrt(1.e0/(P*P)+1.e0/(R*R));
                    double SQQS=sqrt(1.e0/(Q*Q)+1.e0/(S*S));
                    double EPR=exp(X*SQPR);
                    double EQS=exp(X*SQQS);

                    //DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
                    //                                  AND N=2 IS FOR THE SECOND ONE
                    double DX=0, DY=0, DZ=0;
                    for (int N=1; N<=2; N++) {

                        L=L+1;
                        //IF (M.EQ.1) THEN
                        if (M==1) {
                            //IF (N.EQ.1) THEN
                            if (N==1) {
                                DX=-SQPR*EPR*CYPI*SZRK;
                                DY=EPR/P*SYPI*SZRK;
                                DZ=-EPR/R*CYPI*CZRK;
                                *HX=*HX+A[L]*DX;
                                *HY=*HY+A[L]*DY;
                                *HZ=*HZ+A[L]*DZ;
                            } else {//ELSE
                                DX=DX*CPS;
                                DY=DY*CPS;
                                DZ=DZ*CPS;
                                *HX=*HX+A[L]*DX;
                                *HY=*HY+A[L]*DY;
                                *HZ=*HZ+A[L]*DZ;
                            }//ENDIF
                        } else {//ELSE
                            //IF (N.EQ.1) THEN
                            if (N==1) {
                                DX=-SPS*SQQS*EQS*CYQI*CZSK;
                                DY=SPS*EQS/Q*SYQI*CZSK;
                                DZ=SPS*EQS/S*CYQI*SZSK;
                                *HX=*HX+A[L]*DX;
                                *HY=*HY+A[L]*DY;
                                *HZ=*HZ+A[L]*DZ;
                            } else {//ELSE
                                DX=DX*S3PS;
                                DY=DY*S3PS;
                                DZ=DZ*S3PS;
                                *HX=*HX+A[L]*DX;
                                *HY=*HY+A[L]*DY;
                                *HZ=*HZ+A[L]*DZ;
                            }//ENDIF
                        }//ENDIF

                        //4   CONTINUE
                        //3   CONTINUE
                        //2   CONTINUE
                        //1   CONTINUE
                    }
                }
            }
        }

    }//END

    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     C
     SUBROUTINE TAILDISK(X,Y,Z,BX,BY,BZ)
     C
     c
     c       THIS SUBROUTINE COMPUTES THE COMPONENTS OF THE TAIL CURRENT FIELD,
     C        SIMILAR TO THAT DESCRIBED BY TSYGANENKO AND PEREDO (1994).  THE
     C          DIFFERENCE IS THAT NOW WE USE SPACEWARPING, AS DESCRIBED IN OUR
     C           PAPER ON MODELING BIRKELAND CURRENTS (TSYGANENKO AND STERN, 1996)
     C            INSTEAD OF SHEARING IT IN THE SPIRIT OF T89 TAIL MODEL.
     C
     C          IN  ADDITION, INSTEAD OF 8 TERMS FOR THE TAIL CURRENT MODEL, WE USE
     C           NOW ONLY 4 TERMS
     C
     C             FOR DETAILS, SEE NB #3, PAGES 74-
     */
    static void TAILDISK(double const X,double const Y,double const Z,double *BX,double *BY,double *BZ, struct WARP_COMMON const* WARP_C)
    {
        //DIMENSION F(4),BETA(4)
        //COMMON /WARP/ CPSS,SPSS,DPSRR,XNEXT(3),XS,ZS,DXSX,DXSY,DXSZ,
        //*    OTHER(3),DZETAS,DDZETADX,DDZETADY,DDZETADZ,ZSWW
        double const* CPSS = &WARP_C->CPSS;
        double const* SPSS = &WARP_C->SPSS;
        double const* DPSRR = &WARP_C->DPSRR;
        double const* XS = &WARP_C->XS;
        //double &ZS = WARP_C.ZS;
        double const* DXSX = &WARP_C->DXSX;
        double const* DXSY = &WARP_C->DXSY;
        double const* DXSZ = &WARP_C->DXSZ;
        double const* DZETAS = &WARP_C->DZETAS;
        double const* DDZETADX = &WARP_C->DDZETADX;
        double const* DDZETADY = &WARP_C->DDZETADY;
        double const* DDZETADZ = &WARP_C->DDZETADZ;
        double const* ZSWW = &WARP_C->ZSWW;

        //DATA XSHIFT /4.5/
        static double const XSHIFT = 4.5;

        //DATA F,BETA
        //* / -745796.7338D0,1176470.141D0,-444610.529D0,-57508.01028D0,
        //*   7.9250000D0,8.0850000D0,8.4712500D0,27.89500D0/
        static double const F[5] = {0./*dummy*/, -745796.7338e0,1176470.141e0,-444610.529e0,-57508.01028e0};
        static double const BETA[5] = {0./*dummy*/, 7.9250000e0,8.0850000e0,8.4712500e0,27.89500e0};
        /*
         c  here original F(I) are multiplied by BETA(I), to economize
         c    calculations
         */
        double RHOS=sqrt(pow((*XS-XSHIFT), 2)+Y*Y);
        double DRHOSDX, DRHOSDY, DRHOSDZ;
        //IF (RHOS.LT.1.D-5) THEN
        if (RHOS<1.e-5) {
            DRHOSDX=0.e0;
            DRHOSDY= (Y>=0. ? 1. : -1.);//DRHOSDY=DSIGN(1.D0,Y)
            DRHOSDZ=0.e0;
        } else {//ELSE
            DRHOSDX=(*XS-XSHIFT)**DXSX/RHOS;
            DRHOSDY=((*XS-XSHIFT)**DXSY+Y)/RHOS;
            DRHOSDZ=(*XS-XSHIFT)**DXSZ/RHOS;
        }//ENDIF

        *BX=0.e0;
        *BY=0.e0;
        *BZ=0.e0;

        //DO 1 I=1,4
        for (int I=1; I<=4; I++) {
            double BI=BETA[I];

            double S1=sqrt(pow((*DZETAS+BI), 2)+pow((RHOS+BI), 2));
            double S2=sqrt(pow((*DZETAS+BI), 2)+pow((RHOS-BI), 2));
            double DS1DDZ=(*DZETAS+BI)/S1;
            double DS2DDZ=(*DZETAS+BI)/S2;
            double DS1DRHOS=(RHOS+BI)/S1;
            double DS2DRHOS=(RHOS-BI)/S2;

            double DS1DX=DS1DDZ**DDZETADX+DS1DRHOS*DRHOSDX;
            double DS1DY=DS1DDZ**DDZETADY+DS1DRHOS*DRHOSDY;
            double DS1DZ=DS1DDZ**DDZETADZ+DS1DRHOS*DRHOSDZ;

            double DS2DX=DS2DDZ**DDZETADX+DS2DRHOS*DRHOSDX;
            double DS2DY=DS2DDZ**DDZETADY+DS2DRHOS*DRHOSDY;
            double DS2DZ=DS2DDZ**DDZETADZ+DS2DRHOS*DRHOSDZ;

            double S1TS2=S1*S2;
            double S1PS2=S1+S2;
            double S1PS2SQ=S1PS2*S1PS2;
            double FAC1=sqrt(S1PS2SQ-pow((2.e0*BI), 2));
            double AS=FAC1/(S1TS2*S1PS2SQ);
            double TERM1=1.e0/(S1TS2*S1PS2*FAC1);
            double FAC2=AS/S1PS2SQ;
            double DASDS1=TERM1-FAC2/S1*(S2*S2+S1*(3.e0*S1+4.e0*S2));
            double DASDS2=TERM1-FAC2/S2*(S1*S1+S2*(3.e0*S2+4.e0*S1));

            double DASDX=DASDS1*DS1DX+DASDS2*DS2DX;
            double DASDY=DASDS1*DS1DY+DASDS2*DS2DY;
            double DASDZ=DASDS1*DS1DZ+DASDS2*DS2DZ;

            *BX=*BX+F[I]*((2.e0*AS+Y*DASDY)**SPSS-(*XS-XSHIFT)*DASDZ +AS**DPSRR*(Y*Y**CPSS+Z**ZSWW));
            *BY=*BY-F[I]*Y*(AS**DPSRR**XS+DASDZ**CPSS+DASDX**SPSS);
            *BZ=*BZ+F[I]*((2.e0*AS+Y*DASDY)**CPSS+(*XS-XSHIFT)*DASDX -AS**DPSRR*(X**ZSWW+Y*Y**SPSS)); //           1
        }

    }//END

    /*-------------------------------------------------------------------------
     C
     SUBROUTINE TAIL87(X,Z,BX,BZ)

     IMPLICIT REAL*8 (A-H,O-Z)

     COMMON /WARP/ FIRST(3), RPS,WARP,D, OTHER(13)
     C
     C      'LONG' VERSION OF THE 1987 TAIL MAGNETIC FIELD MODEL
     C              (N.A.TSYGANENKO, PLANET. SPACE SCI., V.35, P.1347, 1987)
     C
     C      D   IS THE Y-DEPENDENT SHEET HALF-THICKNESS (INCREASING TOWARDS FLANKS)
     C      RPS  IS THE TILT-DEPENDENT SHIFT OF THE SHEET IN THE Z-DIRECTION,
     C           CORRESPONDING TO THE ASYMPTOTIC HINGING DISTANCE, DEFINED IN THE
     C           MAIN SUBROUTINE (TAILRC96) FROM THE PARAMETERS RH AND DR OF THE
     C           T96-TYPE MODULE, AND
     C      WARP  IS THE BENDING OF THE SHEET FLANKS IN THE Z-DIRECTION, DIRECTED
     C           OPPOSITE TO RPS, AND INCREASING WITH DIPOLE TILT AND |Y|
     */
    static void TAIL87(double const X,double const Z,double *BX,double *BZ, struct WARP_COMMON const* WARP_C)
    {
        //COMMON /WARP/ FIRST(3), RPS,WARP,D, OTHER(13)
        double const* RPS = &WARP_C->RPS;
        double const* WARP = &WARP_C->WARP;
        //double &D = WARP_C.D;

        //DATA DD/3./
        static double const DD = 3.;

        //DATA HPI,RT,XN,X1,X2,B0,B1,B2,XN21,XNR,ADLN
        //* /1.5707963,40.,-10.,
        //* -1.261,-0.663,0.391734,5.89715,24.6833,76.37,-0.1071,0.13238005/
        static double const HPI = 1.5707963,RT = 40.,XN = -10.,X1 = -1.261,X2 = -0.663,B0 = 0.391734,B1 = 5.89715,B2 = 24.6833,XN21 = 76.37,XNR = -0.1071,ADLN = 0.13238005;
        /*
         C                !!!   THESE ARE NEW VALUES OF  X1, X2, B0, B1, B2,
         C                       CORRESPONDING TO TSCALE=1, INSTEAD OF TSCALE=0.6
         C
         C  THE ABOVE QUANTITIES WERE DEFINED AS FOLLOWS:------------------------
         C       HPI=PI/2
         C       RT=40.      !  Z-POSITION OF UPPER AND LOWER ADDITIONAL SHEETS
         C       XN=-10.     !  INNER EDGE POSITION
         C
         C       TSCALE=1  !  SCALING FACTOR, DEFINING THE RATE OF INCREASE OF THE
         C                       CURRENT DENSITY TAILWARDS
         C
         c  ATTENTION !  NOW I HAVE CHANGED TSCALE TO:  TSCALE=1.0, INSTEAD OF 0.6
         c                  OF THE PREVIOUS VERSION
         c
         C       B0=0.391734
         C       B1=5.89715 *TSCALE
         C       B2=24.6833 *TSCALE**2
         C
         C    HERE ORIGINAL VALUES OF THE MODE AMPLITUDES (P.77, NB#3) WERE NORMALIZED
         C      SO THAT ASYMPTOTIC  BX=1  AT X=-200RE
         C
         C      X1=(4.589  -5.85) *TSCALE -(TSCALE-1.)*XN ! NONLINEAR PARAMETERS OF THE
         C                                                         CURRENT FUNCTION
         C      X2=(5.187  -5.85) *TSCALE -(TSCALE-1.)*XN
         c
         c
         C      XN21=(XN-X1)**2
         C      XNR=1./(XN-X2)
         C      ADLN=-DLOG(XNR**2*XN21)
         C
         C---------------------------------------------------------------
         */
        double ZS=Z - *RPS + *WARP;
        double ZP=Z-RT;
        double ZM=Z+RT;

        double XNX=XN-X;
        double XNX2=XNX*XNX;
        double XC1=X-X1;
        double XC2=X-X2;
        double XC22=XC2*XC2;
        double XR2=XC2*XNR;
        double XC12=XC1*XC1;
        double D2=DD*DD;    //  SQUARE OF THE TOTAL HALFTHICKNESS (DD=3Re for this mode)
        double B20=ZS*ZS+D2;
        double B2P=ZP*ZP+D2;
        double B2M=ZM*ZM+D2;
        double B=sqrt(B20);
        double BP=sqrt(B2P);
        double BM=sqrt(B2M);
        double XA1=XC12+B20;
        double XAP1=XC12+B2P;
        double XAM1=XC12+B2M;
        double XA2=1./(XC22+B20);
        double XAP2=1./(XC22+B2P);
        double XAM2=1./(XC22+B2M);
        double XNA=XNX2+B20;
        double XNAP=XNX2+B2P;
        double XNAM=XNX2+B2M;
        double F=B20-XC22;
        double FP=B2P-XC22;
        double FM=B2M-XC22;
        double XLN1=log(XN21/XNA);
        double XLNP1=log(XN21/XNAP);
        double XLNM1=log(XN21/XNAM);
        double XLN2=XLN1+ADLN;
        double XLNP2=XLNP1+ADLN;
        double XLNM2=XLNM1+ADLN;
        double ALN=0.25*(XLNP1+XLNM1-2.*XLN1);
        double S0=(atan(XNX/B)+HPI)/B;
        double S0P=(atan(XNX/BP)+HPI)/BP;
        double S0M=(atan(XNX/BM)+HPI)/BM;
        double S1=(XLN1*.5+XC1*S0)/XA1;
        double S1P=(XLNP1*.5+XC1*S0P)/XAP1;
        double S1M=(XLNM1*.5+XC1*S0M)/XAM1;
        double S2=(XC2*XA2*XLN2-XNR-F*XA2*S0)*XA2;
        double S2P=(XC2*XAP2*XLNP2-XNR-FP*XAP2*S0P)*XAP2;
        double S2M=(XC2*XAM2*XLNM2-XNR-FM*XAM2*S0M)*XAM2;
        double G1=(B20*S0-0.5*XC1*XLN1)/XA1;
        double G1P=(B2P*S0P-0.5*XC1*XLNP1)/XAP1;
        double G1M=(B2M*S0M-0.5*XC1*XLNM1)/XAM1;
        double G2=((0.5*F*XLN2+2.*S0*B20*XC2)*XA2+XR2)*XA2;
        double G2P=((0.5*FP*XLNP2+2.*S0P*B2P*XC2)*XAP2+XR2)*XAP2;
        double G2M=((0.5*FM*XLNM2+2.*S0M*B2M*XC2)*XAM2+XR2)*XAM2;

        *BX=B0*(ZS*S0-0.5*(ZP*S0P+ZM*S0M)) +B1*(ZS*S1-0.5*(ZP*S1P+ZM*S1M))+B2*(ZS*S2-0.5*(ZP*S2P+ZM*S2M));
        *BZ=B0*ALN+B1*(G1-0.5*(G1P+G1M))+B2*(G2-0.5*(G2P+G2M));

    }//END

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     SUBROUTINE TAILRC96(SPS,X,Y,Z,BXRC,BYRC,BZRC,BXT2,BYT2,BZT2,
     *   BXT3,BYT3,BZT3)
     c
     c  COMPUTES THE COMPONENTS OF THE FIELD OF THE MODEL RING CURRENT AND THREE
     c                   TAIL MODES WITH UNIT AMPLITUDES
     C      (FOR THE RING CURRENT, IT MEANS THE DISTURBANCE OF Bz=-1nT AT ORIGIN,
     C   AND FOR THE TAIL MODES IT MEANS MAXIMAL BX JUST ABOVE THE SHEET EQUAL 1 nT.
     */
    static void TAILRC96(double const SPS,double const X,double const Y,double const Z,double *BXRC,double *BYRC,double *BZRC,double *BXT2,double *BYT2,double *BZT2,double *BXT3,double *BYT3,double *BZT3)
    {
        //DIMENSION ARC(48),ATAIL2(48),ATAIL3(48)
        //COMMON /WARP/ CPSS,SPSS,DPSRR,RPS,WARP,D,XS,ZS,DXSX,DXSY,DXSZ,
        //*   DZSX,DZSY,DZSZ,DZETAS,DDZETADX,DDZETADY,DDZETADZ,ZSWW
        struct WARP_COMMON WARP_C;
        double *CPSS = &WARP_C.CPSS;
        double *SPSS = &WARP_C.SPSS;
        double *DPSRR = &WARP_C.DPSRR;
        double *RPS = &WARP_C.RPS;
        double *WARP = &WARP_C.WARP;
        double *D = &WARP_C.D;
        double *XS = &WARP_C.XS;
        double *ZS = &WARP_C.ZS;
        double *DXSX = &WARP_C.DXSX;
        double *DXSY = &WARP_C.DXSY;
        double *DXSZ = &WARP_C.DXSZ;
        double *DZSX = &WARP_C.DZSX;
        double *DZSY = &WARP_C.DZSY;
        double *DZSZ = &WARP_C.DZSZ;
        double *DZETAS = &WARP_C.DZETAS;
        double *DDZETADX = &WARP_C.DDZETADX;
        double *DDZETADY = &WARP_C.DDZETADY;
        double *DDZETADZ = &WARP_C.DDZETADZ;
        double *ZSWW = &WARP_C.ZSWW;

        //DATA ARC/-3.087699646,3.516259114,18.81380577,-13.95772338,
        //*  -5.497076303,0.1712890838,2.392629189,-2.728020808,-14.79349936,
        //*  11.08738083,4.388174084,0.2492163197E-01,0.7030375685,
        //*-.7966023165,-3.835041334,2.642228681,-0.2405352424,-0.7297705678,
        //* -0.3680255045,0.1333685557,2.795140897,-1.078379954,0.8014028630,
        //* 0.1245825565,0.6149982835,-0.2207267314,-4.424578723,1.730471572,
        //* -1.716313926,-0.2306302941,-0.2450342688,0.8617173961E-01,
        //*  1.54697858,-0.6569391113,-0.6537525353,0.2079417515,12.75434981,
        //*  11.37659788,636.4346279,1.752483754,3.604231143,12.83078674,
        //* 7.412066636,9.434625736,676.7557193,1.701162737,3.580307144,
        //*  14.64298662/
        static double const ARC[49] = {0./*dummy*/, -3.087699646,3.516259114,18.81380577,-13.95772338,
            -5.497076303,0.1712890838,2.392629189,-2.728020808,-14.79349936,
            11.08738083,4.388174084,0.2492163197E-01,0.7030375685,
            -.7966023165,-3.835041334,2.642228681,-0.2405352424,-0.7297705678,
            -0.3680255045,0.1333685557,2.795140897,-1.078379954,0.8014028630,
            0.1245825565,0.6149982835,-0.2207267314,-4.424578723,1.730471572,
            -1.716313926,-0.2306302941,-0.2450342688,0.8617173961E-01,
            1.54697858,-0.6569391113,-0.6537525353,0.2079417515,12.75434981,
            11.37659788,636.4346279,1.752483754,3.604231143,12.83078674,
            7.412066636,9.434625736,676.7557193,1.701162737,3.580307144,
            14.64298662};

        //DATA ATAIL2/.8747515218,-.9116821411,2.209365387,-2.159059518,
        //* -7.059828867,5.924671028,-1.916935691,1.996707344,-3.877101873,
        //* 3.947666061,11.38715899,-8.343210833,1.194109867,-1.244316975,
        //* 3.73895491,-4.406522465,-20.66884863,3.020952989,.2189908481,
        //* -.09942543549,-.927225562,.1555224669,.6994137909,-.08111721003,
        //* -.7565493881,.4686588792,4.266058082,-.3717470262,-3.920787807,
        //* .02298569870,.7039506341,-.5498352719,-6.675140817,.8279283559,
        //* -2.234773608,-1.622656137,5.187666221,6.802472048,39.13543412,
        //*  2.784722096,6.979576616,25.71716760,4.495005873,8.068408272,
        //* 93.47887103,4.158030104,9.313492566,57.18240483/
        static double const ATAIL2[49] = {0./*dummy*/, .8747515218,-.9116821411,2.209365387,-2.159059518,
            -7.059828867,5.924671028,-1.916935691,1.996707344,-3.877101873,
            3.947666061,11.38715899,-8.343210833,1.194109867,-1.244316975,
            3.73895491,-4.406522465,-20.66884863,3.020952989,.2189908481,
            -.09942543549,-.927225562,.1555224669,.6994137909,-.08111721003,
            -.7565493881,.4686588792,4.266058082,-.3717470262,-3.920787807,
            .02298569870,.7039506341,-.5498352719,-6.675140817,.8279283559,
            -2.234773608,-1.622656137,5.187666221,6.802472048,39.13543412,
            2.784722096,6.979576616,25.71716760,4.495005873,8.068408272,
            93.47887103,4.158030104,9.313492566,57.18240483};

        //DATA ATAIL3/-19091.95061,-3011.613928,20582.16203,4242.918430,
        //* -2377.091102,-1504.820043,19884.04650,2725.150544,-21389.04845,
        //* -3990.475093,2401.610097,1548.171792,-946.5493963,490.1528941,
        //* 986.9156625,-489.3265930,-67.99278499,8.711175710,-45.15734260,
        //* -10.76106500,210.7927312,11.41764141,-178.0262808,.7558830028,
        //*  339.3806753,9.904695974,69.50583193,-118.0271581,22.85935896,
        //* 45.91014857,-425.6607164,15.47250738,118.2988915,65.58594397,
        //* -201.4478068,-14.57062940,19.69877970,20.30095680,86.45407420,
        //* 22.50403727,23.41617329,48.48140573,24.61031329,123.5395974,
        //* 223.5367692,39.50824342,65.83385762,266.2948657/
        static double const ATAIL3[49] = {0./*dummy*/, -19091.95061,-3011.613928,20582.16203,4242.918430,
            -2377.091102,-1504.820043,19884.04650,2725.150544,-21389.04845,
            -3990.475093,2401.610097,1548.171792,-946.5493963,490.1528941,
            986.9156625,-489.3265930,-67.99278499,8.711175710,-45.15734260,
            -10.76106500,210.7927312,11.41764141,-178.0262808,.7558830028,
            339.3806753,9.904695974,69.50583193,-118.0271581,22.85935896,
            45.91014857,-425.6607164,15.47250738,118.2988915,65.58594397,
            -201.4478068,-14.57062940,19.69877970,20.30095680,86.45407420,
            22.50403727,23.41617329,48.48140573,24.61031329,123.5395974,
            223.5367692,39.50824342,65.83385762,266.2948657};

        //DATA RH,DR,G,D0,DELTADY/9.,4.,10.,2.,10./
        static double const RH = 9.,DR = 4.,G = 10.,D0 = 2.,DELTADY = 10.;
        /*
         C   TO ECONOMIZE THE CODE, WE FIRST CALCULATE COMMON VARIABLES, WHICH ARE
         C      THE SAME FOR ALL MODES, AND PUT THEM IN THE COMMON-BLOCK /WARP/
         */
        double DR2=DR*DR;
        double C11=sqrt(pow((1.e0+RH), 2)+DR2);
        double C12=sqrt(pow((1.e0-RH), 2)+DR2);
        double C1=C11-C12;
        double SPSC1=SPS/C1;
        *RPS=0.5*(C11+C12)*SPS;  //  THIS IS THE SHIFT OF OF THE SHEET WITH RESPECT
        //                            TO GSM EQ.PLANE FOR THE 3RD (ASYMPTOTIC) TAIL MODE

        double R=sqrt(X*X+Y*Y+Z*Z);
        double SQ1=sqrt(pow((R+RH), 2)+DR2);
        double SQ2=sqrt(pow((R-RH), 2)+DR2);
        double C=SQ1-SQ2;
        double CS=(R+RH)/SQ1-(R-RH)/SQ2;
        *SPSS=SPSC1/R*C;
        *CPSS=sqrt(1.e0-*SPSS**SPSS);
        *DPSRR=SPS/(R*R)*(CS*R-C)/sqrt(pow((R*C1), 2)-pow((C*SPS), 2));

        double WFAC=Y/(pow(Y, 4)+1.e4);   //   WARPING
        double W=WFAC * Y*Y*Y;
        double WS=4.e4*Y*WFAC*WFAC;
        *WARP=G*SPS*W;
        *XS=X**CPSS-Z**SPSS;
        *ZSWW=Z**CPSS+X**SPSS;  // "WW" MEANS "WITHOUT Y-Z WARPING" (IN X-Z ONLY)
        *ZS=*ZSWW + *WARP;

        *DXSX=*CPSS-X**ZSWW**DPSRR;
        *DXSY=-Y**ZSWW**DPSRR;
        *DXSZ=-*SPSS-Z**ZSWW**DPSRR;
        *DZSX=*SPSS+X**XS**DPSRR;
        *DZSY=*XS*Y**DPSRR  +G*SPS*WS;  //  THE LAST TERM IS FOR THE Y-Z WARP
        *DZSZ=*CPSS+*XS*Z**DPSRR;        //      (TAIL MODES ONLY)

        *D=D0+DELTADY*pow((Y/20.e0), 2);   //  SHEET HALF-THICKNESS FOR THE TAIL MODES
        double DDDY=DELTADY*Y*0.005e0;      //  (THICKENS TO FLANKS, BUT NO VARIATION
        //                                         ALONG X, IN CONTRAST TO RING CURRENT)

        *DZETAS=sqrt(*ZS**ZS+*D**D);  //  THIS IS THE SAME SIMPLE WAY TO SPREAD
        //                                        OUT THE SHEET, AS THAT USED IN T89
        *DDZETADX=*ZS**DZSX/ *DZETAS;
        *DDZETADY=(*ZS**DZSY+*D*DDDY)/ *DZETAS;
        *DDZETADZ=*ZS**DZSZ/ *DZETAS;

        double WX, WY, WZ;
        SHLCAR3X3(ARC,X,Y,Z,SPS,&WX,&WY,&WZ);
        double HX, HY, HZ;
        RINGCURR96(X,Y,Z,&HX,&HY,&HZ, &WARP_C);
        *BXRC=WX+HX;
        *BYRC=WY+HY;
        *BZRC=WZ+HZ;

        SHLCAR3X3(ATAIL2,X,Y,Z,SPS,&WX,&WY,&WZ);
        TAILDISK(X,Y,Z,&HX,&HY,&HZ, &WARP_C);
        *BXT2=WX+HX;
        *BYT2=WY+HY;
        *BZT2=WZ+HZ;

        SHLCAR3X3(ATAIL3,X,Y,Z,SPS,&WX,&WY,&WZ);
        TAIL87(X,Z,&HX,&HZ, &WARP_C);
        *BXT3=WX+HX;
        *BYT3=WY;
        *BZT3=WZ+HZ;

    }//END

    /********************************************************************

     SUBROUTINE DIPXYZ(X,Y,Z,BXX,BYX,BZX,BXY,BYY,BZY,BXZ,BYZ,BZZ)
     C
     C       RETURNS THE FIELD COMPONENTS PRODUCED BY THREE DIPOLES, EACH
     C        HAVING M=Me AND ORIENTED PARALLEL TO X,Y, and Z AXIS, RESP.
     */
    static void DIPXYZ(double const X,double const Y,double const Z,double *BXX,double *BYX,double *BZX,double *BXY,double *BYY,double *BZY,double *BXZ,double *BYZ,double *BZZ)
    {
        double X2=X*X;
        double Y2=Y*Y;
        double Z2=Z*Z;
        double R2=X2+Y2+Z2;

        double XMR5=30574.e0/(R2*R2*sqrt(R2));
        double XMR53=3.e0*XMR5;
        *BXX=XMR5*(3.e0*X2-R2);
        *BYX=XMR53*X*Y;
        *BZX=XMR53*X*Z;

        *BXY=*BYX;
        *BYY=XMR5*(3.e0*Y2-R2);
        *BZY=XMR53*Y*Z;

        *BXZ=*BZX;
        *BYZ=*BZY;
        *BZZ=XMR5*(3.e0*Z2-R2);

    }//END

    /*-------------------------------------------------------------------------
     C
     SUBROUTINE CIRCLE(X,Y,Z,RL,BX,BY,BZ)
     C
     C  RETURNS COMPONENTS OF THE FIELD FROM A CIRCULAR CURRENT LOOP OF RADIUS RL
     C  USES THE SECOND (MORE ACCURATE) APPROXIMATION GIVEN IN ABRAMOWITZ AND STEGUN
     */
    static void CIRCLE(double const X,double const Y,double const Z,double const RL,double *BX,double *BY,double *BZ)
    {
        double K, E;
        //DATA PI/3.141592654D0/
        static double const PI = 3.141592654e0;

        double RHO2=X*X+Y*Y;
        double RHO=sqrt(RHO2);
        double R22=Z*Z+pow((RHO+RL), 2);
        double R2=sqrt(R22);
        double R12=R22-4.e0*RHO*RL;
        double R32=0.5e0*(R12+R22);
        double XK2=1.e0-R12/R22;
        double XK2S=1.e0-XK2;
        double DL=log(1.e0/XK2S);
        K=1.38629436112e0+XK2S*(0.09666344259e0+XK2S*(0.03590092383+ XK2S*(0.03742563713+XK2S*0.01451196212))) +DL* (0.5e0+XK2S*(0.12498593597e0+XK2S*(0.06880248576e0+ XK2S*(0.03328355346e0+XK2S*0.00441787012e0))));
        E=1.e0+XK2S*(0.44325141463e0+XK2S*(0.0626060122e0+XK2S* (0.04757383546e0+XK2S*0.01736506451e0))) +DL* XK2S*(0.2499836831e0+XK2S*(0.09200180037e0+XK2S* (0.04069697526e0+XK2S*0.00526449639e0)));

        double BRHO;
        //IF (RHO.GT.1.D-6) THEN
        if (RHO>1.e-6) {
            BRHO=Z/(RHO2*R2)*(R32/R12*E-K); //  THIS IS NOT EXACTLY THE B-RHO COM-
        } else {//ELSE           !   PONENT - NOTE THE ADDITIONAL
            BRHO=PI*RL/R2*(RL-RHO)/R12*Z/(R32-RHO2);  //      DIVISION BY RHO
        }//ENDIF

        *BX=BRHO*X;
        *BY=BRHO*Y;
        *BZ=(K-E*(R32-2.e0*RL*RL)/R12)/R2;

    }//END

    /*-------------------------------------------------------------
     C
     SUBROUTINE CROSSLP(X,Y,Z,BX,BY,BZ,XC,RL,AL)
     C
     c   RETURNS FIELD COMPONENTS OF A PAIR OF LOOPS WITH A COMMON CENTER AND
     C    DIAMETER,  COINCIDING WITH THE X AXIS. THE LOOPS ARE INCLINED TO THE
     C    EQUATORIAL PLANE BY THE ANGLE AL (RADIANS) AND SHIFTED IN THE POSITIVE
     C     X-DIRECTION BY THE DISTANCE  XC.
     */
    static void CROSSLP(double const X,double const Y,double const Z,double *BX,double *BY,double *BZ,double const XC,double const RL,double const AL)
    {
        double CAL=cos(AL);
        double SAL=sin(AL);

        double Y1=Y*CAL-Z*SAL;
        double Z1=Y*SAL+Z*CAL;
        double Y2=Y*CAL+Z*SAL;
        double Z2=-Y*SAL+Z*CAL;
        double BX1,BY1,BZ1;
        CIRCLE(X-XC,Y1,Z1,RL,&BX1,&BY1,&BZ1);
        double BX2,BY2,BZ2;
        CIRCLE(X-XC,Y2,Z2,RL,&BX2,&BY2,&BZ2);
        *BX=BX1+BX2;
        *BY= (BY1+BY2)*CAL+(BZ1-BZ2)*SAL;
        *BZ=-(BY1-BY2)*SAL+(BZ1+BZ2)*CAL;

    }//END

    /**------------------------------------------------------------------------------
     C
     C
     SUBROUTINE  DIPLOOP1(XI,D)
     C
     C
     C      Calculates dependent model variables and their deriva-
     C  tives for given independent variables and model parame-
     C  ters.  Specifies model functions with free parameters which
     C  must be determined by means of least squares fits (RMS
     C  minimization procedure).
     C
     C      Description of parameters:
     C
     C  XI  - input vector containing independent variables;
     C  D   - output double precision vector containing
     C        calculated values for derivatives of dependent
     C        variables with respect to LINEAR model parameters;
     C
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C
     c  The  26 coefficients are moments (Z- and X-components) of 12 dipoles placed
     C    inside the  R1-shell,  PLUS amplitudes of two octagonal double loops.
     C     The dipoles with nonzero  Yi appear in pairs with equal moments.
     c                  (see the notebook #2, pp.102-103, for details)
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     */
    static void  DIPLOOP1(double const XI[],double Dx[], double Dy[], double Dz[], struct COORD11_COMMON const* COORD11_C, struct LOOPDIP1_COMMON const* LOOPDIP1_C, struct RHDR_COMMON const* RHDR_C)
    {
        //COMMON /COORD11/ XX(12),YY(12)
        double const* XX = COORD11_C->XX1;// [13]
        double const* YY = COORD11_C->YY1;// [13]
        //COMMON /LOOPDIP1/ TILT,XCENTRE(2),RADIUS(2),  DIPX,DIPY
        double const* TILT = &LOOPDIP1_C->TILT;
        double const* XCENTRE = LOOPDIP1_C->XCENTRE;// [3]
        double const* RADIUS = LOOPDIP1_C->RADIUS; // [3]
        double const* DIPX = &LOOPDIP1_C->DIPX;
        double const* DIPY = &LOOPDIP1_C->DIPY;
        //COMMON /RHDR/RH,DR
        double const* RH = &RHDR_C->RH;
        double const* DR = &RHDR_C->DR;

        //DIMENSION XI(4),D(3,26)

        double X = XI[1];
        double Y = XI[2];
        double Z = XI[3];
        double PS= XI[4];
        double SPS=sin(PS);

        //DO 1 I=1,12
        for(int I=1; I<=12; I++) {
            double R2=pow((XX[I]**DIPX), 2)+pow((YY[I]**DIPY), 2);
            double R=sqrt(R2);
            double RMRH=R-*RH;
            double RPRH=R+*RH;
            double DR2=*DR**DR;
            double SQM=sqrt(RMRH*RMRH+DR2);
            double SQP=sqrt(RPRH*RPRH+DR2);
            double C=SQP-SQM;
            double Q=sqrt(pow((*RH+1.e0), 2)+DR2)-sqrt(pow((*RH-1.e0), 2)+DR2);
            double SPSAS=SPS/R*C/Q;
            double CPSAS=sqrt(1.e0-SPSAS*SPSAS);
            double XD= (XX[I]**DIPX)*CPSAS;
            double YD= (YY[I]**DIPY);
            double ZD=-(XX[I]**DIPX)*SPSAS;

            double BX1X,BY1X,BZ1X,BX1Y,BY1Y,BZ1Y,BX1Z,BY1Z,BZ1Z;
            DIPXYZ(X-XD,Y-YD,Z-ZD,&BX1X,&BY1X,&BZ1X,&BX1Y,&BY1Y,&BZ1Y,&BX1Z,&BY1Z,&BZ1Z);

            double BX2X,BY2X,BZ2X,BX2Y,BY2Y,BZ2Y,BX2Z,BY2Z,BZ2Z;
            //IF (DABS(YD).GT.1.D-10) THEN
            if (fabs(YD)>1.e-10) {
                DIPXYZ(X-XD,Y+YD,Z-ZD,&BX2X,&BY2X,&BZ2X,&BX2Y,&BY2Y,&BZ2Y,&BX2Z,&BY2Z,&BZ2Z);
            } else {//ELSE
                BX2X=0.e0;
                BY2X=0.e0;
                BZ2X=0.e0;

                BX2Z=0.e0;
                BY2Z=0.e0;
                BZ2Z=0.e0;
            }//ENDIF

            Dx[I]=BX1Z+BX2Z;
            Dy[I]=BY1Z+BY2Z;
            Dz[I]=BZ1Z+BZ2Z;
            Dx[I+12]=(BX1X+BX2X)*SPS;
            Dy[I+12]=(BY1X+BY2X)*SPS;
            Dz[I+12]=(BZ1X+BZ2X)*SPS;
        }//1   CONTINUE

        double R2=pow((XCENTRE[1]+RADIUS[1]), 2);
        double R=sqrt(R2);
        double RMRH=R-*RH;
        double RPRH=R+*RH;
        double DR2=*DR**DR;
        double SQM=sqrt(RMRH*RMRH+DR2);
        double SQP=sqrt(RPRH*RPRH+DR2);
        double C=SQP-SQM;
        double Q=sqrt(pow((*RH+1.e0), 2)+DR2)-sqrt(pow((*RH-1.e0), 2)+DR2);
        double SPSAS=SPS/R*C/Q;
        double CPSAS=sqrt(1.e0-SPSAS*SPSAS);
        double XOCT1= X*CPSAS-Z*SPSAS;
        double YOCT1= Y;
        double ZOCT1= X*SPSAS+Z*CPSAS;

        double BXOCT1,BYOCT1,BZOCT1;
        CROSSLP(XOCT1,YOCT1,ZOCT1,&BXOCT1,&BYOCT1,&BZOCT1,XCENTRE[1],RADIUS[1],*TILT);
        Dx[25]=BXOCT1*CPSAS+BZOCT1*SPSAS;
        Dy[25]=BYOCT1;
        Dz[25]=-BXOCT1*SPSAS+BZOCT1*CPSAS;

        R2=pow((RADIUS[2]-XCENTRE[2]), 2);
        R=sqrt(R2);
        RMRH=R-*RH;
        RPRH=R+*RH;
        DR2=*DR**DR;
        SQM=sqrt(RMRH*RMRH+DR2);
        SQP=sqrt(RPRH*RPRH+DR2);
        C=SQP-SQM;
        Q=sqrt(pow((*RH+1.e0), 2)+DR2)-sqrt(pow((*RH-1.e0), 2)+DR2);
        SPSAS=SPS/R*C/Q;
        CPSAS=sqrt(1.e0-SPSAS*SPSAS);
        double XOCT2= X*CPSAS-Z*SPSAS -XCENTRE[2];
        double YOCT2= Y;
        double ZOCT2= X*SPSAS+Z*CPSAS;

        double BX,BY,BZ;
        CIRCLE(XOCT2,YOCT2,ZOCT2,RADIUS[2],&BX,&BY,&BZ);
        Dx[26] =  BX*CPSAS+BZ*SPSAS;
        Dy[26] =  BY;
        Dz[26] = -BX*SPSAS+BZ*CPSAS;

    }//END

    /*------------------------------------------------------------------------------
     SUBROUTINE  CONDIP1(XI,D)
     C
     C      Calculates dependent model variables and their derivatives for given
     C  independent variables and model parameters.  Specifies model functions with
     C  free parameters which must be determined by means of least squares fits
     C  (RMS minimization procedure).
     C
     C      Description of parameters:
     C
     C  XI  - input vector containing independent variables;
     C  D   - output double precision vector containing
     C        calculated values for derivatives of dependent
     C        variables with respect to LINEAR model parameters;
     C
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C
     c  The  79 coefficients are (1) 5 amplitudes of the conical harmonics, plus
     c                           (2) (9x3+5x2)x2=74 components of the dipole moments
     c              (see the notebook #2, pp.113-..., for details)
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     */
    static void CONDIP1(double const XI[],double Dx[],double Dy[],double Dz[], struct DX1_COMMON const* DX1_C, struct COORD21_COMMON const* COORD21_C)
    {
        //COMMON /DX1/ DX,SCALEIN,SCALEOUT
        double const* DX = &DX1_C->DX;
        double const* SCALEIN = &DX1_C->SCALEIN;
        double const* SCALEOUT = &DX1_C->SCALEOUT;
        //COMMON /COORD21/ XX(14),YY(14),ZZ(14)
        double const* XX = COORD21_C->XX2; // [15]
        double const* YY = COORD21_C->YY2; // [15]
        double const* ZZ = COORD21_C->ZZ2; // [15]

        //DIMENSION XI(4),D(3,79),CF(5),SF(5)
        double CF[6], SF[6];

        double X = XI[1];
        double Y = XI[2];
        double Z = XI[3];
        double PS= XI[4];
        double SPS=sin(PS);
        double CPS=cos(PS);

        double XSM=X*CPS-Z*SPS - *DX;
        double ZSM=Z*CPS+X*SPS;
        double RO2=XSM*XSM+Y*Y;
        double RO=sqrt(RO2);

        CF[1]=XSM/RO;
        SF[1]=Y/RO;

        CF[2]=CF[1]*CF[1]-SF[1]*SF[1];
        SF[2]=2.*SF[1]*CF[1];
        CF[3]=CF[2]*CF[1]-SF[2]*SF[1];
        SF[3]=SF[2]*CF[1]+CF[2]*SF[1];
        CF[4]=CF[3]*CF[1]-SF[3]*SF[1];
        SF[4]=SF[3]*CF[1]+CF[3]*SF[1];
        CF[5]=CF[4]*CF[1]-SF[4]*SF[1];
        SF[5]=SF[4]*CF[1]+CF[4]*SF[1];

        double R2=RO2+ZSM*ZSM;
        double R=sqrt(R2);
        double C=ZSM/R;
        double S=RO/R;
        double CH=sqrt(0.5e0*(1.e0+C));
        double SH=sqrt(0.5e0*(1.e0-C));
        double TNH=SH/CH;
        double CNH=1.e0/TNH;

        //DO 1 M=1,5
        for (int M=1; M<=5; M++) {
            double BT=M*CF[M]/(R*S)*(pow(TNH, M)+pow(CNH, M));
            double BF=-0.5e0*M*SF[M]/R*(pow(TNH, (M-1))/(CH*CH)-pow(CNH, (M-1))/(SH*SH));
            double BXSM=BT*C*CF[1]-BF*SF[1];
            double BY=BT*C*SF[1]+BF*CF[1];
            double BZSM=-BT*S;

            Dx[M]=BXSM*CPS+BZSM*SPS;
            Dy[M]=BY;
            Dz[M]=-BXSM*SPS+BZSM*CPS; //         1
        }

        XSM = X*CPS-Z*SPS;
        ZSM = Z*CPS+X*SPS;

        //DO 2 I=1,9
        for (int I=1; I<=9; I++) {
            double XD, YD;
            //IF (I.EQ.3.OR.I.EQ.5.OR.I.EQ.6) THEN
            if (I==3 || I==5 || I==6) {
                XD =  XX[I]**SCALEIN;
                YD =  YY[I]**SCALEIN;
            } else {//ELSE
                XD =  XX[I]**SCALEOUT;
                YD =  YY[I]**SCALEOUT;
            }//ENDIF

            double ZD =  ZZ[I];

            double BX1X,BY1X,BZ1X,BX1Y,BY1Y,BZ1Y,BX1Z,BY1Z,BZ1Z;
            DIPXYZ(XSM-XD,Y-YD,ZSM-ZD,&BX1X,&BY1X,&BZ1X,&BX1Y,&BY1Y,&BZ1Y,&BX1Z,&BY1Z,&BZ1Z);
            double BX2X,BY2X,BZ2X,BX2Y,BY2Y,BZ2Y,BX2Z,BY2Z,BZ2Z;
            DIPXYZ(XSM-XD,Y+YD,ZSM-ZD,&BX2X,&BY2X,&BZ2X,&BX2Y,&BY2Y,&BZ2Y,&BX2Z,&BY2Z,&BZ2Z);
            double BX3X,BY3X,BZ3X,BX3Y,BY3Y,BZ3Y,BX3Z,BY3Z,BZ3Z;
            DIPXYZ(XSM-XD,Y-YD,ZSM+ZD,&BX3X,&BY3X,&BZ3X,&BX3Y,&BY3Y,&BZ3Y,&BX3Z,&BY3Z,&BZ3Z);
            double BX4X,BY4X,BZ4X,BX4Y,BY4Y,BZ4Y,BX4Z,BY4Z,BZ4Z;
            DIPXYZ(XSM-XD,Y+YD,ZSM+ZD,&BX4X,&BY4X,&BZ4X,&BX4Y,&BY4Y,&BZ4Y,&BX4Z,&BY4Z,&BZ4Z);

            int IX=I*3+3;
            int IY=IX+1;
            int IZ=IY+1;

            Dx[IX]=(BX1X+BX2X-BX3X-BX4X)*CPS+(BZ1X+BZ2X-BZ3X-BZ4X)*SPS;
            Dy[IX]= BY1X+BY2X-BY3X-BY4X;
            Dz[IX]=(BZ1X+BZ2X-BZ3X-BZ4X)*CPS-(BX1X+BX2X-BX3X-BX4X)*SPS;

            Dx[IY]=(BX1Y-BX2Y-BX3Y+BX4Y)*CPS+(BZ1Y-BZ2Y-BZ3Y+BZ4Y)*SPS;
            Dy[IY]= BY1Y-BY2Y-BY3Y+BY4Y;
            Dz[IY]=(BZ1Y-BZ2Y-BZ3Y+BZ4Y)*CPS-(BX1Y-BX2Y-BX3Y+BX4Y)*SPS;

            Dx[IZ]=(BX1Z+BX2Z+BX3Z+BX4Z)*CPS+(BZ1Z+BZ2Z+BZ3Z+BZ4Z)*SPS;
            Dy[IZ]= BY1Z+BY2Z+BY3Z+BY4Z;
            Dz[IZ]=(BZ1Z+BZ2Z+BZ3Z+BZ4Z)*CPS-(BX1Z+BX2Z+BX3Z+BX4Z)*SPS;

            IX=IX+27;
            IY=IY+27;
            IZ=IZ+27;

            Dx[IX]=SPS*((BX1X+BX2X+BX3X+BX4X)*CPS+(BZ1X+BZ2X+BZ3X+BZ4X)*SPS);
            Dy[IX]=SPS*(BY1X+BY2X+BY3X+BY4X);
            Dz[IX]=SPS*((BZ1X+BZ2X+BZ3X+BZ4X)*CPS-(BX1X+BX2X+BX3X+BX4X)*SPS);

            Dx[IY]=SPS*((BX1Y-BX2Y+BX3Y-BX4Y)*CPS+(BZ1Y-BZ2Y+BZ3Y-BZ4Y)*SPS);
            Dy[IY]=SPS*(BY1Y-BY2Y+BY3Y-BY4Y);
            Dz[IY]=SPS*((BZ1Y-BZ2Y+BZ3Y-BZ4Y)*CPS-(BX1Y-BX2Y+BX3Y-BX4Y)*SPS);

            Dx[IZ]=SPS*((BX1Z+BX2Z-BX3Z-BX4Z)*CPS+(BZ1Z+BZ2Z-BZ3Z-BZ4Z)*SPS);
            Dy[IZ]=SPS*(BY1Z+BY2Z-BY3Z-BY4Z);
            Dz[IZ]=SPS*((BZ1Z+BZ2Z-BZ3Z-BZ4Z)*CPS-(BX1Z+BX2Z-BX3Z-BX4Z)*SPS);
        }//2   CONTINUE

        //DO 3 I=1,5
        for (int I=1; I<=5; I++) {
            double ZD=ZZ[I+9];
            double BX1X,BY1X,BZ1X,BX1Y,BY1Y,BZ1Y,BX1Z,BY1Z,BZ1Z;
            DIPXYZ(XSM,Y,ZSM-ZD,&BX1X,&BY1X,&BZ1X,&BX1Y,&BY1Y,&BZ1Y,&BX1Z,&BY1Z,&BZ1Z);
            double BX2X,BY2X,BZ2X,BX2Y,BY2Y,BZ2Y,BX2Z,BY2Z,BZ2Z;
            DIPXYZ(XSM,Y,ZSM+ZD,&BX2X,&BY2X,&BZ2X,&BX2Y,&BY2Y,&BZ2Y,&BX2Z,&BY2Z,&BZ2Z);
            int IX=58+I*2;
            int IZ=IX+1;
            Dx[IX]=(BX1X-BX2X)*CPS+(BZ1X-BZ2X)*SPS;
            Dy[IX]=BY1X-BY2X;
            Dz[IX]=(BZ1X-BZ2X)*CPS-(BX1X-BX2X)*SPS;

            Dx[IZ]=(BX1Z+BX2Z)*CPS+(BZ1Z+BZ2Z)*SPS;
            Dy[IZ]=BY1Z+BY2Z;
            Dz[IZ]=(BZ1Z+BZ2Z)*CPS-(BX1Z+BX2Z)*SPS;

            IX=IX+10;
            IZ=IZ+10;
            Dx[IX]=SPS*((BX1X+BX2X)*CPS+(BZ1X+BZ2X)*SPS);
            Dy[IX]=SPS*(BY1X+BY2X);
            Dz[IX]=SPS*((BZ1X+BZ2X)*CPS-(BX1X+BX2X)*SPS);

            Dx[IZ]=SPS*((BX1Z-BX2Z)*CPS+(BZ1Z-BZ2Z)*SPS);
            Dy[IZ]=SPS*(BY1Z-BY2Z);
            Dz[IZ]=SPS*((BZ1Z-BZ2Z)*CPS-(BX1Z-BX2Z)*SPS); //      3
        }

    }//END

    /**$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     C
     SUBROUTINE  BIRK1SHLD(PS,X,Y,Z,BX,BY,BZ)
     C
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C
     C  The 64 linear parameters are amplitudes of the "box" harmonics.
     c The 16 nonlinear parameters are the scales Pi, and Qk entering the arguments
     C  of sines/cosines and exponents in each of  32 cartesian harmonics
     c  N.A. Tsyganenko, Spring 1994, adjusted for the Birkeland field Aug.22, 1995
     c    Revised  June 12, 1996.
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     */
    static void  BIRK1SHLD(double const PS,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ)
    {
        //DIMENSION A(80)
        //DIMENSION P1(4),R1(4),Q1(4),S1(4),RP(4),RR(4),RQ(4),RS(4)
        double RP[5],RR[5],RQ[5],RS[5];
        double const *P1, *R1, *Q1, *S1;

        //EQUIVALENCE (P1(1),A(65)),(R1(1),A(69)),(Q1(1),A(73)),(S1(1),A(77))

        /*DATA A/1.174198045,-1.463820502,4.840161537,-3.674506864,
         * 82.18368896,-94.94071588,-4122.331796,4670.278676,-21.54975037,
         * 26.72661293,-72.81365728,44.09887902,40.08073706,-51.23563510,
         * 1955.348537,-1940.971550,794.0496433,-982.2441344,1889.837171,
         * -558.9779727,-1260.543238,1260.063802,-293.5942373,344.7250789,
         * -773.7002492,957.0094135,-1824.143669,520.7994379,1192.484774,
         * -1192.184565,89.15537624,-98.52042999,-0.8168777675E-01,
         * 0.4255969908E-01,0.3155237661,-0.3841755213,2.494553332,
         * -0.6571440817E-01,-2.765661310,0.4331001908,0.1099181537,
         * -0.6154126980E-01,-0.3258649260,0.6698439193,-5.542735524,
         * 0.1604203535,5.854456934,-0.8323632049,3.732608869,-3.130002153,
         * 107.0972607,-32.28483411,-115.2389298,54.45064360,-0.5826853320,
         * -3.582482231,-4.046544561,3.311978102,-104.0839563,30.26401293,
         * 97.29109008,-50.62370872,-296.3734955,127.7872523,5.303648988,
         * 10.40368955,69.65230348,466.5099509,1.645049286,3.825838190,
         * 11.66675599,558.9781177,1.826531343,2.066018073,25.40971369,
         * 990.2795225,2.319489258,4.555148484,9.691185703,591.8280358/ */
        static double const A[81] = {0./*dummy*/, 1.174198045,-1.463820502,4.840161537,-3.674506864,
            82.18368896,-94.94071588,-4122.331796,4670.278676,-21.54975037,
            26.72661293,-72.81365728,44.09887902,40.08073706,-51.23563510,
            1955.348537,-1940.971550,794.0496433,-982.2441344,1889.837171,
            -558.9779727,-1260.543238,1260.063802,-293.5942373,344.7250789,
            -773.7002492,957.0094135,-1824.143669,520.7994379,1192.484774,
            -1192.184565,89.15537624,-98.52042999,-0.8168777675E-01,
            0.4255969908E-01,0.3155237661,-0.3841755213,2.494553332,
            -0.6571440817E-01,-2.765661310,0.4331001908,0.1099181537,
            -0.6154126980E-01,-0.3258649260,0.6698439193,-5.542735524,
            0.1604203535,5.854456934,-0.8323632049,3.732608869,-3.130002153,
            107.0972607,-32.28483411,-115.2389298,54.45064360,-0.5826853320,
            -3.582482231,-4.046544561,3.311978102,-104.0839563,30.26401293,
            97.29109008,-50.62370872,-296.3734955,127.7872523,5.303648988,
            10.40368955,69.65230348,466.5099509,1.645049286,3.825838190,
            11.66675599,558.9781177,1.826531343,2.066018073,25.40971369,
            990.2795225,2.319489258,4.555148484,9.691185703,591.8280358};

        P1 = A + 64;
        R1 = A + 68;
        Q1 = A + 72;
        S1 = A + 76;

        *BX=0.e0;
        *BY=0.e0;
        *BZ=0.e0;
        double CPS=cos(PS);
        double SPS=sin(PS);
        double S3PS=4.e0*CPS*CPS-1.e0;

        //DO 11 I=1,4
        for (int I=1; I<=4; I++) {
            RP[I]=1.e0/P1[I];
            RR[I]=1.e0/R1[I];
            RQ[I]=1.e0/Q1[I];
            RS[I]=1.e0/S1[I]; //      11
        }

        int L=0;

        /*DO 1 M=1,2     !    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
         C                           AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY) */
        for (int M=1; M<=2; M++) {
            //DO 2 I=1,4
            for (int I=1; I<=4; I++) {
                double CYPI=cos(Y*RP[I]);
                double CYQI=cos(Y*RQ[I]);
                double SYPI=sin(Y*RP[I]);
                double SYQI=sin(Y*RQ[I]);

                //DO 3 K=1,4
                for (int K=1; K<=4; K++) {
                    double SZRK=sin(Z*RR[K]);
                    double CZSK=cos(Z*RS[K]);
                    double CZRK=cos(Z*RR[K]);
                    double SZSK=sin(Z*RS[K]);
                    double SQPR=sqrt(RP[I]*RP[I]+RR[K]*RR[K]);
                    double SQQS=sqrt(RQ[I]*RQ[I]+RS[K]*RS[K]);
                    double EPR=exp(X*SQPR);
                    double EQS=exp(X*SQQS);

                    /*DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
                     C                                  AND N=2 IS FOR THE SECOND ONE */
                    double HX = 0.0, HY = 0.0, HZ = 0.0;
                    for (int N=1; N<=2; N++) {
                        //IF (M.EQ.1) THEN
                        if (M==1) {
                            //IF (N.EQ.1) THEN
                            if (N==1) {
                                HX=-SQPR*EPR*CYPI*SZRK;
                                HY=RP[I]*EPR*SYPI*SZRK;
                                HZ=-RR[K]*EPR*CYPI*CZRK;
                            } else {//ELSE
                                HX=HX*CPS;
                                HY=HY*CPS;
                                HZ=HZ*CPS;
                            }//ENDIF
                        } else {//ELSE
                            //IF (N.EQ.1) THEN
                            if (N==1) {
                                HX=-SPS*SQQS*EQS*CYQI*CZSK;
                                HY=SPS*RQ[I]*EQS*SYQI*CZSK;
                                HZ=SPS*RS[K]*EQS*CYQI*SZSK;
                            } else {//ELSE
                                HX=HX*S3PS;
                                HY=HY*S3PS;
                                HZ=HZ*S3PS;
                            }//ENDIF
                        }//ENDIF
                        L=L+1;

                        *BX=*BX+A[L]*HX;
                        *BY=*BY+A[L]*HY;
                        *BZ=*BZ+A[L]*HZ; //        4
                    }
                }//3   CONTINUE
            }//2   CONTINUE
        }//1   CONTINUE

    }//END

    /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     C
     SUBROUTINE BIRK1TOT_02(PS,X,Y,Z,BX,BY,BZ)
     C
     C  THIS IS THE SECOND VERSION OF THE ANALYTICAL MODEL OF THE REGION 1 FIELD
     C   BASED ON A SEPARATE REPRESENTATION OF THE POTENTIAL FIELD IN THE INNER AND
     C   OUTER SPACE, MAPPED BY MEANS OF A SPHERO-DIPOLAR COORDINATE SYSTEM (NB #3,
     C   P.91).   THE DIFFERENCE FROM THE FIRST ONE IS THAT INSTEAD OF OCTAGONAL
     C   CURRENT LOOPS, CIRCULAR ONES ARE USED IN THIS VERSION FOR APPROXIMATING THE
     C   FIELD IN THE OUTER REGION, WHICH IS FASTER.
     */
    static void BIRK1TOT_02(double const PS,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ)
    {
        //DIMENSION D1(3,26),D2(3,79),XI(4),C1(26),C2(79)
        double XI[5];

        //COMMON /COORD11/ XX1(12),YY1(12)
        struct COORD11_COMMON COORD11_C;
        double *XX1 = COORD11_C.XX1;// [13]
        double *YY1 = COORD11_C.YY1;// [13]
        //COMMON /RHDR/ RH,DR
        struct RHDR_COMMON RHDR_C;
        double *RH = &RHDR_C.RH;
        double *DR = &RHDR_C.DR;
        //COMMON /LOOPDIP1/ TILT,XCENTRE(2),RADIUS(2), DIPX,DIPY
        struct LOOPDIP1_COMMON LOOPDIP1_C;
        double *TILT = &LOOPDIP1_C.TILT;
        double *XCENTRE = LOOPDIP1_C.XCENTRE;// [3]
        double *RADIUS = LOOPDIP1_C.RADIUS; // [3]
        double *DIPX = &LOOPDIP1_C.DIPX;
        double *DIPY = &LOOPDIP1_C.DIPY;

        //COMMON /COORD21/ XX2(14),YY2(14),ZZ2(14)
        struct COORD21_COMMON COORD21_C;
        double *XX2 = COORD21_C.XX2; // [15]
        double *YY2 = COORD21_C.YY2; // [15]
        double *ZZ2 = COORD21_C.ZZ2; // [15]
        //COMMON /DX1/ DX,SCALEIN,SCALEOUT
        struct DX1_COMMON DX1_C;
        double *DX = &DX1_C.DX;
        double *SCALEIN = &DX1_C.SCALEIN;
        double *SCALEOUT = &DX1_C.SCALEOUT;

        //DATA C1/-0.911582E-03,-0.376654E-02,-0.727423E-02,-0.270084E-02,
        //* -0.123899E-02,-0.154387E-02,-0.340040E-02,-0.191858E-01,
        //* -0.518979E-01,0.635061E-01,0.440680,-0.396570,0.561238E-02,
        //*  0.160938E-02,-0.451229E-02,-0.251810E-02,-0.151599E-02,
        //* -0.133665E-02,-0.962089E-03,-0.272085E-01,-0.524319E-01,
        //*  0.717024E-01,0.523439,-0.405015,-89.5587,23.2806/
        static double const C1[27] = {0./*dummy*/, -0.911582e-03,-0.376654e-02,-0.727423e-02,-0.270084e-02,
            -0.123899e-02,-0.154387e-02,-0.340040e-02,-0.191858e-01,
            -0.518979e-01,0.635061e-01,0.440680,-0.396570,0.561238e-02,
            0.160938e-02,-0.451229e-02,-0.251810e-02,-0.151599e-02,
            -0.133665e-02,-0.962089e-03,-0.272085e-01,-0.524319e-01,
            0.717024e-01,0.523439,-0.405015,-89.5587,23.2806};

        /* DATA C2/6.04133,.305415,.606066E-02,.128379E-03,-.179406E-04,
         * 1.41714,-27.2586,-4.28833,-1.30675,35.5607,8.95792,.961617E-03,
         * -.801477E-03,-.782795E-03,-1.65242,-16.5242,-5.33798,.424878E-03,
         * .331787E-03,-.704305E-03,.844342E-03,.953682E-04,.886271E-03,
         * 25.1120,20.9299,5.14569,-44.1670,-51.0672,-1.87725,20.2998,
         * 48.7505,-2.97415,3.35184,-54.2921,-.838712,-10.5123,70.7594,
         * -4.94104,.106166E-03,.465791E-03,-.193719E-03,10.8439,-29.7968,
         *  8.08068,.463507E-03,-.224475E-04,.177035E-03,-.317581E-03,
         * -.264487E-03,.102075E-03,7.71390,10.1915,-4.99797,-23.1114,
         *-29.2043,12.2928,10.9542,33.6671,-9.3851,.174615E-03,-.789777E-06,
         * .686047E-03,.460104E-04,-.345216E-02,.221871E-02,.110078E-01,
         * -.661373E-02,.249201E-02,.343978E-01,-.193145E-05,.493963E-05,
         * -.535748E-04,.191833E-04,-.100496E-03,-.210103E-03,-.232195E-02,
         * .315335E-02,-.134320E-01,-.263222E-01/ */
        static double const C2[80] = {0./*dummy*/, 6.04133,.305415,.606066e-02,.128379e-03,-.179406e-04,
            1.41714,-27.2586,-4.28833,-1.30675,35.5607,8.95792,.961617e-03,
            -.801477e-03,-.782795e-03,-1.65242,-16.5242,-5.33798,.424878e-03,
            .331787e-03,-.704305e-03,.844342e-03,.953682e-04,.886271e-03,
            25.1120,20.9299,5.14569,-44.1670,-51.0672,-1.87725,20.2998,
            48.7505,-2.97415,3.35184,-54.2921,-.838712,-10.5123,70.7594,
            -4.94104,.106166e-03,.465791e-03,-.193719e-03,10.8439,-29.7968,
            8.08068,.463507e-03,-.224475e-04,.177035e-03,-.317581e-03,
            -.264487e-03,.102075e-03,7.71390,10.1915,-4.99797,-23.1114,
            -29.2043,12.2928,10.9542,33.6671,-9.3851,.174615e-03,-.789777e-06,
            .686047e-03,.460104e-04,-.345216e-02,.221871e-02,.110078e-01,
            -.661373e-02,.249201e-02,.343978e-01,-.193145e-05,.493963e-05,
            -.535748e-04,.191833e-04,-.100496e-03,-.210103e-03,-.232195e-02,
            .315335e-02,-.134320e-01,-.263222e-01};

        //DATA TILT,XCENTRE,RADIUS,DIPX,DIPY /1.00891,2.28397,-5.60831,
        //* 1.86106,7.83281,1.12541,0.945719/
        *TILT = 1.00891;
        XCENTRE[1] = 2.28397;
        XCENTRE[2] = -5.60831;
        RADIUS[1] = 1.86106;
        RADIUS[2] = 7.83281;
        *DIPX = 1.12541;
        *DIPY = 0.945719;

        //DATA DX,SCALEIN,SCALEOUT /-0.16D0,0.08D0,0.4D0/
        *DX = -0.16e0, *SCALEIN = 0.08e0, *SCALEOUT = 0.4e0;
        //DATA XX1/-11.D0,2*-7.D0,2*-3.D0,3*1.D0,2*5.D0,2*9.D0/
        //DATA YY1/2.D0,0.D0,4.D0,2.D0,6.D0,0.D0,4.D0,8.D0,2.D0,6.D0,0.D0,
        static double const XX1tmp[13] = {0./*dummy*/, -11.e0,-7.e0, -7.e0,-3.e0, -3.e0,1.e0, 1.e0, 1.e0,5.e0, 5.e0,9.e0, 9.e0};
        static double const YY1tmp[13] = {0./*dummy*/, 2.e0,0.e0,4.e0,2.e0,6.e0,0.e0,4.e0,8.e0,2.e0,6.e0,0.e0,4.e0};
        memcpy(XX1, XX1tmp, sizeof(XX1tmp));
        memcpy(YY1, YY1tmp, sizeof(YY1tmp));

        //DATA XX2/-10.D0,-7.D0,2*-4.D0,0.D0,2*4.D0,7.D0,10.D0,5*0.D0/
        //DATA YY2/3.D0,6.D0,3.D0,9.D0,6.D0,3.D0,9.D0,6.D0,3.D0,5*0.D0/
        //DATA ZZ2/2*20.D0,4.D0,20.D0,2*4.D0,3*20.D0,2.D0,3.D0,4.5D0,
        //*  7.D0,10.D0/
        static double const XX2tmp[15] = {0./*dummy*/, -10.e0,-7.e0,-4.e0, -4.e0,0.e0,4.e0, 4.e0,7.e0,10.e0,0.e0, 0.e0, 0.e0, 0.e0, 0.e0};
        static double const YY2tmp[15] = {0./*dummy*/, 3.e0,6.e0,3.e0,9.e0,6.e0,3.e0,9.e0,6.e0,3.e0, 0.e0, 0.e0, 0.e0, 0.e0, 0.e0};
        static double const ZZ2tmp[15] = {0./*dummy*/, 20.e0, 20.e0,4.e0,20.e0, 4.e0, 4.e0, 20.e0, 20.e0, 20.e0,2.e0,3.e0,4.5e0,7.e0,10.e0};
        memcpy(XX2, XX2tmp, sizeof(XX2tmp));
        memcpy(YY2, YY2tmp, sizeof(YY2tmp));
        memcpy(ZZ2, ZZ2tmp, sizeof(ZZ2tmp));

        /*DATA RH,DR /9.D0,4.D0/   !  RH IS THE "HINGING DISTANCE" AND DR IS THE
         C                                TRANSITION SCALE LENGTH, DEFINING THE
         C                                CURVATURE  OF THE WARPING (SEE P.89, NB #2)*/
        *RH = 9.e0;
        *DR = 4.e0;

        /*DATA XLTDAY,XLTNGHT /78.D0,70.D0/  !  THESE ARE LATITUDES OF THE R-1 OVAL
         C                                             AT NOON AND AT MIDNIGHT*/
        static double const XLTDAY = 78.e0,XLTNGHT = 70.e0;

        /*DATA DTET0 /0.034906/   !   THIS IS THE LATITUDINAL HALF-THICKNESS OF THE
         C                                  R-1 OVAL (THE INTERPOLATION REGION BETWEEN
         C                                    THE HIGH-LAT. AND THE PLASMA SHEET)*/
        static double const DTET0 = 0.034906;

        double TNOONN=(90.e0-XLTDAY)*0.01745329e0;
        double TNOONS=3.141592654e0-TNOONN;     /* HERE WE ASSUME THAT THE POSITIONS OF
                                                 C                                          THE NORTHERN AND SOUTHERN R-1 OVALS
                                                 C                                          ARE SYMMETRIC IN THE SM-COORDINATES*/
        double DTETDN=(XLTDAY-XLTNGHT)*0.01745329e0;
        double DR2=*DR**DR;

        double SPS=sin(PS);
        double R2=X*X+Y*Y+Z*Z;
        double R=sqrt(R2);
        double R3=R*R2;

        double RMRH=R-*RH;
        double RPRH=R+*RH;
        double SQM=sqrt(RMRH*RMRH+DR2);
        double SQP=sqrt(RPRH*RPRH+DR2);
        double C=SQP-SQM;
        double Q=sqrt(pow((*RH+1.e0), 2)+DR2)-sqrt(pow((*RH-1.e0), 2)+DR2);
        double SPSAS=SPS/R*C/Q;
        double CPSAS=sqrt(1.e0-SPSAS*SPSAS);
        double XAS = X*CPSAS-Z*SPSAS;
        double ZAS = X*SPSAS+Z*CPSAS;
        double PAS;
        //IF (XAS.NE.0.D0.OR.Y.NE.0.D0) THEN
        if (XAS!=0.e0 || Y!=0.e0) {
            PAS = atan2(Y,XAS);
        } else {//ELSE
            PAS=0.e0;
        }//ENDIF

        double TAS=atan2(sqrt(XAS*XAS+Y*Y),ZAS);
        double STAS=sin(TAS);
        double F=STAS/pow((pow(STAS, 6)*(1.e0-R3)+R3), 0.1666666667e0);

        double TET0=asin(F);
        //IF (TAS.GT.1.5707963D0) TET0=3.141592654D0-TET0
        if (TAS>1.5707963e0) TET0=3.141592654e0-TET0;
        double DTET=DTETDN*pow(sin(PAS*0.5e0), 2);
        double TETR1N=TNOONN+DTET;
        double TETR1S=TNOONS-DTET;
        /*
         C NOW LET'S DEFINE WHICH OF THE FOUR REGIONS (HIGH-LAT., NORTHERN PSBL,
         C   PLASMA SHEET, SOUTHERN PSBL) DOES THE POINT (X,Y,Z) BELONG TO:
         */
        int LOC;
        //IF (TET0.LT.TETR1N-DTET0.OR.TET0.GT.TETR1S+DTET0)  LOC=1 ! HIGH-LAT.
        if (TET0<TETR1N-DTET0 || TET0>TETR1S+DTET0)  LOC=1;
        //IF (TET0.GT.TETR1N+DTET0.AND.TET0.LT.TETR1S-DTET0) LOC=2 ! PL.SHEET
        if (TET0>TETR1N+DTET0 && TET0<TETR1S-DTET0) LOC=2;
        //IF (TET0.GE.TETR1N-DTET0.AND.TET0.LE.TETR1N+DTET0) LOC=3 ! NORTH PSBL
        if (TET0>=TETR1N-DTET0 && TET0<=TETR1N+DTET0) LOC=3;
        //IF (TET0.GE.TETR1S-DTET0.AND.TET0.LE.TETR1S+DTET0) LOC=4 ! SOUTH PSBL
        if (TET0>=TETR1S-DTET0 && TET0<=TETR1S+DTET0) LOC=4;

        //IF (LOC.EQ.1) THEN   ! IN THE HIGH-LAT. REGION USE THE SUBROUTINE DIPOCT
        if (LOC==1) {
            //      print *, '  LOC=1 (HIGH-LAT)'    !  (test printout; disabled now)
            XI[1]=X;
            XI[2]=Y;
            XI[3]=Z;
            XI[4]=PS;
            double D1x[27], D1y[27], D1z[27]; // D1(3, 26)
            DIPLOOP1(XI,D1x, D1y, D1z, &COORD11_C, &LOOPDIP1_C, &RHDR_C);
            *BX=0.0;
            *BY=0.0;
            *BZ=0.0;
            //DO 1 I=1,26
            for(int I=1; I<=26; I++) {
                *BX=*BX+C1[I]*D1x[I];
                *BY=*BY+C1[I]*D1y[I];
                *BZ=*BZ+C1[I]*D1z[I]; //       1
            }
        }//ENDIF                                           !  END OF THE CASE 1

        //IF (LOC.EQ.2) THEN
        if (LOC==2) {
            //           print *, '  LOC=2 (PLASMA SHEET)'  !  (test printout; disabled now)
            XI[1]=X;
            XI[2]=Y;
            XI[3]=Z;
            XI[4]=PS;
            double D2x[80], D2y[80], D2z[80]; // D1(3, 26)
            CONDIP1(XI,D2x, D2y, D2z, &DX1_C, &COORD21_C);
            *BX=0.e0;
            *BY=0.e0;
            *BZ=0.e0;
            //DO 2 I=1,79
            for (int I=1; I<=79; I++) {
                *BX=*BX+C2[I]*D2x[I];
                *BY=*BY+C2[I]*D2y[I];
                *BZ=*BZ+C2[I]*D2z[I]; //       2
            }
        }//ENDIF                                           !   END OF THE CASE 2

        //IF (LOC.EQ.3) THEN
        if (LOC==3) {
            //       print *, '  LOC=3 (north PSBL)'  !  (test printout; disabled now)
            double T01=TETR1N-DTET0;
            double T02=TETR1N+DTET0;
            double SQR=sqrt(R);
            double ST01AS=SQR/pow((R3+1.e0/pow(sin(T01), 6)-1.e0), 0.1666666667);
            double ST02AS=SQR/pow((R3+1.e0/pow(sin(T02), 6)-1.e0), 0.1666666667);
            double CT01AS=sqrt(1.e0-ST01AS*ST01AS);
            double CT02AS=sqrt(1.e0-ST02AS*ST02AS);
            double XAS1=R*ST01AS*cos(PAS);
            double Y1=  R*ST01AS*sin(PAS);
            double ZAS1=R*CT01AS;
            double X1=XAS1*CPSAS+ZAS1*SPSAS;
            double Z1=-XAS1*SPSAS+ZAS1*CPSAS; // X1,Y1,Z1 ARE COORDS OF THE NORTHERN
            //                                                      BOUNDARY POINT
            XI[1]=X1;
            XI[2]=Y1;
            XI[3]=Z1;
            XI[4]=PS;
            double D1x[27], D1y[27], D1z[27]; // D1(3, 26)
            DIPLOOP1(XI,D1x, D1y, D1z, &COORD11_C, &LOOPDIP1_C, &RHDR_C);
            double BX1=0.e0;
            double BY1=0.e0;
            double BZ1=0.e0;
            //DO 11 I=1,26
            for (int I=1; I<=26; I++) {
                BX1=BX1+C1[I]*D1x[I]; //   BX1,BY1,BZ1  ARE FIELD COMPONENTS
                BY1=BY1+C1[I]*D1y[I];  //  IN THE NORTHERN BOUNDARY POINT
                BZ1=BZ1+C1[I]*D1z[I];  //     11
            }

            double XAS2=R*ST02AS*cos(PAS);
            double Y2=  R*ST02AS*sin(PAS);
            double ZAS2=R*CT02AS;
            double X2=XAS2*CPSAS+ZAS2*SPSAS;
            double Z2=-XAS2*SPSAS+ZAS2*CPSAS; // X2,Y2,Z2 ARE COORDS OF THE SOUTHERN
            //BOUNDARY POINT
            XI[1]=X2;
            XI[2]=Y2;
            XI[3]=Z2;
            XI[4]=PS;
            double D2x[80], D2y[80], D2z[80]; // D1(3, 26)
            CONDIP1(XI,D2x, D2y, D2z, &DX1_C, &COORD21_C);
            double BX2=0.e0;
            double BY2=0.e0;
            double BZ2=0.e0;
            //DO 12 I=1,79
            for (int I=1; I<=79; I++) {
                BX2=BX2+C2[I]*D2x[I]; //  BX2,BY2,BZ2  ARE FIELD COMPONENTS
                BY2=BY2+C2[I]*D2y[I]; //  IN THE SOUTHERN BOUNDARY POINT
                BZ2=BZ2+C2[I]*D2z[I]; //      12
            }

            //  NOW INTERPOLATE:

            double SS=sqrt(pow((X2-X1), 2)+pow((Y2-Y1), 2)+pow((Z2-Z1), 2));
            double DS=sqrt(pow((X-X1), 2)+pow((Y-Y1), 2)+pow((Z-Z1), 2));
            double FRAC=DS/SS;
            *BX=BX1*(1.e0-FRAC)+BX2*FRAC;
            *BY=BY1*(1.e0-FRAC)+BY2*FRAC;
            *BZ=BZ1*(1.e0-FRAC)+BZ2*FRAC;
        }//ENDIF                                              ! END OF THE CASE 3

        //IF (LOC.EQ.4) THEN
        if (LOC==4) {
            //       print *, '  LOC=4 (south PSBL)'  !  (test printout; disabled now)
            double T01=TETR1S-DTET0;
            double T02=TETR1S+DTET0;
            double SQR=sqrt(R);
            double ST01AS=SQR/pow((R3+1.e0/pow(sin(T01), 6)-1.e0), 0.1666666667);
            double ST02AS=SQR/pow((R3+1.e0/pow(sin(T02), 6)-1.e0), 0.1666666667);
            double CT01AS=-sqrt(1.e0-ST01AS*ST01AS);
            double CT02AS=-sqrt(1.e0-ST02AS*ST02AS);
            double XAS1=R*ST01AS*cos(PAS);
            double Y1=  R*ST01AS*sin(PAS);
            double ZAS1=R*CT01AS;
            double X1=XAS1*CPSAS+ZAS1*SPSAS;
            double Z1=-XAS1*SPSAS+ZAS1*CPSAS; // X1,Y1,Z1 ARE COORDS OF THE NORTHERN
            //                                               BOUNDARY POINT
            XI[1]=X1;
            XI[2]=Y1;
            XI[3]=Z1;
            XI[4]=PS;
            double D2x[80], D2y[80], D2z[80]; // D1(3, 26)
            CONDIP1(XI,D2x, D2y, D2z, &DX1_C, &COORD21_C);
            double BX1=0.e0;
            double BY1=0.e0;
            double BZ1=0.e0;
            //DO 21 I=1,79
            for (int I=1; I<=79; I++) {
                BX1=BX1+C2[I]*D2x[I]; //  BX1,BY1,BZ1  ARE FIELD COMPONENTS
                BY1=BY1+C2[I]*D2y[I]; //  IN THE NORTHERN BOUNDARY POINT
                BZ1=BZ1+C2[I]*D2z[I]; //      21
            }

            double XAS2=R*ST02AS*cos(PAS);
            double Y2=  R*ST02AS*sin(PAS);
            double ZAS2=R*CT02AS;
            double X2=XAS2*CPSAS+ZAS2*SPSAS;
            double Z2=-XAS2*SPSAS+ZAS2*CPSAS; // X2,Y2,Z2 ARE COORDS OF THE SOUTHERN
            //                                          BOUNDARY POINT
            XI[1]=X2;
            XI[2]=Y2;
            XI[3]=Z2;
            XI[4]=PS;
            double D1x[27], D1y[27], D1z[27]; // D1(3, 26)
            DIPLOOP1(XI,D1x, D1y, D1z, &COORD11_C, &LOOPDIP1_C, &RHDR_C);
            double BX2=0.e0;
            double BY2=0.e0;
            double BZ2=0.e0;
            //DO 22 I=1,26
            for (int I=1; I<=26; I++) {
                BX2=BX2+C1[I]*D1x[I]; //  BX2,BY2,BZ2  ARE FIELD COMPONENTS
                BY2=BY2+C1[I]*D1y[I]; //     IN THE SOUTHERN BOUNDARY POINT
                BZ2=BZ2+C1[I]*D1z[I]; //         22
            }

            //  NOW INTERPOLATE:

            double SS=sqrt(pow((X2-X1), 2)+pow((Y2-Y1), 2)+pow((Z2-Z1), 2));
            double DS=sqrt(pow((X-X1), 2)+pow((Y-Y1), 2)+pow((Z-Z1), 2));
            double FRAC=DS/SS;
            *BX=BX1*(1.e0-FRAC)+BX2*FRAC;
            *BY=BY1*(1.e0-FRAC)+BY2*FRAC;
            *BZ=BZ1*(1.e0-FRAC)+BZ2*FRAC;

        }//ENDIF                                        ! END OF THE CASE 4

        //   NOW, LET US ADD THE SHIELDING FIELD

        double BSX, BSY, BSZ;
        BIRK1SHLD(PS,X,Y,Z,&BSX,&BSY,&BSZ);
        *BX=*BX+BSX;
        *BY=*BY+BSY;
        *BZ=*BZ+BSZ;

    }//END

    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     C
     C THIS CODE IS FOR THE FIELD FROM  2x2x2=8 "CARTESIAN" HARMONICS
     C
     SUBROUTINE  BIRK2SHL(X,Y,Z,PS,HX,HY,HZ)
     C
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C    The model parameters are provided to this module via common-block /A/.
     C  The 16 linear parameters enter in pairs in the amplitudes of the
     c       "cartesian" harmonics.
     c    The 8 nonlinear parameters are the scales Pi,Ri,Qi,and Si entering the
     c  arguments of exponents, sines, and cosines in each of the 8 "Cartesian"
     c   harmonics
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     */
    static void BIRK2SHL(double const X,double const Y,double const Z,double const PS,double *HX,double *HY,double *HZ)
    {
        //DIMENSION P(2),R(2),Q(2),S(2)
        double const *P, *R, *Q, *S;
        //DIMENSION A(24)

        //EQUIVALENCE(P(1),A(17)),(R(1),A(19)),(Q(1),A(21)),(S(1),A(23))
        /*DATA A/-111.6371348,124.5402702,110.3735178,-122.0095905,
         * 111.9448247,-129.1957743,-110.7586562,126.5649012,-0.7865034384,
         * -0.2483462721,0.8026023894,0.2531397188,10.72890902,0.8483902118,
         * -10.96884315,-0.8583297219,13.85650567,14.90554500,10.21914434,
         * 10.09021632,6.340382460,14.40432686,12.71023437,12.83966657/ */
        static double const A[25] = {0./*dummy*/, -111.6371348,124.5402702,110.3735178,-122.0095905,
            111.9448247,-129.1957743,-110.7586562,126.5649012,-0.7865034384,
            -0.2483462721,0.8026023894,0.2531397188,10.72890902,0.8483902118,
            -10.96884315,-0.8583297219,13.85650567,14.90554500,10.21914434,
            10.09021632,6.340382460,14.40432686,12.71023437,12.83966657};

        P = A + 16;
        R = A + 18;
        Q = A + 20;
        S = A + 22;

        double CPS=cos(PS);
        double SPS=sin(PS);
        double S3PS=4.e0*CPS*CPS-1.e0;   //  THIS IS SIN(3*PS)/SIN(PS)

        *HX=0.e0;
        *HY=0.e0;
        *HZ=0.e0;
        int L=0;

        /*DO 1 M=1,2     !    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
         C                           AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY) */
        for (int M=1; M<=2; M++) {
            //DO 2 I=1,2
            for (int I=1; I<=2; I++) {
                double CYPI=cos(Y/P[I]);
                double CYQI=cos(Y/Q[I]);
                double SYPI=sin(Y/P[I]);
                double SYQI=sin(Y/Q[I]);

                //DO 3 K=1,2
                for (int K=1; K<=2; K++) {
                    double SZRK=sin(Z/R[K]);
                    double CZSK=cos(Z/S[K]);
                    double CZRK=cos(Z/R[K]);
                    double SZSK=sin(Z/S[K]);
                    double SQPR=sqrt(1.e0/pow(P[I], 2)+1.e0/pow(R[K], 2));
                    double SQQS=sqrt(1.e0/pow(Q[I], 2)+1.e0/pow(S[K], 2));
                    double EPR=exp(X*SQPR);
                    double EQS=exp(X*SQQS);

                    /*DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
                     C                                  AND N=2 IS FOR THE SECOND ONE */
                    double DX = 0.0, DY = 0.0, DZ = 0.0;
                    for (int N=1; N<=2; N++) {
                        L=L+1;
                        //IF (M.EQ.1) THEN
                        if (M==1) {
                            //IF (N.EQ.1) THEN
                            if (N==1) {
                                DX=-SQPR*EPR*CYPI*SZRK;
                                DY=EPR/P[I]*SYPI*SZRK;
                                DZ=-EPR/R[K]*CYPI*CZRK;
                                *HX=*HX+A[L]*DX;
                                *HY=*HY+A[L]*DY;
                                *HZ=*HZ+A[L]*DZ;
                            } else {//ELSE
                                DX=DX*CPS;
                                DY=DY*CPS;
                                DZ=DZ*CPS;
                                *HX=*HX+A[L]*DX;
                                *HY=*HY+A[L]*DY;
                                *HZ=*HZ+A[L]*DZ;
                            }//ENDIF
                        } else {//ELSE
                            //IF (N.EQ.1) THEN
                            if (N==1) {
                                DX=-SPS*SQQS*EQS*CYQI*CZSK;
                                DY=SPS*EQS/Q[I]*SYQI*CZSK;
                                DZ=SPS*EQS/S[K]*CYQI*SZSK;
                                *HX=*HX+A[L]*DX;
                                *HY=*HY+A[L]*DY;
                                *HZ=*HZ+A[L]*DZ;
                            } else {//ELSE
                                DX=DX*S3PS;
                                DY=DY*S3PS;
                                DZ=DZ*S3PS;
                                *HX=*HX+A[L]*DX;
                                *HY=*HY+A[L]*DY;
                                *HZ=*HZ+A[L]*DZ;
                            }//ENDIF
                        }//ENDIF

                    }//4   CONTINUE
                }//3   CONTINUE
            }//2   CONTINUE
        }//1   CONTINUE

    }//END

    /*^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

     DOUBLE PRECISION FUNCTION XKSI(X,Y,Z)
     */
    static double XKSI(double const X,double const Y,double const Z)
    {
        double XKSI = 0.;
        //   A11 - C72, R0, and DR below  ARE STRETCH PARAMETERS (P.26-27, NB# 3),
        /*DATA A11A12,A21A22,A41A42,A51A52,A61A62,B11B12,B21B22,C61C62,
         *  C71C72,R0,DR /0.305662,-0.383593,0.2677733,-0.097656,-0.636034,
         *  -0.359862,0.424706,-0.126366,0.292578,1.21563,7.50937/ */
        static double const A11A12 = 0.305662;
        static double const A21A22 = -0.383593;
        static double const A41A42 = 0.2677733;
        static double const A51A52 = -0.097656;
        static double const A61A62 = -0.636034;
        static double const B11B12 = -0.359862;
        static double const B21B22 = 0.424706;
        static double const C61C62 = -0.126366;
        static double const C71C72 = 0.292578;
        static double const R0 = 1.21563;
        static double const DR = 7.50937;

        /*DATA TNOON,DTETA/0.3665191,0.09599309/ ! Correspond to noon and midnight
         C                                         latitudes 69 and 63.5 degs, resp. */
        static double const TNOON = 0.3665191,DTETA = 0.09599309; // Correspond to noon and midnight

        double DR2=DR*DR;

        double X2=X*X;
        double Y2=Y*Y;
        double Z2=Z*Z;
        //double XY=X*Y;
        //double XYZ=XY*Z;
        double R2=X2+Y2+Z2;
        double R=sqrt(R2);
        //double R3=R2*R;
        //double R4=R2*R2;
        double XR=X/R;
        double YR=Y/R;
        double ZR=Z/R;

        double PR;
        //IF (R.LT.R0) THEN
        if (R<R0) {
            PR=0.e0;
        } else {//ELSE
            PR=sqrt(pow((R-R0), 2)+DR2)-DR;
        }//ENDIF

        double F=X+PR*(A11A12+A21A22*XR+A41A42*XR*XR+A51A52*YR*YR+ A61A62*ZR*ZR);
        double G=Y+PR*(B11B12*YR+B21B22*XR*YR);
        double H=Z+PR*(C61C62*ZR+C71C72*XR*ZR);
        double G2=G*G;

        double FGH=F*F+G2+H*H;
        double FGH32=pow(sqrt(FGH), 3);
        double FCHSG2=F*F+G2;

        //IF (FCHSG2.LT.1.D-5) THEN
        if (FCHSG2<1.e-5) {
            XKSI=-1.e0;               //  THIS IS JUST FOR ELIMINATING PROBLEMS
            return XKSI;                    //  ON THE Z-AXIS
        }//ENDIF

        double SQFCHSG2=sqrt(FCHSG2);
        double ALPHA=FCHSG2/FGH32;
        double THETA=TNOON+0.5e0*DTETA*(1.e0-F/SQFCHSG2);
        double PHI=pow(sin(THETA), 2);

        XKSI=ALPHA-PHI;

        return XKSI;
    }//END

    /*--------------------------------------------------------------------
     C
     SUBROUTINE LOOPS4(X,Y,Z,BX,BY,BZ,XC,YC,ZC,R,THETA,PHI)
     C
     C   RETURNS FIELD COMPONENTS FROM A SYSTEM OF 4 CURRENT LOOPS, POSITIONED
     C     SYMMETRICALLY WITH RESPECT TO NOON-MIDNIGHT MERIDIAN AND EQUATORIAL
     C      PLANES.
     C  INPUT: X,Y,Z OF A POINT OF SPACE
     C        XC,YC,ZC (YC > 0 AND ZC > 0) - POSITION OF THE CENTER OF THE
     C                                         1ST-QUADRANT LOOP
     C        R - LOOP RADIUS (THE SAME FOR ALL FOUR)
     C        THETA, PHI  -  SPECIFY THE ORIENTATION OF THE NORMAL OF THE 1ST LOOP
     c      -----------------------------------------------------------
     */
    static void LOOPS4(double const X,double const Y,double const Z,double *BX,double *BY,double *BZ,double const XC,double const YC,double const ZC,double const R,double const THETA,double const PHI)
    {
        double CT=cos(THETA);
        double ST=sin(THETA);
        double CP=cos(PHI);
        double SP=sin(PHI);
        //------------------------------------1ST QUADRANT:
        double XS=(X-XC)*CP+(Y-YC)*SP;
        double YSS=(Y-YC)*CP-(X-XC)*SP;
        double ZS=Z-ZC;
        double XSS=XS*CT-ZS*ST;
        double ZSS=ZS*CT+XS*ST;

        double BXSS,BYS,BZSS;
        CIRCLE(XSS,YSS,ZSS,R,&BXSS,&BYS,&BZSS);
        double BXS=BXSS*CT+BZSS*ST;
        double BZ1=BZSS*CT-BXSS*ST;
        double BX1=BXS*CP-BYS*SP;
        double BY1=BXS*SP+BYS*CP;
        //-------------------------------------2nd QUADRANT:
        XS=(X-XC)*CP-(Y+YC)*SP;
        YSS=(Y+YC)*CP+(X-XC)*SP;
        ZS=Z-ZC;
        XSS=XS*CT-ZS*ST;
        ZSS=ZS*CT+XS*ST;

        CIRCLE(XSS,YSS,ZSS,R,&BXSS,&BYS,&BZSS);
        BXS=BXSS*CT+BZSS*ST;
        double BZ2=BZSS*CT-BXSS*ST;
        double BX2=BXS*CP+BYS*SP;
        double BY2=-BXS*SP+BYS*CP;
        //-------------------------------------3RD QUADRANT:
        XS=-(X-XC)*CP+(Y+YC)*SP;
        YSS=-(Y+YC)*CP-(X-XC)*SP;
        ZS=Z+ZC;
        XSS=XS*CT-ZS*ST;
        ZSS=ZS*CT+XS*ST;

        CIRCLE(XSS,YSS,ZSS,R,&BXSS,&BYS,&BZSS);
        BXS=BXSS*CT+BZSS*ST;
        double BZ3=BZSS*CT-BXSS*ST;
        double BX3=-BXS*CP-BYS*SP;
        double BY3=BXS*SP-BYS*CP;
        //-------------------------------------4TH QUADRANT:
        XS=-(X-XC)*CP-(Y-YC)*SP;
        YSS=-(Y-YC)*CP+(X-XC)*SP;
        ZS=Z+ZC;
        XSS=XS*CT-ZS*ST;
        ZSS=ZS*CT+XS*ST;

        CIRCLE(XSS,YSS,ZSS,R,&BXSS,&BYS,&BZSS);
        BXS=BXSS*CT+BZSS*ST;
        double BZ4=BZSS*CT-BXSS*ST;
        double BX4=-BXS*CP+BYS*SP;
        double BY4=-BXS*SP-BYS*CP;

        *BX=BX1+BX2+BX3+BX4;
        *BY=BY1+BY2+BY3+BY4;
        *BZ=BZ1+BZ2+BZ3+BZ4;

    }//END

    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     C
     SUBROUTINE R2OUTER (X,Y,Z,BX,BY,BZ)
     */
    static void R2OUTER (double const X,double const Y,double const Z,double *BX,double *BY,double *BZ)
    {
        //DATA PL1,PL2,PL3,PL4,PL5/-34.105,-2.00019,628.639,73.4847,12.5162/
        static double const PL1 = -34.105,PL2 = -2.00019,PL3 = 628.639,PL4 = 73.4847,PL5 = 12.5162;
        /*DATA PN1,PN2,PN3,PN4,PN5,PN6,PN7,PN8,PN9,PN10,PN11,PN12,PN13,PN14,
         *  PN15,PN16,PN17 /.55,.694,.0031,1.55,2.8,.1375,-.7,.2,.9625,
         * -2.994,2.925,-1.775,4.3,-.275,2.7,.4312,1.55/ */
        static double const PN1 = .55,PN2 = .694,PN3 = .0031,PN4 = 1.55,PN5 = 2.8,PN6 = .1375,PN7 = -.7,PN8 = .2,PN9 = .9625,PN10 = -2.994,PN11 = 2.925,PN12 = -1.775,PN13 = 4.3,PN14 = -.275, PN15 = 2.7,PN16 = .4312,PN17 = 1.55;

        //    THREE PAIRS OF CROSSED LOOPS:

        double DBX1,DBY1,DBZ1;
        CROSSLP(X,Y,Z,&DBX1,&DBY1,&DBZ1,PN1,PN2,PN3);
        double DBX2,DBY2,DBZ2;
        CROSSLP(X,Y,Z,&DBX2,&DBY2,&DBZ2,PN4,PN5,PN6);
        double DBX3,DBY3,DBZ3;
        CROSSLP(X,Y,Z,&DBX3,&DBY3,&DBZ3,PN7,PN8,PN9);

        //    NOW AN EQUATORIAL LOOP ON THE NIGHTSIDE

        double DBX4,DBY4,DBZ4;
        CIRCLE(X-PN10,Y,Z,PN11,&DBX4,&DBY4,&DBZ4);

        //   NOW A 4-LOOP SYSTEM ON THE NIGHTSIDE

        double DBX5,DBY5,DBZ5;
        LOOPS4(X,Y,Z,&DBX5,&DBY5,&DBZ5,PN12,PN13,PN14,PN15,PN16,PN17);

        //---------------------------------------------------------------------

        //                          NOW COMPUTE THE FIELD COMPONENTS:

        *BX=PL1*DBX1+PL2*DBX2+PL3*DBX3+PL4*DBX4+PL5*DBX5;
        *BY=PL1*DBY1+PL2*DBY2+PL3*DBY3+PL4*DBY4+PL5*DBY5;
        *BZ=PL1*DBZ1+PL2*DBZ2+PL3*DBZ3+PL4*DBZ4+PL5*DBZ5;

    }//END

    //--------------------------------------------------------------------
    // FUNCTION FEXP(S,A)
    static double FEXP(double const S,double const A)
    {
        double FEXP = 0.;
        //DATA E/2.718281828459D0/
        static double const E = 2.718281828459e0;

        //IF (A.LT.0.D0) FEXP=DSQRT(-2.D0*A*E)*S*DEXP(A*S*S)
        if (A<0.e0) FEXP=sqrt(-2.e0*A*E)*S*exp(A*S*S);
        //IF (A.GE.0.D0) FEXP=S*DEXP(A*(S*S-1.D0))
        if (A>=0.e0) FEXP=S*exp(A*(S*S-1.e0));

        return FEXP;
    }//END

    //-----------------------------------------------------------------------
    // FUNCTION FEXP1(S,A)
    static double FEXP1(double const S,double const A)
    {
        //IF (A.LE.0.D0) FEXP1=DEXP(A*S*S)
        if (A<=0.e0) return exp(A*S*S);
        //IF (A.GT.0.D0) FEXP1=DEXP(A*(S*S-1.D0))
        else return exp(A*(S*S-1.e0));
    }//END

    //************************************************************************
    // DOUBLE PRECISION FUNCTION TKSI(XKSI,XKS0,DXKSI)
    static double TKSI(double const XKSI,double const XKS0,double const DXKSI)
    {
        double TKSI = 0.;
        //SAVE M,TDZ3
        //DATA M/0/

        //IF (M.EQ.0) THEN
        double TDZ3=2.*pow(DXKSI, 3);
        //M=1
        //ENDIF

        //IF (XKSI-XKS0.LT.-DXKSI) TKSII=0.
        //IF (XKSI-XKS0.GE.DXKSI)  TKSII=1.
        double TKSII = 1.;
        if (XKSI-XKS0<-DXKSI) TKSII=0.;

        //IF (XKSI.GE.XKS0-DXKSI.AND.XKSI.LT.XKS0) THEN
        if (XKSI>=XKS0-DXKSI && XKSI<XKS0) {
            double BR3=pow((XKSI-XKS0+DXKSI), 3);
            TKSII=1.5*BR3/(TDZ3+BR3);
        }//ENDIF

        //IF (XKSI.GE.XKS0.AND.XKSI.LT.XKS0+DXKSI) THEN
        if (XKSI>=XKS0 && XKSI<XKS0+DXKSI) {
            double BR3=pow((XKSI-XKS0-DXKSI), 3);
            TKSII=1.+1.5*BR3/(TDZ3-BR3);
        }//ENDIF
        TKSI=TKSII;

        return TKSI;
    }//END


    /*@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
     C
     SUBROUTINE R2SHEET(X,Y,Z,BX,BY,BZ)
     */
    static void R2SHEET(double const X,double const Y,double const Z,double *BX,double *BY,double *BZ)
    {
        /*DATA PNONX1,PNONX2,PNONX3,PNONX4,PNONX5,PNONX6,PNONX7,PNONX8,
         *     PNONY1,PNONY2,PNONY3,PNONY4,PNONY5,PNONY6,PNONY7,PNONY8,
         *     PNONZ1,PNONZ2,PNONZ3,PNONZ4,PNONZ5,PNONZ6,PNONZ7,PNONZ8
         * /-19.0969D0,-9.28828D0,-0.129687D0,5.58594D0,22.5055D0,
         *  0.483750D-01,0.396953D-01,0.579023D-01,-13.6750D0,-6.70625D0,
         *  2.31875D0,11.4062D0,20.4562D0,0.478750D-01,0.363750D-01,
         * 0.567500D-01,-16.7125D0,-16.4625D0,-0.1625D0,5.1D0,23.7125D0,
         * 0.355625D-01,0.318750D-01,0.538750D-01/ */
        static double const PNONX1 = -19.0969e0;
        static double const PNONX2 = -9.28828e0;
        static double const PNONX3 = -0.129687e0;
        static double const PNONX4 = 5.58594e0;
        static double const PNONX5 = 22.5055e0;
        static double const PNONX6 = 0.483750e-01;
        static double const PNONX7 = 0.396953e-01;
        static double const PNONX8 = 0.579023e-01;
        static double const PNONY1 = -13.6750e0;
        static double const PNONY2 = -6.70625e0;
        static double const PNONY3 = 2.31875e0;
        static double const PNONY4 = 11.4062e0;
        static double const PNONY5 = 20.4562e0;
        static double const PNONY6 = 0.478750e-01;
        static double const PNONY7 = 0.363750e-01;
        static double const PNONY8 = 0.567500e-01;
        static double const PNONZ1 = -16.7125e0;
        static double const PNONZ2 = -16.4625e0;
        static double const PNONZ3 = -0.1625e0;
        static double const PNONZ4 = 5.1e0;
        static double const PNONZ5 = 23.7125e0;
        static double const PNONZ6 = 0.355625e-01;
        static double const PNONZ7 = 0.318750e-01;
        static double const PNONZ8 = 0.538750e-01;

        /*DATA A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,
         *  A18,A19,A20,A21,A22,A23,A24,A25,A26,A27,A28,A29,A30,A31,A32,A33,
         *  A34,A35,A36,A37,A38,A39,A40,A41,A42,A43,A44,A45,A46,A47,A48,A49,
         *  A50,A51,A52,A53,A54,A55,A56,A57,A58,A59,A60,A61,A62,A63,A64,A65,
         *  A66,A67,A68,A69,A70,A71,A72,A73,A74,A75,A76,A77,A78,A79,A80
         * /8.07190D0,-7.39582D0,-7.62341D0,0.684671D0,-13.5672D0,11.6681D0,
         * 13.1154,-0.890217D0,7.78726D0,-5.38346D0,-8.08738D0,0.609385D0,
         * -2.70410D0, 3.53741D0,3.15549D0,-1.11069D0,-8.47555D0,0.278122D0,
         *  2.73514D0,4.55625D0,13.1134D0,1.15848D0,-3.52648D0,-8.24698D0,
         * -6.85710D0,-2.81369D0, 2.03795D0, 4.64383D0,2.49309D0,-1.22041D0,
         * -1.67432D0,-0.422526D0,-5.39796D0,7.10326D0,5.53730D0,-13.1918D0,
         *  4.67853D0,-7.60329D0,-2.53066D0, 7.76338D0, 5.60165D0,5.34816D0,
         * -4.56441D0,7.05976D0,-2.62723D0,-0.529078D0,1.42019D0,-2.93919D0,
         *  55.6338D0,-1.55181D0,39.8311D0,-80.6561D0,-46.9655D0,32.8925D0,
         * -6.32296D0,19.7841D0,124.731D0,10.4347D0,-30.7581D0,102.680D0,
         * -47.4037D0,-3.31278D0,9.37141D0,-50.0268D0,-533.319D0,110.426D0,
         *  1000.20D0,-1051.40D0, 1619.48D0,589.855D0,-1462.73D0,1087.10D0,
         *  -1994.73D0,-1654.12D0,1263.33D0,-260.210D0,1424.84D0,1255.71D0,
         *  -956.733D0, 219.946D0/ */
        static double const A1 = 8.07190e0;
        static double const A2 = -7.39582e0;
        static double const A3 = -7.62341e0;
        static double const A4 = 0.684671e0;
        static double const A5 = -13.5672e0;
        static double const A6 = 11.6681e0;
        static double const A7 = 13.1154;
        static double const A8 = -0.890217e0;
        static double const A9 = 7.78726e0;
        static double const A10 = -5.38346e0;
        static double const A11 = -8.08738e0;
        static double const A12 = 0.609385e0;
        static double const A13 = -2.70410e0;
        static double const A14 = 3.53741e0;
        static double const A15 = 3.15549e0;
        static double const A16 = -1.11069e0;
        static double const A17 = -8.47555e0;
        static double const A18 = 0.278122e0;
        static double const A19 = 2.73514e0;
        static double const A20 = 4.55625e0;
        static double const A21 = 13.1134e0;
        static double const A22 = 1.15848e0;
        static double const A23 = -3.52648e0;
        static double const A24 = -8.24698e0;
        static double const A25 = -6.85710e0;
        static double const A26 = -2.81369e0;
        static double const A27 = 2.03795e0;
        static double const A28 =  4.64383e0;
        static double const A29 = 2.49309e0;
        static double const A30 = -1.22041e0;
        static double const A31 = -1.67432e0;
        static double const A32 = -0.422526e0;
        static double const A33 = -5.39796e0;
        static double const A34 = 7.10326e0;
        static double const A35 = 5.53730e0;
        static double const A36 = -13.1918e0;
        static double const A37 = 4.67853e0;
        static double const A38 = -7.60329e0;
        static double const A39 = -2.53066e0;
        static double const A40 =  7.76338e0;
        static double const A41 =  5.60165e0;
        static double const A42 = 5.34816e0;
        static double const A43 = -4.56441e0;
        static double const A44 = 7.05976e0;
        static double const A45 = -2.62723e0;
        static double const A46 = -0.529078e0;
        static double const A47 = 1.42019e0;
        static double const A48 = -2.93919e0;
        static double const A49 = 55.6338e0;
        static double const A50 = -1.55181e0;
        static double const A51 = 39.8311e0;
        static double const A52 = -80.6561e0;
        static double const A53 = -46.9655e0;
        static double const A54 = 32.8925e0;
        static double const A55 = -6.32296e0;
        static double const A56 = 19.7841e0;
        static double const A57 = 124.731e0;
        static double const A58 = 10.4347e0;
        static double const A59 = -30.7581e0;
        static double const A60 = 102.680e0;
        static double const A61 = -47.4037e0;
        static double const A62 = -3.31278e0;
        static double const A63 = 9.37141e0;
        static double const A64 = -50.0268e0;
        static double const A65 = -533.319e0;
        static double const A66 = 110.426e0;
        static double const A67 = 1000.20e0;
        static double const A68 = -1051.40e0;
        static double const A69 =  1619.48e0;
        static double const A70 = 589.855e0;
        static double const A71 = -1462.73e0;
        static double const A72 = 1087.10e0;
        static double const A73 = -1994.73e0;
        static double const A74 = -1654.12e0;
        static double const A75 = 1263.33e0;
        static double const A76 = -260.210e0;
        static double const A77 = 1424.84e0;
        static double const A78 = 1255.71e0;
        static double const A79 = -956.733e0;
        static double const A80 =  219.946e0;

        /*DATA B1,B2,B3,B4,B5,B6,B7,B8,B9,B10,B11,B12,B13,B14,B15,B16,B17,
         *  B18,B19,B20,B21,B22,B23,B24,B25,B26,B27,B28,B29,B30,B31,B32,B33,
         *  B34,B35,B36,B37,B38,B39,B40,B41,B42,B43,B44,B45,B46,B47,B48,B49,
         *  B50,B51,B52,B53,B54,B55,B56,B57,B58,B59,B60,B61,B62,B63,B64,B65,
         *  B66,B67,B68,B69,B70,B71,B72,B73,B74,B75,B76,B77,B78,B79,B80
         * /-9.08427D0,10.6777D0,10.3288D0,-0.969987D0,6.45257D0,-8.42508D0,
         * -7.97464D0,1.41996D0,-1.92490D0,3.93575D0,2.83283D0,-1.48621D0,
         *0.244033D0,-0.757941D0,-0.386557D0,0.344566D0,9.56674D0,-2.5365D0,
         * -3.32916D0,-5.86712D0,-6.19625D0,1.83879D0,2.52772D0,4.34417D0,
         * 1.87268D0,-2.13213D0,-1.69134D0,-.176379D0,-.261359D0,.566419D0,
         * 0.3138D0,-0.134699D0,-3.83086D0,-8.4154D0,4.77005D0,-9.31479D0,
         * 37.5715D0,19.3992D0,-17.9582D0,36.4604D0,-14.9993D0,-3.1442D0,
         * 6.17409D0,-15.5519D0,2.28621D0,-0.891549D-2,-.462912D0,2.47314D0,
         * 41.7555D0,208.614D0,-45.7861D0,-77.8687D0,239.357D0,-67.9226D0,
         * 66.8743D0,238.534D0,-112.136D0,16.2069D0,-40.4706D0,-134.328D0,
         * 21.56D0,-0.201725D0,2.21D0,32.5855D0,-108.217D0,-1005.98D0,
         * 585.753D0,323.668D0,-817.056D0,235.750D0,-560.965D0,-576.892D0,
         * 684.193D0,85.0275D0,168.394D0,477.776D0,-289.253D0,-123.216D0,
         * 75.6501D0,-178.605D0/ */
        static double const B1 = -9.08427e0;
        static double const B2 = 10.6777e0;
        static double const B3 = 10.3288e0;
        static double const B4 = -0.969987e0;
        static double const B5 = 6.45257e0;
        static double const B6 = -8.42508e0;
        static double const B7 = -7.97464e0;
        static double const B8 = 1.41996e0;
        static double const B9 = -1.92490e0;
        static double const B10 = 3.93575e0;
        static double const B11 = 2.83283e0;
        static double const B12 = -1.48621e0;
        static double const B13 = 0.244033e0;
        static double const B14 = -0.757941e0;
        static double const B15 = -0.386557e0;
        static double const B16 = 0.344566e0;
        static double const B17 = 9.56674e0;
        static double const B18 = -2.5365e0;
        static double const B19 = -3.32916e0;
        static double const B20 = -5.86712e0;
        static double const B21 = -6.19625e0;
        static double const B22 = 1.83879e0;
        static double const B23 = 2.52772e0;
        static double const B24 = 4.34417e0;
        static double const B25 = 1.87268e0;
        static double const B26 = -2.13213e0;
        static double const B27 = -1.69134e0;
        static double const B28 = -.176379e0;
        static double const B29 = -.261359e0;
        static double const B30 = .566419e0;
        static double const B31 = 0.3138e0;
        static double const B32 = -0.134699e0;
        static double const B33 = -3.83086e0;
        static double const B34 = -8.4154e0;
        static double const B35 = 4.77005e0;
        static double const B36 = -9.31479e0;
        static double const B37 = 37.5715e0;
        static double const B38 = 19.3992e0;
        static double const B39 = -17.9582e0;
        static double const B40 = 36.4604e0;
        static double const B41 = -14.9993e0;
        static double const B42 = -3.1442e0;
        static double const B43 = 6.17409e0;
        static double const B44 = -15.5519e0;
        static double const B45 = 2.28621e0;
        static double const B46 = -0.891549e-2;
        static double const B47 = -.462912e0;
        static double const B48 = 2.47314e0;
        static double const B49 = 41.7555e0;
        static double const B50 = 208.614e0;
        static double const B51 = -45.7861e0;
        static double const B52 = -77.8687e0;
        static double const B53 = 239.357e0;
        static double const B54 = -67.9226e0;
        static double const B55 = 66.8743e0;
        static double const B56 = 238.534e0;
        static double const B57 = -112.136e0;
        static double const B58 = 16.2069e0;
        static double const B59 = -40.4706e0;
        static double const B60 = -134.328e0;
        static double const B61 = 21.56e0;
        static double const B62 = -0.201725e0;
        static double const B63 = 2.21e0;
        static double const B64 = 32.5855e0;
        static double const B65 = -108.217e0;
        static double const B66 = -1005.98e0;
        static double const B67 = 585.753e0;
        static double const B68 = 323.668e0;
        static double const B69 = -817.056e0;
        static double const B70 = 235.750e0;
        static double const B71 = -560.965e0;
        static double const B72 = -576.892e0;
        static double const B73 = 684.193e0;
        static double const B74 = 85.0275e0;
        static double const B75 = 168.394e0;
        static double const B76 = 477.776e0;
        static double const B77 = -289.253e0;
        static double const B78 = -123.216e0;
        static double const B79 = 75.6501e0;
        static double const B80 = -178.605e0;

        /*DATA C1,C2,C3,C4,C5,C6,C7,C8,C9,C10,C11,C12,C13,C14,C15,C16,C17,
         *  C18,C19,C20,C21,C22,C23,C24,C25,C26,C27,C28,C29,C30,C31,C32,C33,
         *  C34,C35,C36,C37,C38,C39,C40,C41,C42,C43,C44,C45,C46,C47,C48,C49,
         *  C50,C51,C52,C53,C54,C55,C56,C57,C58,C59,C60,C61,C62,C63,C64,C65,
         *  C66,C67,C68,C69,C70,C71,C72,C73,C74,C75,C76,C77,C78,C79,C80
         * / 1167.61D0,-917.782D0,-1253.2D0,-274.128D0,-1538.75D0,1257.62D0,
         * 1745.07D0,113.479D0,393.326D0,-426.858D0,-641.1D0,190.833D0,
         * -29.9435D0,-1.04881D0,117.125D0,-25.7663D0,-1168.16D0,910.247D0,
         * 1239.31D0,289.515D0,1540.56D0,-1248.29D0,-1727.61D0,-131.785D0,
         * -394.577D0,426.163D0,637.422D0,-187.965D0,30.0348D0,0.221898D0,
         * -116.68D0,26.0291D0,12.6804D0,4.84091D0,1.18166D0,-2.75946D0,
         * -17.9822D0,-6.80357D0,-1.47134D0,3.02266D0,4.79648D0,0.665255D0,
         * -0.256229D0,-0.857282D-1,-0.588997D0,0.634812D-1,0.164303D0,
         * -0.15285D0,22.2524D0,-22.4376D0,-3.85595D0,6.07625D0,-105.959D0,
         * -41.6698D0,0.378615D0,1.55958D0,44.3981D0,18.8521D0,3.19466D0,
         *  5.89142D0,-8.63227D0,-2.36418D0,-1.027D0,-2.31515D0,1035.38D0,
         *  2040.66D0,-131.881D0,-744.533D0,-3274.93D0,-4845.61D0,482.438D0,
         * 1567.43D0,1354.02D0,2040.47D0,-151.653D0,-845.012D0,-111.723D0,
         * -265.343D0,-26.1171D0,216.632D0/ */
        static double const C1 = 1167.61e0;
        static double const C2 = -917.782e0;
        static double const C3 = -1253.2e0;
        static double const C4 = -274.128e0;
        static double const C5 = -1538.75e0;
        static double const C6 = 1257.62e0;
        static double const C7 = 1745.07e0;
        static double const C8 = 113.479e0;
        static double const C9 = 393.326e0;
        static double const C10 = -426.858e0;
        static double const C11 = -641.1e0;
        static double const C12 = 190.833e0;
        static double const C13 = -29.9435e0;
        static double const C14 = -1.04881e0;
        static double const C15 = 117.125e0;
        static double const C16 = -25.7663e0;
        static double const C17 = -1168.16e0;
        static double const C18 = 910.247e0;
        static double const C19 = 1239.31e0;
        static double const C20 = 289.515e0;
        static double const C21 = 1540.56e0;
        static double const C22 = -1248.29e0;
        static double const C23 = -1727.61e0;
        static double const C24 = -131.785e0;
        static double const C25 = -394.577e0;
        static double const C26 = 426.163e0;
        static double const C27 = 637.422e0;
        static double const C28 = -187.965e0;
        static double const C29 = 30.0348e0;
        static double const C30 = 0.221898e0;
        static double const C31 = -116.68e0;
        static double const C32 = 26.0291e0;
        static double const C33 = 12.6804e0;
        static double const C34 = 4.84091e0;
        static double const C35 = 1.18166e0;
        static double const C36 = -2.75946e0;
        static double const C37 = -17.9822e0;
        static double const C38 = -6.80357e0;
        static double const C39 = -1.47134e0;
        static double const C40 = 3.02266e0;
        static double const C41 = 4.79648e0;
        static double const C42 = 0.665255e0;
        static double const C43 = -0.256229e0;
        static double const C44 = -0.857282e-1;
        static double const C45 = -0.588997e0;
        static double const C46 = 0.634812e-1;
        static double const C47 = 0.164303e0;
        static double const C48 = -0.15285e0;
        static double const C49 = 22.2524e0;
        static double const C50 = -22.4376e0;
        static double const C51 = -3.85595e0;
        static double const C52 = 6.07625e0;
        static double const C53 = -105.959e0;
        static double const C54 = -41.6698e0;
        static double const C55 = 0.378615e0;
        static double const C56 = 1.55958e0;
        static double const C57 = 44.3981e0;
        static double const C58 = 18.8521e0;
        static double const C59 = 3.19466e0;
        static double const C60 = 5.89142e0;
        static double const C61 = -8.63227e0;
        static double const C62 = -2.36418e0;
        static double const C63 = -1.027e0;
        static double const C64 = -2.31515e0;
        static double const C65 = 1035.38e0;
        static double const C66 = 2040.66e0;
        static double const C67 = -131.881e0;
        static double const C68 = -744.533e0;
        static double const C69 = -3274.93e0;
        static double const C70 = -4845.61e0;
        static double const C71 = 482.438e0;
        static double const C72 = 1567.43e0;
        static double const C73 = 1354.02e0;
        static double const C74 = 2040.47e0;
        static double const C75 = -151.653e0;
        static double const C76 = -845.012e0;
        static double const C77 = -111.723e0;
        static double const C78 = -265.343e0;
        static double const C79 = -26.1171e0;
        static double const C80 = 216.632e0;

        //------------------------------------------------------------------

        double XKS=XKSI(X,Y,Z);    //  variation across the current sheet
        double T1X=XKS/sqrt(XKS*XKS+PNONX6*PNONX6);
        double T2X=pow(PNONX7, 3)/pow(sqrt(XKS*XKS+PNONX7*PNONX7), 3);
        double T3X=XKS/pow(sqrt(XKS*XKS+PNONX8*PNONX8), 5) *3.493856e0*pow(PNONX8, 4);

        double T1Y=XKS/sqrt(XKS*XKS+PNONY6*PNONY6);
        double T2Y=pow(PNONY7, 3)/pow(sqrt(XKS*XKS+PNONY7*PNONY7), 3);
        double T3Y=XKS/pow(sqrt(XKS*XKS+PNONY8*PNONY8), 5) *3.493856e0*pow(PNONY8, 4);

        double T1Z=XKS/sqrt(XKS*XKS+PNONZ6*PNONZ6);
        double T2Z=pow(PNONZ7, 3)/pow(sqrt(XKS*XKS+PNONZ7*PNONZ7), 3);
        double T3Z=XKS/pow(sqrt(XKS*XKS+PNONZ8*PNONZ8), 5) *3.493856e0*pow(PNONZ8, 4);

        double RHO2=X*X+Y*Y;
        double R=sqrt(RHO2+Z*Z);
        double RHO=sqrt(RHO2);

        double C1P=X/RHO;
        double S1P=Y/RHO;
        double S2P=2.e0*S1P*C1P;
        double C2P=C1P*C1P-S1P*S1P;
        double S3P=S2P*C1P+C2P*S1P;
        double C3P=C2P*C1P-S2P*S1P;
        double S4P=S3P*C1P+C3P*S1P;
        double CT=Z/R;
        //double ST=RHO/R;

        double S1=FEXP(CT,PNONX1);
        double S2=FEXP(CT,PNONX2);
        double S3=FEXP(CT,PNONX3);
        double S4=FEXP(CT,PNONX4);
        double S5=FEXP(CT,PNONX5);

        //                   NOW COMPUTE THE GSM FIELD COMPONENTS:
        *BX=S1*((A1+A2*T1X+A3*T2X+A4*T3X)
                +C1P*(A5+A6*T1X+A7*T2X+A8*T3X)
                +C2P*(A9+A10*T1X+A11*T2X+A12*T3X)
                +C3P*(A13+A14*T1X+A15*T2X+A16*T3X))
        +S2*((A17+A18*T1X+A19*T2X+A20*T3X)
             +C1P*(A21+A22*T1X+A23*T2X+A24*T3X)
             +C2P*(A25+A26*T1X+A27*T2X+A28*T3X)
             +C3P*(A29+A30*T1X+A31*T2X+A32*T3X))
        +S3*((A33+A34*T1X+A35*T2X+A36*T3X)
             +C1P*(A37+A38*T1X+A39*T2X+A40*T3X)
             +C2P*(A41+A42*T1X+A43*T2X+A44*T3X)
             +C3P*(A45+A46*T1X+A47*T2X+A48*T3X))
        +S4*((A49+A50*T1X+A51*T2X+A52*T3X)
             +C1P*(A53+A54*T1X+A55*T2X+A56*T3X)
             +C2P*(A57+A58*T1X+A59*T2X+A60*T3X)
             +C3P*(A61+A62*T1X+A63*T2X+A64*T3X))
        +S5*((A65+A66*T1X+A67*T2X+A68*T3X)
             +C1P*(A69+A70*T1X+A71*T2X+A72*T3X)
             +C2P*(A73+A74*T1X+A75*T2X+A76*T3X)
             +C3P*(A77+A78*T1X+A79*T2X+A80*T3X));

        S1=FEXP(CT,PNONY1);
        S2=FEXP(CT,PNONY2);
        S3=FEXP(CT,PNONY3);
        S4=FEXP(CT,PNONY4);
        S5=FEXP(CT,PNONY5);

        *BY=S1*(S1P*(B1+B2*T1Y+B3*T2Y+B4*T3Y)
                +S2P*(B5+B6*T1Y+B7*T2Y+B8*T3Y)
                +S3P*(B9+B10*T1Y+B11*T2Y+B12*T3Y)
                +S4P*(B13+B14*T1Y+B15*T2Y+B16*T3Y))
        +S2*(S1P*(B17+B18*T1Y+B19*T2Y+B20*T3Y)
             +S2P*(B21+B22*T1Y+B23*T2Y+B24*T3Y)
             +S3P*(B25+B26*T1Y+B27*T2Y+B28*T3Y)
             +S4P*(B29+B30*T1Y+B31*T2Y+B32*T3Y))
        +S3*(S1P*(B33+B34*T1Y+B35*T2Y+B36*T3Y)
             +S2P*(B37+B38*T1Y+B39*T2Y+B40*T3Y)
             +S3P*(B41+B42*T1Y+B43*T2Y+B44*T3Y)
             +S4P*(B45+B46*T1Y+B47*T2Y+B48*T3Y))
        +S4*(S1P*(B49+B50*T1Y+B51*T2Y+B52*T3Y)
             +S2P*(B53+B54*T1Y+B55*T2Y+B56*T3Y)
             +S3P*(B57+B58*T1Y+B59*T2Y+B60*T3Y)
             +S4P*(B61+B62*T1Y+B63*T2Y+B64*T3Y))
        +S5*(S1P*(B65+B66*T1Y+B67*T2Y+B68*T3Y)
             +S2P*(B69+B70*T1Y+B71*T2Y+B72*T3Y)
             +S3P*(B73+B74*T1Y+B75*T2Y+B76*T3Y)
             +S4P*(B77+B78*T1Y+B79*T2Y+B80*T3Y));

        S1=FEXP1(CT,PNONZ1);
        S2=FEXP1(CT,PNONZ2);
        S3=FEXP1(CT,PNONZ3);
        S4=FEXP1(CT,PNONZ4);
        S5=FEXP1(CT,PNONZ5);

        *BZ=S1*((C1+C2*T1Z+C3*T2Z+C4*T3Z)
                +C1P*(C5+C6*T1Z+C7*T2Z+C8*T3Z)
                +C2P*(C9+C10*T1Z+C11*T2Z+C12*T3Z)
                +C3P*(C13+C14*T1Z+C15*T2Z+C16*T3Z))
        +S2*((C17+C18*T1Z+C19*T2Z+C20*T3Z)
             +C1P*(C21+C22*T1Z+C23*T2Z+C24*T3Z)
             +C2P*(C25+C26*T1Z+C27*T2Z+C28*T3Z)
             +C3P*(C29+C30*T1Z+C31*T2Z+C32*T3Z))
        +S3*((C33+C34*T1Z+C35*T2Z+C36*T3Z)
             +C1P*(C37+C38*T1Z+C39*T2Z+C40*T3Z)
             +C2P*(C41+C42*T1Z+C43*T2Z+C44*T3Z)
             +C3P*(C45+C46*T1Z+C47*T2Z+C48*T3Z))
        +S4*((C49+C50*T1Z+C51*T2Z+C52*T3Z)
             +C1P*(C53+C54*T1Z+C55*T2Z+C56*T3Z)
             +C2P*(C57+C58*T1Z+C59*T2Z+C60*T3Z)
             +C3P*(C61+C62*T1Z+C63*T2Z+C64*T3Z))
        +S5*((C65+C66*T1Z+C67*T2Z+C68*T3Z)
             +C1P*(C69+C70*T1Z+C71*T2Z+C72*T3Z)
             +C2P*(C73+C74*T1Z+C75*T2Z+C76*T3Z)
             +C3P*(C77+C78*T1Z+C79*T2Z+C80*T3Z));

    }//END

    /*-----------------------------------------------------------------------

     SUBROUTINE BCONIC(X,Y,Z,CBX,CBY,CBZ,NMAX)
     C
     c   "CONICAL" HARMONICS
     */
    static void BCONIC(double const X,double const Y,double const Z,double CBX[],double CBY[],double CBZ[],int const NMAX)
    {
        //DIMENSION CBX(NMAX),CBY(NMAX),CBZ(NMAX)

        double RO2=X*X+Y*Y;
        double RO=sqrt(RO2);

        double CF=X/RO;
        double SF=Y/RO;
        double CFM1=1.e0;
        double SFM1=0.e0;

        double R2=RO2+Z*Z;
        double R=sqrt(R2);
        double C=Z/R;
        double S=RO/R;
        double CH=sqrt(0.5e0*(1.e0+C));
        double SH=sqrt(0.5e0*(1.e0-C));
        double TNHM1=1.e0;
        double CNHM1=1.e0;
        double TNH=SH/CH;
        double CNH=1.e0/TNH;

        //DO 1 M=1,NMAX
        for (int M=1; M<=NMAX; M++) {
            double CFM=CFM1*CF-SFM1*SF;
            double SFM=CFM1*SF+SFM1*CF;
            CFM1=CFM;
            SFM1=SFM;
            double TNHM=TNHM1*TNH;
            double CNHM=CNHM1*CNH;
            double BT=M*CFM/(R*S)*(TNHM+CNHM);
            double BF=-0.5e0*M*SFM/R*(TNHM1/(CH*CH)-CNHM1/(SH*SH));
            TNHM1=TNHM;
            CNHM1=CNHM;
            CBX[M]=BT*C*CF-BF*SF;
            CBY[M]=BT*C*SF+BF*CF;
            CBZ[M]=-BT*S; //         1
        }

    }//END

    /*-------------------------------------------------------------------
     C
     SUBROUTINE DIPDISTR(X,Y,Z,BX,BY,BZ,MODE)
     C
     C   RETURNS FIELD COMPONENTS FROM A LINEAR DISTRIBUTION OF DIPOLAR SOURCES
     C     ON THE Z-AXIS.  THE PARAMETER MODE DEFINES HOW THE DIPOLE STRENGTH
     C     VARIES ALONG THE Z-AXIS:  MODE=0 IS FOR A STEP-FUNCTION (Mx=const > 0
     c         FOR Z > 0, AND Mx=-const < 0 FOR Z < 0)
     C      WHILE MODE=1 IS FOR A LINEAR VARIATION OF THE DIPOLE MOMENT DENSITY
     C       SEE NB#3, PAGE 53 FOR DETAILS.
     C
     C
     C INPUT: X,Y,Z OF A POINT OF SPACE, AND MODE
     */
    static void DIPDISTR(double const X,double const Y,double const Z,double *BX,double *BY,double *BZ,int MODE)
    {
        double X2=X*X;
        double RHO2=X2+Y*Y;
        double R2=RHO2+Z*Z;
        double R3=R2*sqrt(R2);

        //IF (MODE.EQ.0) THEN
        if (MODE==0) {
            *BX=Z/(RHO2*RHO2)*(R2*(Y*Y-X2)-RHO2*X2)/R3;
            *BY=-X*Y*Z/(RHO2*RHO2)*(2.e0*R2+RHO2)/R3;
            *BZ=X/R3;
        } else {//ELSE
            *BX=Z/(RHO2*RHO2)*(Y*Y-X2);
            *BY=-2.e0*X*Y*Z/(RHO2*RHO2);
            *BZ=X/RHO2;
        }//ENDIF

    }//END

    /*****************************************************************

     c
     SUBROUTINE R2INNER (X,Y,Z,BX,BY,BZ)
     C
     */
    static void R2INNER (double const X,double const Y,double const Z,double *BX,double *BY,double *BZ)
    {
        //DIMENSION CBX(5),CBY(5),CBZ(5)
        double CBX[6],CBY[6],CBZ[6];

        /*DATA PL1,PL2,PL3,PL4,PL5,PL6,PL7,PL8/154.185,-2.12446,.601735E-01,
         * -.153954E-02,.355077E-04,29.9996,262.886,99.9132/ */
        static double const PL1 = 154.185,PL2 = -2.12446,PL3 = .601735e-01,PL4 = -.153954e-02,PL5 = .355077e-04,PL6 = 29.9996,PL7 = 262.886,PL8 = 99.9132;
        /*DATA PN1,PN2,PN3,PN4,PN5,PN6,PN7,PN8/-8.1902,6.5239,5.504,7.7815,
         * .8573,3.0986,.0774,-.038/ */
        static double const PN1 = -8.1902,PN2 = 6.5239,PN3 = 5.504,PN4 = 7.7815,PN5 = .8573,PN6 = 3.0986,PN7 = .0774,PN8 = -.038;

        BCONIC(X,Y,Z,CBX,CBY,CBZ,5);

        //   NOW INTRODUCE  ONE  4-LOOP SYSTEM:

        double DBX8,DBY8,DBZ8;
        LOOPS4(X,Y,Z,&DBX8,&DBY8,&DBZ8,PN1,PN2,PN3,PN4,PN5,PN6);

        double DBX6,DBY6,DBZ6;
        DIPDISTR(X-PN7,Y,Z,&DBX6,&DBY6,&DBZ6,0);
        double DBX7,DBY7,DBZ7;
        DIPDISTR(X-PN8,Y,Z,&DBX7,&DBY7,&DBZ7,1);

        //                           NOW COMPUTE THE FIELD COMPONENTS:

        *BX=PL1*CBX[1]+PL2*CBX[2]+PL3*CBX[3]+PL4*CBX[4]+PL5*CBX[5] +PL6*DBX6+PL7*DBX7+PL8*DBX8;
        *BY=PL1*CBY[1]+PL2*CBY[2]+PL3*CBY[3]+PL4*CBY[4]+PL5*CBY[5] +PL6*DBY6+PL7*DBY7+PL8*DBY8;
        *BZ=PL1*CBZ[1]+PL2*CBZ[2]+PL3*CBZ[3]+PL4*CBZ[4]+PL5*CBZ[5] +PL6*DBZ6+PL7*DBZ7+PL8*DBZ8;

    }//END

    /********************************************************************
     C
     SUBROUTINE R2_BIRK(X,Y,Z,PS,BX,BY,BZ)
     C
     C  RETURNS THE MODEL FIELD FOR THE REGION 2 BIRKELAND CURRENT/PARTIAL RC
     C    (WITHOUT SHIELDING FIELD)
     */
    static void R2_BIRK(double const X,double const Y,double const Z,double const PS,double *BX,double *BY,double *BZ)
    {
        //SAVE PSI,CPS,SPS
        //DATA DELARG/0.030D0/,DELARG1/0.015D0/,PSI/10.D0/
        static double const DELARG = 0.030e0,DELARG1 = 0.015e0;

        //IF (DABS(PSI-PS).GT.1.D-10) THEN
        //PSI=PS
        double CPS=cos(PS);
        double SPS=sin(PS);
        //ENDIF

        double XSM=X*CPS-Z*SPS;
        double ZSM=Z*CPS+X*SPS;

        double XKS=XKSI(XSM,Y,ZSM);
        double BXSM,BZSM;
        //IF (XKS.LT.-(DELARG+DELARG1)) THEN
        if (XKS<-(DELARG+DELARG1)) {
            R2OUTER(XSM,Y,ZSM,&BXSM,BY,&BZSM);
            BXSM=-BXSM*0.02;      //  ALL COMPONENTS ARE MULTIPLIED BY THE
            *BY=-*BY*0.02;          //  FACTOR -0.02, IN ORDER TO NORMALIZE THE
            BZSM=-BZSM*0.02;      //  FIELD (SO THAT Bz=-1 nT at X=-5.3 RE, Y=Z=0)
        }//ENDIF

        //IF (XKS.GE.-(DELARG+DELARG1).AND.XKS.LT.-DELARG+DELARG1) THEN
        if (XKS>=-(DELARG+DELARG1) && XKS<-DELARG+DELARG1) {
            double BXSM1,BY1,BZSM1;
            R2OUTER(XSM,Y,ZSM,&BXSM1,&BY1,&BZSM1);
            double BXSM2,BY2,BZSM2;
            R2SHEET(XSM,Y,ZSM,&BXSM2,&BY2,&BZSM2);
            double F2=-0.02*TKSI(XKS,-DELARG,DELARG1);
            double F1=-0.02-F2;
            BXSM=BXSM1*F1+BXSM2*F2;
            *BY=BY1*F1+BY2*F2;
            BZSM=BZSM1*F1+BZSM2*F2;
        }//ENDIF

        //IF (XKS.GE.-DELARG+DELARG1.AND.XKS.LT.DELARG-DELARG1) THEN
        if (XKS>=-DELARG+DELARG1 && XKS<DELARG-DELARG1) {
            R2SHEET(XSM,Y,ZSM,&BXSM,BY,&BZSM);
            BXSM=-BXSM*0.02;
            *BY=-*BY*0.02;
            BZSM=-BZSM*0.02;
        }//ENDIF
        //IF (XKS.GE.DELARG-DELARG1.AND.XKS.LT.DELARG+DELARG1) THEN
        if (XKS>=DELARG-DELARG1 && XKS<DELARG+DELARG1) {
            double BXSM1,BY1,BZSM1;
            R2INNER(XSM,Y,ZSM,&BXSM1,&BY1,&BZSM1);
            double BXSM2,BY2,BZSM2;
            R2SHEET(XSM,Y,ZSM,&BXSM2,&BY2,&BZSM2);
            double F1=-0.02*TKSI(XKS,DELARG,DELARG1);
            double F2=-0.02-F1;
            BXSM=BXSM1*F1+BXSM2*F2;
            *BY=BY1*F1+BY2*F2;
            BZSM=BZSM1*F1+BZSM2*F2;
        }//ENDIF
        //IF (XKS.GE.DELARG+DELARG1) THEN
        if (XKS>=DELARG+DELARG1) {
            R2INNER(XSM,Y,ZSM,&BXSM,BY,&BZSM);
            BXSM=-BXSM*0.02;
            *BY=-*BY*0.02;
            BZSM=-BZSM*0.02;
        }//ENDIF

        *BX=BXSM*CPS+BZSM*SPS;
        *BZ=BZSM*CPS-BXSM*SPS;

    }//END

    /*##########################################################################
     C
     SUBROUTINE BIRK2TOT_02(PS,X,Y,Z,BX,BY,BZ)
     */
    static void BIRK2TOT_02(double const PS,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ)
    {
        double WX,WY,WZ;
        BIRK2SHL(X,Y,Z,PS,&WX,&WY,&WZ);
        double HX,HY,HZ;
        R2_BIRK(X,Y,Z,PS,&HX,&HY,&HZ);
        *BX=WX+HX;
        *BY=WY+HY;
        *BZ=WZ+HZ;

    }//END

    /*~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
     c
     SUBROUTINE DIPOLE(PS,X,Y,Z,BX,BY,BZ)
     C
     C  CALCULATES GSM COMPONENTS OF GEODIPOLE FIELD WITH THE DIPOLE MOMENT
     C  CORRESPONDING TO THE EPOCH OF 1980.
     C------------INPUT PARAMETERS:
     C   PS - GEODIPOLE TILT ANGLE IN RADIANS, X,Y,Z - GSM COORDINATES IN RE
     C------------OUTPUT PARAMETERS:
     C   BX,BY,BZ - FIELD COMPONENTS IN GSM SYSTEM, IN NANOTESLA.
     C
     C
     C     WRITEN BY: N. A. TSYGANENKO
     */
    static void DIPOLE(double const PS,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ)
    {
        //DATA M,PSI/0,5./
        //double const PSI = 5.;
        //SAVE M,PSI,SPS,CPS
        //IF(M.EQ.1.AND.ABS(PS-PSI).LT.1.E-5) GOTO 1
        double SPS=sin(PS);
        double CPS=cos(PS);
        //double PSI=PS;
        //M=1

        double P=X*X; //        1
        double U=Z*Z;
        double V=3.*Z*X;
        double T=Y*Y;
        double Q=30574./pow(sqrt(P+T+U), 5);
        *BX=Q*((T+U-2.*P)*SPS-V*CPS);
        *BY=-3.*Y*Q*(X*SPS+Z*CPS);
        *BZ=Q*((P+T-2.*U)*CPS-V*SPS);

    }//END

    /*%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
     C
     C  THIS CODE YIELDS THE SHIELDING FIELD FOR THE PARALLEL DIPOLE
     C
     SUBROUTINE  CYLHAR1(A, X,Y,Z, BX,BY,BZ)
     C
     C
     C   ***  N.A. Tsyganenko ***  Sept. 14-18, 1993; revised March 16, 1994 ***
     C
     C   An approximation of the Chapman-Ferraro field by a sum of 6 cylin-
     c   drical harmonics (see pages 97-113 in the brown GSFC notebook #1)
     c
     C      Description of parameters:
     C
     C  A   - input vector containing model parameters;
     C  X,Y,Z - input GSM coordinates,
     C  BX,BY,BZ - output GSM components of the shielding field
     C
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C
     C      The 6 linear parameters A(1)-A(6) are amplitudes of the cylindrical
     c  harmonic terms.
     c      The 6 nonlinear parameters A(7)-A(12) are the corresponding scale
     c  lengths for each term (see GSFC brown notebook).
     c
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     */
    static void CYLHAR1(double const A[], double const X, double const Y, double const Z, double *BX, double *BY,double *BZ)
    {
        double RHO=sqrt(Y*Y+Z*Z);
        double SINFI, COSFI;
        //IF (RHO.LT.1.D-10) THEN
        if (RHO<1.e-10) {
            SINFI=1.e0;
            COSFI=0.e0;
            //GOTO 1
        } else { //ENDIF
            SINFI=Z/RHO;
            COSFI=Y/RHO;
        }
        
        *BX=0.e0;   //  1
        *BY=0.e0;
        *BZ=0.e0;
        
        //DO 11 I=1,3
        for (int I=1; I<=3; I++) {
            double DZETA=RHO/A[I+6];
            double XKSI=X/A[I+6];
            double XJ0=BES(DZETA,0);
            double XJ1=BES(DZETA,1);
            double XEXP=exp(XKSI);
            double BRHO=XJ1*XEXP;
            *BX=*BX-A[I]*XJ0*XEXP;
            *BY=*BY+A[I]*BRHO*COSFI;
            *BZ=*BZ+A[I]*BRHO*SINFI;
            //11        CONTINUE
        }
        
        //DO 12 I=4,6
        for (int I=4; I<=6; I++) {
            double DZETA=RHO/A[I+6];
            double XKSI=X/A[I+6];
            double XJ0=BES(DZETA,0);
            double XJ1=BES(DZETA,1);
            double XEXP=exp(XKSI);
            double BRHO=(DZETA*XJ0+XKSI*XJ1)*XEXP;
            *BX=*BX+A[I]*(DZETA*XJ1-XJ0*(XKSI+1.e0))*XEXP;
            *BY=*BY+A[I]*BRHO*COSFI;
            *BZ=*BZ+A[I]*BRHO*SINFI;
            //12        CONTINUE
        }
        
    }//END
    
    /*+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
     C
     C  THIS CODE YIELDS THE SHIELDING FIELD FOR THE PERPENDICULAR DIPOLE
     C
     SUBROUTINE  CYLHARM( A, X,Y,Z, BX,BY,BZ)
     C
     C
     C   ***  N.A. Tsyganenko ***  Sept. 14-18, 1993; revised March 16, 1994 ***
     C
     C   An approximation for the Chapman-Ferraro field by a sum of 6 cylin-
     c   drical harmonics (see pp. 97-113 in the brown GSFC notebook #1)
     c
     C      Description of parameters:
     C
     C  A   - input vector containing model parameters;
     C  X,Y,Z   -  input GSM coordinates
     C  BX,BY,BZ - output GSM components of the shielding field
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C  The 6 linear parameters A(1)-A(6) are amplitudes of the cylindrical harmonic
     c       terms.
     c  The 6 nonlinear parameters A(7)-A(12) are the corresponding scale lengths
     C       for each term (see GSFC brown notebook).
     c
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     */
    static void CYLHARM(double const A[], double const X, double const Y, double const Z, double *BX,double *BY,double *BZ)
    {
        double RHO=sqrt(Y*Y+Z*Z);
        double SINFI, COSFI;
        //IF (RHO.LT.1.D-8) THEN
        if (RHO<1.e-8) {
            SINFI=1.e0;
            COSFI=0.e0;
            RHO=1.e-8;
            //GOTO 1
        } else { //ENDIF
            SINFI=Z/RHO;
            COSFI=Y/RHO;
        }
        double SINFI2=SINFI*SINFI; //    1
        double SI2CO2=SINFI2-COSFI*COSFI;
        
        *BX=0.e0;
        *BY=0.e0;
        *BZ=0.e0;
        
        //DO 11 I=1,3
        for (int I=1; I<=3; I++) {
            double DZETA=RHO/A[I+6];
            double XJ0=BES(DZETA,0);
            double XJ1=BES(DZETA,1);
            double XEXP=exp(X/A[I+6]);
            *BX=*BX-A[I]*XJ1*XEXP*SINFI;
            *BY=*BY+A[I]*(2.e0*XJ1/DZETA-XJ0)*XEXP*SINFI*COSFI;
            *BZ=*BZ+A[I]*(XJ1/DZETA*SI2CO2-XJ0*SINFI2)*XEXP;
            //11        CONTINUE
        }
        
        //DO 12 I=4,6
        for (int I=4; I<=6; I++) {
            double DZETA=RHO/A[I+6];
            double XKSI=X/A[I+6];
            double XJ0=BES(DZETA,0);
            double XJ1=BES(DZETA,1);
            double XEXP=exp(XKSI);
            double BRHO=(XKSI*XJ0-(DZETA*DZETA+XKSI-1.e0)*XJ1/DZETA)*XEXP*SINFI;
            double BPHI=(XJ0+XJ1/DZETA*(XKSI-1.e0))*XEXP*COSFI;
            *BX=*BX+A[I]*(DZETA*XJ0+XKSI*XJ1)*XEXP*SINFI;
            *BY=*BY+A[I]*(BRHO*COSFI-BPHI*SINFI);
            *BZ=*BZ+A[I]*(BRHO*SINFI+BPHI*COSFI);
            //12        CONTINUE
        }
        
    }//END
    
    /*=====================================================================
     
     SUBROUTINE DIPSHLD(PS,X,Y,Z,BX,BY,BZ)
     C
     C   CALCULATES GSM COMPONENTS OF THE EXTERNAL MAGNETIC FIELD DUE TO
     C    SHIELDING OF THE EARTH'S DIPOLE ONLY
     */
    static void DIPSHLD(double const PS, double const X, double const Y, double const Z, double *BX, double *BY, double *BZ)
    {
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION A1(12),A2(12)
        //DATA A1 /.24777,-27.003,-.46815,7.0637,-1.5918,-.90317E-01,57.522,
        //* 13.757,2.0100,10.458,4.5798,2.1695/
        static double const A1[13] = {0./*dummy*/, .24777,-27.003,-.46815,7.0637,-1.5918,-.90317E-01,57.522, 13.757,2.0100,10.458,4.5798,2.1695};
        //DATA A2/-.65385,-18.061,-.40457,-5.0995,1.2846,.78231E-01,39.592,
        //* 13.291,1.9970,10.062,4.5140,2.1558/
        static double const A2[13] = {0./*dummy*/, -.65385,-18.061,-.40457,-5.0995,1.2846,.78231E-01,39.592,13.291,1.9970,10.062,4.5140,2.1558};
        
        double CPS=cos(PS);
        double SPS=sin(PS);
        double HX, HY, HZ;
        CYLHARM(A1,X,Y,Z,&HX,&HY,&HZ);
        double FX, FY, FZ;
        CYLHAR1(A2,X,Y,Z,&FX,&FY,&FZ);
        
        *BX=HX*CPS+FX*SPS;
        *BY=HY*CPS+FY*SPS;
        *BZ=HZ*CPS+FZ*SPS;
        
    }//END
    
    /*
     C----------------------------------------------------------------------
     c
     SUBROUTINE T96_01 (IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
     C
     c     RELEASE DATE OF THIS VERSION:   JUNE 22, 1996.
     C     LAST UPDATE: MAY 01, 2006:  IN THE S/R DIPOLE, SPS AND CPS WERE ADDED IN THE SAVE STATEMENT
     
     C----------------------------------------------------------------------
     C
     C  WITH TWO CORRECTIONS, SUGGESTED BY T.SOTIRELIS' COMMENTS (APR.7, 1997)
     C
     C  (1) A "STRAY "  CLOSING PARENTHESIS WAS REMOVED IN THE S/R   R2_BIRK
     C  (2) A 0/0 PROBLEM ON THE Z-AXIS WAS SIDESTEPPED (LINES 44-46 OF THE
     c       DOUBLE PRECISION FUNCTION XKSI)
     c--------------------------------------------------------------------
     C DATA-BASED MODEL CALIBRATED BY (1) SOLAR WIND PRESSURE PDYN (NANOPASCALS),
     C           (2) DST (NANOTESLA),  (3) BYIMF, AND (4) BZIMF (NANOTESLA).
     c THESE INPUT PARAMETERS SHOULD BE PLACED IN THE FIRST 4 ELEMENTS
     c OF THE ARRAY PARMOD(10).
     C
     C   THE REST OF THE INPUT VARIABLES ARE: THE GEODIPOLE TILT ANGLE PS (RADIANS),
     C AND   X,Y,Z -  GSM POSITION (RE)
     C
     c   IOPT  IS JUST A DUMMY INPUT PARAMETER, NECESSARY TO MAKE THIS SUBROUTINE
     C COMPATIBLE WITH THE NEW RELEASE (APRIL 1996) OF THE TRACING SOFTWARE
     C PACKAGE (GEOPACK). IOPT VALUE DOES NOT AFFECT THE OUTPUT FIELD.
     c
     C
     c OUTPUT:  GSM COMPONENTS OF THE EXTERNAL MAGNETIC FIELD (BX,BY,BZ, nanotesla)
     C            COMPUTED AS A SUM OF CONTRIBUTIONS FROM PRINCIPAL FIELD SOURCES
     C
     c  (C) Copr. 1995, 1996, Nikolai A. Tsyganenko, Hughes STX, Code 695, NASA GSFC
     c      Greenbelt, MD 20771, USA
     c
     C                            REFERENCES:
     C
     C               (1) N.A. TSYGANENKO AND D.P. STERN, A NEW-GENERATION GLOBAL
     C           MAGNETOSPHERE FIELD MODEL  , BASED ON SPACECRAFT MAGNETOMETER DATA,
     C           ISTP NEWSLETTER, V.6, NO.1, P.21, FEB.1996.
     C
     c              (2) N.A.TSYGANENKO,  MODELING THE EARTH'S MAGNETOSPHERIC
     C           MAGNETIC FIELD CONFINED WITHIN A REALISTIC MAGNETOPAUSE,
     C           J.GEOPHYS.RES., V.100, P. 5599, 1995.
     C
     C              (3) N.A. TSYGANENKO AND M.PEREDO, ANALYTICAL MODELS OF THE
     C           MAGNETIC FIELD OF DISK-SHAPED CURRENT SHEETS, J.GEOPHYS.RES.,
     C           V.99, P. 199, 1994.
     C
     c----------------------------------------------------------------------
     */
    /*
     */
    static void T96_01 (int const IOPT,double const PARMOD[],double const PS, double const X, double const Y, double const Z, double *BX,double *BY,double *BZ)
    {
        //IMPLICIT REAL*8 (A-H,O-Z)
        //REAL PDYN,DST,BYIMF,BZIMF,PS,X,Y,Z,BX,BY,BZ,QX,QY,QZ,PARMOD(10),
        //*   A(9)
        double PDYN, DST, BYIMF, BZIMF, QX=0., QY=0., QZ=0.;
        //DATA PDYN0,EPS10 /2.,3630.7/
        static double const PDYN0 = 2., EPS10 = 3630.7;
        //DATA A/1.162,22.344,18.50,2.602,6.903,5.287,0.5790,0.4462,0.7850/
        static double const A[10] = {0./*dummy*/, 1.162,22.344,18.50,2.602,6.903,5.287,0.5790,0.4462,0.7850};
        
        //DATA  AM0,S0,X00,DSIG/70.,1.08,5.48,0.005/
        static double const AM0 = 70., S0 = 1.08, X00 = 5.48, DSIG = 0.005;
        //DATA  DELIMFX,DELIMFY /20.,10./
        static double const DELIMFX = 20., DELIMFY = 10.;
        
        PDYN=PARMOD[0]; // Start from 0!
        DST=PARMOD[1];
        BYIMF=PARMOD[2];
        BZIMF=PARMOD[3];
        
        double SPS=sin(PS);
        double PPS=PS;
        
        double DEPR=0.8*DST-13.*sqrt(PDYN);  //  DEPR is an estimate of total near-Earth
        //                                  depression, based on DST and Pdyn
        //                                      (usually, DEPR < 0 )
        
        //  CALCULATE THE IMF-RELATED QUANTITIES:
        
        double BT=sqrt(BYIMF*BYIMF+BZIMF*BZIMF);
        
        //IF (BYIMF.EQ.0..AND.BZIMF.EQ.0.) THEN
        double THETA;
        if (BYIMF==0. && BZIMF==0.) {
            THETA=0.;
            //GOTO 1
        } else {
            THETA=atan2(BYIMF,BZIMF);
            //IF (THETA.LE.0.D0) THETA=THETA+6.2831853
            if (THETA<=0.e0) THETA=THETA+6.2831853;
        }//ENDIF
        double CT=cos(THETA);        //  1
        double ST=sin(THETA);
        double EPS=718.5*sqrt(PDYN)*BT*sin(THETA/2.);
        
        double FACTEPS=EPS/EPS10-1.;
        double FACTPD=sqrt(PDYN/PDYN0)-1.;
        
        double RCAMPL=-A[1]*DEPR;     //   RCAMPL is the amplitude of the ring current
        //                  (positive and equal to abs.value of RC depression at origin
        
        double TAMPL2=A[2]+A[3]*FACTPD+A[4]*FACTEPS;
        double TAMPL3=A[5]+A[6]*FACTPD;
        double B1AMPL=A[7]+A[8]*FACTEPS;
        double B2AMPL=20.*B1AMPL;  // IT IS EQUIVALENT TO ASSUMING THAT THE TOTAL CURRENT
        //                           IN THE REGION 2 SYSTEM IS 40% OF THAT IN REGION 1
        double RECONN=A[9];
        
        double XAPPA=pow(PDYN/PDYN0, 0.14);
        double XAPPA3=XAPPA*XAPPA*XAPPA;
        double YS=Y*CT-Z*ST;
        double ZS=Z*CT+Y*ST;
        
        double FACTIMF=exp(X/DELIMFX-(YS/DELIMFY)*(YS/DELIMFY));
        
        //  CALCULATE THE "IMF" COMPONENTS OUTSIDE THE LAYER  (HENCE BEGIN WITH "O")
        
        double OIMFX=0.;
        double OIMFY=RECONN*BYIMF*FACTIMF;
        double OIMFZ=RECONN*BZIMF*FACTIMF;
        
        double RIMFAMPL=RECONN*BT;
        
        PPS=PS;
        double XX=X*XAPPA;
        double YY=Y*XAPPA;
        double ZZ=Z*XAPPA;
        /*
         C  SCALE AND CALCULATE THE MAGNETOPAUSE PARAMETERS FOR THE INTERPOLATION ACROSS
         C   THE BOUNDARY LAYER (THE COORDINATES XX,YY,ZZ  ARE ALREADY SCALED)
         */
        double X0=X00/XAPPA;
        double AM=AM0/XAPPA;
        double RHO2=Y*Y+Z*Z;
        double ASQ=AM*AM;
        double XMXM=AM+X-X0;
        //IF (XMXM.LT.0.) XMXM=0. ! THE BOUNDARY IS A CYLINDER TAILWARD OF X=X0-AM
        if (XMXM<0.) XMXM=0.; // THE BOUNDARY IS A CYLINDER TAILWARD OF X=X0-AM
        double AXX0=XMXM*XMXM;
        double ARO=ASQ+RHO2;
        double SIGMA=sqrt((ARO+AXX0+sqrt((ARO+AXX0)*(ARO+AXX0)-4.*ASQ*AXX0))/(2.*ASQ));
        /*
         C   NOW, THERE ARE THREE POSSIBLE CASES:
         C    (1) INSIDE THE MAGNETOSPHERE
         C    (2) IN THE BOUNDARY LAYER
         C    (3) OUTSIDE THE MAGNETOSPHERE AND B.LAYER
         C       FIRST OF ALL, CONSIDER THE CASES (1) AND (2):
         */
        //IF (SIGMA.LT.S0+DSIG) THEN  !  CALCULATE THE T95_06 FIELD (WITH THE
        if (SIGMA<S0+DSIG) {  //  CALCULATE THE T95_06 FIELD (WITH THE
            //                       POTENTIAL "PENETRATED" INTERCONNECTION FIELD):
            
            double CFX, CFY, CFZ;
            DIPSHLD(PPS,XX,YY,ZZ,&CFX,&CFY,&CFZ);
            double BXRC,BYRC,BZRC,BXT2,BYT2,BZT2,BXT3,BYT3,BZT3;
            TAILRC96(SPS,XX,YY,ZZ,&BXRC,&BYRC,&BZRC,&BXT2,&BYT2,&BZT2,&BXT3,&BYT3,&BZT3);
            double R1X, R1Y, R1Z;
            BIRK1TOT_02(PPS,XX,YY,ZZ,&R1X,&R1Y,&R1Z);
            double R2X, R2Y, R2Z;
            BIRK2TOT_02(PPS,XX,YY,ZZ,&R2X,&R2Y,&R2Z);
            double RIMFX, RIMFYS, RIMFZS;
            INTERCON(XX,YS*XAPPA,ZS*XAPPA,&RIMFX,&RIMFYS,&RIMFZS);
            double RIMFY=RIMFYS*CT+RIMFZS*ST;
            double RIMFZ=RIMFZS*CT-RIMFYS*ST;
            
            double FX=CFX*XAPPA3+RCAMPL*BXRC +TAMPL2*BXT2+TAMPL3*BXT3 +B1AMPL*R1X +B2AMPL*R2X +RIMFAMPL*RIMFX;
            double FY=CFY*XAPPA3+RCAMPL*BYRC +TAMPL2*BYT2+TAMPL3*BYT3 +B1AMPL*R1Y +B2AMPL*R2Y +RIMFAMPL*RIMFY;
            double FZ=CFZ*XAPPA3+RCAMPL*BZRC +TAMPL2*BZT2+TAMPL3*BZT3 +B1AMPL*R1Z +B2AMPL*R2Z +RIMFAMPL*RIMFZ;
            
            //  NOW, LET US CHECK WHETHER WE HAVE THE CASE (1). IF YES - WE ARE DONE:
            
            //IF (SIGMA.LT.S0-DSIG) THEN
            if (SIGMA<S0-DSIG) {
                *BX=FX;
                *BY=FY;
                *BZ=FZ;
                //ELSE  !  THIS IS THE MOST COMPLEX CASE: WE ARE INSIDE
            } else {  //  THIS IS THE MOST COMPLEX CASE: WE ARE INSIDE
                //                  THE INTERPOLATION REGION
                double FINT=0.5*(1.-(SIGMA-S0)/DSIG);
                double FEXT=0.5*(1.+(SIGMA-S0)/DSIG);
                
                DIPOLE(PS,X,Y,Z,&QX,&QY,&QZ);
                *BX=(FX+QX)*FINT+OIMFX*FEXT -QX;
                *BY=(FY+QY)*FINT+OIMFY*FEXT -QY;
                *BZ=(FZ+QZ)*FINT+OIMFZ*FEXT -QZ;
                
            }// ENDIF    THE CASES (1) AND (2) ARE EXHAUSTED; THE ONLY REMAINING
            //                      POSSIBILITY IS NOW THE CASE (3):
            //ELSE
        } else {//ELSE
            DIPOLE(PS,X,Y,Z,&QX,&QY,&QZ);
            *BX=OIMFX-QX;
            *BY=OIMFY-QY;
            *BZ=OIMFZ-QZ;
        }//ENDIF
        
    }//END

}