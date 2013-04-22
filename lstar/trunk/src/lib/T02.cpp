//
//  T02.cpp
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

#include "T02.h"
#include "Geopack.h"
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace UBK {
    using namespace std;

    T02::T02 (Geopack const* geopack, double const parmod[]) : TSExternalField(geopack, 1, parmod)
    {
        if (NULL == parmod) {
            throw invalid_argument("PARMOD must be 10-element array.");
        }
    }

    static void T01_01 (int const IOPT,double const PARMOD[],double const PS,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ);
    void T02::getFieldInGSW_atPoint(Point *bOut, const Point ptgsw) const
    {
        T01_01(this->iopt(), this->parmod(), this->geopack()->psi(), ptgsw.x, ptgsw.y, ptgsw.z, &bOut->x, &bOut->y, &bOut->z);
    }

    /*-----------------------------Model---------------------*/
    /* Common blocks
     */
    struct TAIL_COMMON {
        double DXSHIFT1,DXSHIFT2,D,DELTADY;
    };
    struct BIRKPAR_COMMON {
        double XKAPPA1,XKAPPA2;
    };
    struct RCPAR_COMMON {
        double SC_SY,SC_AS,PHI;
    };
    struct G_COMMON {
        double G;
    };
    struct RH0_COMMON {
        double RH0;
    };
    struct DPHI_B_RHO0_COMMON {
        double DPHI, B, RHO_0, XKAPPA;
    };

    //======================================================================================
    static double R_S(double const A[],double const R,double const THETA) {
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION A(31)

        return R+A[2]/R+A[3]*R/sqrt(R*R+pow(A[11], 2))+A[4]*R/(R*R+pow(A[12], 2))
        +(A[5]+A[6]/R+A[7]*R/sqrt(R*R+pow(A[13], 2))+A[8]*R/(R*R+pow(A[14], 2)))*
        cos(THETA)
        +(A[9]*R/sqrt(R*R+pow(A[15], 2))+A[10]*R/pow((R*R+pow(A[16], 2)), 2))
        *cos(2.e0*THETA);

        //RETURN
    } // END

    //-----------------------------------------------------------------------------

    static double THETA_S(double const A[],double const R,double const THETA) {
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION A(31)

        return THETA+(A[17]+A[18]/R+A[19]/(R*R)
                      +A[20]*R/sqrt(R*R+pow(A[27], 2)))*sin(THETA)
        +(A[21]+A[22]*R/sqrt(R*R+pow(A[28], 2))
          +A[23]*R/(R*R+pow(A[29], 2)))*sin(2.e0*THETA)
        +(A[24]+A[25]/R+A[26]*R/(R*R+pow(A[30], 2)))*sin(3.e0*THETA);

        //RETURN
    } // END

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    static double AP(double const R,double const SINT,double const COST) {
        /*
         C      Calculates azimuthal component of the vector potential of the symmetric
         c  part of the model ring current.
         */
        //IMPLICIT  REAL * 8  (A - H, O - Z)
        bool PROX;//LOGICAL PROX   !  INDICATES WHETHER WE ARE TOO CLOSE TO THE AXIS OF SYMMETRY, WHERE THE INVERSION
        //                                                            OF DIPOLAR COORDINATES BECOMES INACCURATE
        //DATA A1,A2,RRC1,DD1,RRC2,DD2,P1,R1,DR1,DLA1,P2,R2,DR2,DLA2,P3,
        //*R3,DR3/
        static double const A1=-456.5289941,A2=375.9055332,RRC1=4.274684950,DD1=2.439528329,RRC2=3.367557287,          //    CORRECTED VALUES
        DD2=3.146382545,P1=-0.2291904607,R1=3.746064740,DR1=1.508802177,DLA1=0.5873525737,        //   (UPDATED 04/20/06 (SEE NB#9, P.37))
        P2=0.1556236119,R2=4.993638842,DR2=3.324180497,DLA2=0.4368407663,P3=0.1855957207,
        R3=2.969226745,DR3=2.243367377;

        PROX=false;
        double SINT1=SINT;
        double COST1=COST;
        if (SINT1<1.e-2) {  //  TOO CLOSE TO Z-AXIS;  USE LINEAR INTERPOLATION BETWEEN SINT=0 & SINT=0.01
            SINT1=1.e-2;
            COST1=.99994999875;
            PROX=true;
        }//ENDIF

        double ALPHA=SINT1*SINT1/R;         //  R,THETA -> ALPHA,GAMMA
        double GAMMA=COST1/(R*R);

        double ARG1=-pow(((R-R1)/DR1), 2)-pow((COST1/DLA1), 2);
        double ARG2=-pow(((R-R2)/DR2), 2)-pow((COST1/DLA2), 2);
        double ARG3=-pow(((R-R3)/DR3), 2);
        double DEXP1;
        if (ARG1<-500.e0) {        //   TO PREVENT "FLOATING UNDERFLOW" CRASHES
            DEXP1=0.e0;
        } else {//ELSE
            DEXP1=exp(ARG1);
        }//ENDIF
        double DEXP2;
        if (ARG2<-500.e0) {
            DEXP2=0.e0;
        } else {//ELSE
            DEXP2=exp(ARG2);
        }//ENDIF
        double DEXP3;
        if (ARG3<-500.e0) {
            DEXP3=0.e0;
        } else {//ELSE
            DEXP3=exp(ARG3);
        }//ENDIF


        double ALPHA_S=ALPHA*(1.e0+P1*DEXP1+P2*DEXP2+P3*DEXP3);     //  ALPHA -> ALPHA_S  (DEFORMED)

        double GAMMA_S=GAMMA;
        double GAMMAS2=GAMMA_S*GAMMA_S;


        double ALSQH=ALPHA_S*ALPHA_S/2.e0;            //  ALPHA_S,GAMMA_S -> RS,SINTS,COSTS
        double F=64.e0/27.e0*GAMMAS2+ALSQH*ALSQH;
        double Q=pow((sqrt(F)+ALSQH), (1.e0/3.e0));
        double C=Q-4.e0*pow(GAMMAS2, (1.e0/3.e0))/(3.e0*Q);
        if (C<0.e0) C=0.e0;
        double G=sqrt(C*C+4.e0*pow(GAMMAS2, (1.e0/3.e0)));
        double RS=4.e0/((sqrt(2.e0*G-C)+sqrt(C))*(G+C));
        double COSTS=GAMMA_S*RS*RS;
        double SINTS=sqrt(1.e0-COSTS*COSTS);
        double RHOS=RS*SINTS;
        //double RHOS2=RHOS*RHOS;
        double ZS=RS*COSTS;

        //  1st loop:

        double P=pow((RRC1+RHOS), 2)+ZS*ZS+DD1*DD1;
        double XK2=4.e0*RRC1*RHOS/P;
        double XK=sqrt(XK2);
        double XKRHO12=XK*sqrt(RHOS);     //   SEE NB#4, P.3

        double XK2S=1.e0-XK2;
        double DL=log(1.e0/XK2S);
        double ELK=1.38629436112e0+XK2S*(0.09666344259e0+XK2S*(0.03590092383+
                                                               XK2S*(0.03742563713+XK2S*0.01451196212))) +DL*
        (0.5e0+XK2S*(0.12498593597e0+XK2S*(0.06880248576e0+
                                           XK2S*(0.03328355346e0+XK2S*0.00441787012e0))));
        double ELE=1.e0+XK2S*(0.44325141463e0+XK2S*(0.0626060122e0+XK2S*
                                                    (0.04757383546e0+XK2S*0.01736506451e0))) +DL*
        XK2S*(0.2499836831e0+XK2S*(0.09200180037e0+XK2S*
                                   (0.04069697526e0+XK2S*0.00526449639e0)));

        double APHI1=((1.e0-XK2*0.5e0)*ELK-ELE)/XKRHO12;

        //  2nd loop:

        P=pow((RRC2+RHOS), 2)+ZS*ZS+DD2*DD2;
        XK2=4.e0*RRC2*RHOS/P;
        XK=sqrt(XK2);
        XKRHO12=XK*sqrt(RHOS);     //   SEE NB#4, P.3

        XK2S=1.e0-XK2;
        DL=log(1.e0/XK2S);
        ELK=1.38629436112e0+XK2S*(0.09666344259e0+XK2S*(0.03590092383+
                                                        XK2S*(0.03742563713+XK2S*0.01451196212))) +DL*
        (0.5e0+XK2S*(0.12498593597e0+XK2S*(0.06880248576e0+
                                           XK2S*(0.03328355346e0+XK2S*0.00441787012e0))));
        ELE=1.e0+XK2S*(0.44325141463e0+XK2S*(0.0626060122e0+XK2S*
                                             (0.04757383546e0+XK2S*0.01736506451e0))) +DL*
        XK2S*(0.2499836831e0+XK2S*(0.09200180037e0+XK2S*
                                   (0.04069697526e0+XK2S*0.00526449639e0)));

        double APHI2=((1.e0-XK2*0.5e0)*ELK-ELE)/XKRHO12;

        double AP=A1*APHI1+A2*APHI2;
        if (PROX) AP=AP*SINT/SINT1;   //   LINEAR INTERPOLATION, IF TOO CLOSE TO THE Z-AXIS

        return AP;
    }//END

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    static double APPRC(double const R,double const SINT,double const COST) {
        /*
         C      Calculates azimuthal component of the vector potential of the symmetric
         c  part of the model PARTIAL ring current.
         */
        //IMPLICIT  REAL * 8  (A - H, O - Z)
        bool PROX;//LOGICAL PROX
        //DATA A1,A2,RRC1,DD1,RRC2,DD2,P1,ALPHA1,DAL1,BETA1,DG1,P2,ALPHA2,
        //* DAL2,BETA2,DG2,BETA3,P3,ALPHA3,DAL3,BETA4,DG3,BETA5,Q0,Q1,ALPHA4,
        //* DAL4,DG4,Q2,ALPHA5,DAL5,DG5,BETA6,BETA7
        static double const A1=-80.11202281,A2=12.58246758,RRC1=6.560486035,DD1=1.930711037,RRC2=3.827208119,
        DD2=.7789990504,P1=.3058309043,ALPHA1=.1817139853,DAL1=.1257532909,BETA1=3.422509402,
        DG1=.04742939676,P2=-4.800458958,ALPHA2=-.02845643596,DAL2=.2188114228,BETA2=2.545944574,
        DG2=.00813272793,BETA3=.35868244,P3=103.1601001,ALPHA3=-.00764731187,DAL3=.1046487459,
        BETA4=2.958863546,DG3=.01172314188,BETA5=.4382872938,Q0=.01134908150,Q1=14.51339943,
        ALPHA4=.2647095287,DAL4=.07091230197,DG4=.01512963586,Q2=6.861329631,ALPHA5=.1677400816,
        DAL5=.04433648846,DG5=.05553741389,BETA6=.7665599464,BETA7=.7277854652;

        PROX=false;
        double SINT1=SINT;
        double COST1=COST;
        if (SINT1<1.e-2) {//  !  TOO CLOSE TO Z-AXIS;  USE LINEAR INTERPOLATION BETWEEN SINT=0 & SINT=0.01
            SINT1=1.e-2;
            COST1=.99994999875;
            PROX=true;
        }//ENDIF

        double ALPHA=SINT1*SINT1/R;         //  R,THETA -> ALPHA,GAMMA
        double GAMMA=COST1/(R*R);

        double ARG1=-pow((GAMMA/DG1), 2);
        double ARG2=-pow(((ALPHA-ALPHA4)/DAL4), 2)-pow((GAMMA/DG4), 2);
        double DEXP1;
        if (ARG1<-500.e0) {        //   TO PREVENT "FLOATING UNDERFLOW" CRASHES
            DEXP1=0.e0;
        } else {//ELSE
            DEXP1=exp(ARG1);
        }//ENDIF
        double DEXP2;
        if (ARG2<-500.e0) {        //   TO PREVENT "FLOATING UNDERFLOW" CRASHES
            DEXP2=0.e0;
        } else {//ELSE
            DEXP2=exp(ARG2);
        }//ENDIF

        double ALPHA_S=ALPHA*(1.e0+P1/pow((1.e0+pow(((ALPHA-ALPHA1)/DAL1), 2)), BETA1)
                              *DEXP1+P2*(ALPHA-ALPHA2)/pow((1.e0+pow(((ALPHA-ALPHA2)/DAL2), 2)), BETA2)
                              /pow((1.e0+pow((GAMMA/DG2), 2)), BETA3)
                              +P3*pow((ALPHA-ALPHA3), 2)/pow((1.e0+pow(((ALPHA-ALPHA3)/DAL3), 2)), BETA4)
                              /pow((1.e0+pow((GAMMA/DG3), 2)), BETA5));     //  ALPHA -> ALPHA_S  (DEFORMED)

        double GAMMA_S=GAMMA*(1.e0+Q0+Q1*(ALPHA-ALPHA4)*DEXP2              //  GAMMA -> GAMMA_  (DEFORMED)
                              +Q2*(ALPHA-ALPHA5)/pow((1.e0+pow(((ALPHA-ALPHA5)/DAL5), 2)), BETA6)
                              /pow((1.e0+pow((GAMMA/DG5), 2)), BETA7));

        double GAMMAS2=GAMMA_S*GAMMA_S;

        double ALSQH=ALPHA_S*ALPHA_S/2.e0;                            //  ALPHA_S,GAMMA_S -> RS,SINTS,COSTS
        double F=64.e0/27.e0*GAMMAS2+ALSQH*ALSQH;
        double Q=pow((sqrt(F)+ALSQH), (1.e0/3.e0));
        double C=Q-4.e0*pow(GAMMAS2, (1.e0/3.e0))/(3.e0*Q);
        if (C<0.e0) C=0.e0;
        double G=sqrt(C*C+4.e0*pow(GAMMAS2, (1.e0/3.e0)));
        double RS=4.e0/((sqrt(2.e0*G-C)+sqrt(C))*(G+C));
        double COSTS=GAMMA_S*RS*RS;
        double SINTS=sqrt(1.e0-COSTS*COSTS);
        double RHOS=RS*SINTS;
        //double RHOS2=RHOS*RHOS;
        double ZS=RS*COSTS;

        //  1st loop:

        double P=pow((RRC1+RHOS), 2)+ZS*ZS+DD1*DD1;
        double XK2=4.e0*RRC1*RHOS/P;
        double XK=sqrt(XK2);
        double XKRHO12=XK*sqrt(RHOS);     //   SEE NB#4, P.3

        double XK2S=1.e0-XK2;
        double DL=log(1.e0/XK2S);
        double ELK=1.38629436112e0+XK2S*(0.09666344259e0+XK2S*(0.03590092383+
                                                               XK2S*(0.03742563713+XK2S*0.01451196212))) +DL*
        (0.5e0+XK2S*(0.12498593597e0+XK2S*(0.06880248576e0+
                                           XK2S*(0.03328355346e0+XK2S*0.00441787012e0))));
        double ELE=1.e0+XK2S*(0.44325141463e0+XK2S*(0.0626060122e0+XK2S*
                                                    (0.04757383546e0+XK2S*0.01736506451e0))) +DL*
        XK2S*(0.2499836831e0+XK2S*(0.09200180037e0+XK2S*
                                   (0.04069697526e0+XK2S*0.00526449639e0)));

        double APHI1=((1.e0-XK2*0.5e0)*ELK-ELE)/XKRHO12;

        //  2nd loop:

        P=pow((RRC2+RHOS), 2)+ZS*ZS+DD2*DD2;
        XK2=4.e0*RRC2*RHOS/P;
        XK=sqrt(XK2);
        XKRHO12=XK*sqrt(RHOS);     //   SEE NB#4, P.3

        XK2S=1.e0-XK2;
        DL=log(1.e0/XK2S);
        ELK=1.38629436112e0+XK2S*(0.09666344259e0+XK2S*(0.03590092383+
                                                        XK2S*(0.03742563713+XK2S*0.01451196212))) +DL*
        (0.5e0+XK2S*(0.12498593597e0+XK2S*(0.06880248576e0+
                                           XK2S*(0.03328355346e0+XK2S*0.00441787012e0))));
        ELE=1.e0+XK2S*(0.44325141463e0+XK2S*(0.0626060122e0+XK2S*
                                             (0.04757383546e0+XK2S*0.01736506451e0))) +DL*
        XK2S*(0.2499836831e0+XK2S*(0.09200180037e0+XK2S*
                                   (0.04069697526e0+XK2S*0.00526449639e0)));

        double APHI2=((1.e0-XK2*0.5e0)*ELK-ELE)/XKRHO12;

        double APPRC=A1*APHI1+A2*APHI2;
        if (PROX) APPRC=APPRC*SINT/SINT1;   //   LINEAR INTERPOLATION, IF TOO CLOSE TO THE Z-AXIS

        return APPRC;
    }//END

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    static void FFS(double const A,double const A0,double const DA,double *F,double *FA,double *FS) {
        //IMPLICIT  REAL * 8  (A - H, O - Z)
        double SQ1=sqrt(pow((A+A0), 2)+DA*DA);
        double SQ2=sqrt(pow((A-A0), 2)+DA*DA);
        *FA=2.e0/(SQ1+SQ2);
        *F=*FA*A;
        *FS=0.5e0*(SQ1+SQ2)/(SQ1*SQ2)*(1.e0-*F**F);
        //RETURN
    }//END

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    static double BR_PRC_Q (double const R,double const SINT,double const COST) {

        //alculates the radial component of the "quadrupole" part of the model partial ring current.

        //IMPLICIT  REAL * 8  (A - H, O - Z)

        //DATA A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,   ! ALL LINEAR PARAMETERS HERE
        //* A18,XK1,AL1,DAL1,B1,BE1,XK2,AL2,DAL2,B2,BE2,XK3,XK4,AL3,DAL3,B3,  ! WERE MULTIPLIED BY 0.1,
        //* BE3,AL4,DAL4,DG1,AL5,DAL5,DG2,C1,C2,C3,AL6,DAL6,DRM/-21.2666329,  ! SO THAT THEY CORRESPOND TO P_0=1 nPa,
        static double const A1 = -21.2666329,A2 = 32.24527521,A3 = -6.062894078,A4 = 7.515660734, //! RATHER THAN THE ORIGINAL VALUE OF 10 nPa
        A5 = 233.7341288,A6 = -227.1195714,A7 = 8.483233889,A8 = 16.80642754,A9 = -24.63534184, // ! ASSUMED IN THE BIOT-SAVART INTEGRAL
        A10 = 9.067120578,A11 = -1.052686913,A12 = -12.08384538,A13 = 18.61969572,
        A14 = -12.71686069,A15 = 47017.35679,A16 = -50646.71204,A17 = 7746.058231,
        A18 = 1.531069371,XK1 = 2.318824273,AL1 = .1417519429,DAL1 = .6388013110E-02,
        B1 = 5.303934488,BE1 = 4.213397467,XK2 = .7955534018,AL2 = .1401142771,
        DAL2 = .2306094179E-01,B2 = 3.462235072,BE2 = 2.568743010,XK3 = 3.477425908,
        XK4 = 1.922155110,AL3 = .1485233485,DAL3 = .2319676273E-01,B3 = 7.830223587,
        BE3 = 8.492933868,AL4 = .1295221828,DAL4 = .01753008801,DG1 = .01125504083,
        AL5 = .1811846095,DAL5 = .04841237481,DG2 = .01981805097,C1 = 6.557801891,
        C2 = 6.348576071,C3 = 5.744436687,AL6 = .2265212965,DAL6 = .1301957209,
        DRM = .5654023158;

        double SINT2=SINT*SINT;
        double COST2=COST*COST;
        double SC=SINT*COST;
        double ALPHA=SINT2/R;
        double GAMMA=COST/(R*R);
        double F,FA,FS;
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

        return A1*D1+A2*D2+A3*D3+A4*D4+A5*D5+A6*D6+A7*D7+A8*D8+A9*D9+
        A10*D10+A11*D11+A12*D12+A13*D13+A14*D14+A15*D15+A16*D16+A17*D17+
        A18*D18;

        //RETURN
    }//END

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    static double BT_PRC_Q (double const R,double const SINT,double const COST) {

        //alculates the Theta component of the "quadrupole" part of the model partial ring current.

        //IMPLICIT  REAL * 8  (A - H, O - Z)

        //DATA A1,A2,A3,A4,A5,A6,A7,A8,A9,A10,A11,A12,A13,A14,A15,A16,A17,  ! ALL LINEAR PARAMETERS HERE
        //*XK1,AL1,DAL1,B1,BE1,XK2,AL2,DAL2,BE2,XK3,XK4,AL3,DAL3,B3,BE3,AL4, ! WERE MULTIPLIED BY 0.1,
        //*DAL4,DG1,AL5,DAL5,DG2,C1,C2,C3/12.74640393,-7.516393516,          ! SO THAT THEY CORRESPOND TO P_0=1 nPa,
        static double const A1=12.74640393, A2=-7.516393516, A3=-5.476233865, //  ! RATHER THAN THE ORIGINAL VALUE OF 10 nPa
        A4=3.212704645, A5=-59.10926169, A6=46.62198189, A7=-.01644280062, //   ! ASSUMED IN THE BIOT-SAVART INTEGRAL
        A8=.1234229112, A9=-.08579198697, A10=.01321366966, A11=.8970494003,
        A12=9.136186247, A13=-38.19301215, A14=21.73775846, A15=-410.0783424,
        A16=-69.90832690, A17=-848.8543440, XK1=1.243288286, AL1=.2071721360,
        DAL1=.05030555417, B1=7.471332374, BE1=3.180533613, XK2=1.376743507,
        AL2=.1568504222, DAL2=.02092910682, BE2=1.985148197, XK3=.3157139940,
        XK4=1.056309517, AL3=.1701395257, DAL3=.1019870070, B3=6.293740981,
        BE3=5.671824276, AL4=.1280772299, DAL4=.02189060799, DG1=.01040696080,
        AL5=.1648265607, DAL5=.04701592613, DG2=.01526400086, C1=12.88384229,
        C2=3.361775101, C3=23.44173897;

        double SINT2=SINT*SINT;
        double COST2=COST*COST;
        //double SC=SINT*COST;
        double ALPHA=SINT2/R;
        double GAMMA=COST/(R*R);
        double F,FA,FS;
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

        return A1*D1+A2*D2+A3*D3+A4*D4+A5*D5+A6*D6+A7*D7+A8*D8+A9*D9+
        A10*D10+A11*D11+A12*D12+A13*D13+A14*D14+A15*D15+A16*D16+A17*D17;

        //RETURN
    }//END

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    static void SHLCAR3X3(double const X,double const Y,double const Z,double const PS,double *BX,double *BY,double *BZ) {
        /*
         C   THIS S/R RETURNS THE SHIELDING FIELD FOR THE EARTH'S DIPOLE,
         C   REPRESENTED BY  2x3x3=18 "CARTESIAN" HARMONICS, tilted with respect
         C   to the z=0 plane  (NB#4, p.74)
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
        //IMPLICIT  REAL * 8  (A - H, O - Z)

        //DIMENSION A(50)
        static double const A[51] = {0./*dummy*/, -901.2327248,895.8011176,817.6208321,-845.5880889,
            -83.73539535,86.58542841,336.8781402,-329.3619944,-311.2947120,
            308.6011161,31.94469304,-31.30824526,125.8739681,-372.3384278,
            -235.4720434,286.7594095,21.86305585,-27.42344605,-150.4874688,
            2.669338538,1.395023949,-.5540427503,-56.85224007,3.681827033,
            -43.48705106,5.103131905,1.073551279,-.6673083508,12.21404266,
            4.177465543,5.799964188,-.3977802319,-1.044652977,.5703560010,
            3.536082962,-3.222069852,9.620648151,6.082014949,27.75216226,
            12.44199571,5.122226936,6.982039615,20.12149582,6.150973118,
            4.663639687,15.73319647,2.303504968,5.840511214,.8385953499E-01,
            .3477844929};

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

        double T1 =A[49];
        double T2 =A[50];

        double CPS=cos(PS);
        double SPS=sin(PS);
        double S2PS=2.e0*CPS;      //   MODIFIED HERE (SIN(2*PS) INSTEAD OF SIN(3*PS))

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
         C       I=1:
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
         C      I=3:
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

        /*  MAKE THE TERMS IN THE 2ND SUM ("PARALLEL" SYMMETRY):
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
         C       I=2:
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
         C       I=3:
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

        //RETURN
    }//END

    //$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$

    static void TAILDISK(double const D0,double const DELTADX,double const DELTADY,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {
        /*
         c       THIS SUBROUTINE COMPUTES THE COMPONENTS OF THE TAIL CURRENT FIELD,
         C       SIMILAR TO THAT DESCRIBED BY TSYGANENKO AND PEREDO (1994).  THE
         C       DIFFERENCE IS THAT NOW WE USE SPACEWARPING, AS DESCRIBED IN OUR
         C       PAPER ON MODELING BIRKELAND CURRENTS (TSYGANENKO AND STERN, 1996)
         C       INSTEAD OF SHEARING IT IN THE SPIRIT OF THE T89 TAIL MODEL.
         */
        //IMPLICIT REAL*8 (A-H,O-Z)

        //DIMENSION F(5),B(5),C(5)

        static double const F[6] = {0./*dummy*/, -71.09346626e0,-1014.308601e0,-1272.939359e0,
            -3224.935936e0,-44546.86232e0};
        static double const B[6] = {0./*dummy*/, 10.90101242e0,12.68393898e0,13.51791954e0,14.86775017e0,
            15.12306404e0};
        static double const C[6] = {0./*dummy*/, .7954069972e0,.6716601849e0,1.174866319e0,2.565249920e0,
            10.01986790e0};

        double RHO = sqrt(X*X+Y*Y);
        double DRHODX = X/RHO;
        double DRHODY = Y/RHO;

        double DEX = exp(X/7.e0);
        double D = D0+DELTADY*pow((Y/20.e0), 2)  +DELTADX*DEX;// THE LAST TERM (INTRODUCED 10/11/2000) MAKES THE SHEET
        double DDDY = DELTADY*Y*0.005e0;// THICKEN SUNWARD, TO AVOID PROBLEMS IN THE SUBSOLAR REGION
        double DDDX = DELTADX/7.e0*DEX;

        double DZETA = sqrt(Z*Z+D*D);// THIS IS THE SAME SIMPLE WAY TO SPREAD
        //                                        OUT THE SHEET, AS THAT USED IN T89
        double DDZETADX = D*DDDX/DZETA;
        double DDZETADY = D*DDDY/DZETA;
        double DDZETADZ = Z/DZETA;


        double DBX = 0.e0;
        double DBY = 0.e0;
        double DBZ = 0.e0;

        for (int I = 1; I<= 5; I++) {// DO 1 I = 1,5

            double BI = B[I];
            double CI = C[I];

            double S1 = sqrt(pow((RHO+BI), 2)+pow((DZETA+CI), 2));
            double S2 = sqrt(pow((RHO-BI), 2)+pow((DZETA+CI), 2));

            double DS1DRHO = (RHO+BI)/S1;
            double DS2DRHO = (RHO-BI)/S2;
            double DS1DDZ = (DZETA+CI)/S1;
            double DS2DDZ = (DZETA+CI)/S2;

            double DS1DX = DS1DRHO*DRHODX  +DS1DDZ*DDZETADX;
            double DS1DY = DS1DRHO*DRHODY  +   DS1DDZ*DDZETADY;
            double DS1DZ = DS1DDZ*DDZETADZ;

            double DS2DX = DS2DRHO*DRHODX  +DS2DDZ*DDZETADX;
            double DS2DY = DS2DRHO*DRHODY  +   DS2DDZ*DDZETADY;
            double DS2DZ = DS2DDZ*DDZETADZ;

            double S1TS2 = S1*S2;
            double S1PS2 = S1+S2;
            double S1PS2SQ = S1PS2*S1PS2;

            double FAC1 = sqrt(S1PS2SQ-pow((2.e0*BI), 2));
            double AS = FAC1/(S1TS2*S1PS2SQ);
            double DASDS1 = (1.e0/(FAC1*S2)-AS/S1PS2*(S2*S2+S1*(3.e0*S1+4.e0*S2)))
            /(S1*S1PS2);
            double DASDS2 = (1.e0/(FAC1*S1)-AS/S1PS2*(S1*S1+S2*(3.e0*S2+4.e0*S1)))
            /(S2*S1PS2);

            double DASDX = DASDS1*DS1DX+DASDS2*DS2DX;
            double DASDY = DASDS1*DS1DY+DASDS2*DS2DY;
            double DASDZ = DASDS1*DS1DZ+DASDS2*DS2DZ;

            DBX = DBX-F[I]*X*DASDZ;
            DBY = DBY-F[I]*Y*DASDZ;
            DBZ = DBZ+F[I]*(2.e0*AS+X*DASDX+Y*DASDY);// 1
        }
        *BX = DBX;
        *BY = DBY;
        *BZ = DBZ;

        //RETURN
    }// END

    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     C
     C THIS CODE RETURNS THE SHIELDING FIELD REPRESENTED BY  5x5 = 25 "CARTESIAN"
     C    HARMONICS
     */
    static void  SHLCAR5X5(double const A[],double const X,double const Y,double const Z,double const DSHIFT,double *HX,double *HY,double *HZ) {
        /*
         C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         C  The NLIN coefficients are the amplitudes of the "cartesian"
         c    harmonics (A(1)-A(NLIN).
         c  The NNP nonlinear parameters (A(NLIN+1)-A(NTOT) are the scales Pi and Ri
         C   entering the arguments of exponents, sines, and cosines in each of the
         C   NLIN "Cartesian" harmonics
         C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
         */
        //IMPLICIT  REAL * 8  (A - H, O - Z)

        //DIMENSION A(60)

        double DHX = 0.e0;
        double DHY = 0.e0;
        double DHZ = 0.e0;

        int L = 0;

        for (int I = 1; I<= 5; I++) {// DO 2 I = 1,5
            double RP = 1.e0/A[50+I];
            double CYPI = cos(Y*RP);
            double SYPI = sin(Y*RP);

            for (int K = 1; K<= 5; K++) {// DO 2 K = 1,5
                double RR = 1.e0/A[55+K];
                double SZRK = sin(Z*RR);
                double CZRK = cos(Z*RR);
                double SQPR = sqrt(RP*RP+RR*RR);
                double EPR = exp(X*SQPR);

                double DBX = -SQPR*EPR*CYPI*SZRK;
                double DBY = RP*EPR*SYPI*SZRK;
                double DBZ = -RR*EPR*CYPI*CZRK;

                L = L+2;
                double COEF = A[L-1]+A[L]*DSHIFT;

                DHX = DHX+COEF*DBX;
                DHY = DHY+COEF*DBY;
                DHZ = DHZ+COEF*DBZ;
            }
        } // 2      CONTINUE

        *HX = DHX;
        *HY = DHY;
        *HZ = DHZ;

        //RETURN
    } // END

    //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    static void UNWARPED (int const IOPT,double const X,double const Y,double const Z,double *BX1,double *BY1,double *BZ1,double *BX2,double *BY2,double *BZ2, struct TAIL_COMMON const* TAIL_C) {
        /*
         C   IOPT - TAIL FIELD MODE FLAG:   IOPT = 0 - THE TWO TAIL MODES ARE ADDED UP
         C                                  IOPT = 1 - MODE 1 ONLY
         C                                  IOPT = 2 - MODE 2 ONLY
         C
         C    CALCULATES GSM COMPONENTS OF THE SHIELDED FIELD OF TWO TAIL MODES WITH UNIT
         C    AMPLITUDES,  WITHOUT ANY WARPING OR BENDING.  NONLINEAR PARAMETERS OF THE MODES
         C    ARE FORWARDED HERE VIA A COMMON BLOCK /TAIL/.
         */
        //IMPLICIT REAL*8 (A-H,O-Z)

        //DIMENSION A1(60),A2(60)  !   TAIL SHIELDING FIELD PARAMETERS FOR THE MODES #1 & #2

        double const* DXSHIFT1 = &TAIL_C->DXSHIFT1, *DXSHIFT2 = &TAIL_C->DXSHIFT2, *D0 = &TAIL_C->D, *DELTADY = &TAIL_C->DELTADY;// COMMON /TAIL/ DXSHIFT1,DXSHIFT2,D0,DELTADY

        static double const DELTADX1 = 1.e0,ALPHA1 = 1.1e0,XSHIFT1 = 6.e0;
        static double const DELTADX2 = 0.e0,ALPHA2 = .25e0,XSHIFT2 = 4.e0;

        static double const A1[61] = {0./*dummy*/, -25.45869857,57.35899080,317.5501869,-2.626756717,
            -93.38053698,-199.6467926,-858.8129729,34.09192395,845.4214929,
            -29.07463068,47.10678547,-128.9797943,-781.7512093,6.165038619,
            167.8905046,492.0680410,1654.724031,-46.77337920,-1635.922669,
            40.86186772,-.1349775602,-.9661991179E-01,-.1662302354,
            .002810467517,.2487355077,.1025565237,-14.41750229,-.8185333989,
            11.07693629,.7569503173,-9.655264745,112.2446542,777.5948964,
            -5.745008536,-83.03921993,-490.2278695,-1155.004209,39.08023320,
            1172.780574,-39.44349797,-14.07211198,-40.41201127,-313.2277343,
            2.203920979,8.232835341,197.7065115,391.2733948,-18.57424451,
            -437.2779053,23.04976898,11.75673963,13.60497313,4.691927060,
            18.20923547,27.59044809,6.677425469,1.398283308,2.839005878,
            31.24817706,24.53577264};

        static double const A2[61] = {0./*dummy*/, -287187.1962,4970.499233,410490.1952,-1347.839052,
            -386370.3240,3317.983750,-143462.3895,5706.513767,171176.2904,
            250.8882750,-506570.8891,5733.592632,397975.5842,9771.762168,
            -941834.2436,7990.975260,54313.10318,447.5388060,528046.3449,
            12751.04453,-21920.98301,-21.05075617,31971.07875,3012.641612,
            -301822.9103,-3601.107387,1797.577552,-6.315855803,142578.8406,
            13161.93640,804184.8410,-14168.99698,-851926.6360,-1890.885671,
            972475.6869,-8571.862853,26432.49197,-2554.752298,-482308.3431,
            -4391.473324,105155.9160,-1134.622050,-74353.53091,-5382.670711,
            695055.0788,-916.3365144,-12111.06667,67.20923358,-367200.9285,
            -21414.14421,14.75567902,20.75638190,59.78601609,16.86431444,
            32.58482365,23.69472951,17.24977936,13.64902647,68.40989058,
            11.67828167};

        static double const XM1 = -12.e0,XM2 = -12.e0;

        if (IOPT!=2) {// IF (IOPT.EQ.2) GOTO 1

            double XSC1 = (X-XSHIFT1-*DXSHIFT1)*ALPHA1-XM1*(ALPHA1-1.e0);
            double YSC1 = Y*ALPHA1;
            double ZSC1 = Z*ALPHA1;
            double D0SC1 = *D0*ALPHA1;// HERE WE USE A SINGLE VALUE D0 OF THE THICKNESS FOR BOTH MODES
            double FX1,FY1,FZ1, HX1,HY1,HZ1;
            TAILDISK(D0SC1,DELTADX1,*DELTADY,XSC1,YSC1,ZSC1,&FX1,&FY1,&FZ1);
            SHLCAR5X5(A1,X,Y,Z,*DXSHIFT1,&HX1,&HY1,&HZ1);

            *BX1 = FX1+HX1;
            *BY1 = FY1+HY1;
            *BZ1 = FZ1+HZ1;

            if (IOPT==1) {
                *BX2 = 0.e0;
                *BY2 = 0.e0;
                *BZ2 = 0.e0;
                return;
            }// ENDIF
        }
        double XSC2 = (X-XSHIFT2-*DXSHIFT2)*ALPHA2-XM2*(ALPHA2-1.e0);// 1
        double YSC2 = Y*ALPHA2;
        double ZSC2 = Z*ALPHA2;
        double D0SC2 = *D0*ALPHA2;// HERE WE USE A SINGLE VALUE D0 OF THE THICKNESS FOR BOTH MODES
        double FX2,FY2,FZ2,HX2,HY2,HZ2;
        TAILDISK(D0SC2,DELTADX2,*DELTADY,XSC2,YSC2,ZSC2,&FX2,&FY2,&FZ2);
        SHLCAR5X5(A2,X,Y,Z,*DXSHIFT2,&HX2,&HY2,&HZ2);

        *BX2 = FX2+HX2;
        *BY2 = FY2+HY2;
        *BZ2 = FZ2+HZ2;

        if (IOPT==2) {
            *BX1 = 0.e0;
            *BY1 = 0.e0;
            *BZ1 = 0.e0;
            return;
        }// ENDIF

        //RETURN
    }// END

    //------------------------------------------------------------------

    static void WARPED (int const IOPT,double const PS,double const X,double const Y,double const Z,double *BX1,double *BY1,double *BZ1,double *BX2,double *BY2,double *BZ2, struct G_COMMON const* G_C, struct TAIL_COMMON const* TAIL_C) {
        /*
         C   CALCULATES GSM COMPONENTS OF THE WARPED FIELD FOR TWO TAIL UNIT MODES.
         C   THE WARPING DEFORMATION IS IMPOSED ON THE UNWARPED FIELD, COMPUTED
         C   BY THE S/R "UNWARPED".  THE WARPING PARAMETER G WAS OBTAINED BY LEAST
         C   SQUARES FITTING TO THE ENTIRE DATASET.
         C
         C   IOPT - TAIL FIELD MODE FLAG:   IOPT=0 - THE TWO TAIL MODES ARE ADDED UP
         C                                  IOPT=1 - MODE 1 ONLY
         C                                  IOPT=2 - MODE 2 ONLY
         */
        //IMPLICIT REAL*8 (A-H,O-Z)

        double const* G = &G_C->G;//COMMON /G/ G
        double DGDX=0.e0;
        double XL=20.e0;
        double DXLDX=0.e0;

        double SPS=sin(PS);
        double RHO2=Y*Y+Z*Z;
        double RHO=sqrt(RHO2);
        double PHI, CPHI, SPHI;
        if (Y==0.e0 && Z==0.e0) {
            PHI=0.e0;
            CPHI=1.e0;
            SPHI=0.e0;
        } else {//ELSE
            PHI = atan2(Z,Y);
            CPHI = Y/RHO;
            SPHI = Z/RHO;
        }//ENDIF

        double RR4L4 = RHO/(RHO2*RHO2+pow(XL, 4));

        double F = PHI+*G*RHO2*RR4L4*CPHI*SPS;
        double DFDPHI = 1.e0-*G*RHO2*RR4L4*SPHI*SPS;
        double DFDRHO = *G*RR4L4*RR4L4*(3.e0*pow(XL, 4)-RHO2*RHO2)*CPHI*SPS;
        double DFDX = RR4L4*CPHI*SPS*(DGDX*RHO2-*G*RHO*RR4L4*4.e0*pow(XL, 3)*DXLDX);

        double CF = cos(F);
        double SF = sin(F);
        double YAS = RHO*CF;
        double ZAS = RHO*SF;
        double BX_AS1,BY_AS1,BZ_AS1, BX_AS2,BY_AS2,BZ_AS2;
        UNWARPED (IOPT,X,YAS,ZAS,&BX_AS1,&BY_AS1,&BZ_AS1,
                  &BX_AS2,&BY_AS2,&BZ_AS2, TAIL_C);

        double BRHO_AS = BY_AS1*CF+BZ_AS1*SF; // DEFORM THE 1ST MODE
        double BPHI_AS = -BY_AS1*SF+BZ_AS1*CF;

        double BRHO_S = BRHO_AS*DFDPHI;
        double BPHI_S = BPHI_AS-RHO*(BX_AS1*DFDX+BRHO_AS*DFDRHO);
        *BX1 = BX_AS1*DFDPHI;

        *BY1 = BRHO_S*CPHI-BPHI_S*SPHI;
        *BZ1 = BRHO_S*SPHI+BPHI_S*CPHI; // DONE

        BRHO_AS = BY_AS2*CF+BZ_AS2*SF; // DEFORM THE 2ND MODE
        BPHI_AS = -BY_AS2*SF+BZ_AS2*CF;

        BRHO_S = BRHO_AS*DFDPHI;
        BPHI_S = BPHI_AS-RHO*(BX_AS2*DFDX+BRHO_AS*DFDRHO);
        *BX2 = BX_AS2*DFDPHI;

        *BY2 = BRHO_S*CPHI-BPHI_S*SPHI;
        *BZ2 = BRHO_S*SPHI+BPHI_S*CPHI; // DONE

        //RETURN
    }//END

    //############################################################################

    static void DEFORMED (int const IOPT,double const PS,double const X,double const Y,double const Z,double *BX1,double *BY1,double *BZ1,double *BX2,double *BY2,double *BZ2, struct RH0_COMMON const* RH0_C, struct G_COMMON const* G_C, struct TAIL_COMMON const* TAIL_C) {
        /*
         C   IOPT - TAIL FIELD MODE FLAG:   IOPT=0 - THE TWO TAIL MODES ARE ADDED UP
         C                                  IOPT=1 - MODE 1 ONLY
         C                                  IOPT=2 - MODE 2 ONLY
         C
         C   CALCULATES GSM COMPONENTS OF TWO UNIT-AMPLITUDE TAIL FIELD MODES,
         C    TAKING INTO ACCOUNT BOTH EFFECTS OF DIPOLE TILT:
         C    WARPING IN Y-Z (DONE BY THE S/R WARPED) AND BENDING IN X-Z (DONE BY THIS SUBROUTINE)
         */
        //IMPLICIT REAL*8 (A-H,O-Z)
        double const* RH0 = &RH0_C->RH0;//COMMON /RH0/ RH0
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

        //     DEFORM:
        double BXAS1,BYAS1,BZAS1,BXAS2,BYAS2,BZAS2;
        WARPED(IOPT,PS,XAS,Y,ZAS,&BXAS1,&BYAS1,&BZAS1,&BXAS2,&BYAS2,&BZAS2, G_C, TAIL_C);

        *BX1=BXAS1*DZASDZ-BZAS1*DXASDZ +BYAS1*FAC1;
        *BY1=BYAS1*FAC2;
        *BZ1=BZAS1*DXASDX-BXAS1*DZASDX +BYAS1*FAC3;

        *BX2=BXAS2*DZASDZ-BZAS2*DXASDZ +BYAS2*FAC1;
        *BY2=BYAS2*FAC2;
        *BZ2=BZAS2*DXASDX-BXAS2*DZASDX +BYAS2*FAC3;

        //RETURN
    }//END

    //!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    static void FIALCOS(double const R,double const THETA,double const PHI,double *BTHETA,double *BPHI,int const N,double const THETA0,double const DT) {
        /*
         C  CONICAL MODEL OF BIRKELAND CURRENT FIELD; BASED ON THE OLD S/R FIALCO (OF 1990-91)
         C  NB OF 1985-86-88, NOTE OF MARCH 5, BUT HERE BOTH INPUT AND OUTPUT ARE IN SPHERICAL CDS.

         C  BTN, AND BPN ARE THE ARRAYS OF BTHETA AND BPHI (BTN(i), BPN(i) CORRESPOND TO i-th MODE).
         C   ONLY FIRST  N  MODE AMPLITUDES ARE COMPUTED (N< = 10).
         C    THETA0 IS THE ANGULAR HALF-WIDTH OF THE CONE, DT IS THE ANGULAR H.-W. OF THE CURRENT LAYER

         C   NOTE:  BR = 0  (BECAUSE ONLY RADIAL CURRENTS ARE PRESENT IN THIS MODEL)
         */
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION  BTN(10),BPN(10),CCOS(10),SSIN(10)
        double BTN[11],BPN[11],CCOS[11],SSIN[11];
        double SINTE=sin(THETA);
        double RO=R*SINTE;
        double COSTE=cos(THETA);
        double SINFI=sin(PHI);
        double COSFI=cos(PHI);
        double TG=SINTE/(1.e0+COSTE);   //        TAN(THETA/2)
        double CTG=SINTE/(1.e0-COSTE);  //        COT(THETA/2)

        double TGP,TGM,TGM2,TGP2;
        double TETANP=THETA0+DT;
        double TETANM=THETA0-DT;
        if (THETA>=TETANM) {//IF(THETA.LT.TETANM) GOTO 1
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
        double T,DTT,DTT0,FC,FC1;
        for (int M=1; M<=N; M++) {//DO 2 M=1,N
            TM=TM*TG;
            CCOS[M]=COSM1*COSFI-SINM1*SINFI;
            SSIN[M]=SINM1*COSFI+COSM1*SINFI;
            COSM1=CCOS[M];
            SINM1=SSIN[M];
            if (THETA<TETANM) {
                T=TM;
                DTT=0.5e0*M*TM*(TG+CTG);
                DTT0=0.e0;
            } else if (THETA<TETANP) {//ELSE IF(THETA.LT.TETANP) THEN
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

            BTN[M]=M*T*CCOS[M]/RO;
            BPN[M]=-DTT*SSIN[M]/R; //         2
        }
        *BTHETA=BTN[N] *800.;
        *BPHI  =BPN[N] *800.;

        //RETURN
    }//END

    //-------------------------------------------------------------------------

    static void ONE_CONE(double const A[],double const X,double const Y,double const Z,double *BX,double *BY,double *BZ, int const* MODENUM_C, double const* DTHETA_C) {
        /*
         c  RETURNS FIELD COMPONENTS FOR A DEFORMED CONICAL CURRENT SYSTEM, FITTED TO A BIOSAVART FIELD
         c  HERE ONLY THE NORTHERN CONE IS TAKEN INTO ACCOUNT.
         */
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION A(31)

        double const* DTHETA = DTHETA_C;// COMMON /DTHETA/ DTHETA
        int const* M = MODENUM_C;// COMMON /MODENUM/ M

        static double const DR = 1.e-6,DT = 1.e-6;// JUST FOR NUMERICAL DIFFERENTIATION

        double THETA0 = A[31];

        double RHO2 = X*X+Y*Y;
        double RHO = sqrt(RHO2);
        double R = sqrt(RHO2+Z*Z);
        double THETA = atan2(RHO,Z);
        double PHI = atan2(Y,X);
        /*
         C   MAKE THE DEFORMATION OF COORDINATES:
         */
        double RS = R_S(A,R,THETA);
        double THETAS = THETA_S(A,R,THETA);
        double PHIS = PHI;

        //   CALCULATE FIELD COMPONENTS AT THE NEW POSITION (ASTERISKED):
        double BTAST,BFAST;
        FIALCOS (RS,THETAS,PHIS,&BTAST,&BFAST,*M,THETA0,*DTHETA); // MODE #M

        /*   NOW TRANSFORM B{R,T,F}_AST BY THE DEFORMATION TENSOR:
         C
         C      FIRST OF ALL, FIND THE DERIVATIVES: */

        double DRSDR = (R_S(A,R+DR,THETA)-R_S(A,R-DR,THETA))/(2.e0*DR);
        double DRSDT = (R_S(A,R,THETA+DT)-R_S(A,R,THETA-DT))/(2.e0*DT);
        double DTSDR = (THETA_S(A,R+DR,THETA)-THETA_S(A,R-DR,THETA))/(2.e0*DR);
        double DTSDT = (THETA_S(A,R,THETA+DT)-THETA_S(A,R,THETA-DT))/(2.e0*DT);

        double STSST = sin(THETAS)/sin(THETA);
        double RSR = RS/R;

        double BR = -RSR/R*STSST*BTAST*DRSDT;// NB#6, P.43    BRAST DOES NOT ENTER HERE
        double BTHETA = RSR*STSST*BTAST*DRSDR;// (IT IS IDENTICALLY ZERO IN OUR CASE)
        double BPHI = RSR*BFAST*(DRSDR*DTSDT-DRSDT*DTSDR);

        double S = RHO/R;
        double C = Z/R;
        double SF = Y/RHO;
        double CF = X/RHO;

        double BE = BR*S+BTHETA*C;

        *BX = A[1]*(BE*CF-BPHI*SF);
        *BY = A[1]*(BE*SF+BPHI*CF);
        *BZ = A[1]*(BR*C-BTHETA*S);

        //RETURN
    }//END

    //==================================================================================

    static void  TWOCONES (double const A[],double const X,double const Y,double const Z,double *BX,double *BY,double *BZ, int const* MODENUM_C, double const* DTHETA_C) {
        /*
         C  ADDS FIELDS FROM TWO CONES (NORTHERN AND SOUTHERN), WITH A PROPER SYMMETRY
         C  OF THE CURRENT AND FIELD, CORRESPONDING TO THE REGION 1 BIRKELAND CURRENTS. (NB #6, P.58).
         */
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION A(31)
        double BXN,BYN,BZN, BXS,BYS,BZS;
        ONE_CONE (A,X,Y,Z,&BXN,&BYN,&BZN, MODENUM_C, DTHETA_C);
        ONE_CONE (A,X,-Y,-Z,&BXS,&BYS,&BZS, MODENUM_C, DTHETA_C);
        *BX = BXN-BXS;
        *BY = BYN+BYS;
        *BZ = BZN+BZS;

        //RETURN
    } // END

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    static void BIRK_1N2 (int const NUMB,int const MODE,double const PS,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ, struct DPHI_B_RHO0_COMMON * DPHI_B_RHO0_C) {
        /*
         C  CALCULATES COMPONENTS  OF REGION 1/2 FIELD IN SPHERICAL COORDS.  DERIVED FROM THE S/R DIPDEF2C (WHICH
         C    DOES THE SAME JOB, BUT INPUT/OUTPUT THERE WAS IN SPHERICAL COORDS, WHILE HERE WE USE CARTESIAN ONES)
         C
         C   INPUT:  NUMB = 1 (2) FOR REGION 1 (2) CURRENTS
         C           MODE = 1 YIELDS SIMPLE SINUSOIDAL MLT VARIATION, WITH MAXIMUM CURRENT AT DAWN/DUSK MERIDIAN
         C     WHILE MODE = 2 YIELDS THE SECOND HARMONIC.
         C
         */
        int MODENUM_C; // IMPLICIT REAL*8 (A-H,O-Z)
        double DTHETA_C; // DIMENSION A11(31),A12(31),A21(31),A22(31)
        int *M = &MODENUM_C;// COMMON /MODENUM/ M
        double *DTHETA = &DTHETA_C;// COMMON /DTHETA/ DTHETA

        double * DPHI = &DPHI_B_RHO0_C->DPHI, *B = &DPHI_B_RHO0_C->B, *RHO_0 = &DPHI_B_RHO0_C->RHO_0, *XKAPPA = &DPHI_B_RHO0_C->XKAPPA;// COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! THESE PARAMETERS CONTROL DAY-NIGHT ASYMMETRY OF F.A.C., AS FOLLOWS:
        /*
         C  (1) DPHI:   HALF-DIFFERENCE (IN RADIANS) BETWEEN DAY AND NIGHT LATITUDE OF FAC OVAL AT IONOSPHERIC ALTITUDE;
         C              TYPICAL VALUE: 0.06
         C  (2) B:      AN ASYMMETRY FACTOR AT HIGH-ALTITUDES;  FOR B = 0, THE ONLY ASYMMETRY IS THAT FROM DPHI
         C              TYPICAL VALUES: 0.35-0.70
         C  (3) RHO_0:  A FIXED PARAMETER, DEFINING THE DISTANCE RHO, AT WHICH THE LATITUDE SHIFT GRADUALLY SATURATES AND
         C              STOPS INCREASING
         C              ITS VALUE WAS ASSUMED FIXED, EQUAL TO 7.0.
         C  (4) XKAPPA: AN OVERALL SCALING FACTOR, WHICH CAN BE USED FOR CHANGING THE SIZE OF THE F.A.C. OVAL
         */
        static double const BETA = 0.9e0,RH = 10.e0,EPS = 3.e0;// ! parameters of the tilt-dependent deformation of the untilted F.A.C. field

        static double const A11[32] = {0./*dummy*/, .1618068350,-.1797957553,2.999642482,-.9322708978,
            -.6811059760,.2099057262,-8.358815746,-14.86033550,.3838362986,
            -16.30945494,4.537022847,2.685836007,27.97833029,6.330871059,
            1.876532361,18.95619213,.9651528100,.4217195118,-.08957770020,
            -1.823555887,.7457045438,-.5785916524,-1.010200918,.01112389357,
            .09572927448,-.3599292276,8.713700514,.9763932955,3.834602998,
            2.492118385,.7113544659};
        static double const A12[32] = {0./*dummy*/, .7058026940,-.2845938535,5.715471266,-2.472820880,
            -.7738802408,.3478293930,-11.37653694,-38.64768867,.6932927651,
            -212.4017288,4.944204937,3.071270411,33.05882281,7.387533799,
            2.366769108,79.22572682,.6154290178,.5592050551,-.1796585105,
            -1.654932210,.7309108776,-.4926292779,-1.130266095,-.009613974555,
            .1484586169,-.2215347198,7.883592948,.02768251655,2.950280953,
            1.212634762,.5567714182};
        static double const A21[32] = {0./*dummy*/, .1278764024,-.2320034273,1.805623266,-32.37241440,
            -.9931490648,.3175085630,-2.492465814,-16.21600096,.2695393416,
            -6.752691265,3.971794901,14.54477563,41.10158386,7.912889730,
            1.258297372,9.583547721,1.014141963,.5104134759,-.1790430468,
            -1.756358428,.7561986717,-.6775248254,-.04014016420,.01446794851,
            .1200521731,-.2203584559,4.508963850,.8221623576,1.779933730,
            1.102649543,.8867880020};
        static double const A22[32] = {0./*dummy*/, .4036015198,-.3302974212,2.827730930,-45.44405830,
            -1.611103927,.4927112073,-.003258457559,-49.59014949,.3796217108,
            -233.7884098,4.312666980,18.05051709,28.95320323,11.09948019,
            .7471649558,67.10246193,.5667096597,.6468519751,-.1560665317,
            -1.460805289,.7719653528,-.6658988668,.2515179349E-05,
            .02426021891,.1195003324,-.2625739255,4.377172556,.2421190547,
            2.503482679,1.071587299,.7247997430};

        *B = 0.5;
        *RHO_0 = 7.0;

        *M = MODE;
        if (NUMB==1) {
            *DPHI = 0.055e0;
            *DTHETA = 0.06e0;
        } // ENDIF

        if (NUMB==2) {
            *DPHI = 0.030e0;
            *DTHETA = 0.09e0;
        } // ENDIF

        double XSC = X**XKAPPA;
        double YSC = Y**XKAPPA;
        double ZSC = Z**XKAPPA;
        double RHO = sqrt(XSC*XSC+ZSC*ZSC);

        double RSC = sqrt(XSC*XSC+YSC*YSC+ZSC*ZSC);// SCALED
        double RHO2 = *RHO_0**RHO_0;
        double PHI;
        if (XSC==0.e0 && ZSC==0.e0) {
            PHI = 0.e0;
        } else { // ELSE
            PHI = atan2(-ZSC,XSC);// FROM CARTESIAN TO CYLINDRICAL (RHO,PHI,Y)
        } // ENDIF

        double SPHIC = sin(PHI);
        double CPHIC = cos(PHI);// "C" means "CYLINDRICAL", TO DISTINGUISH FROM SPHERICAL PHI

        double BRACK = *DPHI+*B*RHO2/(RHO2+1.e0)*(RHO*RHO-1.e0)/(RHO2+RHO*RHO);
        double R1RH = (RSC-1.e0)/RH;
        double PSIAS = BETA*PS/pow((1.e0+pow(R1RH, EPS)), (1.e0/EPS));

        double PHIS = PHI-BRACK*sin(PHI) -PSIAS;
        double DPHISPHI = 1.e0-BRACK*cos(PHI);
        double DPHISRHO = -2.e0**B*RHO2*RHO/pow((RHO2+RHO*RHO), 2) *sin(PHI)
        +BETA*PS*pow(R1RH, (EPS-1.e0))*RHO/(RH*RSC*
                                            pow((1.e0+pow(R1RH,EPS)), (1.e0/EPS+1.e0)));
        double DPHISDY = BETA*PS*pow(R1RH, (EPS-1.e0))*YSC/(RH*RSC*
                                                            pow((1.e0+pow(R1RH, EPS)), (1.e0/EPS+1.e0)));

        double SPHICS = sin(PHIS);
        double CPHICS = cos(PHIS);

        double XS = RHO*CPHICS;
        double ZS = -RHO*SPHICS;
        double BXS,BYAS,BZS;
        if (NUMB==1) {
            if (MODE==1) TWOCONES (A11,XS,YSC,ZS,&BXS,&BYAS,&BZS, &MODENUM_C, &DTHETA_C);
            if (MODE==2) TWOCONES (A12,XS,YSC,ZS,&BXS,&BYAS,&BZS, &MODENUM_C, &DTHETA_C);
        } else { // ELSE
            if (MODE==1) TWOCONES (A21,XS,YSC,ZS,&BXS,&BYAS,&BZS, &MODENUM_C, &DTHETA_C);
            if (MODE==2) TWOCONES (A22,XS,YSC,ZS,&BXS,&BYAS,&BZS, &MODENUM_C, &DTHETA_C);
        } // ENDIF

        double BRHOAS = BXS*CPHICS-BZS*SPHICS;
        double BPHIAS = -BXS*SPHICS-BZS*CPHICS;

        double BRHO_S = BRHOAS*DPHISPHI                             **XKAPPA;// SCALING
        double BPHI_S = (BPHIAS-RHO*(BYAS*DPHISDY+BRHOAS*DPHISRHO)) **XKAPPA;
        double BY_S = BYAS*DPHISPHI                                 **XKAPPA;

        *BX = BRHO_S*CPHIC-BPHI_S*SPHIC;
        *BY = BY_S;
        *BZ = -BRHO_S*SPHIC-BPHI_S*CPHIC;

        //RETURN
    } // END

    //-------------------------------------------------------------------------


    static void BIRK_SHL (double const A[],double const PS,double const X_SC,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {

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

        for (int M=1; M<=2; M++) {//DO 1 M=1,2     !    M=1 IS FOR THE 1ST SUM ("PERP." SYMMETRY)
            //                          AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
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

                    for (int N=1; N<=2; N++) {//DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
                        //                                AND N=2 IS FOR THE SECOND ONE

                        for (int NN=1; NN<=2; NN++) {//DO 5 NN=1,2 !   NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
                            //                                         TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE
                            double HX,HY,HZ,FX,FY,FZ;
                            if (M==1) {
                                FX=-SQPR*EPR*CYPI*SZRK;
                                FY=EPR*SYPI*SZRK/P;
                                FZ=-EPR*CYPI*CZRK/R;
                                if (N==1) {
                                    if (NN==1) {
                                        HX=FX;
                                        HY=FY;
                                        HZ=FZ;
                                    } else {//ELSE
                                        HX=FX*X_SC;
                                        HY=FY*X_SC;
                                        HZ=FZ*X_SC;
                                    }//ENDIF
                                } else {//ELSE
                                    if (NN==1) {
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
                                if (N==1) {
                                    if (NN==1) {
                                        HX=FX;
                                        HY=FY;
                                        HZ=FZ;
                                    } else {//ELSE
                                        HX=FX*X_SC;
                                        HY=FY*X_SC;
                                        HZ=FZ*X_SC;
                                    }//ENDIF
                                } else {//ELSE
                                    if (NN==1) {
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
                            if (M==1) {
                                HXR=HX*CT1+HZ*ST1;
                                HZR=-HX*ST1+HZ*CT1;
                            } else {//ELSE
                                HXR=HX*CT2+HZ*ST2;
                                HZR=-HX*ST2+HZ*CT2;
                            }//ENDIF

                            GX=GX+HXR*A[L];
                            GY=GY+HY *A[L];
                            GZ=GZ+HZR*A[L]; //       5
                        }
                    }//4   CONTINUE
                }//3   CONTINUE
            }//2   CONTINUE
        }//1   CONTINUE

        *BX=GX;
        *BY=GY;
        *BZ=GZ;

        //RETURN
    }//END

    // %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    static void BIRK_TOT (int const IOPB,double const PS,double const X,double const Y,double const Z,double *BX11,double *BY11,double *BZ11,double *BX12,double *BY12,double *BZ12,double *BX21,double *BY21,double *BZ21,double *BX22,double *BY22,double *BZ22, struct BIRKPAR_COMMON const* BIRKPAR_C) {
        /*
         C      IOPB -  BIRKELAND FIELD MODE FLAG:
         C         IOPB = 0 - ALL COMPONENTS
         C         IOPB = 1 - REGION 1, MODES 1 & 2
         C         IOPB = 2 - REGION 2, MODES 1 & 2
         */
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION SH11(86),SH12(86),SH21(86),SH22(86)
        double const* XKAPPA1 = &BIRKPAR_C->XKAPPA1, *XKAPPA2 = &BIRKPAR_C->XKAPPA2;// COMMON /BIRKPAR/ XKAPPA1,XKAPPA2   !  INPUT PARAMETERS, SPECIFIED FROM S/R EXTALL
        struct DPHI_B_RHO0_COMMON DPHI_B_RHO0_C; // COMMON /DPHI_B_RHO0/ DPHI,B,RHO_0,XKAPPA ! PARAMETERS, CONTROLLING THE DAY-NIGHT ASYMMETRY OF F.A.C.
        //double &DPHI = DPHI_B_RHO0_C.DPHI, &B = DPHI_B_RHO0_C.B, &RHO_0 = DPHI_B_RHO0_C.RHO_0;
        double *XKAPPA = &DPHI_B_RHO0_C.XKAPPA;

        static double const SH11[87] = {0./*dummy*/, 46488.84663,-15541.95244,-23210.09824,-32625.03856,
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

        static double const SH12[87] = {0./*dummy*/, 210260.4816,-1443587.401,-1468919.281,281939.2993,
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

        static double const SH21[87] = {0./*dummy*/, 162294.6224,503885.1125,-27057.67122,-531450.1339,
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

        static double const SH22[87] = {0./*dummy*/, -131287.8986,-631927.6885,-318797.4173,616785.8782,
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

        *XKAPPA = *XKAPPA1;// FORWARDED IN BIRK_1N2
        double X_SC = *XKAPPA1-1.1e0;// FORWARDED IN BIRK_SHL

        if (IOPB==0 || IOPB==1) {
            double FX11,FY11,FZ11, HX11,HY11,HZ11;
            BIRK_1N2 (1,1,PS,X,Y,Z,&FX11,&FY11,&FZ11, &DPHI_B_RHO0_C);// REGION 1, MODE 1
            BIRK_SHL (SH11,PS,X_SC,X,Y,Z,&HX11,&HY11,&HZ11);
            *BX11 = FX11+HX11;
            *BY11 = FY11+HY11;
            *BZ11 = FZ11+HZ11;
            double FX12,FY12,FZ12,HX12,HY12,HZ12;
            BIRK_1N2 (1,2,PS,X,Y,Z,&FX12,&FY12,&FZ12, &DPHI_B_RHO0_C); // REGION 1, MODE 2
            BIRK_SHL (SH12,PS,X_SC,X,Y,Z,&HX12,&HY12,&HZ12);
            *BX12 = FX12+HX12;
            *BY12 = FY12+HY12;
            *BZ12 = FZ12+HZ12;

        } // ENDIF

        *XKAPPA = *XKAPPA2;// FORWARDED IN BIRK_1N2
        X_SC = *XKAPPA2-1.0e0;// FORWARDED IN BIRK_SHL

        if (IOPB==0 || IOPB==2) {
            double FX21,FY21,FZ21,HX21,HY21,HZ21;
            BIRK_1N2 (2,1,PS,X,Y,Z,&FX21,&FY21,&FZ21, &DPHI_B_RHO0_C); // REGION 2, MODE 1
            BIRK_SHL (SH21,PS,X_SC,X,Y,Z,&HX21,&HY21,&HZ21);
            *BX21 = FX21+HX21;
            *BY21 = FY21+HY21;
            *BZ21 = FZ21+HZ21;
            double FX22,FY22,FZ22,HX22,HY22,HZ22;
            BIRK_1N2 (2,2,PS,X,Y,Z,&FX22,&FY22,&FZ22, &DPHI_B_RHO0_C); // REGION 2, MODE 2
            BIRK_SHL (SH22,PS,X_SC,X,Y,Z,&HX22,&HY22,&HZ22);
            *BX22 = FX22+HX22;
            *BY22 = FY22+HY22;
            *BZ22 = FZ22+HZ22;

        } // ENDIF

        //RETURN
    } // END

    //&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&

    static void RC_SYMM (double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {
        //IMPLICIT  REAL * 8  (A - H, O - Z)
        static double const DS=1.e-2,DC=0.99994999875e0, D=1.e-4,DRD=5.e3;  // DS=SIN(THETA) AT THE BOUNDARY OF THE LINEARITY
        //     REGION; DC=sqrt(1-DS**2);  DRD=1/(2*D)
        double RHO2=X*X+Y*Y;
        double R2=RHO2+Z*Z;
        double R=sqrt(R2);
        double RP=R+D;
        double RM=R-D;
        double SINT=sqrt(RHO2)/R;
        double COST=Z/R;

        if (SINT<DS) {  //  TOO CLOSE TO THE Z-AXIS; USING A LINEAR APPROXIMATION A_PHI~SINT,
            //                                    TO AVOID THE SINGULARITY PROBLEM
            double A=AP(R,DS,DC)/DS;
            double DARDR=(RP*AP(RP,DS,DC)-RM*AP(RM,DS,DC))*DRD;
            double FXY=Z*(2.e0*A-DARDR)/(R*R2);
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
            double FXY=(BR+BT*COST/SINT)/R;
            *BX=FXY*X;
            *BY=FXY*Y;
            *BZ=BR*COST-BT*SINT;

        }//ENDIF

        //RETURN
    }//END

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    static void PRC_SYMM (double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {
        //IMPLICIT  REAL * 8  (A - H, O - Z)
        static double const DS=1.e-2,DC=0.99994999875e0, D=1.e-4,DRD=5.e3;  // DS=SIN(THETA) AT THE BOUNDARY OF THE LINEARITY
        //    REGION; DC=sqrt(1-DS**2);  DRD=1/(2*D)
        double RHO2=X*X+Y*Y;
        double R2=RHO2+Z*Z;
        double R=sqrt(R2);
        double RP=R+D;
        double RM=R-D;
        double SINT=sqrt(RHO2)/R;
        double COST=Z/R;

        if (SINT<DS) {  //  TOO CLOSE TO THE Z-AXIS; USING A LINEAR APPROXIMATION A_PHI~SINT,
            //                                    TO AVOID THE SINGULARITY PROBLEM
            double A=APPRC(R,DS,DC)/DS;
            double DARDR=(RP*APPRC(RP,DS,DC)-RM*APPRC(RM,DS,DC))*DRD;
            double FXY=Z*(2.e0*A-DARDR)/(R*R2);
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
            double BR=(SINTP*APPRC(R,SINTP,COSTP)-SINTM*APPRC(R,SINTM,COSTM))
            /(R*SINT)*DRD;
            double BT=(RM*APPRC(RM,SINT,COST)-RP*APPRC(RP,SINT,COST))/R*DRD;
            double FXY=(BR+BT*COST/SINT)/R;
            *BX=FXY*X;
            *BY=FXY*Y;
            *BZ=BR*COST-BT*SINT;

        }//ENDIF

        //RETURN
    }//END

    //@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@

    static void PRC_QUAD (double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {

        //  CALCULATES COMPONENTS OF THE FIELD FROM THE "QUADRUPOLE" COMPONENT OF THE PRC

        //IMPLICIT  REAL * 8  (A - H, O - Z)

        static double const D=1.e-4,DD=2.e-4, DS=1.e-2,DC=0.99994999875e0;

        double RHO2=X*X+Y*Y;
        double R=sqrt(RHO2+Z*Z);
        double RHO=sqrt(RHO2);
        double SINT=RHO/R;
        double COST=Z/R;
        double RP=R+D;
        double RM=R-D;

        if (SINT>DS) {
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
            if (Z<0.e0) CT=-DC;
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

        //RETURN
    }//END

    //---------------------------------------------------------------------------------------

    static void SRC_PRC (int const IOPR,double const SC_SY,double const SC_PR,double const PHI,double const PS,double const X,double const Y,double const Z,double *BXSRC,double *BYSRC,double *BZSRC,double *BXPRC,double *BYPRC,double *BZPRC) {
        /*
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
         IMPLICIT REAL*8 (A-H,O-Z)
         c
         c   1.  TRANSFORM TO TILTED COORDINATES (i.e., SM coordinates):
         */
        double CPS=cos(PS);
        double SPS=sin(PS);

        double XT=X*CPS-Z*SPS;
        double ZT=Z*CPS+X*SPS;
        /*
         C   2.  SCALE THE COORDINATES FOR THE SYMMETRIC AND PARTIAL RC COMPONENTS:
         */
        double XTS=XT/SC_SY;    //  SYMMETRIC
        double YTS=Y /SC_SY;
        double ZTS=ZT/SC_SY;

        double XTA=XT/SC_PR;   //  PARTIAL
        double YTA=Y /SC_PR;
        double ZTA=ZT/SC_PR;
        /*
         C   3.  CALCULATE COMPONENTS OF THE TOTAL FIELD IN THE TILTED (SOLAR-MAGNETIC) COORDINATE SYSTEM:
         C
         C==========   ONLY FOR LEAST SQUARES FITTING: */
        double BXS=0.e0;
        double BYS=0.e0;
        double BZS=0.e0;
        double BXA_S=0.e0;
        double BYA_S=0.e0;
        double BZA_S=0.e0;
        double BXA_QR=0.e0;
        double BYA_QR=0.e0;
        double BZA_Q=0.e0;
        /*============================================
         C
         C    3a. SYMMETRIC FIELD:
         */
        if (IOPR<=1) RC_SYMM(XTS,YTS,ZTS,&BXS,&BYS,&BZS);
        if (IOPR==0 || IOPR==2) PRC_SYMM(XTA,YTA,ZTA,&BXA_S,&BYA_S,&BZA_S);

        /*    3b. ROTATE THE SCALED SM COORDINATES BY PHI AROUND ZSM AXIS AND CALCULATE QUADRUPOLE PRC FIELD
         C         IN THOSE COORDS:*/

        double CP=cos(PHI);
        double SP=sin(PHI);
        double XR=XTA*CP-YTA*SP;
        double YR=XTA*SP+YTA*CP;

        if (IOPR==0 || IOPR==2)
            PRC_QUAD(XR,YR,ZTA,&BXA_QR,&BYA_QR,&BZA_Q);

        //    3c. TRANSFORM THE QUADRUPOLE FIELD COMPONENTS BACK TO THE SM COORDS:

        double BXA_Q= BXA_QR*CP+BYA_QR*SP;
        double BYA_Q=-BXA_QR*SP+BYA_QR*CP;

        //    3d. FIND THE TOTAL FIELD OF PRC (SYMM.+QUADR.) IN THE SM COORDS:

        double BXP=BXA_S+BXA_Q;
        double BYP=BYA_S+BYA_Q;
        double BZP=BZA_S+BZA_Q;

        //   4.  TRANSFORM THE FIELDS OF BOTH PARTS OF THE RING CURRENT BACK TO THE GSM SYSTEM:

        *BXSRC=BXS*CPS+BZS*SPS;   //    SYMMETRIC RC
        *BYSRC=BYS;
        *BZSRC=BZS*CPS-BXS*SPS;

        *BXPRC=BXP*CPS+BZP*SPS;   //    PARTIAL RC
        *BYPRC=BYP;
        *BZPRC=BZP*CPS-BXP*SPS;

        //RETURN
    }//END

    //||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||||

    static void RC_SHIELD (double const A[],double const PS,double const X_SC,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {
        /*
         C   COMPUTES THE COMPONENTS OF THE SHIELDING FIELD FOR THE RING CURRENT
         C       (EITHER PARTIAL OR AXISYMMETRICAL)
         C   INPUT:   A - AN ARRAY CONTAINING THE HARMONIC COEFFICIENTS AND NONLINEAR PARAMETERS
         C            PS - GEODIPOLE TILT ANGLE IN RADIANS
         C            X_SC - SCALING FACTOR ( X_SC>1 AND X_SC<1 CORRESPOND TO LARGER/SMALLER
         C                  RING CURRENT, RESP.)
         C            X,Y,Z - POSITION IN RE (GSM COORDS)
         C   OUTPUT:  BX,BY,BZ - SHIELDING FIELD COMPONENTS (GSM)
         */
        //IMPLICIT  REAL * 8  (A - H, O - Z)
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
            //                           AND M=2 IS FOR THE SECOND SUM ("PARALL." SYMMETRY)
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

                    for (int N=1; N<=2; N++) {//DO 4 N=1,2  ! N=1 IS FOR THE FIRST PART OF EACH COEFFICIENT
                        //                                AND N=2 IS FOR THE SECOND ONE

                        for (int NN=1; NN<=2; NN++) {//DO 5 NN=1,2 !   NN = 1,2 FURTHER SPLITS THE COEFFICIENTS INTO 2 PARTS,
                            //                                         TO TAKE INTO ACCOUNT THE SCALE FACTOR DEPENDENCE
                            double FX, FY, FZ, HX,HY,HZ;
                            if (M==1) {
                                FX=-SQPR*EPR*CYPI*SZRK  *FAC_SC;
                                FY=EPR*SYPI*SZRK/P   *FAC_SC;
                                FZ=-EPR*CYPI*CZRK/R  *FAC_SC;
                                if (N==1) {
                                    if (NN==1) {
                                        HX=FX;
                                        HY=FY;
                                        HZ=FZ;
                                    } else {//ELSE
                                        HX=FX*X_SC;
                                        HY=FY*X_SC;
                                        HZ=FZ*X_SC;
                                    }//ENDIF
                                } else {//ELSE
                                    if (NN==1) {
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
                                if (N==1) {
                                    if (NN==1) {
                                        HX=FX;
                                        HY=FY;
                                        HZ=FZ;
                                    } else {//ELSE
                                        HX=FX*X_SC;
                                        HY=FY*X_SC;
                                        HZ=FZ*X_SC;
                                    }//ENDIF
                                } else {//ELSE
                                    if (NN==1) {
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
                            double HXR,HZR;
                            if (M==1) {
                                HXR=HX*CT1+HZ*ST1;
                                HZR=-HX*ST1+HZ*CT1;
                            } else {//ELSE
                                HXR=HX*CT2+HZ*ST2;
                                HZR=-HX*ST2+HZ*CT2;
                            }//ENDIF

                            GX=GX+HXR*A[L];
                            GY=GY+HY *A[L];
                            GZ=GZ+HZR*A[L]; //       5
                        }
                    }//4   CONTINUE
                }//3   CONTINUE
            }//2   CONTINUE
        }//1   CONTINUE

        *BX=GX;
        *BY=GY;
        *BZ=GZ;

        //RETURN
    }//END

    //************************************************************************************

    static void FULL_RC (int const IOPR,double const PS,double const X,double const Y,double const Z,double *BXSRC,double *BYSRC,double *BZSRC,double *BXPRC,double *BYPRC,
                         double *BZPRC, struct RCPAR_COMMON const* RCPAR_C) {
        /*
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
         */
        //IMPLICIT REAL*8 (A-H,O-Z)
        //DIMENSION C_SY(86),C_PR(86)
        double const* SC_SY=&RCPAR_C->SC_SY, *SC_PR=&RCPAR_C->SC_AS, *PHI=&RCPAR_C->PHI;//COMMON /RCPAR/ SC_SY,SC_PR,PHI

        static double const C_SY[87] = {0./*dummy*/, -957.2534900,-817.5450246,583.2991249,758.8568270,     //   CORRECTED VALUES (AS OF MAY 2006)
            13.17029064,68.94173502,-15.29764089,-53.43151590,27.34311724,
            149.5252826,-11.00696044,-179.7031814,953.0914774,817.2340042,
            -581.0791366,-757.5387665,-13.10602697,-68.58155678,15.22447386,
            53.15535633,-27.07982637,-149.1413391,10.91433279,179.3251739,
            -6.028703251,1.303196101,-1.345909343,-1.138296330,-0.06642634348,
            -0.3795246458,.07487833559,.2891156371,-.5506314391,-.4443105812,
            0.2273682152,0.01086886655,-9.130025352,1.118684840,1.110838825,
            .1219761512,-.06263009645,-.1896093743,.03434321042,.01523060688,
            -.4913171541,-.2264814165,-.04791374574,.1981955976,-68.32678140,
            -48.72036263,14.03247808,16.56233733,2.369921099,6.200577111,
            -1.415841250,-0.8184867835,-3.401307527,-8.490692287,3.217860767,
            -9.037752107,66.09298105,48.23198578,-13.67277141,-16.27028909,
            -2.309299411,-6.016572391,1.381468849,0.7935312553,3.436934845,
            8.260038635,-3.136213782,8.833214943,8.041075485,8.024818618,
            35.54861873,12.55415215,1.738167799,3.721685353,23.06768025,
            6.871230562,6.806229878,21.35990364,1.687412298,3.500885177,
            0.3498952546,0.6595919814};

        static double const C_PR[87] = {0./*dummy*/, -64820.58481,-63965.62048,66267.93413,135049.7504,
            -36.56316878,124.6614669,56.75637955,-87.56841077,5848.631425,
            4981.097722,-6233.712207,-10986.40188,68716.52057,65682.69473,
            -69673.32198,-138829.3568,43.45817708,-117.9565488,-62.14836263,
            79.83651604,-6211.451069,-5151.633113,6544.481271,11353.03491,
            23.72352603,-256.4846331,25.77629189,145.2377187,-4.472639098,
            -3.554312754,2.936973114,2.682302576,2.728979958,26.43396781,
            -9.312348296,-29.65427726,-247.5855336,-206.9111326,74.25277664,
            106.4069993,15.45391072,16.35943569,-5.965177750,-6.079451700,
            115.6748385,-35.27377307,-32.28763497,-32.53122151,93.74409310,
            84.25677504,-29.23010465,-43.79485175,-6.434679514,-6.620247951,
            2.443524317,2.266538956,-43.82903825,6.904117876,12.24289401,
            17.62014361,152.3078796,124.5505289,-44.58690290,-63.02382410,
            -8.999368955,-9.693774119,3.510930306,3.770949738,-77.96705716,
            22.07730961,20.46491655,18.67728847,9.451290614,9.313661792,
            644.7620970,418.2515954,7.183754387,35.62128817,19.43180682,
            39.57218411,15.69384715,7.123215241,2.300635346,21.90881131,
            -.01775839370,.3996346710};

        double HXSRC,HYSRC,HZSRC, HXPRC,HYPRC,HZPRC;
        SRC_PRC (IOPR,*SC_SY,*SC_PR,*PHI,PS,X,Y,Z,&HXSRC,&HYSRC,&HZSRC,
                 &HXPRC,&HYPRC,&HZPRC);
        double FSX,FSY,FSZ;
        double X_SC=*SC_SY-1.e0;
        if (IOPR==0 || IOPR==1) {
            RC_SHIELD (C_SY,PS,X_SC,X,Y,Z,&FSX,&FSY,&FSZ);
        } else {//ELSE
            FSX=0.e0;
            FSY=0.e0;
            FSZ=0.e0;
        }//ENDIF
        double FPX,FPY,FPZ;
        X_SC=*SC_PR-1.e0;
        if (IOPR==0 || IOPR==2) {
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

        //RETURN
    }//END

    //===========================================================================

    static void DIPOLE (double const PS,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {
        /*
         C     THIS IS A DOUBLE PRECISION ROUTINE, OTHERWISE IDENTICAL TO THE S/R DIP OF GEOPACK
         C
         C     CALCULATES GSM COMPONENTS OF A GEODIPOLE FIELD WITH THE DIPOLE MOMENT
         C     CORRESPONDING TO THE EPOCH OF 2000.
         C
         C------INPUT PARAMETERS:
         C       PS - GEODIPOLE TILT ANGLE IN RADIANS,
         C       X,Y,Z - GSM COORDINATES IN RE (1 RE = 6371.2 km)
         C
         C----OUTPUT PARAMETERS:
         C     BX,BY,BZ - FIELD COMPONENTS IN GSM SYSTEM, IN NANOTESLA.
         C
         C    LAST MODIFICATION: JAN. 5, 2001. THE VALUE OF THE DIPOLE MOMENT WAS UPDATED TO 2000.
         C      AND A "SAVE" STATEMENT HAS BEEN ADDED, TO AVOID POTENTIAL PROBLEMS WITH SOME
         C      FORTRAN COMPILERS
         C
         C    WRITTEN BY: N. A. TSYGANENKO
         */
        //IMPLICIT REAL*8 (A-H,O-Z)
        //SAVE M,PSI
        //DATA M,PSI/0,5.D0/
        //IF(M.EQ.1.AND.fabs(PS-PSI).LT.1.D-5) GOTO 1   !   THIS IS TO AVOID MULTIPLE CALCULATIONS
        double SPS=sin(PS);                                  //   OF SIN(PS) AND COS(PS), IF THE ANGLE PS
        double CPS=cos(PS);                                  //   REMAINS UNCHANGED
        //PSI=PS
        //M=1
        double P=X*X; //        1
        double U=Z*Z;
        double V=3.e0*Z*X;
        double T=Y*Y;
        double Q=30115.e0/pow(sqrt(P+T+U), 5);
        *BX=Q*((T+U-2.e0*P)*SPS-V*CPS);
        *BY=-3.e0*Y*Q*(X*SPS+Z*CPS);
        *BZ=Q*((P+T-2.e0*U)*CPS-V*SPS);
        //RETURN
    }//END

    //================================================================
    static void EXTALL (int const IOPGEN,int const IOPT,int const IOPB,int const IOPR,double const A[],int const NTOT,
                        double const PDYN,double const DST,double const BYIMF,double const BZIMF,double const VBIMF1,double const VBIMF2,double const PS,double const X,double const Y,double const Z,
                        double *BXCF,double *BYCF,double *BZCF,double *BXT1,double *BYT1,double *BZT1,double *BXT2,double *BYT2,double *BZT2,
                        double *BXSRC,double *BYSRC,double *BZSRC,double *BXPRC,double *BYPRC,double *BZPRC, double *BXR11,double *BYR11,double *BZR11,
                        double *BXR12,double *BYR12,double *BZR12,double *BXR21,double *BYR21,double *BZR21,double *BXR22,double *BYR22,double *BZR22,double *HXIMF,
                        double *HYIMF,double *HZIMF,double *BX,double *BY,double *BZ) {
        /*
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
         */
        //IMPLICIT  REAL * 8  (A - H, O - Z)

        //DIMENSION A(NTOT)
        struct TAIL_COMMON TAIL_C;
        double *DXSHIFT1=&TAIL_C.DXSHIFT1, *DXSHIFT2=&TAIL_C.DXSHIFT2, *D=&TAIL_C.D, *DELTADY=&TAIL_C.DELTADY;  // THE COMMON BLOCKS FORWARD NONLINEAR PARAMETERS
        struct BIRKPAR_COMMON BIRKPAR_C;
        double *XKAPPA1=&BIRKPAR_C.XKAPPA1, *XKAPPA2=&BIRKPAR_C.XKAPPA2;
        struct RCPAR_COMMON RCPAR_C;
        double *SC_SY=&RCPAR_C.SC_SY, *SC_AS=&RCPAR_C.SC_AS, *PHI=&RCPAR_C.PHI;
        struct G_COMMON G_C;
        double *G=&G_C.G;
        struct RH0_COMMON RH0_C;
        double *RH0=&RH0_C.RH0;

        static double const A0_A=34.586e0,A0_S0=1.1960e0,A0_X0=3.4397e0;   //   SHUE ET AL. PARAMETERS
        static double const DSIG=0.003e0,RH2=-5.2e0;

        *RH0=8.0e0;

        double XAPPA=pow((PDYN/2.), A[39]);   //  NOW THIS IS A VARIABLE PARAMETER
        *RH0=A[40];
        *G=A[41];

        double XAPPA3=XAPPA*XAPPA*XAPPA;

        double XX=X*XAPPA;
        double YY=Y*XAPPA;
        double ZZ=Z*XAPPA;

        double SPS=sin(PS);

        double X0=A0_X0/XAPPA;
        double AM=A0_A/XAPPA;
        double S0=A0_S0;

        //double BPERP=sqrt(BYIMF*BYIMF+BZIMF*BZIMF);

        //   CALCULATE THE IMF CLOCK ANGLE:
        double THETA;
        if (BYIMF==0.e0 && BZIMF==0.e0) {
            THETA=0.e0;
        } else {//ELSE
            THETA=atan2(BYIMF,BZIMF);
            if (THETA<=0.e0) THETA=THETA+6.283185307e0;
        }//ENDIF

        //double CT=cos(THETA);
        //double ST=sin(THETA);
        //double YS=Y*CT-Z*ST;
        //double ZS=Z*CT+Y*ST;

        double STHETAH=pow(sin(THETA/2.), 2);
        /*
         C  CALCULATE "IMF" COMPONENTS OUTSIDE THE MAGNETOPAUSE LAYER (HENCE BEGIN WITH "O")
         C  THEY ARE NEEDED ONLY IF THE POINT (X,Y,Z) IS WITHIN THE TRANSITION MAGNETOPAUSE LAYER
         C  OR OUTSIDE THE MAGNETOSPHERE:
         */
        double FACTIMF=A[24]+A[25]*STHETAH;

        double OIMFX=0.e0;
        double OIMFY=BYIMF*FACTIMF;
        double OIMFZ=BZIMF*FACTIMF;

        double R=sqrt(X*X+Y*Y+Z*Z);
        double XSS=X;
        double ZSS=Z;
        double DD = 0.;
        do {
            double XSOLD=XSS;      //  1   BEGIN ITERATIVE SEARCH OF UNWARPED COORDS (TO FIND SIGMA)
            double ZSOLD=ZSS;

            double RH=*RH0+RH2*pow((ZSS/R), 2);
            double SINPSAS=SPS/pow((1.e0+pow((R/RH), 3)), 0.33333333e0);
            double COSPSAS=sqrt(1.e0-SINPSAS*SINPSAS);
            ZSS=X*SINPSAS+Z*COSPSAS;
            XSS=X*COSPSAS-Z*SINPSAS;
            DD=fabs(XSS-XSOLD)+fabs(ZSS-ZSOLD);
        } while (DD>1.e-6);//IF (DD.GT.1.D-6) GOTO 1
        //                                END OF ITERATIVE SEARCH
        double RHO2=Y*Y+ZSS*ZSS;
        double ASQ=AM*AM;
        double XMXM=AM+XSS-X0;
        if (XMXM<0.) XMXM=0.; // THE BOUNDARY IS A CYLINDER TAILWARD OF X=X0-AM
        double AXX0=XMXM*XMXM;
        double ARO=ASQ+RHO2;
        double SIGMA=sqrt((ARO+AXX0+sqrt(pow((ARO+AXX0), 2)-4.*ASQ*AXX0))/(2.*ASQ));
        /*
         C   NOW, THERE ARE THREE POSSIBLE CASES:
         C    (1) INSIDE THE MAGNETOSPHERE   (SIGMA
         C    (2) IN THE BOUNDARY LAYER
         C    (3) OUTSIDE THE MAGNETOSPHERE AND B.LAYER
         C       FIRST OF ALL, CONSIDER THE CASES (1) AND (2):
         C
         C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
        if (SIGMA<S0+DSIG) {  //  CASES (1) OR (2); CALCULATE THE MODEL FIELD
            /*                              (WITH THE POTENTIAL "PENETRATED" INTERCONNECTION FIELD):
             C+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
             */
            if (IOPGEN<=1) {
                double CFX,CFY,CFZ;
                SHLCAR3X3(XX,YY,ZZ,PS,&CFX,&CFY,&CFZ);         //  DIPOLE SHIELDING FIELD
                *BXCF=CFX*XAPPA3;
                *BYCF=CFY*XAPPA3;
                *BZCF=CFZ*XAPPA3;
            } else {//ELSE
                *BXCF=0.e0;
                *BYCF=0.e0;
                *BZCF=0.e0;
            }//ENDIF
            
            if (IOPGEN==0 || IOPGEN==2) {
                *DXSHIFT1=A[26]+A[27]*VBIMF2;
                *DXSHIFT2=0.e0;
                *D=A[28];
                *DELTADY=A[29];
                DEFORMED (IOPT,PS,XX,YY,ZZ,                //  TAIL FIELD (THREE MODES)
                          BXT1,BYT1,BZT1,BXT2,BYT2,BZT2, &RH0_C, &G_C, &TAIL_C);
            } else {//ELSE
                *BXT1=0.e0;
                *BYT1=0.e0;
                *BZT1=0.e0;
                *BXT2=0.e0;
                *BYT2=0.e0;
                *BZT2=0.e0;
            }//ENDIF
            
            if (IOPGEN==0 || IOPGEN==3) {
                *XKAPPA1=A[35]+A[36]*VBIMF2;
                *XKAPPA2=A[37]+A[38]*VBIMF2;
                BIRK_TOT (IOPB,PS,XX,YY,ZZ,BXR11,BYR11,BZR11,BXR12,BYR12,
                          BZR12,BXR21,BYR21,BZR21,BXR22,BYR22,BZR22, &BIRKPAR_C);    //   BIRKELAND FIELD (TWO MODES FOR R1 AND TWO MODES FOR R2)
            } else {//ELSE
                *BXR11=0.e0;
                *BYR11=0.e0;
                *BZR11=0.e0;
                *BXR12=0.e0;
                *BYR12=0.e0;
                *BZR12=0.e0;
                *BXR21=0.e0;
                *BYR21=0.e0;
                *BZR21=0.e0;
                *BXR22=0.e0;
                *BYR22=0.e0;
                *BZR22=0.e0;
            }//ENDIF
            
            if (IOPGEN==0 || IOPGEN==4) {
                *PHI=1.5707963e0*tanh(fabs(DST)/A[34]);
                double ZNAM=fabs(DST);
                if (ZNAM<20.e0) ZNAM=20.e0;
                *SC_SY=A[30]*pow((20.e0/ZNAM), A[31]) *XAPPA;
                *SC_AS=A[32]*pow((20.e0/ZNAM), A[33]) *XAPPA;
                FULL_RC(IOPR,PS,XX,YY,ZZ,BXSRC,BYSRC,BZSRC,BXPRC,BYPRC,
                        BZPRC, &RCPAR_C);  //  SHIELDED RING CURRENT (SRC AND PRC)
            } else {//ELSE
                *BXSRC=0.e0;
                *BYSRC=0.e0;
                *BZSRC=0.e0;
                *BXPRC=0.e0;
                *BYPRC=0.e0;
                *BZPRC=0.e0;
            }//ENDIF
            
            if (IOPGEN==0 || IOPGEN==5) {
                *HXIMF=0.e0;
                *HYIMF=BYIMF;
                *HZIMF=BZIMF;   /* THESE ARE COMPONENTS OF THE PENETRATED FIELD PER UNIT OF THE PENETRATION COEFFICIENT.
                                 C                        IN OTHER WORDS, THESE ARE DERIVATIVES OF THE PENETRATION FIELD COMPONENTS WITH RESPECT
                                 C                        TO THE PENETRATION COEFFICIENT.   WE ASSUME THAT ONLY THE TRANSVERSE COMPONENT OF THE
                                 C                        FIELD PENETRATES INSIDE.*/
            } else {//ELSE
                *HXIMF=0.e0;
                *HYIMF=0.e0;
                *HZIMF=0.e0;
            }//ENDIF
            /*
             C-----------------------------------------------------------
             C
             C    NOW, ADD UP ALL THE COMPONENTS:*/
            
            double DLP1=pow((PDYN/2.e0), A[42]);
            double DLP2=pow((PDYN/2.e0), A[43]);
            
            double TAMP1=A[2]+A[3]*DLP1+A[4]*VBIMF1+A[5]*DST;
            double TAMP2=A[6]+A[7]*DLP2+A[8]*VBIMF1+A[9]*DST;
            double A_SRC=A[10]+A[11]*DST+A[12]*sqrt(PDYN);
            double A_PRC=A[13]+A[14]*DST+A[15]*sqrt(PDYN);
            double A_R11=A[16]+A[17]*VBIMF2;
            double A_R12=A[18]+A[19]*VBIMF2;
            double A_R21=A[20]+A[21]*VBIMF2;
            double A_R22=A[22]+A[23]*VBIMF2;
            
            double BBX=A[1]**BXCF+TAMP1**BXT1+TAMP2**BXT2+A_SRC**BXSRC+A_PRC**BXPRC
            +A_R11**BXR11+A_R12**BXR12+A_R21**BXR21+A_R22**BXR22
            +A[24]**HXIMF+A[25]**HXIMF*STHETAH;
            
            double BBY=A[1]**BYCF+TAMP1**BYT1+TAMP2**BYT2+A_SRC**BYSRC+A_PRC**BYPRC
            +A_R11**BYR11+A_R12**BYR12+A_R21**BYR21+A_R22**BYR22
            +A[24]**HYIMF+A[25]**HYIMF*STHETAH;
            
            double BBZ=A[1]**BZCF+TAMP1**BZT1+TAMP2**BZT2+A_SRC**BZSRC+A_PRC**BZPRC
            +A_R11**BZR11+A_R12**BZR12+A_R21**BZR21+A_R22**BZR22
            +A[24]**HZIMF+A[25]**HZIMF*STHETAH;
            
            
            /*        print *, '   from inside T01: Bz (T+RC)=',
             c     +  TAMP1*BZT1+TAMP2*BZT2+A_SRC*BZSRC
             
             c        print *, '   from inside T01: Bz (CF)=',BZCF*/
            
            /*
             C   AND WE HAVE THE TOTAL EXTERNAL FIELD.
             C
             C  NOW, LET US CHECK WHETHER WE HAVE THE CASE (1). IF YES - WE ARE DONE:
             */
            if (SIGMA<S0-DSIG) {    //  (X,Y,Z) IS INSIDE THE MAGNETOSPHERE
                //-------------------------------------------------------------------------
                *BX=BBX;
                *BY=BBY;
                *BZ=BBZ;
                //-------------------------------------------------------------------------
            } else {//ELSE           !  THIS IS THE MOST COMPLEX CASE: WE ARE INSIDE
                //                                             THE INTERPOLATION REGION
                double FINT=0.5*(1.-(SIGMA-S0)/DSIG);
                double FEXT=0.5*(1.+(SIGMA-S0)/DSIG);
                double  QX,QY,QZ;
                DIPOLE (PS,X,Y,Z,&QX,&QY,&QZ);
                *BX=(BBX+QX)*FINT+OIMFX*FEXT -QX;
                *BY=(BBY+QY)*FINT+OIMFY*FEXT -QY;
                *BZ=(BBZ+QZ)*FINT+OIMFZ*FEXT -QZ;
                
            }//ENDIF  !   THE CASES (1) AND (2) ARE EXHAUSTED; THE ONLY REMAINING
            //                      POSSIBILITY IS NOW THE CASE (3):
            //--------------------------------------------------------------------------
        } else {//ELSE
            double  QX,QY,QZ;
            DIPOLE (PS,X,Y,Z,&QX,&QY,&QZ);
            *BX=OIMFX-QX;
            *BY=OIMFY-QY;
            *BZ=OIMFZ-QZ;
        }//ENDIF
        
    }//END
    
    /*$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$
     C
     c
     SUBROUTINE T01_01 (IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
     c
     c     RELEASE DATE OF THIS VERSION:   AUGUST 8, 2001.
     c
     c    LATEST MODIFICATIONS/BUGS REMOVED:  JUNE 24, 2006:  REPLACED COEFFICIENTS IN:
     c        (i) DATA statement in FUNCTION AP,
     C        (ii)  DATA C_SY statement in SUBROUTINE FULL_RC, and
     c        (iii) DATA A statement in SUBROUTINE T01_01.
     C    This correction was needed because of a bug found in the symmetric ring current module.
     c    Its impact is a minor (a few percent) change of the model field in the inner magnetosphere.
     C
     c--------------------------------------------------------------------
     C   A DATA-BASED MODEL OF THE EXTERNAL (I.E., WITHOUT EARTH'S CONTRIBUTION) PART OF THE
     C   MAGNETOSPHERIC MAGNETIC FIELD, CALIBRATED BY
     C    (1) SOLAR WIND PRESSURE PDYN (NANOPASCALS),
     C    (2) DST (NANOTESLA),
     C    (3) BYIMF,
     C    (4) BZIMF (NANOTESLA)
     C    (5) G1-INDEX
     C    (6) G2-INDEX  (SEE TSYGANENKO [2001] FOR AN EXACT DEFINITION OF THESE TWO INDICES)
     
     c   THESE INPUT PARAMETERS SHOULD BE PLACED IN THE FIRST 6 ELEMENTS
     c   OF THE ARRAY PARMOD(10).
     C
     C   THE REST OF THE INPUT VARIABLES ARE: THE GEODIPOLE TILT ANGLE PS (RADIANS),
     C     AND   X,Y,Z -  GSM POSITION (RE)
     C
     c   IOPT  IS JUST A DUMMY INPUT PARAMETER, NECESSARY TO MAKE THIS SUBROUTINE
     C   COMPATIBLE WITH THE TRACING SOFTWARE PACKAGE (GEOPACK). IN THIS MODEL
     C   IT DOES NOT AFFECT THE OUTPUT FIELD.
     c
     C*******************************************************************************************
     c** ATTENTION:  THE MODEL IS BASED ON DATA TAKEN SUNWARD FROM X=-15Re, AND HENCE BECOMES   *
     C**              INVALID AT LARGER TAILWARD DISTANCES !!!                                  *
     C*******************************************************************************************
     C
     c   OUTPUT:  GSM COMPONENTS OF THE EXTERNAL MAGNETIC FIELD (BX,BY,BZ, nanotesla)
     C            COMPUTED AS A SUM OF CONTRIBUTIONS FROM PRINCIPAL FIELD SOURCES
     C
     c  (C) Copr. 2001, Nikolai A. Tsyganenko, USRA, Code 690.2, NASA GSFC
     c      Greenbelt, MD 20771, USA
     c
     C                            REFERENCE:
     C
     C    N. A. Tsyganenko, A new data-based model of the near magnetosphere magnetic field:
     c       1. Mathematical structure.
     c       2. Parameterization and fitting to observations.
     c
     c             (submitted to JGR, July 2001)
     C
     C
     c----------------------------------------------------------------------
     */
    static void T01_01 (int const IOPT,double const PARMOD[],double const PS,double const X,double const Y,double const Z,double *BX,double *BY,double *BZ) {
        //REAL PARMOD(10),PS,X,Y,Z,BX,BY,BZ
        //REAL*8 A(43),PDYN,DST_AST,BYIMF,BZIMF,G1,G2,PSS,XX,YY,ZZ,
        //*  BXCF,BYCF,BZCF,BXT1,BYT1,BZT1,BXT2,BYT2,BZT2,
        //*  BXSRC,BYSRC,BZSRC,BXPRC,BYPRC,BZPRC, BXR11,BYR11,BZR11,
        //*  BXR12,BYR12,BZR12,BXR21,BYR21,BZR21,BXR22,BYR22,BZR22,HXIMF,
        //*  HYIMF,HZIMF,BBX,BBY,BBZ
        
        static double const A[44] = {0./*dummy*/,1.00000,2.47341,0.40791,0.30429,-0.10637,-0.89108,3.29350,
            -0.05413,-0.00696,1.07869,-0.02314,-0.66173,-0.68018,-0.03246,
            0.02681,0.28062,0.16535,-0.02939,0.02639,-0.24891,-0.08063,
            0.08900,-0.02475,0.05887,0.57691,0.65256,-0.03230,2.24733,
            4.10546,1.13665,0.05506,0.97669,0.21164,0.64594,1.12556,0.01389,
            1.02978,0.02968,0.15821,9.00519,28.17582,1.35285,0.42279};
        /*
         c   the disclaimer below is temporarily disabled:
         C
         c      IF (X.LT.-20.) THEN
         c      PRINT *,
         c     * '  ATTENTION:  THE MODEL IS VALID SUNWARD FROM X=-15 Re ONLY,'
         c      PRINT *,'              WHILE YOU ARE TRYING TO USE IT AT X=', X
         c      PAUSE
         c      ENDIF
         */
        double PDYN=PARMOD[0];
        double DST_AST=PARMOD[1]*0.8-13.*sqrt(PDYN);
        double BYIMF=PARMOD[2];
        double BZIMF=PARMOD[3];
        
        double G1=PARMOD[4];
        double G2=PARMOD[5];
        double PSS=PS;
        double XX=X;
        double YY=Y;
        double ZZ=Z;
        double BXCF,BYCF,BZCF,BXT1,BYT1,BZT1,BXT2,BYT2,BZT2, BXSRC,BYSRC,BZSRC,BXPRC,BYPRC,BZPRC, BXR11,BYR11,BZR11, BXR12,BYR12,BZR12,BXR21,BYR21,BZR21,BXR22,BYR22,BZR22,HXIMF, HYIMF,HZIMF,BBX,BBY,BBZ;
        EXTALL (0,0,0,0,A,43,PDYN,DST_AST,BYIMF,BZIMF,G1,G2,
                PSS,XX,YY,ZZ,&BXCF,&BYCF,&BZCF,&BXT1,&BYT1,&BZT1,&BXT2,&BYT2,&BZT2,
                &BXSRC,&BYSRC,&BZSRC,&BXPRC,&BYPRC,&BZPRC, &BXR11,&BYR11,&BZR11,
                &BXR12,&BYR12,&BZR12,&BXR21,&BYR21,&BZR21,&BXR22,&BYR22,&BZR22,&HXIMF,
                &HYIMF,&HZIMF,&BBX,&BBY,&BBZ);
        
        *BX=BBX;
        *BY=BBY;
        *BZ=BBZ;
        
        //RETURN
    }//END

}