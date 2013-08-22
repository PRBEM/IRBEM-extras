//
//  T89.cpp
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

#include "T89.h"
#include "Geopack.h"
#include <cassert>
#include <cmath>
#include <stdexcept>

namespace UBK {
    using namespace std;

    T89::T89 (Geopack const* geopack, int iopt) : TSExternalField(geopack, iopt, NULL)
    {
        if (iopt<1 || iopt>7) {
            throw out_of_range("IOPT for T89 model should be in the range 1<=IOPT<=7.");
        }
    }

    static void T89C (int const IOPT, double const PARMOD[], double const PS, double const X, double const Y, double const Z, double *BX, double *BY, double *BZ);
    void T89::getFieldInGSW_atPoint (Point *bOut, Point const ptgsw) const
    {
        T89C(this->iopt(), this->parmod(), this->geopack()->psi(), ptgsw.x, ptgsw.y, ptgsw.z, &bOut->x, &bOut->y, &bOut->z);
    }

    /*-----------------------------Model---------------------*/
    /*
     c  The small main program below is an example of how to compute field
     c   components with T89C.
     c    See GEOPACK.DOC for an example of field line tracing.
     c
     c     dimension parmod(10)
     c 1    print *, '  enter x,y,z,ps,iopt'
     c     read*, x,y,z,ps,iopt
     c     call t89c(iopt,parmod,ps,x,y,z,bx,by,bz)
     c     print *, bx,by,bz
     c     goto 1
     c     end

     C
     SUBROUTINE T89C(IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
     C
     C
     C   COMPUTES GSM COMPONENTS OF THE MAGNETIC FIELD PRODUCED BY EXTRA-
     C   TERRESTRIAL CURRENT SYSTEMS IN THE GEOMAGNETOSPHERE. THE MODEL IS
     C   VALID UP TO GEOCENTRIC DISTANCES OF 70 RE AND IS BASED ON THE MER-
     C   GED IMP-A,C,D,E,F,G,H,I,J (1966-1974), HEOS-1 AND -2 (1969-1974),
     C   AND ISEE-1 AND -2  SPACECRAFT DATA SET.
     C
     C   THIS IS A MODIFIED VERSION (T89c), WHICH REPLACED THE ORIGINAL ONE
     C     IN 1992 AND DIFFERS FROM IT IN THE FOLLOWING:
     C
     C   (1)  ISEE-1,2 DATA WERE ADDED TO THE ORIGINAL IMP-HEOS DATASET
     C   (2)  TWO TERMS WERE ADDED TO THE ORIGINAL TAIL FIELD MODES, ALLOWING
     C          A MODULATION OF THE CURRENT BY THE GEODIPOLE TILT ANGLE
     C
     C
     C  REFERENCE FOR THE ORIGINAL MODEL: N.A. TSYGANENKO, A MAGNETOSPHERIC MAGNETIC
     C       FIELD MODEL WITH A WARPED TAIL CURRENT SHEET: PLANET.SPACE SCI., V.37,
     C         PP.5-20, 1989.
     C
     C----INPUT PARAMETERS: IOPT - SPECIFIES THE GROUND DISTURBANCE LEVEL:
     C
     C   IOPT= 1       2        3        4        5        6      7
     C                  CORRESPOND TO:
     C    KP= 0,0+  1-,1,1+  2-,2,2+  3-,3,3+  4-,4,4+  5-,5,5+  > =6-
     C
     C    PS - GEODIPOLE TILT ANGLE IN RADIANS
     C    X, Y, Z  - GSM COORDINATES OF THE POINT IN EARTH RADII
     C
     C----OUTPUT PARAMETERS: BX,BY,BZ - GSM COMPONENTS OF THE MODEL MAGNETIC
     C                        FIELD IN NANOTESLAS
     c
     c   THE PARAMETER PARMOD(10) IS A DUMMY ARRAY.  IT IS NOT USED IN THIS
     C        SUBROUTINE AND IS PROVIDED JUST FOR MAKING IT COMPATIBLE WITH THE
     C           NEW VERSION (4/16/96) OF THE GEOPACK SOFTWARE.
     C
     C   THIS RELEASE OF T89C IS DATED  FEB 12, 1996;
     C--------------------------------------------------------------------------
     C
     C
     C              AUTHOR:     NIKOLAI A. TSYGANENKO
     C                          HSTX CORP./NASA GSFC
     C
     */
    /*
     C-------------------------------------------------------------------
     C
     SUBROUTINE  T89 (ID, A, XI, F, DER)
     C
     C        ***  N.A. Tsyganenko ***  8-10.12.1991  ***
     C
     C      Calculates dependent model variables and their deriva-
     C  tives for given independent variables and model parame-
     C  ters.  Specifies model functions with free parameters which
     C  must be determined by means of least squares fits (RMS
     C  minimization procedure).
     C
     C      Description of parameters:
     C
     C  ID  - number of the data point in a set (initial assignments are performed
     c        only for ID=1, saving thus CPU time)
     C  A   - input vector containing model parameters;
     C  XI  - input vector containing independent variables;
     C  F   - output double precision vector containing
     C        calculated values of dependent variables;
     C  DER   - output double precision vector containing
     C        calculated values for derivatives of dependent
     C        variables with respect to model parameters;
     C
     C - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C
     C      T89 represents external magnetospheric magnetic field
     C  in Cartesian SOLAR MAGNETOSPHERIC coordinates (Tsyganenko N.A.,
     C  Planet. Space Sci., 1989, v.37, p.5-20; the "T89 model" with the warped
     c  tail current sheet) + A MODIFICATION ADDED IN APRIL 1992 (SEE BELOW)
     C
     C      Model formulas for the magnetic field components contain in total
     c  30 free parameters (17 linear and 13 nonlinear parameters).
     C      First 2 independent linear parameters A(1)-A(2) correspond to contribu-
     c  tion from the tail current system, then follow A(3) and A(4) which are the
     c  amplitudes of symmetric and antisymmetric terms in the contribution from
     c  the closure currents; A(5) is the ring current amplitude. Then follow the
     c coefficients A(6)-A(15) which define Chapman-Ferraro+Birkeland current field.
     c    The coefficients c16-c19  (see Formula 20 in the original paper),
     c   due to DivB=0 condition, are expressed through A(6)-A(15) and hence are not
     c    independent ones.
     c  A(16) AND A(17) CORRESPOND TO THE TERMS WHICH YIELD THE TILT ANGLE DEPEN-
     C    DENCE OF THE TAIL CURRENT INTENSITY (ADDED ON APRIL 9, 1992)
     C
     C      Nonlinear parameters:
     C
     C    A(18) : DX - Characteristic scale of the Chapman-Ferraro field along the
     c        X-axis
     C    A(19) : ADR (aRC) - Characteristic radius of the ring current
     c    A(20) : D0 - Basic half-thickness of the tail current sheet
     C    A(21) : DD (GamRC)- defines rate of thickening of the ring current, as
     c             we go from night- to dayside
     C    A(22) : Rc - an analog of "hinging distance" entering formula (11)
     C    A(23) : G - amplitude of tail current warping in the Y-direction
     C    A(24) : aT - Characteristic radius of the tail current
     c    A(25) : Dy - characteristic scale distance in the Y direction entering
     c                 in W(x,y) in (13)
     c    A(26) : Delta - defines the rate of thickening of the tail current sheet
     c                 in the Y-direction (in T89 it was fixed at 0.01)
     c    A(27) : Q - this parameter was fixed at 0 in the final version of T89;
     c              initially it was introduced for making Dy to depend on X
     c    A(28) : Sx (Xo) - enters in W(x,y) ; see (13)
     c    A(29) : Gam (GamT) - enters in DT in (13) and defines rate of tail sheet
     c              thickening on going from night to dayside; in T89 fixed at 4.0
     c    A(30) : Dyc - the Dy parameter for closure current system; in T89 fixed
     c               at 20.0
     c  - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
     C
     */
    static void T89 (int const ID,double const A[],double const XI[],double F[]/*,double *DER*/)
    {
        static double const A02 = 25.;
        static double const XLW2 = 170.;
        //double const YN = 30.;
        //double const RPI = 0.31830989;
        static double const RT = 30.;
        //DATA A02,XLW2,YN,RPI,RT/25.D0,170.D0,30.D0,0.31830989D0,30.D0/
        static double const XD = 0.;
        static double const XLD2 = 40.;
        //DATA XD,XLD2/0.D0,40.D0/
        /*
         C   The last four quantities define variation of tail sheet thickness along X
         */
        static double const SXC = 4.;
        static double const XLWC2 = 50.;
        //DATA SXC,XLWC2/4.D0,50.D0/
        /*
         C   The two quantities belong to the function WC which confines tail closure
         c    current in X- and Y- direction
         */
        //double const DXL = 20.;
        //DATA DXL/20.D0/

        //IF (ID.NE.1)  GOTO  3
        //DO  2  I = 1, 30
        //DO  1  L = 1, 3
        //1            DER(L,I) = 0.0D0
        //2        CONTINUE
        double DER[30+1][3+1] = {0,}; // DER[0][0] are dummy. Fortran index start from 1.

        double DYC=A[29];
        double DYC2=DYC*DYC;
        double DX=A[17];
        double HA02=0.5*A02;
        //double RDX2M=-1./(DX*DX);
        //double RDX2=-RDX2M;
        double RDYC2=1./DYC2;
        double HLWC2M=-0.5*XLWC2;
        double DRDYC2=-2.*RDYC2;
        //double DRDYC3=2.*RDYC2*sqrt(RDYC2);
        double HXLW2M=-0.5*XLW2;
        double ADR=A[18];
        double D0=A[19];
        double DD=A[20];
        double RC=A[21];
        double G=A[22];
        double AT=A[23];
        double DT=D0;
        double DEL=A[25];
        double P=A[24];
        double Q=A[26];
        double SX=A[27];
        double GAM=A[28];
        double HXLD2M=-0.5*XLD2;
        double ADSL=0.;
        double XGHS=0.;
        double H=0.;
        double HS=0.;
        double GAMH=0.;
        double W1=-0.5/DX;
        double DBLDEL=2.*DEL;
        double W2=W1*2.;
        double W4=-1./3.;
        double W3=W4/DX;
        double W5=-0.5;
        double W6=-3.;
        double AK1=A[0];
        double AK2=A[1];
        double AK3=A[2];
        double AK4=A[3];
        double AK5=A[4];
        double AK6=A[5];
        double AK7=A[6];
        double AK8=A[7];
        double AK9=A[8];
        double AK10=A[9];
        double AK11=A[10];
        double AK12=A[11];
        double AK13=A[12];
        double AK14=A[13];
        double AK15=A[14];
        //double AK16=A[15];
        //double AK17=A[16];
        double SXA=0.;
        double SYA=0.;
        double SZA=0.;
        double AK610=AK6*W1+AK10*W5;
        double AK711=AK7*W2-AK11;
        double AK812=AK8*W2+AK12*W6;
        double AK913=AK9*W3+AK13*W4;
        //double RDXL=1./DXL;
        //double HRDXL=0.5*RDXL;
        //double A6H=AK6*0.5;
        //double A9T=AK9/3.;
        //double YNP=RPI/YN*0.5;
        //double YND=2.*YN;

        //3      CONTINUE

        double X  = XI[0];
        double Y  = XI[1];
        double Z  = XI[2];
        double TILT=XI[3];
        double TLT2=TILT*TILT;
        double SPS = sin(TILT);
        double CPS = sqrt(1.0 - SPS*SPS);

        double X2=X*X;
        double Y2=Y*Y;
        double Z2=Z*Z;
        double TPS=SPS/CPS;
        double HTP=TPS*0.5;
        //double GSP=G*SPS;
        double XSM=X*CPS-Z*SPS;
        double ZSM=X*SPS+Z*CPS;
        /*
         C   CALCULATE THE FUNCTION ZS DEFINING THE SHAPE OF THE TAIL CURRENT SHEET
         C    AND ITS SPATIAL DERIVATIVES:
         */
        double XRC=XSM+RC;
        double XRC16=XRC*XRC+16.;
        double SXRC=sqrt(XRC16);
        double Y4=Y2*Y2;
        double Y410=Y4+1.e4;
        double SY4=SPS/Y410;
        double GSY4=G*SY4;
        double ZS1=HTP*(XRC-SXRC);
        double DZSX=-ZS1/SXRC;
        double ZS=ZS1-GSY4*Y4;
        double D2ZSGY=-SY4/Y410*4.e4*Y2*Y;
        double DZSY=G*D2ZSGY;
        /*
         C   CALCULATE THE COMPONENTS OF THE RING CURRENT CONTRIBUTION:
         */
        double XSM2=XSM*XSM;
        double DSQT=sqrt(XSM2+A02);
        double FA0=0.5e0*(1.e0+XSM/DSQT);
        double DDR=D0+DD*FA0;
        double DFA0=HA02/(DSQT*DSQT*DSQT);
        double ZR=ZSM-ZS;
        double TR=sqrt(ZR*ZR + DDR*DDR);
        double RTR=1.e0/TR;
        double RO2=XSM2+Y2;
        double ADRT=ADR+TR;
        double ADRT2=ADRT*ADRT;
        double FK=1.e0/(ADRT2+RO2);
        double DSFC=sqrt(FK);
        double FC=FK*FK*DSFC;
        double FACXY=3.0e0*ADRT*FC*RTR;
        double XZR=XSM*ZR;
        double YZR=Y*ZR;
        double DBXDP=FACXY*XZR;
        DER[5][2]=FACXY*YZR;
        double XZYZ=XSM*DZSX+Y*DZSY;
        double FAQ=ZR*XZYZ-DDR*DD*DFA0*XSM;
        double DBZDP=FC*(2.e0*ADRT2-RO2)+FACXY*FAQ;
        DER[5][1]=DBXDP*CPS+DBZDP*SPS;
        DER[5][3]=DBZDP*CPS-DBXDP*SPS;
        /*
         C  CALCULATE THE TAIL CURRENT SHEET CONTRIBUTION:
         */
        double DELY2=DEL*Y2;
        double D=DT+DELY2;
        if (fabs(GAM) >= 1e-6) {
            //IF (DABS(GAM).LT.1.D-6) GOTO 8
            double XXD=XSM-XD;
            double RQD=1.e0/(XXD*XXD+XLD2);
            double RQDS=sqrt(RQD);
            H=0.5e0*(1.e0+XXD*RQDS);
            HS=-HXLD2M*RQD*RQDS;
            GAMH=GAM*H;
            D=D+GAMH;
            XGHS=XSM*GAM*HS;
            ADSL=-D*XGHS;
        }
        //8   D2=D**2
        double D2=D*D;
        double T=sqrt(ZR*ZR+D2);
        double XSMX=XSM-SX;
        double RDSQ2=1.e0/(XSMX*XSMX+XLW2);
        double RDSQ=sqrt(RDSQ2);
        double V=0.5e0*(1.e0-XSMX*RDSQ);
        double DVX=HXLW2M*RDSQ*RDSQ2;
        double OM=sqrt(sqrt(XSM2+16.e0)-XSM);
        double OMS=-OM/(OM*OM+XSM)*0.5e0;
        double RDY=1.e0/(P+Q*OM);
        double OMSV=OMS*V;
        double RDY2=RDY*RDY;
        double FY=1.e0/(1.e0+Y2*RDY2);
        double W=V*FY;
        double YFY1=2.e0*FY*Y2*RDY2;
        double FYPR=YFY1*RDY;
        double FYDY=FYPR*FY;
        double DWX=DVX*FY+FYDY*Q*OMSV;
        double YDWY=-V*YFY1*FY;
        double DDY=DBLDEL*Y;
        double ATT=AT+T;
        double S1=sqrt(ATT*ATT+RO2);
        double F5=1.e0/S1;
        double F7=1.e0/(S1+ATT);
        double F1=F5*F7;
        double F3=F5*F5*F5;
        double F9=ATT*F3;
        double FS=ZR*XZYZ-D*Y*DDY+ADSL;
        double XDWX=XSM*DWX+YDWY;
        double RTT=1.e0/T;
        double WT=W*RTT;
        double BRRZ1=WT*F1;
        double BRRZ2=WT*F3;
        double DBXC1=BRRZ1*XZR;
        double DBXC2=BRRZ2*XZR;
        DER[1][2]=BRRZ1*YZR;
        DER[2][2]=BRRZ2*YZR;
        DER[16][2]=DER[1][2]*TLT2;
        DER[17][2]=DER[2][2]*TLT2;
        double WTFS=WT*FS;
        double DBZC1=W*F5+XDWX*F7+WTFS*F1;
        double DBZC2=W*F9+XDWX*F1+WTFS*F3;
        DER[1][1]=DBXC1*CPS+DBZC1*SPS;
        DER[2][1]=DBXC2*CPS+DBZC2*SPS;
        DER[1][3]=DBZC1*CPS-DBXC1*SPS;
        DER[2][3]=DBZC2*CPS-DBXC2*SPS;
        DER[16][1]=DER[1][1]*TLT2;
        DER[17][1]=DER[2][1]*TLT2;
        DER[16][3]=DER[1][3]*TLT2;
        DER[17][3]=DER[2][3]*TLT2;
        /*
         C  CALCULATE CONTRIBUTION FROM THE CLOSURE CURRENTS
         */
        double ZPL=Z+RT;
        double ZMN=Z-RT;
        double ROGSM2=X2+Y2;
        double SPL=sqrt(ZPL*ZPL+ROGSM2);
        double SMN=sqrt(ZMN*ZMN+ROGSM2);
        double XSXC=X-SXC;
        double RQC2=1.e0/(XSXC*XSXC+XLWC2);
        double RQC=sqrt(RQC2);
        double FYC=1.e0/(1.e0+Y2*RDYC2);
        double WC=0.5e0*(1.e0-XSXC*RQC)*FYC;
        double DWCX=HLWC2M*RQC2*RQC*FYC;
        double DWCY=DRDYC2*WC*FYC*Y;
        double SZRP=1.e0/(SPL+ZPL);
        double SZRM=1.e0/(SMN-ZMN);
        double XYWC=X*DWCX+Y*DWCY;
        double WCSP=WC/SPL;
        double WCSM=WC/SMN;
        double FXYP=WCSP*SZRP;
        double FXYM=WCSM*SZRM;
        double FXPL=X*FXYP;
        double FXMN=-X*FXYM;
        double FYPL=Y*FXYP;
        double FYMN=-Y*FXYM;
        double FZPL=WCSP+XYWC*SZRP;
        double FZMN=WCSM+XYWC*SZRM;
        DER[3][1]=FXPL+FXMN;
        DER[4][1]=(FXPL-FXMN)*SPS;
        DER[3][2]=FYPL+FYMN;
        DER[4][2]=(FYPL-FYMN)*SPS;
        DER[3][3]=FZPL+FZMN;
        DER[4][3]=(FZPL-FZMN)*SPS;
        /*
         C   NOW CALCULATE CONTRIBUTION FROM CHAPMAN-FERRARO SOURCES + ALL OTHER
         */
        double EX=exp(X/DX);
        double EC=EX*CPS;
        double ES=EX*SPS;
        double ECZ=EC*Z;
        double ESZ=ES*Z;
        double ESZY2=ESZ*Y2;
        double ESZZ2=ESZ*Z2;
        double ECZ2=ECZ*Z;
        double ESY=ES*Y;

        DER[6][1]=ECZ;
        DER[7][1]=ES;
        DER[8][1]=ESY*Y;
        DER[9][1]=ESZ*Z;
        DER[10][2]=ECZ*Y;
        DER[11][2]=ESY;
        DER[12][2]=ESY*Y2;
        DER[13][2]=ESY*Z2;
        DER[14][3]=EC;
        DER[15][3]=EC*Y2;
        DER[6][3]=ECZ2*W1;
        DER[10][3]=ECZ2*W5;
        DER[7][3]=ESZ*W2;
        DER[11][3]=-ESZ;
        DER[8][3]=ESZY2*W2;
        DER[12][3]=ESZY2*W6;
        DER[9][3]=ESZZ2*W3;
        DER[13][3]=ESZZ2*W4;
        /*
         C  FINALLY, CALCULATE NET EXTERNAL MAGNETIC FIELD COMPONENTS,
         C    BUT FIRST OF ALL THOSE FOR C.-F. FIELD:
         */
        double SX1=AK6*DER[6][1]+AK7*DER[7][1]+AK8*DER[8][1]+AK9*DER[9][1];
        double SY1=AK10*DER[10][2]+AK11*DER[11][2]+AK12*DER[12][2]+AK13*DER[13][2];
        double SZ1=AK14*DER[14][3]+AK15*DER[15][3]+AK610*ECZ2+AK711*ESZ+AK812*ESZY2+AK913*ESZZ2;
        //c     SX1=0D0;SZ1=0D0;SY1=0D0;
        double BXCL=AK3*DER[3][1]+AK4*DER[4][1];
        double BYCL=AK3*DER[3][2]+AK4*DER[4][2];
        double BZCL=AK3*DER[3][3]+AK4*DER[4][3];
        //c      BXCL=0D0;BYCL=0D0;BZCL=0D0
        //c      BXT=AK1*DER(1,1)+AK2*DER(1,2)+BXCL +AK16*DER(1,16)+AK17*DER(1,17)
        //c      BYT=AK1*DER(2,1)+AK2*DER(2,2)+BYCL +AK16*DER(2,16)+AK17*DER(2,17)
        //c      BZT=AK1*DER(3,1)+AK2*DER(3,2)+BZCL +AK16*DER(3,16)+AK17*DER(3,17)
        double BXT=AK1*DER[1][1]+AK2*DER[2][1]+BXCL;
        double BYT=AK1*DER[1][2]+AK2*DER[2][2]+BYCL;
        double BZT=AK1*DER[1][3]+AK2*DER[2][3]+BZCL;
        F[0]=BXT+AK5*DER[5][1]+SX1+SXA;
        F[1]=BYT+AK5*DER[5][2]+SY1+SYA;
        F[2]=BZT+AK5*DER[5][3]+SZ1+SZA;
    }
    /*
     C-------------------------------------------------------------------
     c  The small main program below is an example of how to compute field
     c   components with T89C.
     c    See GEOPACK.DOC for an example of field line tracing.
     c
     c     dimension parmod(10)
     c 1    print *, '  enter x,y,z,ps,iopt'
     c     read*, x,y,z,ps,iopt
     c     call t89c(iopt,parmod,ps,x,y,z,bx,by,bz)
     c     print *, bx,by,bz
     c     goto 1
     c     end

     C
     SUBROUTINE T89C(IOPT,PARMOD,PS,X,Y,Z,BX,BY,BZ)
     C
     C
     C   COMPUTES GSM COMPONENTS OF THE MAGNETIC FIELD PRODUCED BY EXTRA-
     C   TERRESTRIAL CURRENT SYSTEMS IN THE GEOMAGNETOSPHERE. THE MODEL IS
     C   VALID UP TO GEOCENTRIC DISTANCES OF 70 RE AND IS BASED ON THE MER-
     C   GED IMP-A,C,D,E,F,G,H,I,J (1966-1974), HEOS-1 AND -2 (1969-1974),
     C   AND ISEE-1 AND -2  SPACECRAFT DATA SET.
     C
     C   THIS IS A MODIFIED VERSION (T89c), WHICH REPLACED THE ORIGINAL ONE
     C     IN 1992 AND DIFFERS FROM IT IN THE FOLLOWING:
     C
     C   (1)  ISEE-1,2 DATA WERE ADDED TO THE ORIGINAL IMP-HEOS DATASET
     C   (2)  TWO TERMS WERE ADDED TO THE ORIGINAL TAIL FIELD MODES, ALLOWING
     C          A MODULATION OF THE CURRENT BY THE GEODIPOLE TILT ANGLE
     C
     C
     C  REFERENCE FOR THE ORIGINAL MODEL: N.A. TSYGANENKO, A MAGNETOSPHERIC MAGNETIC
     C       FIELD MODEL WITH A WARPED TAIL CURRENT SHEET: PLANET.SPACE SCI., V.37,
     C         PP.5-20, 1989.
     C
     C----INPUT PARAMETERS: IOPT - SPECIFIES THE GROUND DISTURBANCE LEVEL:
     C
     C   IOPT= 1       2        3        4        5        6      7
     C                  CORRESPOND TO:
     C    KP= 0,0+  1-,1,1+  2-,2,2+  3-,3,3+  4-,4,4+  5-,5,5+  > =6-
     C
     C    PS - GEODIPOLE TILT ANGLE IN RADIANS
     C    X, Y, Z  - GSM COORDINATES OF THE POINT IN EARTH RADII
     C
     C----OUTPUT PARAMETERS: BX,BY,BZ - GSM COMPONENTS OF THE MODEL MAGNETIC
     C                        FIELD IN NANOTESLAS
     c
     c   THE PARAMETER PARMOD(10) IS A DUMMY ARRAY.  IT IS NOT USED IN THIS
     C        SUBROUTINE AND IS PROVIDED JUST FOR MAKING IT COMPATIBLE WITH THE
     C           NEW VERSION (4/16/96) OF THE GEOPACK SOFTWARE.
     C
     C   THIS RELEASE OF T89C IS DATED  FEB 12, 1996;
     C--------------------------------------------------------------------------
     C
     C
     C              AUTHOR:     NIKOLAI A. TSYGANENKO
     C                          HSTX CORP./NASA GSFC
     C
     */
    static void T89C (int const IOPT, double const PARMOD[], double const PS, double const X, double const Y, double const Z, double *BX, double *BY, double *BZ)
    {
        double XI[4],F[3],/*DER[3*30],*/A[30]/*,PARMOD[10]*/;
        static double const PARAM[7][30] = {
            -116.53,-10719.,42.375,59.753,-11363.,1.7844,30.268,
            -0.35372E-01,-0.66832E-01,0.16456E-01,-1.3024,0.16529E-02,
            0.20293E-02,20.289,0.25203E-01,224.91,-9234.8,22.788,7.8813,
            1.8362,-0.27228,8.8184,2.8714,14.468,32.177,0.01,0.0,
            7.0459,4.0,20.0,-55.553,-13198.,60.647,61.072,-16064.,
            2.2534,34.407,-0.38887E-01,-0.94571E-01,0.27154E-01,-1.3901,
            0.13460E-02,0.13238E-02,23.005,-0.30565E-01,55.047,-3875.7,
            20.178,7.9693,1.4575,0.89471,9.4039,3.5215,14.474,36.555,
            0.01,0.0,7.0787,4.0,20.0,-101.34,-13480.,111.35,12.386,-24699.,
            2.6459,38.948,-0.34080E-01,-0.12404,0.29702E-01,-1.4052,
            0.12103E-02,0.16381E-02,24.49,-0.37705E-01,-298.32,4400.9,18.692,
            7.9064,1.3047,2.4541,9.7012,7.1624,14.288,33.822,0.01,0.0,6.7442,
            4.0,20.0,-181.69,-12320.,173.79,-96.664,-39051.,3.2633,44.968,
            -0.46377E-01,-0.16686,0.048298,-1.5473,0.10277E-02,0.31632E-02,
            27.341,-0.50655E-01,-514.10,12482.,16.257,8.5834,1.0194,3.6148,
            8.6042,5.5057,13.778,32.373,0.01,0.0,7.3195,4.0,20.0,-436.54,
            -9001.0,323.66,-410.08,-50340.,3.9932,58.524,-0.38519E-01,
            -0.26822,0.74528E-01,-1.4268,-0.10985E-02,0.96613E-02,27.557,
            -0.56522E-01,-867.03,20652.,14.101,8.3501,0.72996,3.8149,9.2908,
            6.4674,13.729,28.353,0.01,0.0,7.4237,4.0,20.0,-707.77,-4471.9,
            432.81,-435.51,-60400.,4.6229,68.178,-0.88245E-01,-0.21002,
            0.11846,-2.6711,0.22305E-02,0.10910E-01,27.547,-0.54080E-01,
            -424.23,1100.2,13.954,7.5337,0.89714,3.7813,8.2945,5.174,14.213,
            25.237,0.01,0.0,7.0037,4.0,20.0,-1190.4,2749.9,742.56,-1110.3,
            -77193.,7.6727,102.05,-0.96015E-01,-0.74507,0.11214,-1.3614,
            0.15157E-02,0.22283E-01,23.164,-0.74146E-01,-2219.1,48253.,
            12.714,7.6777,0.57138,2.9633,9.3909,9.7263,11.123,21.558,0.01,
            0.0,4.4518,4.0,20.0};
        
        int ID=1;
        //int IOP=IOPT;
        for (int I=0; I<30; I++) {
            A[I]=PARAM[IOPT-1][I];
        }
        
        XI[0]=X;
        XI[1]=Y;
        XI[2]=Z;
        XI[3]=PS;
        T89(ID,A,XI,F/*,DER*/);
        *BX=F[0];
        *BY=F[1];
        *BZ=F[2];
    }

}