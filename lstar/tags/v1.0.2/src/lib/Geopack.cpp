//
//  Geopack.cpp
//  UBJDevelopment
//
//  Translated by Kyungguk Min on 4/5/13.
//  Original code by N. Tsyganenko.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#include "Geopack.h"
#include <cmath>
#include <cstdio>
#include <cassert>
#include <stdexcept>

namespace UBK {
    using namespace std;

    //
    // Unit test
    //
#ifdef DEBUG
    void Geopack::test()
    {
        Geopack g(2007,1,12,00,00,400,kGeopackDipoleField);
        //g.setPsi(0.);
        printf("PSI = %f[deg]\n", g.psi()*180/M_PI);

        Point pt(6.);
        Point b;
        g.getFieldInGSW_atPoint(&b, pt);

        printf("B(%s) = %s\n", pt.desc().c_str(), b.desc().c_str());

        assert(0);
    }
#endif

    //////////////////////////////////////////////////////////////////////////////
    // GEOPACK common block definitions. Instead of using global storage, the common blocks are stored in user program space so that the code is more thread safe. Instances of GeopackContext should not be modified. One should call RECALC_08 to initialize contents.
    typedef struct _GeopackContext {
        // GEOPACK 1
        double G[106];
        double H[106];
        double REC[106];
        double DIPMOM;
        // GEOPACK 2. PSI, is the calculated Earth's magnetic dipole tilt angle.
        double ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33,E11,E21,E31,E12,E22,E32,E13,E23,E33;
    } GeopackContext;

    static void RECALC_08 (int const IYEAR,int const IDAY,int const IHOUR,int const MIN,int const ISEC,double const VGSEX,double const VGSEY,double const VGSEZ, GeopackContext *ctx);
    static void IGRF_GSW_08 (double XGSW,double YGSW,double ZGSW,double *HXGSW,double *HYGSW,double *HZGSW, GeopackContext const* ctx);
    static void DIP_08 (double XGSW,double YGSW,double ZGSW,double *BXGSW,double *BYGSW,double *BZGSW, GeopackContext const* ctx);
    //////////////////////////////////////////////////////////////////////////////

    //
    // Dipole tilt angle
    //
    void Geopack::setPsi(const double psi)
    {
        static_cast<GeopackContext *>(_ctx)->PSI = psi;
        static_cast<GeopackContext *>(_ctx)->SPS = sin(psi);
        static_cast<GeopackContext *>(_ctx)->CPS = cos(psi);
    }
    double const& Geopack::psi () const
    {
        return static_cast<GeopackContext *>(_ctx)->PSI;
    }
    double const& Geopack::sinpsi () const
    {
        return static_cast<GeopackContext *>(_ctx)->SPS;
    }
    double const& Geopack::cospsi () const
    {
        return static_cast<GeopackContext *>(_ctx)->CPS;
    }

    //
    // Constructor
    //
    Geopack::~Geopack()
    {
        delete static_cast<GeopackContext *>(_ctx);
    }
    Geopack::Geopack (int year, int doy, int hour, int min, int sec, Point vgse, GeopackInternalFieldModel flag) : _year(0), _doy(0), _hour(0), _min(0), _sec(0), _vgse()
    {
        _ctx = new GeopackContext;
        this->recaluateWithDate_vGSE(year, doy, hour, min, sec, vgse);

        switch (flag) {
            case kGeopackDipoleField:
                _helper = (GeopackFieldHelper)&DIP_08;
                break;
                case kGeopackIGRFField:
                _helper = (GeopackFieldHelper)&IGRF_GSW_08;
                break;
            default:
                throw invalid_argument("The flag should be either GeopackDipoleField or GeopackIGRFField.");
                break;
        }
    }
    Geopack::Geopack (int year, int doy, int hour, int min, int sec, double vsw, GeopackInternalFieldModel flag) : _year(0), _doy(0), _hour(0), _min(0), _sec(0), _vgse() {
        assert(vsw >= 0.);

        _ctx = new GeopackContext;
        this->recaluateWithDate_vGSE(year, doy, hour, min, sec, Point(-fabs(vsw), 0., 0.));

        switch (flag) {
            case kGeopackDipoleField:
                _helper = (GeopackFieldHelper)&DIP_08;
                break;
            case kGeopackIGRFField:
                _helper = (GeopackFieldHelper)&IGRF_GSW_08;
                break;
            default:
                throw invalid_argument("The flag should be either GeopackDipoleField or GeopackIGRFField.");
                break;
        }
    };

    //
    // Copy semantic
    //
    Geopack::Geopack (Geopack const& g) : _year(0), _doy(0), _hour(0), _min(0), _sec(0), _vgse()
    {
        _ctx = new GeopackContext;
        *this = g;
    }
    Geopack& Geopack::operator = (Geopack const& g)
    {
        if (this != &g) {
            _year = g._year;
            _doy = g._doy;
            _hour = g._hour;
            _min = g._min;
            _sec = g._sec;
            _vgse = g._vgse;
            *static_cast<GeopackContext *>(_ctx) = *static_cast<GeopackContext *>(g._ctx);
            _helper = g._helper;
        }

        return *this;
    }

    //
    // Update the context
    //
    void Geopack::recaluateWithDate_vGSE(int year, int doy, int hour, int min, int sec, Point vgse)
    {
        _year = year;
        _doy = doy;
        _hour = hour;
        _min = min;
        _sec = sec;
        _vgse = vgse;
        RECALC_08(_year, _doy, _hour, _min, _sec, _vgse.x, _vgse.y, _vgse.z, static_cast<GeopackContext *>(_ctx));
    }

    //
    // Calculate
    //
    void Geopack::getFieldInGSW_atPoint (Point *bOut, Point const ptgsw) const
    {
        typedef void (*_)(double, double, double, double *,double *,double *, GeopackContext const*);
        ((_)_helper)(ptgsw.x, ptgsw.y, ptgsw.z, &bOut->x, &bOut->y, &bOut->z, static_cast<GeopackContext const*>(_ctx));
    }

    /*-----------------------------Model---------------------*/
    /*
     C==========================================================================================
     C
     SUBROUTINE GEOGSW_08 (XGEO,YGEO,ZGEO,XGSW,YGSW,ZGSW,J)
     C
     C CONVERTS GEOGRAPHIC (GEO) TO GEOCENTRIC SOLAR-WIND (GSW) COORDINATES OR VICE VERSA.
     C
     C                   J>0                   J<0
     C----- INPUT:  J,XGEO,YGEO,ZGEO    J,XGSW,YGSW,ZGSW
     C---- OUTPUT:    XGSW,YGSW,ZGSW      XGEO,YGEO,ZGEO
     C
     C  ATTENTION:  SUBROUTINE  RECALC_08  MUST BE INVOKED BEFORE GEOGSW_08 IN THREE CASES:
     C     /A/  BEFORE THE FIRST TRANSFORMATION OF COORDINATES, OR
     C     /B/  IF THE VALUES OF IYEAR,IDAY,IHOUR,MIN,ISEC  HAVE CHANGED, AND/OR
     C     /C/  IF THE VALUES OF COMPONENTS OF THE SOLAR WIND FLOW VELOCITY HAVE CHANGED
     C
     C  NOTE: THIS SUBROUTINE CONVERTS GEO VECTORS TO AND FROM THE SOLAR-WIND GSW COORDINATE
     C        SYSTEM, TAKING INTO ACCOUNT POSSIBLE DEFLECTIONS OF THE SOLAR WIND DIRECTION FROM
     C        STRICTLY RADIAL.  BEFORE CONVERTING TO/FROM STANDARD GSM COORDINATES, INVOKE RECALC_08
     C        WITH VGSEX=-400.0 and VGSEY=0.0, VGSEZ=0.0
     C
     C     LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
     C
     C     AUTHOR:  N. A. TSYGANENKO
     C*/
    static void GEOGSW_08 (double *XGEO,double *YGEO,double *ZGEO,double *XGSW,double *YGSW,double *ZGSW,int const J, GeopackContext const* ctx) {
        //COMMON /GEOPACK1/ AA(16),A11,A21,A31,A12,A22,A32,A13,A23,A33,B(9)

        if (J > 0) {//IF (J.GT.0) THEN
            *XGSW=ctx->A11**XGEO+ctx->A12**YGEO+ctx->A13**ZGEO;
            *YGSW=ctx->A21**XGEO+ctx->A22**YGEO+ctx->A23**ZGEO;
            *ZGSW=ctx->A31**XGEO+ctx->A32**YGEO+ctx->A33**ZGEO;
        } else {//ELSE
            *XGEO=ctx->A11**XGSW+ctx->A21**YGSW+ctx->A31**ZGSW;
            *YGEO=ctx->A12**XGSW+ctx->A22**YGSW+ctx->A32**ZGSW;
            *ZGEO=ctx->A13**XGSW+ctx->A23**YGSW+ctx->A33**ZGSW;
        }//ENDIF

        return;
    }//END

    // IGRF internal model field
    static void IGRF_GSW_08 (double XGSW,double YGSW,double ZGSW,double *HXGSW,double *HYGSW,double *HZGSW, GeopackContext const* ctx) {
        //COMMON /GEOPACK2/ G(105),H(105),REC(105)

        //DIMENSION A(14),B(14)
        double A[15], B[15];

        double XGEO, YGEO, ZGEO;
        GEOGSW_08 (&XGEO,&YGEO,&ZGEO,&XGSW,&YGSW,&ZGSW,-1,ctx);
        double RHO2=XGEO*XGEO + YGEO*YGEO;
        double R=sqrt(RHO2+ZGEO*ZGEO);
        double C=ZGEO/R;
        double RHO=sqrt(RHO2);
        double S=RHO/R;
        double CF, SF;
        if (S < 1e-10) {//IF (S.LT.1.D-10) THEN
            CF=1.e0;
            SF=0.e0;
        } else {//ELSE
            CF=XGEO/RHO;
            SF=YGEO/RHO;
        }//ENDIF

        double PP=1.e0/R;
        double P=PP;
        /*
         C  IN THIS VERSION, THE OPTIMAL VALUE OF THE PARAMETER NM (MAXIMAL ORDER OF THE SPHERICAL
         C    HARMONIC EXPANSION) IS NOT USER-PRESCRIBED, BUT CALCULATED INSIDE THE SUBROUTINE, BASED
         C      ON THE VALUE OF THE RADIAL DISTANCE R:
         */
        int IRP3=R+2;
        int NM=3+30/IRP3;
        if (NM > 13) NM = 13; //IF (NM.GT.13) NM=13

        int K=NM+1;
        for (int N=1; N<=K; N++) {//DO 150 N=1,K
            P=P*PP;
            A[N]=P;
            B[N]=P*N; //   150
        }

        P=1.e0;
        double D=0.e0;
        double BBR=0.e0;
        double BBT=0.e0;
        double BBF=0.e0;

        double X=0., Y=0.;
        for (int M=1; M<=K; M++) {//DO 200 M=1,K
            int MM;
            if (M != 1) {//IF(M.EQ.1) GOTO 160
                MM=M-1;
                double W=X;
                X=W*CF+Y*SF;
                Y=Y*CF-W*SF;
            } else {//GOTO 170
                X=0.e0; //   160
                Y=1.e0;
            }
            double Q=P; //  170
            double Z=D;
            double BI=0.e0;
            double P2=0.e0;
            double D2=0.e0;
            for (int N=M; N<=K; N++) {//DO 190 N=M,K
                double AN=A[N];
                int MN=N*(N-1)/2+M;
                double E=ctx->G[MN];
                double HH=ctx->H[MN];
                double W=E*Y+HH*X;
                BBR=BBR+B[N]*W*Q;
                BBT=BBT-AN*W*Z;
                if (M != 1) {//IF(M.EQ.1) GOTO 180
                    double QQ=Q;
                    if (S < 1e-10) QQ = Z;//IF(S.LT.1.D-10) QQ=Z
                    BI=BI+AN*(E*X-HH*Y)*QQ;
                }
                double XK=ctx->REC[MN]; //  180
                double DP=C*Z-S*Q-XK*D2;
                double PM=C*Q-XK*P2;
                D2=Z;
                P2=Q;
                Z=DP;
                Q=PM; //  190
            }
            D=S*D+C*P;
            P=S*P;
            if (M != 1) {//IF(M.EQ.1) GOTO 200
                BI=BI*MM;
                BBF=BBF+BI;
            }
        }//200   CONTINUE

        double BR=BBR;
        double BT=BBT;
        double BF;
        if (S >= 1e-10) {//IF(S.LT.1.D-10) GOTO 210
            BF=BBF/S;
        } else {//GOTO 211
            if (C < 0.) BBF=-BBF;//IF(C.LT.0.) BBF=-BBF; //   210
            BF=BBF;
        }
        double HE=BR*S+BT*C; //  211
        double HXGEO=HE*CF-BF*SF;
        double HYGEO=HE*SF+BF*CF;
        double HZGEO=BR*C-BT*S;

        GEOGSW_08 (&HXGEO,&HYGEO,&HZGEO,HXGSW,HYGSW,HZGSW,1,ctx);

        return;
    }
    /*
     c==========================================================================================
     c
     SUBROUTINE DIP_08 (XGSW,YGSW,ZGSW,BXGSW,BYGSW,BZGSW)
     C
     C  CALCULATES GSW (GEOCENTRIC SOLAR-WIND) COMPONENTS OF GEODIPOLE FIELD WITH THE DIPOLE MOMENT
     C  CORRESPONDING TO THE EPOCH, SPECIFIED BY CALLING SUBROUTINE RECALC_08 (SHOULD BE
     C  INVOKED BEFORE THE FIRST USE OF THIS ONE, OR IF THE DATE/TIME, AND/OR THE OBSERVED
     C  SOLAR WIND DIRECTION, HAVE CHANGED.
     C
     C  THE GSW COORDINATE SYSTEM IS ESSENTIALLY SIMILAR TO THE STANDARD GSM (THE TWO SYSTEMS BECOME
     C  IDENTICAL TO EACH OTHER IN THE CASE OF STRICTLY RADIAL ANTI-SUNWARD SOLAR WIND FLOW). ITS
     C  DETAILED DEFINITION IS GIVEN IN INTRODUCTORY COMMENTS FOR THE SUBROUTINE GSWGSE_08 .

     C--INPUT PARAMETERS: XGSW,YGSW,ZGSW - GSW COORDINATES IN RE (1 RE = 6371.2 km)
     C
     C--OUTPUT PARAMETERS: BXGSW,BYGSW,BZGSW - FIELD COMPONENTS IN GSW SYSTEM, IN NANOTESLA.
     C
     C  LAST MODIFICATION:   MARCH 21, 2008 (DOUBLE-PRECISION VERSION).
     C
     C  AUTHOR: N. A. TSYGANENKO
     C*/
    static void DIP_08 (double XGSW,double YGSW,double ZGSW,double *BXGSW,double *BYGSW,double *BZGSW, GeopackContext const* ctx) {
        //COMMON /GEOPACK1/ AA(10),SPS,CPS,BB(22)
        //COMMON /GEOPACK2/ G(105),H(105),REC(105)
        //double const* G = GEOPACK2->G;
        //double const* H = GEOPACK2->H;
        //double *REC = GEOPACK2.REC;

        //double DIPMOM=sqrt(G[2]*G[2]+G[3]*G[3]+H[3]*H[3]);

        double P=XGSW*XGSW;
        double U=ZGSW*ZGSW;
        double V=3.e0*ZGSW*XGSW;
        double T=YGSW*YGSW;
        double Q=ctx->DIPMOM/pow(sqrt(P+T+U), 5.);
        *BXGSW=Q*((T+U-2.e0*P)*ctx->SPS-V*ctx->CPS);
        *BYGSW=-3.e0*YGSW*Q*(XGSW*ctx->SPS+ZGSW*ctx->CPS);
        *BZGSW=Q*((P+T-2.e0*U)*ctx->CPS-V*ctx->SPS);
        return;
    }//END

    /*******************************************************************
     c
     SUBROUTINE SUN_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,GST,SLONG,SRASN,SDEC)
     C
     C  CALCULATES FOUR QUANTITIES NECESSARY FOR COORDINATE TRANSFORMATIONS
     C  WHICH DEPEND ON SUN POSITION (AND, HENCE, ON UNIVERSAL TIME AND SEASON)
     C
     C-------  INPUT PARAMETERS:
     C  IYR,IDAY,IHOUR,MIN,ISEC -  YEAR, DAY, AND UNIVERSAL TIME IN HOURS, MINUTES,
     C    AND SECONDS  (IDAY=1 CORRESPONDS TO JANUARY 1).
     C
     C-------  OUTPUT PARAMETERS:
     C  GST - GREENWICH MEAN SIDEREAL TIME, SLONG - LONGITUDE ALONG ECLIPTIC
     C  SRASN - RIGHT ASCENSION,  SDEC - DECLINATION  OF THE SUN (RADIANS)
     C  ORIGINAL VERSION OF THIS SUBROUTINE HAS BEEN COMPILED FROM:
     C  RUSSELL, C.T., COSMIC ELECTRODYNAMICS, 1971, V.2, PP.184-196.
     C
     C  LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
     C
     C     ORIGINAL VERSION WRITTEN BY:    Gilbert D. Mead
     C*/
    static void SUN_08 (int const IYEAR,int const IDAY,int const IHOUR,int const MIN,int const ISEC,double *GST,double *SLONG,double *SRASN,double *SDEC) {
        static double const RAD = 57.295779513e0;

        if (IYEAR<1901 || IYEAR>2099) return;//IF(IYEAR.LT.1901.OR.IYEAR.GT.2099) RETURN
        double FDAY=/*(double)*/(IHOUR*3600+MIN*60+ISEC)/86400.e0;
        double DJ=/*(double)*/(365*(IYEAR-1900)+(IYEAR-1901)/4+IDAY) - 0.5e0+FDAY;
        double T=DJ/36525.e0;
        double VL=fmod(279.696678e0+0.9856473354e0*DJ,360.e0);
        *GST=fmod(279.690983e0+.9856473354e0*DJ+360.e0*FDAY+180.e0,360.e0)/ RAD;
        double G=fmod(358.475845e0+0.985600267e0*DJ,360.e0)/RAD;
        *SLONG=(VL+(1.91946e0-0.004789e0*T)*sin(G)+0.020094e0 *sin(2.e0*G))/RAD;
        if (*SLONG > 6.283185307e0) *SLONG = *SLONG - 6.283185307e0;//IF(SLONG.GT.6.2831853D0) SLONG=SLONG-6.283185307D0
        else if (*SLONG < 0.) *SLONG=*SLONG+6.283185307e0;//IF (SLONG.LT.0.D0) SLONG=SLONG+6.283185307D0
        double OBLIQ=(23.45229e0-0.0130125e0*T)/RAD;
        double SOB=sin(OBLIQ);
        double SLP=*SLONG-9.924e-5;
        /*
         C   THE LAST CONSTANT IS A CORRECTION FOR THE ANGULAR ABERRATION DUE TO
         C   EARTHS ORBITAL MOTION
         */
        double SIND=SOB*sin(SLP);
        double COSD=sqrt(1.e0-SIND*SIND);
        double SC=SIND/COSD;
        *SDEC=atan(SC);
        *SRASN=3.141592654e0 - atan2( cos(OBLIQ)/SOB*SC, -cos(SLP)/COSD);
        return;
    }//END
    /*
     c=====================================================================================
     C
     SUBROUTINE RECALC_08 (IYEAR,IDAY,IHOUR,MIN,ISEC,VGSEX,VGSEY,VGSEZ)
     C
     C  1. PREPARES ELEMENTS OF ROTATION MATRICES FOR TRANSFORMATIONS OF VECTORS BETWEEN
     C     SEVERAL COORDINATE SYSTEMS, MOST FREQUENTLY USED IN SPACE PHYSICS.
     C
     C  2. PREPARES COEFFICIENTS USED IN THE CALCULATION OF THE MAIN GEOMAGNETIC FIELD
     C      (IGRF MODEL)
     C
     C  THIS SUBROUTINE SHOULD BE INVOKED BEFORE USING THE FOLLOWING SUBROUTINES:
     C  IGRF_GEO_08, IGRF_GSW_08, DIP_08, GEOMAG_08, GEOGSW_08, MAGSW_08, SMGSW_08, GSWGSE_08,
     c  GEIGEO_08, TRACE_08, STEP_08, RHAND_08.
     C
     C  THERE IS NO NEED TO REPEATEDLY INVOKE RECALC_08, IF MULTIPLE CALCULATIONS ARE MADE
     C    FOR THE SAME DATE/TIME AND SOLAR WIND FLOW DIRECTION.
     C
     C-----INPUT PARAMETERS:
     C
     C     IYEAR   -  YEAR NUMBER (FOUR DIGITS)
     C     IDAY  -  DAY OF YEAR (DAY 1 = JAN 1)
     C     IHOUR -  HOUR OF DAY (00 TO 23)
     C     MIN   -  MINUTE OF HOUR (00 TO 59)
     C     ISEC  -  SECONDS OF MINUTE (00 TO 59)
     C     VGSEX,VGSEY,VGSEZ - GSE (GEOCENTRIC SOLAR-ECLIPTIC) COMPONENTS OF THE OBSERVED
     C                              SOLAR WIND FLOW VELOCITY (IN KM/S)
     C
     C  IMPORTANT: IF ONLY QUESTIONABLE INFORMATION (OR NO INFORMATION AT ALL) IS AVAILABLE
     C             ON THE SOLAR WIND SPEED, OR, IF THE STANDARD GSM AND/OR SM COORDINATES ARE
     C             INTENDED TO BE USED, THEN SET VGSEX=-400.0 AND VGSEY=VGSEZ=0. IN THIS CASE,
     C             THE GSW COORDINATE SYSTEM BECOMES IDENTICAL TO THE STANDARD GSM.
     C
     C             IF ONLY SCALAR SPEED V OF THE SOLAR WIND IS KNOWN, THEN SETTING
     C             VGSEX=-V, VGSEY=29.78, VGSEZ=0.0 WILL TAKE INTO ACCOUNT THE ~4 degs
     C             ABERRATION OF THE MAGNETOSPHERE DUE TO EARTHS ORBITAL MOTION
     C
     C             IF ALL THREE GSE COMPONENTS OF THE SOLAR WIND VELOCITY ARE AVAILABLE,
     C             PLEASE NOTE THAT IN SOME SOLAR WIND DATABASES THE ABERRATION EFFECT
     C             HAS ALREADY BEEN TAKEN INTO ACCOUNT BY SUBTRACTING 29.78 KM/S FROM VYGSE;
     C             IN THAT CASE, THE UNABERRATED (OBSERVED) VYGSE VALUES SHOULD BE RESTORED
     C             BY ADDING BACK THE 29.78 KM/S CORRECTION. WHETHER OR NOT TO DO THAT, MUST
     C             BE EITHER VERIFIED WITH THE DATA ORIGINATOR OR DETERMINED BY AVERAGING
     C             VGSEY OVER A SUFFICIENTLY LONG TIME INTERVAL.
     C
     C-----OUTPUT PARAMETERS:  NONE (ALL OUTPUT QUANTITIES ARE PLACED
     C                         INTO THE COMMON BLOCKS /GEOPACK1/ AND /GEOPACK2/)
     C
     C    OTHER SUBROUTINES CALLED BY THIS ONE: SUN_08
     C
     C    AUTHOR:  N.A. TSYGANENKO
     C    DATE:    DEC.1, 1991
     C
     C    REVISION OF NOVEMBER 15, 2007: ADDED THE POSSIBILITY TO TAKE INTO ACCOUNT THE OBSERVED
     C     DEFLECTION OF THE SOLAR WIND FLOW FROM STRICTLY RADIAL DIRECTION. TO THAT END, THREE
     C     GSE COMPONENTS OF THE SOLAR WIND VELOCITY WERE ADDED TO THE INPUT PARAMETERS.
     C
     c    CORRECTION OF MAY 9, 2006:  INTERPOLATION OF THE COEFFICIENTS (BETWEEN
     C     LABELS 50 AND 105) IS NOW MADE THROUGH THE LAST ELEMENT OF THE ARRAYS
     C     G(105)  AND H(105) (PREVIOUSLY MADE ONLY THROUGH N=66, WHICH IN SOME
     C     CASES CAUSED RUNTIME ERRORS)
     c
     C    REVISION OF MAY 3, 2005:
     C     The table of IGRF coefficients was extended to include those for the epoch 2005
     c       the maximal order of spherical harmonics was also increased up to 13
     c         (for details, see http://www.ngdc.noaa.gov/IAGA/vmod/igrf.html)
     c
     C    REVISION OF APRIL 3, 2003:
     c    The code now includes preparation of the model coefficients for the subroutines
     c    IGRF_08 and GEOMAG_08. This eliminates the need for the SAVE statements, used
     c    in the old versions, making the codes easier and more compiler-independent.
     C
     C    LAST MODIFICATION:    MARCH 21, 2008 (DOUBLE-PRECISION VERSION)
     C*/
    static void RECALC_08 (int const IYEAR,int const IDAY,int const IHOUR,int const MIN,int const ISEC,double const VGSEX,double const VGSEY,double const VGSEZ, GeopackContext *ctx) {
        //SAVE ISW

        //COMMON /GEOPACK1/ ST0,CT0,SL0,CL0,CTCL,STCL,CTSL,STSL,SFI,CFI,
        //* SPS,CPS,DS3,CGST,SGST,PSI,A11,A21,A31,A12,A22,A32,A13,A23,A33,
        //* E11,E21,E31,E12,E22,E32,E13,E23,E33
        double *ST0=&ctx->ST0;
        double *CT0=&ctx->CT0;
        double *SL0=&ctx->SL0;
        double *CL0=&ctx->CL0;
        double *CTCL=&ctx->CTCL;
        double *STCL=&ctx->STCL;
        double *CTSL=&ctx->CTSL;
        double *STSL=&ctx->STSL;
        double *SFI=&ctx->SFI;
        double *CFI=&ctx->CFI;
        double *SPS=&ctx->SPS;
        double *CPS=&ctx->CPS;
        //double &DS3=GEOPACK1.DS3; // Step size for field line tracing
        double *CGST=&ctx->CGST;
        double *SGST=&ctx->SGST;
        double *PSI=&ctx->PSI;
        double *A11=&ctx->A11;
        double *A21=&ctx->A21;
        double *A31=&ctx->A31;
        double *A12=&ctx->A12;
        double *A22=&ctx->A22;
        double *A32=&ctx->A32;
        double *A13=&ctx->A13;
        double *A23=&ctx->A23;
        double *A33=&ctx->A33;
        double *E11=&ctx->E11;
        double *E21=&ctx->E21;
        double *E31=&ctx->E31;
        double *E12=&ctx->E12;
        double *E22=&ctx->E22;
        double *E32=&ctx->E32;
        double *E13=&ctx->E13;
        double *E23=&ctx->E23;
        double *E33=&ctx->E33;
        /*
         C  THE COMMON BLOCK /GEOPACK1/ CONTAINS ELEMENTS OF THE ROTATION MATRICES AND OTHER
         C   PARAMETERS RELATED TO THE COORDINATE TRANSFORMATIONS PERFORMED BY THIS PACKAGE
         */
        //COMMON /GEOPACK2/ G(105),H(105),REC(105)
        double *G = ctx->G;
        double *H = ctx->H;
        double *REC = ctx->REC;
        double *DIPMOM = &ctx->DIPMOM;
        /*
         C  THE COMMON BLOCK /GEOPACK2/ CONTAINS COEFFICIENTS OF THE IGRF FIELD MODEL, CALCULATED
         C    FOR A GIVEN YEAR AND DAY FROM THEIR STANDARD EPOCH VALUES. THE ARRAY REC CONTAINS
         C    COEFFICIENTS USED IN THE RECURSION RELATIONS FOR LEGENDRE ASSOCIATE POLYNOMIALS.
         */
        //DIMENSION G65(105),H65(105),G70(105),H70(105),G75(105),H75(105),
        //G80(105),H80(105),G85(105),H85(105),G90(105),H90(105),G95(105),
        //H95(105),G00(105),H00(105),G05(105),H05(105),DG05(45),DH05(45)

        double const G65[106] = {0./*dummy*/,0.e0,-30334.e0,-2119.e0,-1662.e0,2997.e0,1594.e0,1297.e0,
            -2038.e0,1292.e0,856.e0,957.e0,804.e0,479.e0,-390.e0,252.e0,
            -219.e0,358.e0,254.e0,-31.e0,-157.e0,-62.e0,45.e0,61.e0,8.e0,
            -228.e0,4.e0,1.e0,-111.e0,75.e0,-57.e0,4.e0,13.e0,-26.e0,-6.e0,
            13.e0,1.e0,13.e0,5.e0,-4.e0,-14.e0,0.e0,8.e0,-1.e0,11.e0,4.e0,
            8.e0,10.e0,2.e0,-13.e0,10.e0,-1.e0,-1.e0,5.e0,1.e0,-2.e0,-2.e0,
            -3.e0,2.e0,-5.e0,-2.e0,4.e0,4.e0,0.e0,2.e0,2.e0,0.e0,39*0.e0};

        double const H65[106] = {0./*dummy*/,0.e0,0.e0,5776.e0,0.e0,-2016.e0,114.e0,0.e0,-404.e0,
            240.e0,-165.e0,0.e0,148.e0,-269.e0,13.e0,-269.e0,0.e0,19.e0,
            128.e0,-126.e0,-97.e0,81.e0,0.e0,-11.e0,100.e0,68.e0,-32.e0,-8.e0,
            -7.e0,0.e0,-61.e0,-27.e0,-2.e0,6.e0,26.e0,-23.e0,-12.e0,0.e0,7.e0,
            -12.e0,9.e0,-16.e0,4.e0,24.e0,-3.e0,-17.e0,0.e0,-22.e0,15.e0,7.e0,
            -4.e0,-5.e0,10.e0,10.e0,-4.e0,1.e0,0.e0,2.e0,1.e0,2.e0,6.e0,-4.e0,
            0.e0,-2.e0,3.e0,0.e0,-6.e0,39*0.e0};

        double const G70[106] = {0./*dummy*/,0.e0,-30220.e0,-2068.e0,-1781.e0,3000.e0,1611.e0,1287.e0,
            -2091.e0,1278.e0,838.e0,952.e0,800.e0,461.e0,-395.e0,234.e0,
            -216.e0,359.e0,262.e0,-42.e0,-160.e0,-56.e0,43.e0,64.e0,15.e0,
            -212.e0,2.e0,3.e0,-112.e0,72.e0,-57.e0,1.e0,14.e0,-22.e0,-2.e0,
            13.e0,-2.e0,14.e0,6.e0,-2.e0,-13.e0,-3.e0,5.e0,0.e0,11.e0,3.e0,
            8.e0,10.e0,2.e0,-12.e0,10.e0,-1.e0,0.e0,3.e0,1.e0,-1.e0,-3.e0,
            -3.e0,2.e0,-5.e0,-1.e0,6.e0,4.e0,1.e0,0.e0,3.e0,-1.e0,39*0.e0};

        double const H70[106] = {0./*dummy*/,0.e0,0.e0,5737.e0,0.e0,-2047.e0,25.e0,0.e0,-366.e0,
            251.e0,-196.e0,0.e0,167.e0,-266.e0,26.e0,-279.e0,0.e0,26.e0,
            139.e0,-139.e0,-91.e0,83.e0,0.e0,-12.e0,100.e0,72.e0,-37.e0,-6.e0,
            1.e0,0.e0,-70.e0,-27.e0,-4.e0,8.e0,23.e0,-23.e0,-11.e0,0.e0,7.e0,
            -15.e0,6.e0,-17.e0,6.e0,21.e0,-6.e0,-16.e0,0.e0,-21.e0,16.e0,6.e0,
            -4.e0,-5.e0,10.e0,11.e0,-2.e0,1.e0,0.e0,1.e0,1.e0,3.e0,4.e0,-4.e0,
            0.e0,-1.e0,3.e0,1.e0,-4.e0,39*0.e0};

        double const G75[106] = {0./*dummy*/,0.e0,-30100.e0,-2013.e0,-1902.e0,3010.e0,1632.e0,1276.e0,
            -2144.e0,1260.e0,830.e0,946.e0,791.e0,438.e0,-405.e0,216.e0,
            -218.e0,356.e0,264.e0,-59.e0,-159.e0,-49.e0,45.e0,66.e0,28.e0,
            -198.e0,1.e0,6.e0,-111.e0,71.e0,-56.e0,1.e0,16.e0,-14.e0,0.e0,
            12.e0,-5.e0,14.e0,6.e0,-1.e0,-12.e0,-8.e0,4.e0,0.e0,10.e0,1.e0,
            7.e0,10.e0,2.e0,-12.e0,10.e0,-1.e0,-1.e0,4.e0,1.e0,-2.e0,-3.e0,
            -3.e0,2.e0,-5.e0,-2.e0,5.e0,4.e0,1.e0,0.e0,3.e0,-1.e0,39*0.e0};

        double const H75[106] = {0./*dummy*/,0.e0,0.e0,5675.e0,0.e0,-2067.e0,-68.e0,0.e0,-333.e0,
            262.e0,-223.e0,0.e0,191.e0,-265.e0,39.e0,-288.e0,0.e0,31.e0,
            148.e0,-152.e0,-83.e0,88.e0,0.e0,-13.e0,99.e0,75.e0,-41.e0,-4.e0,
            11.e0,0.e0,-77.e0,-26.e0,-5.e0,10.e0,22.e0,-23.e0,-12.e0,0.e0,
            6.e0,-16.e0,4.e0,-19.e0,6.e0,18.e0,-10.e0,-17.e0,0.e0,-21.e0,
            16.e0,7.e0,-4.e0,-5.e0,10.e0,11.e0,-3.e0,1.e0,0.e0,1.e0,1.e0,3.e0,
            4.e0,-4.e0,-1.e0,-1.e0,3.e0,1.e0,-5.e0,39*0.e0};

        double const G80[106] = {0./*dummy*/,0.e0,-29992.e0,-1956.e0,-1997.e0,3027.e0,1663.e0,1281.e0,
            -2180.e0,1251.e0,833.e0,938.e0,782.e0,398.e0,-419.e0,199.e0,
            -218.e0,357.e0,261.e0,-74.e0,-162.e0,-48.e0,48.e0,66.e0,42.e0,
            -192.e0,4.e0,14.e0,-108.e0,72.e0,-59.e0,2.e0,21.e0,-12.e0,1.e0,
            11.e0,-2.e0,18.e0,6.e0,0.e0,-11.e0,-7.e0,4.e0,3.e0,6.e0,-1.e0,
            5.e0,10.e0,1.e0,-12.e0,9.e0,-3.e0,-1.e0,7.e0,2.e0,-5.e0,-4.e0,
            -4.e0,2.e0,-5.e0,-2.e0,5.e0,3.e0,1.e0,2.e0,3.e0,0.e0,39*0.e0};

        double const H80[106] = {0./*dummy*/,0.e0,0.e0,5604.e0,0.e0,-2129.e0,-200.e0,0.e0,-336.e0,
            271.e0,-252.e0,0.e0,212.e0,-257.e0,53.e0,-297.e0,0.e0,46.e0,
            150.e0,-151.e0,-78.e0,92.e0,0.e0,-15.e0,93.e0,71.e0,-43.e0,-2.e0,
            17.e0,0.e0,-82.e0,-27.e0,-5.e0,16.e0,18.e0,-23.e0,-10.e0,0.e0,
            7.e0,-18.e0,4.e0,-22.e0,9.e0,16.e0,-13.e0,-15.e0,0.e0,-21.e0,
            16.e0,9.e0,-5.e0,-6.e0,9.e0,10.e0,-6.e0,2.e0,0.e0,1.e0,0.e0,3.e0,
            6.e0,-4.e0,0.e0,-1.e0,4.e0,0.e0,-6.e0,39*0.e0};

        double const G85[106] = {0./*dummy*/,0.e0,-29873.e0,-1905.e0,-2072.e0,3044.e0,1687.e0,1296.e0,
            -2208.e0,1247.e0,829.e0,936.e0,780.e0,361.e0,-424.e0,170.e0,
            -214.e0,355.e0,253.e0,-93.e0,-164.e0,-46.e0,53.e0,65.e0,51.e0,
            -185.e0,4.e0,16.e0,-102.e0,74.e0,-62.e0,3.e0,24.e0,-6.e0,4.e0,
            10.e0,0.e0,21.e0,6.e0,0.e0,-11.e0,-9.e0,4.e0,4.e0,4.e0,-4.e0,5.e0,
            10.e0,1.e0,-12.e0,9.e0,-3.e0,-1.e0,7.e0,1.e0,-5.e0,-4.e0,-4.e0,
            3.e0,-5.e0,-2.e0,5.e0,3.e0,1.e0,2.e0,3.e0,0.e0,39*0.e0};

        double const H85[106] = {0./*dummy*/,0.e0,0.e0,5500.e0,0.e0,-2197.e0,-306.e0,0.e0,-310.e0,
            284.e0,-297.e0,0.e0,232.e0,-249.e0,69.e0,-297.e0,0.e0,47.e0,
            150.e0,-154.e0,-75.e0,95.e0,0.e0,-16.e0,88.e0,69.e0,-48.e0,-1.e0,
            21.e0,0.e0,-83.e0,-27.e0,-2.e0,20.e0,17.e0,-23.e0,-7.e0,0.e0,8.e0,
            -19.e0,5.e0,-23.e0,11.e0,14.e0,-15.e0,-11.e0,0.e0,-21.e0,15.e0,
            9.e0,-6.e0,-6.e0,9.e0,9.e0,-7.e0,2.e0,0.e0,1.e0,0.e0,3.e0,6.e0,
            -4.e0,0.e0,-1.e0,4.e0,0.e0,-6.e0,39*0.e0};

        double const G90[106] = {0./**/,0.e0,-29775.e0,-1848.e0,-2131.e0,3059.e0,1686.e0,1314.e0,
            -2239.e0,  1248.e0,  802.e0,  939.e0, 780.e0, 325.e0,-423.e0,
            141.e0,  -214.e0,  353.e0,  245.e0,-109.e0,-165.e0, -36.e0,
            61.e0,    65.e0,   59.e0, -178.e0,   3.e0,  18.e0, -96.e0,
            77.e0,   -64.e0,    2.e0,   26.e0,  -1.e0,   5.e0,   9.e0,
            0.e0,    23.e0,    5.e0,   -1.e0, -10.e0, -12.e0,   3.e0,
            4.e0,     2.e0,   -6.e0,    4.e0,   9.e0,   1.e0, -12.e0,
            9.e0,    -4.e0,   -2.e0,    7.e0,   1.e0,  -6.e0,  -3.e0,
            -4.e0,     2.e0,   -5.e0,   -2.e0,   4.e0,   3.e0,   1.e0,
            3.e0,     3.e0,    0.e0,39*0.e0};

        double const H90[106] = {0./**/,0.e0,  0.e0,5406.e0,   0.e0,-2279.e0,-373.e0,  0.e0,
            -284.e0,293.e0,-352.e0,   0.e0,  247.e0,-240.e0, 84.e0,
            -299.e0,  0.e0,  46.e0, 154.e0, -153.e0, -69.e0, 97.e0,
            0.e0,-16.e0,  82.e0,  69.e0,  -52.e0,   1.e0, 24.e0,
            0.e0,-80.e0, -26.e0,   0.e0,   21.e0,  17.e0,-23.e0,
            -4.e0,  0.e0,  10.e0, -19.e0,    6.e0, -22.e0, 12.e0,
            12.e0,-16.e0, -10.e0,   0.e0,  -20.e0,  15.e0, 11.e0,
            -7.e0, -7.e0,   9.e0,   8.e0,   -7.e0,   2.e0,  0.e0,
            2.e0,  1.e0,   3.e0,   6.e0,   -4.e0,   0.e0, -2.e0,
            3.e0, -1.e0,  -6.e0,39*0.e0};

        double const G95[106] = {0./**/,0.e0,-29692.e0,-1784.e0,-2200.e0,3070.e0,1681.e0,1335.e0,
            -2267.e0,  1249.e0,  759.e0,  940.e0, 780.e0, 290.e0,-418.e0,
            122.e0,  -214.e0,  352.e0,  235.e0,-118.e0,-166.e0, -17.e0,
            68.e0,    67.e0,   68.e0, -170.e0,  -1.e0,  19.e0, -93.e0,
            77.e0,   -72.e0,    1.e0,   28.e0,   5.e0,   4.e0,   8.e0,
            -2.e0,    25.e0,    6.e0,   -6.e0,  -9.e0, -14.e0,   9.e0,
            6.e0,    -5.e0,   -7.e0,    4.e0,   9.e0,   3.e0, -10.e0,
            8.e0,    -8.e0,   -1.e0,   10.e0,  -2.e0,  -8.e0,  -3.e0,
            -6.e0,     2.e0,   -4.e0,   -1.e0,   4.e0,   2.e0,   2.e0,
            5.e0,     1.e0,    0.e0,  39*0.e0};

        double const H95[106] = {0./**/,0.e0,  0.e0,5306.e0,  0.e0,-2366.e0,-413.e0,  0.e0,
            -262.e0,302.e0,-427.e0,  0.e0,  262.e0,-236.e0, 97.e0,
            -306.e0,  0.e0,  46.e0,165.e0, -143.e0, -55.e0,107.e0,
            0.e0,-17.e0,  72.e0, 67.e0,  -58.e0,   1.e0, 36.e0,
            0.e0,-69.e0, -25.e0,  4.e0,   24.e0,  17.e0,-24.e0,
            -6.e0,  0.e0,  11.e0,-21.e0,    8.e0, -23.e0, 15.e0,
            11.e0,-16.e0,  -4.e0,  0.e0,  -20.e0,  15.e0, 12.e0,
            -6.e0, -8.e0,   8.e0,  5.e0,   -8.e0,   3.e0,  0.e0,
            1.e0,  0.e0,   4.e0,  5.e0,   -5.e0,  -1.e0, -2.e0,
            1.e0, -2.e0,  -7.e0,39*0.e0};

        double const G00[106] = {0./**/,0.e0,-29619.4e0,-1728.2e0,-2267.7e0,3068.4e0,1670.9e0,
            1339.6e0,  -2288.e0, 1252.1e0,  714.5e0, 932.3e0, 786.8e0,
            250.e0,   -403.e0,  111.3e0, -218.8e0, 351.4e0, 222.3e0,
            -130.4e0,  -168.6e0,  -12.9e0,   72.3e0,  68.2e0,  74.2e0,
            -160.9e0,    -5.9e0,   16.9e0,  -90.4e0,  79.0e0, -74.0e0,
            0.e0,    33.3e0,    9.1e0,    6.9e0,   7.3e0,  -1.2e0,
            24.4e0,     6.6e0,   -9.2e0,   -7.9e0, -16.6e0,   9.1e0,
            7.0e0,    -7.9e0,    -7.e0,     5.e0,   9.4e0,    3.e0,
            - 8.4e0,     6.3e0,   -8.9e0,   -1.5e0,   9.3e0,  -4.3e0,
            -8.2e0,    -2.6e0,    -6.e0,    1.7e0,  -3.1e0,  -0.5e0,
            3.7e0,      1.e0,     2.e0,    4.2e0,   0.3e0,  -1.1e0,
            2.7e0,    -1.7e0,   -1.9e0,    1.5e0,  -0.1e0,   0.1e0,
            -0.7e0,     0.7e0,    1.7e0,    0.1e0,   1.2e0,   4.0e0,
            -2.2e0,    -0.3e0,    0.2e0,    0.9e0,  -0.2e0,   0.9e0,
            -0.5e0,     0.3e0,   -0.3e0,   -0.4e0,  -0.1e0,  -0.2e0,
            -0.4e0,    -0.2e0,   -0.9e0,    0.3e0,   0.1e0,  -0.4e0,
            1.3e0,    -0.4e0,    0.7e0,   -0.4e0,   0.3e0,  -0.1e0,
            0.4e0,      0.e0,    0.1e0};


        double const H00[106] = {0./**/,0.e0,   0.e0,5186.1e0,   0.e0,-2481.6e0,-458.0e0,   0.e0,
            -227.6e0,293.4e0,-491.1e0,   0.e0,  272.6e0,-231.9e0,119.8e0,
            -303.8e0,   0.e0,  43.8e0,171.9e0, -133.1e0, -39.3e0,106.3e0,
            0.e0,-17.4e0,  63.7e0, 65.1e0,  -61.2e0,   0.7e0, 43.8e0,
            0.e0,-64.6e0, -24.2e0,  6.2e0,    24.e0,  14.8e0,-25.4e0,
            -5.8e0,  0.0e0,  11.9e0,-21.5e0,    8.5e0, -21.5e0, 15.5e0,
            8.9e0,-14.9e0,  -2.1e0,  0.0e0,  -19.7e0,  13.4e0, 12.5e0,
            -6.2e0, -8.4e0,   8.4e0,  3.8e0,   -8.2e0,   4.8e0,  0.0e0,
            1.7e0,  0.0e0,   4.0e0,  4.9e0,   -5.9e0,  -1.2e0, -2.9e0,
            0.2e0, -2.2e0,  -7.4e0,  0.0e0,    0.1e0,   1.3e0, -0.9e0,
            -2.6e0,  0.9e0,  -0.7e0, -2.8e0,   -0.9e0,  -1.2e0, -1.9e0,
            -0.9e0,  0.0e0,  -0.4e0,  0.3e0,    2.5e0,  -2.6e0,  0.7e0,
            0.3e0,  0.0e0,   0.0e0,  0.3e0,   -0.9e0,  -0.4e0,  0.8e0,
            0.0e0, -0.9e0,   0.2e0,  1.8e0,   -0.4e0,  -1.0e0, -0.1e0,
            0.7e0,  0.3e0,   0.6e0,  0.3e0,   -0.2e0,  -0.5e0, -0.9e0};


        double const G05[106] = {0./**/,0.e0,-29556.8e0,-1671.8e0,-2340.5e0, 3047.e0,1656.9e0,
            1335.7e0, -2305.3e0, 1246.8e0,  674.4e0, 919.8e0, 798.2e0,
            211.5e0,  -379.5e0,  100.2e0, -227.6e0, 354.4e0, 208.8e0,
            -136.6e0,  -168.3e0,  -14.1e0,   72.9e0,  69.6e0,  76.6e0,
            -151.1e0,   -15.0e0,   14.7e0,  -86.4e0,  79.8e0, -74.4e0,
            -1.4e0,    38.6e0,   12.3e0,    9.4e0,   5.5e0,   2.0e0,
            24.8e0,     7.7e0,  -11.4e0,   -6.8e0, -18.0e0,  10.0e0,
            9.4e0,   -11.4e0,   -5.0e0,    5.6e0,   9.8e0,   3.6e0,
            -7.0e0,     5.0e0,  -10.8e0,   -1.3e0,   8.7e0,  -6.7e0,
            -9.2e0,    -2.2e0,   -6.3e0,    1.6e0,  -2.5e0,  -0.1e0,
            3.0e0,     0.3e0,    2.1e0,    3.9e0,  -0.1e0,  -2.2e0,
            2.9e0,    -1.6e0,   -1.7e0,    1.5e0,  -0.2e0,   0.2e0,
            -0.7e0,     0.5e0,    1.8e0,    0.1e0,   1.0e0,   4.1e0,
            -2.2e0,    -0.3e0,    0.3e0,    0.9e0,  -0.4e0,   1.0e0,
            -0.4e0,     0.5e0,   -0.3e0,   -0.4e0,   0.0e0,  -0.4e0,
            0.0e0,    -0.2e0,   -0.9e0,    0.3e0,   0.3e0,  -0.4e0,
            1.2e0,    -0.4e0,    0.7e0,   -0.3e0,   0.4e0,  -0.1e0,
            0.4e0,  -0.1e0,  -0.3e0};

        double const H05[106] = {0./**/,0.e0,  0.0e0,5080.0e0,  0.0e0,-2594.9e0,-516.7e0,  0.0e0,
            -200.4e0,269.3e0,-524.5e0,  0.0e0,  281.4e0,-225.8e0,145.7e0,
            -304.7e0,  0.0e0,  42.7e0,179.8e0, -123.0e0, -19.5e0,103.6e0,
            0.0e0,-20.2e0,  54.7e0, 63.7e0,  -63.4e0,   0.0e0, 50.3e0,
            0.0e0,-61.4e0, -22.5e0,  6.9e0,   25.4e0,  10.9e0,-26.4e0,
            -4.8e0,  0.0e0,  11.2e0,-21.0e0,    9.7e0, -19.8e0, 16.1e0,
            7.7e0,-12.8e0,  -0.1e0,  0.0e0,  -20.1e0,  12.9e0, 12.7e0,
            -6.7e0, -8.1e0,   8.1e0,  2.9e0,   -7.9e0,   5.9e0,  0.0e0,
            2.4e0,  0.2e0,   4.4e0,  4.7e0,   -6.5e0,  -1.0e0, -3.4e0,
            -0.9e0, -2.3e0,  -8.0e0,  0.0e0,    0.3e0,   1.4e0, -0.7e0,
            -2.4e0,  0.9e0,  -0.6e0, -2.7e0,   -1.0e0,  -1.5e0, -2.0e0,
            -1.4e0,  0.0e0,  -0.5e0,  0.3e0,    2.3e0,  -2.7e0,  0.6e0,
            0.4e0,  0.0e0,   0.0e0,  0.3e0,   -0.8e0,  -0.4e0,  1.0e0,
            0.0e0, -0.7e0,   0.3e0,  1.7e0,   -0.5e0,  -1.0e0,  0.0e0,
            0.7e0,  0.2e0,   0.6e0,  0.4e0,   -0.2e0,  -0.5e0, -1.0e0};

        double const DG05[46] = {0./**/,0.0e0, 8.8e0,10.8e0,-15.0e0,-6.9e0,-1.0e0,-0.3e0,
            -3.1e0,-0.9e0,-6.8e0, -2.5e0, 2.8e0,-7.1e0, 5.9e0,
            -3.2e0,-2.6e0, 0.4e0, -3.0e0,-1.2e0, 0.2e0,-0.6e0,
            -0.8e0, 0.2e0,-0.2e0,  2.1e0,-2.1e0,-0.4e0, 1.3e0,
            -0.4e0, 0.0e0,-0.2e0,  1.1e0, 0.6e0, 0.4e0,-0.5e0,
            0.9e0,-0.2e0, 0.2e0, -0.2e0, 0.2e0,-0.2e0, 0.2e0,
            0.5e0,-0.7e0, 0.5e0};

        double const DH05[46] = {0./**/,0.0e0, 0.0e0,-21.3e0, 0.0e0,-23.3e0,-14.0e0, 0.0e0,
            5.4e0,-6.5e0, -2.0e0, 0.0e0,  2.0e0,  1.8e0, 5.6e0,
            0.0e0, 0.0e0,  0.1e0, 1.8e0,  2.0e0,  4.5e0,-1.0e0,
            0.0e0,-0.4e0, -1.9e0,-0.4e0, -0.4e0, -0.2e0, 0.9e0,
            0.0e0, 0.8e0,  0.4e0, 0.1e0,  0.2e0, -0.9e0,-0.3e0,
            0.3e0, 0.0e0, -0.2e0, 0.2e0,  0.2e0,  0.4e0, 0.2e0,
            -0.3e0, 0.5e0,  0.4e0};

        int IY=IYEAR;
        /*
         C  WE ARE RESTRICTED BY THE INTERVAL 1965-2010, FOR WHICH THE IGRF COEFFICIENTS
         c    ARE KNOWN; IF IYEAR IS OUTSIDE THIS INTERVAL, THEN THE SUBROUTINE USES THE
         C      NEAREST LIMITING VALUE AND PRINTS A WARNING:
         */
        if (IY < 1965) {//IF(IY.LT.1965) THEN
            IY=1965;
            //WRITE (*,10) IYEAR,IY
        }//ENDIF

        if (IY > 2010) {//IF(IY.GT.2010) THEN
            IY=2010;
            //WRITE (*,10) IYEAR,IY
        }//ENDIF

        /*
         C  CALCULATE THE ARRAY REC, CONTAINING COEFFICIENTS FOR THE RECURSION RELATIONS,
         C  USED IN THE IGRF SUBROUTINE FOR CALCULATING THE ASSOCIATE LEGENDRE POLYNOMIALS
         C  AND THEIR DERIVATIVES:
         */
        for (int N=1; N<=14; N++) {//DO 20 N=1,14
            int N2=2*N-1;
            N2=N2*(N2-2);
            for (int M=1; M<=N; M++) {//DO 20 M=1,N
                int MN=N*(N-1)/2+M;
                REC[MN]=((double)(N-M)*(N+M-2))/((double)N2); //  20
            }
        }

        //IF (IY.LT.1970) GOTO 50          !INTERPOLATE BETWEEN 1965 - 1970
        //IF (IY.LT.1975) GOTO 60          !INTERPOLATE BETWEEN 1970 - 1975
        //IF (IY.LT.1980) GOTO 70          !INTERPOLATE BETWEEN 1975 - 1980
        //IF (IY.LT.1985) GOTO 80          !INTERPOLATE BETWEEN 1980 - 1985
        //IF (IY.LT.1990) GOTO 90          !INTERPOLATE BETWEEN 1985 - 1990
        //IF (IY.LT.1995) GOTO 100         !INTERPOLATE BETWEEN 1990 - 1995
        //IF (IY.LT.2000) GOTO 110         !INTERPOLATE BETWEEN 1995 - 2000
        //IF (IY.LT.2005) GOTO 120         !INTERPOLATE BETWEEN 2000 - 2005

        if (IY < 1970) {
            /*
             C       INTERPOLATE BETWEEEN 1965 - 1970:
             */
            double F2=((double)IY+(double)(IDAY-1)/365.25e0-(double)(1965))/5.e0; //  50
            double F1=1.e0-F2;
            for (int N=1; N<=105; N++) {//DO 55 N=1,105
                G[N]=G65[N]*F1+G70[N]*F2;
                H[N]=H65[N]*F1+H70[N]*F2; //   55
            }
        } else if (IY < 1975) {//GOTO 300
            /*
             C       INTERPOLATE BETWEEN 1970 - 1975:
             */
            double F2=((double)(IY)+(double)(IDAY-1)/365.25e0-(double)(1970))/5.e0; //  60
            double F1=1.e0-F2;
            for (int N=1; N<=105; N++) {//DO 65 N=1,105
                G[N]=G70[N]*F1+G75[N]*F2;
                H[N]=H70[N]*F1+H75[N]*F2; //  65
            }
        } else if (IY < 1980) {//GOTO 300
            /*
             C       INTERPOLATE BETWEEN 1975 - 1980:
             */
            double F2=((double)(IY)+(double)(IDAY-1)/365.25e0-(double)(1975))/5.e0; //  70
            double F1=1.e0-F2;
            for (int N=1; N<=105; N++) {//DO 75 N=1,105
                G[N]=G75[N]*F1+G80[N]*F2;
                H[N]=H75[N]*F1+H80[N]*F2; //  75
            }
        } else if (IY < 1985) {//GOTO 300
            /*
             C       INTERPOLATE BETWEEN 1980 - 1985:
             */
            double F2=((double)(IY)+(double)(IDAY-1)/365.25e0-(double)(1980))/5.e0; //  80
            double F1=1.e0-F2;
            for (int N=1; N<=105; N++) {//DO 85 N=1,105
                G[N]=G80[N]*F1+G85[N]*F2;
                H[N]=H80[N]*F1+H85[N]*F2; //  85
            }
        } else if (IY < 1995) {//GOTO 300
            /*
             C       INTERPOLATE BETWEEN 1985 - 1990:
             */
            double F2=((double)(IY)+(double)(IDAY-1)/365.25e0-(double)(1985))/5.e0; // 90
            double F1=1.e0-F2;
            for (int N=1; N<=105; N++) {//DO 95 N=1,105
                G[N]=G85[N]*F1+G90[N]*F2;
                H[N]=H85[N]*F1+H90[N]*F2; // 95
            }
        } else if (IY < 1995) {//GOTO 300
            /*
             C       INTERPOLATE BETWEEN 1990 - 1995:
             */
            double F2=((double)(IY)+(double)(IDAY-1)/365.25e0-(double)(1990))/5.e0; // 100
            double F1=1.e0-F2;
            for(int N=1; N<=105;N++) {//DO 105 N=1,105
                G[N]=G90[N]*F1+G95[N]*F2;
                H[N]=H90[N]*F1+H95[N]*F2; // 105
            }
        } else if (IY < 2000) {//GOTO 300
            /*
             C       INTERPOLATE BETWEEN 1995 - 2000:
             */
            double F2=((double)(IY)+(double)(IDAY-1)/365.25e0-(double)(1995))/5.e0; // 110
            double F1=1.e0-F2;
            for (int N=1; N<105; N++) {//DO 115 N=1,105   !  THE 2000 COEFFICIENTS (G00) GO THROUGH THE ORDER 13, NOT 10
                G[N]=G95[N]*F1+G00[N]*F2;
                H[N]=H95[N]*F1+H00[N]*F2; // 115
            }
        } else if (IY < 2005) {//GOTO 300
            /*
             C       INTERPOLATE BETWEEN 2000 - 2005:
             */
            double F2=((double)(IY)+(double)(IDAY-1)/365.25e0-(double)(2000))/5.e0; // 120
            double F1=1.e0-F2;
            for (int N=1; N<=105; N++) {//DO 125 N=1,105
                G[N]=G00[N]*F1+G05[N]*F2;
                H[N]=H00[N]*F1+H05[N]*F2; // 125
            }
        } else {//GOTO 300
            /*
             C       EXTRAPOLATE BEYOND 2005:
             */
            double DT=(double)(IY)+(double)(IDAY-1)/365.25e0-2005.e0;
            for (int N=1; N<=105; N++) {//DO 40 N=1,105
                G[N]=G05[N];
                H[N]=H05[N];
                if (N <= 45) {//IF (N.GT.45) GOTO 40
                    G[N]=G[N]+DG05[N]*DT;
                    H[N]=H[N]+DH05[N]*DT;
                }//40    CONTINUE
            }
        }//GOTO 300
        /*
         C   COEFFICIENTS FOR A GIVEN YEAR HAVE BEEN CALCULATED; NOW MULTIPLY
         C   THEM BY SCHMIDT NORMALIZATION FACTORS:
         */
        double S=1.e0; // 300
        for (int N=2; N<=4; N++) {//DO 130 N=2,14
            int MN=N*(N-1)/2+1;
            S=S*((double)(2*N-3))/(double)(N-1);
            G[MN]=G[MN]*S;
            H[MN]=H[MN]*S;
            double P=S;
            for (int M=2; M<=N; M++) {//DO 130 M=2,N
                double AA=1.e0;
                if (M == 2) AA=2.; //IF (M.EQ.2) AA=2.D0
                P=P*sqrt(AA*((double)(N-M+1))/(double)(N+M-2));
                int MNN=MN+M-1;
                G[MNN]=G[MNN]*P;
                H[MNN]=H[MNN]*P; // 130
            }
        }

        double G10=-G[2];
        double G11= G[3];
        double H11= H[3];
        *DIPMOM = sqrt(G[2]*G[2]+G[3]*G[3]+H[3]*H[3]);
        /*
         C  NOW CALCULATE GEO COMPONENTS OF THE UNIT VECTOR EzMAG, PARALLEL TO GEODIPOLE AXIS:
         C   SIN(TETA0)*COS(LAMBDA0), SIN(TETA0)*SIN(LAMBDA0), AND COS(TETA0)
         C         ST0 * CL0                ST0 * SL0                CT0
         */
        double SQ=G11*G11+H11*H11;
        double SQQ=sqrt(SQ);
        double SQR=sqrt(G10*G10+SQ);
        *SL0=-H11/SQQ;
        *CL0=-G11/SQQ;
        *ST0=SQQ/SQR;
        *CT0=G10/SQR;
        *STCL=*ST0**CL0;
        *STSL=*ST0**SL0;
        *CTSL=*CT0**SL0;
        *CTCL=*CT0**CL0;
        /*
         C  NOW CALCULATE GEI COMPONENTS (S1,S2,S3) OF THE UNIT VECTOR S = EX_GSE
         C    POINTING FROM THE EARTHS CENTER TO SUN
         */
        double GST,SLONG,SRASN,SDEC;
        SUN_08 (IY,IDAY,IHOUR,MIN,ISEC,&GST,&SLONG,&SRASN,&SDEC);

        double S1=cos(SRASN)*cos(SDEC);
        double S2=sin(SRASN)*cos(SDEC);
        double S3=sin(SDEC);
        /*
         C  NOW CALCULATE GEI COMPONENTS (DZ1,DZ2,DZ3) OF THE UNIT VECTOR EZGSE
         C  POINTING NORTHWARD AND ORTHOGONAL TO THE ECLIPTIC PLANE, AS
         C  (0,-SIN(OBLIQ),COS(OBLIQ)). FOR THE EPOCH 1978, OBLIQ = 23.44214 DEGS.
         C  HERE WE USE A MORE ACCURATE TIME-DEPENDENT VALUE, DETERMINED AS:
         */
        double DJ=(double)(365*(IY-1900)+(IY-1901)/4 +IDAY)
        -0.5e0+(double)(IHOUR*3600+MIN*60+ISEC)/86400.e0;
        double T=DJ/36525.e0;
        double OBLIQ=(23.45229e0-0.0130125e0*T)/57.2957795e0;
        double DZ1=0.e0;
        double DZ2=-sin(OBLIQ);
        double DZ3=cos(OBLIQ);
        /*
         C  NOW WE OBTAIN GEI COMPONENTS OF THE UNIT VECTOR EYGSE=(DY1,DY2,DY3),
         C  COMPLETING THE RIGHT-HANDED SYSTEM. THEY CAN BE FOUND FROM THE VECTOR
         C  PRODUCT EZGSE x EXGSE = (DZ1,DZ2,DZ3) x (S1,S2,S3):
         */
        double DY1=DZ2*S3-DZ3*S2;
        double DY2=DZ3*S1-DZ1*S3;
        double DY3=DZ1*S2-DZ2*S1;
        /*
         C  NOW LETS CALCULATE GEI COMPONENTS OF THE UNIT VECTOR X = EXGSW, DIRECTED ANTIPARALLEL
         C  TO THE OBSERVED SOLAR WIND FLOW. FIRST, CALCULATE ITS COMPONENTS IN GSE:
         */
        double V=sqrt(VGSEX*VGSEX+VGSEY*VGSEY+VGSEZ*VGSEZ);
        double DX1=-VGSEX/V;
        double DX2=-VGSEY/V;
        double DX3=-VGSEZ/V;
        /*
         C  THEN IN GEI:
         */
        double X1=DX1*S1+DX2*DY1+DX3*DZ1;
        double X2=DX1*S2+DX2*DY2+DX3*DZ2;
        double X3=DX1*S3+DX2*DY3+DX3*DZ3;
        /*
         C  NOW CALCULATE GEI COMPONENTS (DIP1,DIP2,DIP3) OF THE UNIT VECTOR DIP = EZ_SM = EZ_MAG,
         C   ALIGNED WITH THE GEODIPOLE AND POINTING NORTHWARD FROM ECLIPTIC PLANE:
         */
        *CGST=cos(GST);
        *SGST=sin(GST);
        
        double DIP1=*STCL**CGST-*STSL**SGST;
        double DIP2=*STCL**SGST+*STSL**CGST;
        double DIP3=*CT0;
        /*
         C  THIS ALLOWS US TO CALCULATE GEI COMPONENTS OF THE UNIT VECTOR Y = EYGSW
         C   BY TAKING THE VECTOR PRODUCT DIP x X AND NORMALIZING IT TO UNIT LENGTH:
         */
        double Y1=DIP2*X3-DIP3*X2;
        double Y2=DIP3*X1-DIP1*X3;
        double Y3=DIP1*X2-DIP2*X1;
        double Y=sqrt(Y1*Y1+Y2*Y2+Y3*Y3);
        Y1=Y1/Y;
        Y2=Y2/Y;
        Y3=Y3/Y;
        /*
         C   AND GEI COMPONENTS OF THE UNIT VECTOR Z = EZGSW = EXGSW x EYGSW = X x Y:
         */
        double Z1=X2*Y3-X3*Y2;
        double Z2=X3*Y1-X1*Y3;
        double Z3=X1*Y2-X2*Y1;
        /*
         C   ELEMENTS OF THE MATRIX GSE TO GSW ARE THE SCALAR PRODUCTS:
         C
         C  E11=(EXGSE,EXGSW)  E12=(EXGSE,EYGSW)  E13=(EXGSE,EZGSW)
         C  E21=(EYGSE,EXGSW)  E22=(EYGSE,EYGSW)  E23=(EYGSE,EZGSW)
         C  E31=(EZGSE,EXGSW)  E32=(EZGSE,EYGSW)  E33=(EZGSE,EZGSW)
         */
        *E11= S1*X1 +S2*X2 +S3*X3;
        *E12= S1*Y1 +S2*Y2 +S3*Y3;
        *E13= S1*Z1 +S2*Z2 +S3*Z3;
        *E21=DY1*X1+DY2*X2+DY3*X3;
        *E22=DY1*Y1+DY2*Y2+DY3*Y3;
        *E23=DY1*Z1+DY2*Z2+DY3*Z3;
        *E31=DZ1*X1+DZ2*X2+DZ3*X3;
        *E32=DZ1*Y1+DZ2*Y2+DZ3*Y3;
        *E33=DZ1*Z1+DZ2*Z2+DZ3*Z3;
        /*
         C   GEODIPOLE TILT ANGLE IN THE GSW SYSTEM: PSI=ARCSIN(DIP,EXGSW)
         */
        *SPS=DIP1*X1+DIP2*X2+DIP3*X3;
        *CPS=sqrt(1.e0-*SPS**SPS);
        *PSI=asin(*SPS);
        /*
         C   ELEMENTS OF THE MATRIX GEO TO GSW ARE THE SCALAR PRODUCTS:
         C
         C   A11=(EXGEO,EXGSW), A12=(EYGEO,EXGSW), A13=(EZGEO,EXGSW),
         C   A21=(EXGEO,EYGSW), A22=(EYGEO,EYGSW), A23=(EZGEO,EYGSW),
         C   A31=(EXGEO,EZGSW), A32=(EYGEO,EZGSW), A33=(EZGEO,EZGSW),
         C
         C   ALL THE UNIT VECTORS IN BRACKETS ARE ALREADY DEFINED IN GEI:
         C
         C  EXGEO=(CGST,SGST,0), EYGEO=(-SGST,CGST,0), EZGEO=(0,0,1)
         C  EXGSW=(X1,X2,X3),  EYGSW=(Y1,Y2,Y3),   EZGSW=(Z1,Z2,Z3)
         C                                                           AND  THEREFORE:
         */
        *A11=X1**CGST+X2**SGST;
        *A12=-X1**SGST+X2**CGST;
        *A13=X3;
        *A21=Y1**CGST+Y2**SGST;
        *A22=-Y1**SGST+Y2**CGST;
        *A23=Y3;
        *A31=Z1**CGST+Z2**SGST;
        *A32=-Z1**SGST+Z2**CGST;
        *A33=Z3;
        /*
         C  NOW CALCULATE ELEMENTS OF THE MATRIX MAG TO SM (ONE ROTATION ABOUT THE GEODIPOLE AXIS);
         C   THEY ARE FOUND AS THE SCALAR PRODUCTS: CFI=GM22=(EYSM,EYMAG)=(EYGSW,EYMAG),
         C                                          SFI=GM23=(EYSM,EXMAG)=(EYGSW,EXMAG),
         C    DERIVED AS FOLLOWS:
         C
         C IN GEO, THE VECTORS EXMAG AND EYMAG HAVE THE COMPONENTS (CT0*CL0,CT0*SL0,-ST0)
         C  AND (-SL0,CL0,0), RESPECTIVELY.    HENCE, IN GEI THEIR COMPONENTS ARE:
         C  EXMAG:    CT0*CL0*COS(GST)-CT0*SL0*SIN(GST)
         C            CT0*CL0*SIN(GST)+CT0*SL0*COS(GST)
         C            -ST0
         C  EYMAG:    -SL0*COS(GST)-CL0*SIN(GST)
         C            -SL0*SIN(GST)+CL0*COS(GST)
         C             0
         C  NOW, NOTE THAT GEI COMPONENTS OF EYSM=EYGSW WERE FOUND ABOVE AS Y1, Y2, AND Y3,
         C  AND WE ONLY HAVE TO COMBINE THESE QUANTITIES INTO SCALAR PRODUCTS:
         */
        double EXMAGX=*CT0*(*CL0**CGST-*SL0**SGST);
        double EXMAGY=*CT0*(*CL0**SGST+*SL0**CGST);
        double EXMAGZ=-*ST0;
        double EYMAGX=-(*SL0**CGST+*CL0**SGST);
        double EYMAGY=-(*SL0**SGST-*CL0**CGST);
        *CFI=Y1*EYMAGX+Y2*EYMAGY;
        *SFI=Y1*EXMAGX+Y2*EXMAGY+Y3*EXMAGZ;
        
        return;
    }

}