//
//  ubk_dlm_c.cpp
//  UBJDevelopment
//
//  Created by Kyungguk Min on 4/11/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#include "ubk_dlm_c.h"
#include "UBKLstarxx.h"
#include <cstdio>
#include <vector>
#include <cmath>
#include <algorithm>
#include <string>

using namespace std;
using namespace UBK;

//
// Log facility
//
UBK_INLINE void ASSERT(BOOL test, const char *msg) {
    if (!test) {
        IDL_Message(IDL_M_NAMED_GENERIC, IDL_MSG_LONGJMP, msg);
    }
}

//
// Helpers
//
UBK_INLINE void IDL_VectorCopy(IDL_VPTR src, vector<IDL_LONG> &dest) {
    if (IDL_TYP_LONG == src->type) {
        IDL_VarCopy(src, src=IDL_Gettmp());
    } else {
        src = IDL_CvtLng(1, &src);
    }

    IDL_MEMINT n;
    IDL_LONG *pd;
    IDL_VarGetData(src, &n, (char **)&pd, IDL_TRUE);
    dest.resize(n);
    copy(pd, pd + n, dest.begin());

    IDL_Deltmp(src);
}

UBK_INLINE void IDL_VectorCopy(vector<double> const& src, IDL_VPTR dst) {
    {
        IDL_VPTR vartmp;
        IDL_MakeTempVector(IDL_TYP_DOUBLE, src.size(), IDL_ARR_INI_NOP, &vartmp);
        IDL_VarCopy(vartmp, dst);
    }
    copy(src.begin(), src.end(), (double *)dst->value.arr->data);
}
UBK_INLINE void IDL_VectorCopy(IDL_VPTR src, vector<double> &dest) {
    if (IDL_TYP_DOUBLE == src->type) {
        IDL_VarCopy(src, src=IDL_Gettmp());
    } else {
        src = IDL_CvtDbl(1, &src);
    }

    IDL_MEMINT n;
    double *pd;
    IDL_VarGetData(src, &n, (char **)&pd, IDL_TRUE);
    dest.resize(n);
    copy(pd, pd + n, dest.begin());

    IDL_Deltmp(src);
}
UBK_INLINE void IDL_ArrayCopy(long firstDim, vector< vector<double> > const& src, IDL_VPTR dst) {
    {
        IDL_VPTR vartmp;
        IDL_MEMINT dim[2] = {firstDim, static_cast<IDL_MEMINT>(src.size())};
        IDL_MakeTempArray(IDL_TYP_DOUBLE, 2, dim, IDL_ARR_INI_NOP, &vartmp);
        IDL_VarCopy(vartmp, dst);

        double *p = (double *)dst->value.arr->data;
        double *pend = p + dim[0]*dim[1];
        for ( ; p<pend; ) {
            *p++ = NAN;
        }
    }
    double *p = (double *)dst->value.arr->data;
    for (long idx=0, cnt=src.size(); idx<cnt; idx++) {
        copy(src[idx].begin(), src[idx].end(), p + idx*firstDim);
    }
}

////////////////////////////////////////////////////////////////////////////////
// IDL Loader
////////////////////////////////////////////////////////////////////////////////

/*
 * Define message codes and their corresponding printf(3) format
 * strings. Note that message codes start at zero and each one is
 * one less that the previous one. Codes must be monotonic and
 * contiguous.
 */
static IDL_MSG_DEF msg_arr[] =
{
#define M_TM_INPRO                       0
	{  "M_TM_INPRO","%NIDL_UBKL* module has been loaded." },
};

/*
 * The load function fills in this message block handle with the
 * opaque handle to the message block used for this module. The other
 * routines can then use it to throw errors from this block.
 */
static IDL_MSG_BLOCK msg_block;

UBK_C_EXTERN int IDL_Load(void)
{
	int ret;
	/*
	 * These tables contain information on the functions and procedures
	 * that make up the TESTMODULE DLM. The information contained in these
	 * tables must be identical to that contained in testmodule.dlm.
	 */
	/*
	 * typedef struct {
	 *   IDL_SYSRTN_GENERIC funct_addr;
	 *   char *name;
	 *   unsigned short arg_min;
	 *   unsigned short arg_max;
	 *   int flags;
	 *   void *extra;
	 * } IDL_SYSFUN_DEF2;
	 */
	static IDL_SYSFUN_DEF2 procedure_addr[] = {
		{ (IDL_SYSRTN_GENERIC)ubk_cotrans, "UBK_COTRANS_C", 9, 9, 0, 0},
        { (IDL_SYSRTN_GENERIC)ubk_ts_field, "UBK_TS_FIELD_C", 12, 12, 0, 0},
        { (IDL_SYSRTN_GENERIC)ubk_field_line, "UBK_FIELD_LINE_C", 20, 20, 0, 0},
        { (IDL_SYSRTN_GENERIC)ubk_lstar, "UBK_LSTAR_C", 23, 23, 0, 0},
	};

	/*
	 * Create a message block to hold our messages. Save its handle where
	 * the other routines can access it.
	 */
	if (!(msg_block = IDL_MessageDefineBlock("IDL_UBKL*", \
                                             IDL_CARRAY_ELTS(msg_arr), msg_arr)))
		return IDL_FALSE;

	/*
	 * Register our routine. The routines must be specified exactly the same
	 * as in testmodule.dlm.
	 */

	if ((ret=IDL_SysRtnAdd(procedure_addr, IDL_FALSE, IDL_CARRAY_ELTS(procedure_addr))) == IDL_FALSE)
		IDL_Message(IDL_M_NAMED_GENERIC,IDL_MSG_RET, \
                    "Error adding system routine");
	return ret;
}

////////////////////////////////////////////////////////////////////////////////
// IDL Procedures
////////////////////////////////////////////////////////////////////////////////

//
// PRO UBK_COTRANS_C, TILTOUT, XOUT, YOUT, ZOUT, XIN, YIN, ZIN, YEAR, DOY, HOUR, MIN, SEC, TO_CO_SYSTEM, N_THREADS
//
UBK_C_EXTERN void ubk_cotrans(int argc, IDL_VPTR *argv,char *argk)
{
#ifdef DEBUG
    fprintf(stderr, "%%%% DEBUG:%s:%d: %s Enter.\n", __FILE__, __LINE__, __FUNCTION__);
#endif

    //
    // Inputs
    //
    vector<double> xin, yin, zin;
    IDL_VectorCopy(argv[3], xin);
    IDL_VectorCopy(argv[4], yin);
    IDL_VectorCopy(argv[5], zin);

    vector<IDL_LONG> d;
    IDL_VectorCopy(argv[6], d);

    MagneticFieldCoordinateSystem to_co_system = IDL_LongScalar(argv[7]); // 0: To GSM, 1: To SM
    unsigned long n_threads = IDL_ULongScalar(argv[8]);

    //
    // Minimal validation
    //
    IDL_MEMINT M = xin.size();

    ASSERT(M && yin.size()==M && zin.size()==M,
           "The length of the input position array is different.");
    ASSERT(5 == d.size(),
           "The number of element of DATE should be 5 ([year, doy, hour, min, sec]).");
    ASSERT(kMagneticFieldGSM==to_co_system || kMagneticFieldSM==to_co_system,
           "TO_CO_SYSTEM should be either 0 (SM->GSM) or 1 (GSM->SM).");
    ASSERT(n_threads >= 1,
           "N_THREADS parameter should be greater than 0.");

    //
    // Log
    //
#ifdef DEBUG
    fprintf(stderr, "%%%% DEBUG:%s:%d: Options\n", __FUNCTION__, __LINE__);
    fprintf(stderr, "\tM = %lld\n", M);
    fprintf(stderr, "\tDATE = [%d, %d, %d, %d, %d]\n", d[0], d[1], d[2], d[3], d[4]);
    fprintf(stderr, "\tTO_CO_SYSTEM = %ld\n", to_co_system);
    fprintf(stderr, "\tN_THREADS = %ld\n", n_threads);
#endif

    //
    // Output buffer
    //
    vector<double> xout(M), yout(M), zout(M);

    //
    // Calculation start
    //
    {
#ifdef DEBUG
        fprintf(stderr, "%%%% DEBUG:%s:%d: Loop block enter.\n", __FUNCTION__, __LINE__);
#endif

        class _ : public CotransCoordinator {
        public:
            vector<double> const *xin, *yin, *zin;
            vector<double> *xout, *yout, *zout;

            _(FieldModel const& fm) : CotransCoordinator(fm) {};
            unsigned long count() const {
                return xin->size();
            };
            Point from(long idx) const {
                return Point(xin->at(idx), yin->at(idx), zin->at(idx));
            };
            void callback(long idx, Point const& pt) const {
                xout->at(idx) = pt.x;
                yout->at(idx) = pt.y;
                zout->at(idx) = pt.z;
            };
        };

        double const vsw = 400.;
        Date date(d[0], d[1], d[2], d[3], d[4]);
        TSFieldModel fm(date, vsw, kGeopackDipoleField, 1, NULL, kTSNone);
        _ c(fm);

        c.setNThreads(n_threads);
        c.setToCoSystem(to_co_system);

        c.xin = &xin;
        c.yin = &yin;
        c.zin = &zin;
        c.xout = &xout;
        c.yout = &yout;
        c.zout = &zout;

        c.start();

#ifdef DEBUG
        fprintf(stderr, "%%%% DEBUG:%s:%d: Loop block exit.\n", __FUNCTION__, __LINE__);
#endif
    }

    //
    // Push result
    //
    IDL_VectorCopy(xout, argv[0]);
    IDL_VectorCopy(yout, argv[1]);
    IDL_VectorCopy(zout, argv[2]);

#ifdef DEBUG
    fprintf(stderr, "%%%% DEBUG:%s:%d: %s Enter.\n", __FILE__, __LINE__, __FUNCTION__);
#endif
}

//
// PRO UBK_TS_FIELD_C, TILTOUT, XOUT, YOUT, ZOUT, XIN, YIN, ZIN, YEAR, DOY, HOUR, MIN, SEC, TO_CO_SYSTEM, N_THREADS
//
UBK_C_EXTERN void ubk_ts_field(int argc, IDL_VPTR *argv,char *argk)
{
#ifdef DEBUG
    fprintf(stderr, "%%%% DEBUG:%s:%d: %s Enter.\n", __FILE__, __LINE__, __FUNCTION__);
#endif

    //
    // Inputs
    //
    vector<double> x, y, z;
    IDL_VectorCopy(argv[3], x);
    IDL_VectorCopy(argv[4], y);
    IDL_VectorCopy(argv[5], z);

    vector<IDL_LONG> d;
    IDL_VectorCopy(argv[6], d);

    vector<double> ioptparmod;
    IDL_VectorCopy(argv[7], ioptparmod);

    TSExternalFieldModel external = IDL_ULongScalar(argv[8]);
    GeopackInternalFieldModel internal = IDL_ULongScalar(argv[9]);;

    MagneticFieldCoordinateSystem co_system = IDL_LongScalar(argv[10]); // 0: To GSM, 1: To SM
    unsigned long n_threads = IDL_ULongScalar(argv[11]);

    //
    // Minimal validation
    //
    IDL_MEMINT M = x.size();

    ASSERT(M && y.size()==M && z.size()==M,
           "The length of the input position array is different.");
    ASSERT(5 == d.size(),
           "The number of element of DATE should be 5 ([year, doy, hour, min, sec]).");
    ASSERT(kGeopackDipoleField==internal || kGeopackIGRFField==internal,
           "Invalid internal field model flag.");
    ASSERT((kTSNone==external) ||
           (kTS89Model==external) ||
           (kTS96Model==external) ||
           (kTS02Model==external) ||
           (kTS05Model==external),
           "Invalid external field model flag.");
    ASSERT((kTSNone==external) ||
           (kTS89Model==external && 1==ioptparmod.size()) ||
           (10==ioptparmod.size()),
           "Invalid ioptparmod parameter.");
    ASSERT(kMagneticFieldGSM==co_system || kMagneticFieldSM==co_system,
           "CO_SYSTEM should be either 0 (GSM) or 1 (SM).");
    ASSERT(n_threads >= 1,
           "N_THREADS parameter should be greater than 0.");

    //
    // Log
    //
#ifdef DEBUG
    fprintf(stderr, "%%%% DEBUG:%s:%d: Options\n", __FUNCTION__, __LINE__);
    fprintf(stderr, "\tM = %lld\n", M);
    fprintf(stderr, "\tDATE = [%d, %d, %d, %d, %d]\n", d[0], d[1], d[2], d[3], d[4]);
    fprintf(stderr, "\tINTERNAL = %ld\n", internal);
    fprintf(stderr, "\tEXTERNAL = %ld\n", external);
    fprintf(stderr, "\tCO_SYSTEM = %ld\n", co_system);
    fprintf(stderr, "\tN_THREADS = %ld\n", n_threads);
#endif

    //
    // Output buffer
    //
    vector<double> Bx(M), By(M), Bz(M);

    //
    // Calculation start
    //
    {
#ifdef DEBUG
        fprintf(stderr, "%%%% DEBUG:%s:%d: Loop block enter.\n", __FUNCTION__, __LINE__);
#endif

        class _ : public BFieldCoordinator {
        public:
            vector<double> const *x, *y, *z;
            vector<double> *Bx, *By, *Bz;

            _(FieldModel const& fm) : BFieldCoordinator(fm) {};
            unsigned long count() const {
                return x->size();
            };
            Point at(long idx) const {
                return Point(x->at(idx), y->at(idx), z->at(idx));
            };
            void callback(long idx, Point const& pt) const {
                Bx->at(idx) = pt.x;
                By->at(idx) = pt.y;
                Bz->at(idx) = pt.z;
            };
        };

        double const vsw = 400.;
        Date date(d[0], d[1], d[2], d[3], d[4]);
        int iopt = (kTS89Model==external ? ioptparmod[0] : 1);
        double const* parmod = (kTS89Model==external && kTSNone==external ? NULL : &ioptparmod.front());
        TSFieldModel fm(date, vsw, internal, iopt, parmod, external);
        _ c(fm);
        c.setNThreads(n_threads);
        c.setCoSystem(co_system);

        c.x = &x;
        c.y = &y;
        c.z = &z;
        c.Bx = &Bx;
        c.By = &By;
        c.Bz = &Bz;

        c.start();

#ifdef DEBUG
        fprintf(stderr, "%%%% DEBUG:%s:%d: Loop block exit.\n", __FUNCTION__, __LINE__);
#endif
    }

    //
    // Push result
    //
    IDL_VectorCopy(Bx, argv[0]);
    IDL_VectorCopy(By, argv[1]);
    IDL_VectorCopy(Bz, argv[2]);

#ifdef DEBUG
    fprintf(stderr, "%%%% DEBUG:%s:%d: %s Enter.\n", __FILE__, __LINE__, __FUNCTION__);
#endif
}

//
//
//
UBK_C_EXTERN void ubk_field_line(int argc, IDL_VPTR *argv,char *argk)
{
#ifdef DEBUG
    fprintf(stderr, "%%%% DEBUG:%s:%d: %s Enter.\n", __FILE__, __LINE__, __FUNCTION__);
#endif

    //
    // Inputs
    //
    vector<double> x0, y0, z0;
    IDL_VectorCopy(argv[10], x0);
    IDL_VectorCopy(argv[11], y0);
    IDL_VectorCopy(argv[12], z0);

    vector<IDL_LONG> d;
    IDL_VectorCopy(argv[13], d);

    vector<double> ioptparmod;
    IDL_VectorCopy(argv[14], ioptparmod);

    TSExternalFieldModel external = IDL_ULongScalar(argv[15]);
    GeopackInternalFieldModel internal = IDL_ULongScalar(argv[16]);;

    double ionoR = IDL_DoubleScalar(argv[17]);
    double ds = IDL_DoubleScalar(argv[18]);

    unsigned long n_threads = IDL_ULongScalar(argv[19]);

    //
    // Minimal validation
    //
    IDL_MEMINT M = x0.size();

    ASSERT(M && y0.size()==M && z0.size()==M,
           "The length of the input position array is different.");
    ASSERT(5 == d.size(),
           "The number of element of DATE should be 5 ([year, doy, hour, min, sec]).");
    ASSERT(kGeopackDipoleField==internal || kGeopackIGRFField==internal,
           "Invalid internal field model flag.");
    ASSERT((kTSNone==external) ||
           (kTS89Model==external) ||
           (kTS96Model==external) ||
           (kTS02Model==external) ||
           (kTS05Model==external),
           "Invalid external field model flag.");
    ASSERT((kTSNone==external) ||
           (kTS89Model==external && 1==ioptparmod.size()) ||
           (10==ioptparmod.size()),
           "Invalid ioptparmod parameter.");
    ASSERT(ionoR >= 1.,
           "ionoR < 1.");
    ASSERT(ds > 0.,
           "ds <= 0.");
    ASSERT(n_threads >= 1,
           "N_THREADS parameter should be greater than 0.");

    //
    // Log
    //
#ifdef DEBUG
    fprintf(stderr, "%%%% DEBUG:%s:%d: Options\n", __FUNCTION__, __LINE__);
    fprintf(stderr, "\tM = %lld\n", M);
    fprintf(stderr, "\tDATE = [%d, %d, %d, %d, %d]\n", d[0], d[1], d[2], d[3], d[4]);
    fprintf(stderr, "\tINTERNAL = %ld\n", internal);
    fprintf(stderr, "\tEXTERNAL = %ld\n", external);
    fprintf(stderr, "\tIONOR = %f\n", ionoR);
    fprintf(stderr, "\tDS = %f\n", ds);
    fprintf(stderr, "\tN_THREADS = %ld\n", n_threads);
#endif

    //
    // Output buffer
    //
    vector< vector<double> > XYZmeq(M), XYZeq(M), XYZfoot(M);
    vector<double> Bmeq(M), Beq(M);
    vector< vector<double> > K(M), Bm(M), Xfl(M), Yfl(M), Zfl(M);

    //
    // Calculation start
    //
    {
#ifdef DEBUG
        fprintf(stderr, "%%%% DEBUG:%s:%d: Loop block enter.\n", __FUNCTION__, __LINE__);
#endif

        class _ : public FieldLineCoordinator {
        public:
            vector<double> const *x0, *y0, *z0;
            vector<double> *Bmeq, *Beq;
            vector< vector<double> > *XYZmeq, *XYZeq, *XYZfoot;
            vector< vector<double> > *K, *Bm, *Xfl, *Yfl, *Zfl;

            _(FieldModel const& fm) : FieldLineCoordinator(fm) {};
            unsigned long count() const {
                return x0->size();
            };
            Point from(long idx) const {
                return Point(x0->at(idx), y0->at(idx), z0->at(idx));
            };
            void callback(long idx, FieldLine const& fl) const {
                Bmeq->at(idx) = fl.magneticEquatorFieldMagnitude();
                Beq->at(idx) = fl.coordinateEquatorFieldMagnitude();
                XYZmeq->at(idx).insert(XYZmeq->at(idx).begin(), fl.magneticEquator().xyz, fl.magneticEquator().xyz+3);
                XYZeq->at(idx).insert(XYZeq->at(idx).begin(), fl.coordinateEquator().xyz, fl.coordinateEquator().xyz+3);
                XYZfoot->at(idx).insert(XYZfoot->at(idx).begin(), fl.footPoint().xyz, fl.footPoint().xyz+3);
                K->at(idx).insert(K->at(idx).begin(), fl.modifiedInvariants().begin(), fl.modifiedInvariants().end());
                Bm->at(idx).insert(Bm->at(idx).begin(), fl.mirrorMagnitudes().begin(), fl.mirrorMagnitudes().end());
                for (long jdx=0, cnt=fl.coordinates().size(); jdx<cnt; jdx++) {
                    Point const& pt = fl.coordinates()[jdx];
                    Xfl->at(idx).push_back(pt.x);
                    Yfl->at(idx).push_back(pt.y);
                    Zfl->at(idx).push_back(pt.z);
                }
            };
        };

        double const vsw = 400.;
        Date date(d[0], d[1], d[2], d[3], d[4]);
        int iopt = (kTS89Model==external ? ioptparmod[0] : 1);
        double const* parmod = (kTS89Model==external && kTSNone==external ? NULL : &ioptparmod.front());
        TSFieldModel fm(date, vsw, internal, iopt, parmod, external);

        _ c(fm);
        c.setNThreads(n_threads);
        c.setIonoR(ionoR);
        c.setDs(ds);

        c.x0 = &x0;
        c.y0 = &y0;
        c.z0 = &z0;
        c.Bmeq = &Bmeq;
        c.Beq = &Beq;
        c.XYZmeq = &XYZmeq;
        c.XYZeq = &XYZeq;
        c.XYZfoot = &XYZfoot;
        c.K = &K;
        c.Bm = &Bm;
        c.Xfl = &Xfl;
        c.Yfl = &Yfl;
        c.Zfl = &Zfl;

        c.start();

#ifdef DEBUG
        fprintf(stderr, "%%%% DEBUG:%s:%d: Loop block exit.\n", __FUNCTION__, __LINE__);
#endif
    }

    //
    // Push result
    //
    long firstDim = 0;
    for (long idx=0; idx<M; idx++) {
        if (K[idx].size() > firstDim) {
            firstDim = K[idx].size();
        }
    }

    IDL_ArrayCopy(firstDim, K, argv[0]);
    IDL_ArrayCopy(firstDim, Bm, argv[1]);
    IDL_ArrayCopy(3, XYZmeq, argv[2]);
    IDL_VectorCopy(Bmeq, argv[3]);
    IDL_ArrayCopy(3, XYZeq, argv[4]);
    IDL_VectorCopy(Beq, argv[5]);
    IDL_ArrayCopy(3, XYZfoot, argv[6]);
    IDL_ArrayCopy(firstDim, Xfl, argv[7]);
    IDL_ArrayCopy(firstDim, Yfl, argv[8]);
    IDL_ArrayCopy(firstDim, Zfl, argv[9]);

#ifdef DEBUG
    fprintf(stderr, "%%%% DEBUG:%s:%d: %s Enter.\n", __FILE__, __LINE__, __FUNCTION__);
#endif
}

//
//
//
UBK_C_EXTERN void ubk_lstar(int argc, IDL_VPTR *argv,char *argk)
{
#ifdef DEBUG
    fprintf(stderr, "%%%% DEBUG:%s:%d: %s Enter.\n", __FILE__, __LINE__, __FUNCTION__);
#endif

    //
    // Inputs
    //
    vector<double> x0, y0, pa0;
    IDL_VectorCopy(argv[7], x0);
    IDL_VectorCopy(argv[8], y0);
    IDL_VectorCopy(argv[9], pa0);

    vector<IDL_LONG> d;
    IDL_VectorCopy(argv[10], d);

    vector<double> ioptparmod;
    IDL_VectorCopy(argv[11], ioptparmod);

    TSExternalFieldModel external = IDL_ULongScalar(argv[12]);
    GeopackInternalFieldModel internal = IDL_ULongScalar(argv[13]);;

    double ionoR = IDL_DoubleScalar(argv[14]);
    double ds = IDL_DoubleScalar(argv[15]);
    double dx = IDL_DoubleScalar(argv[16]);
    double dy = IDL_DoubleScalar(argv[17]);

    unsigned long n_phi = IDL_ULongScalar(argv[18]);
    unsigned long n_theta = IDL_ULongScalar(argv[19]);

    unsigned long shouldKeepContour = IDL_ULongScalar(argv[20]);
    unsigned long isCartesianGrid = IDL_ULongScalar(argv[21]);

    unsigned long n_threads = IDL_ULongScalar(argv[22]);

    //
    // Minimal validation
    //
    IDL_MEMINT M = x0.size();

    ASSERT(M && y0.size()==M && pa0.size()==M,
           "The length of the input position array is different.");
    ASSERT(5 == d.size(),
           "The number of element of DATE should be 5 ([year, doy, hour, min, sec]).");
    ASSERT(kGeopackDipoleField==internal || kGeopackIGRFField==internal,
           "Invalid internal field model flag.");
    ASSERT((kTSNone==external) ||
           (kTS89Model==external) ||
           (kTS96Model==external) ||
           (kTS02Model==external) ||
           (kTS05Model==external),
           "Invalid external field model flag.");
    ASSERT((kTSNone==external) ||
           (kTS89Model==external && 1==ioptparmod.size()) ||
           (10==ioptparmod.size()),
           "Invalid ioptparmod parameter.");
    ASSERT(ionoR >= 1.,
           "ionoR < 1.");
    ASSERT(ds > 0.,
           "ds <= 0.");
    ASSERT(dx > 0.,
           "dx or dr <= 0.");
    ASSERT(dy > 0.,
           "dy or dphi <= 0.");
    ASSERT(n_phi >= 10,
           "n_phi < 10.");
    ASSERT(n_theta >= 10,
           "n_theta < 10.");
    ASSERT(n_threads >= 1,
           "N_THREADS parameter should be greater than 0.");

    //
    // Log
    //
#ifdef DEBUG
    fprintf(stderr, "%%%% DEBUG:%s:%d: Options\n", __FUNCTION__, __LINE__);
    fprintf(stderr, "\tM = %lld\n", M);
    fprintf(stderr, "\tDATE = [%d, %d, %d, %d, %d]\n", d[0], d[1], d[2], d[3], d[4]);
    fprintf(stderr, "\tINTERNAL = %ld\n", internal);
    fprintf(stderr, "\tEXTERNAL = %ld\n", external);
    fprintf(stderr, "\tIONOR = %f\n", ionoR);
    fprintf(stderr, "\tDS = %f\n", ds);
    fprintf(stderr, "\tDX = %f\n", dx);
    fprintf(stderr, "\tDY = %f\n", dy);
    fprintf(stderr, "\tN_PHI = %ld\n", n_phi);
    fprintf(stderr, "\tN_THETA = %ld\n", n_theta);
    fprintf(stderr, "\tSHOULDKEEPCONTOUR = %d\n", !!shouldKeepContour);
    fprintf(stderr, "\tISCARTESIANGRID = %d\n", !!isCartesianGrid);
    fprintf(stderr, "\tN_THREADS = %ld\n", n_threads);
#endif
    
    //
    // Output buffer
    //
    vector<double> Ls(M), K(M), Phi0(1);
    vector< vector<double> > XorRc(M), YorPhic(M);
    vector< vector<double> > Phif(M), Thetaf(M);

    //
    // Calculation start
    //
    {
#ifdef DEBUG
        fprintf(stderr, "%%%% DEBUG:%s:%d: Loop block enter.\n", __FUNCTION__, __LINE__);
#endif

        class _ : public LstarCoordinator {
        public:
            vector<double> const *x0, *y0, *pa0;
            vector<double> *Ls, *K;
            vector< vector<double> > *XorRc, *YorPhic, *Phif, *Thetaf;
            BOOL shouldKeepContour;

            _(FieldModel const& fm, BOOL isCylindricalGrid) : LstarCoordinator(fm, isCylindricalGrid) {};
            unsigned long count() const {
                return x0->size();
            };
            Point from(long idx) const {
                return Point(x0->at(idx), y0->at(idx), pa0->at(idx));
            };
            void callback(long idx, Particle const& ptl) const {
                Ls->at(idx) = ptl.Lstar();
                K->at(idx) = ptl.K();
                for (long jdx=0, cnt=ptl.coordinates().size(); shouldKeepContour && jdx<cnt; jdx++) {
                    Point const& c = ptl.coordinates()[jdx];
                    XorRc->at(idx).push_back(c.x);
                    YorPhic->at(idx).push_back(c.y);
                    Point const& f = ptl.footPoints()[jdx];
                    Phif->at(idx).push_back(f.phi);
                    Thetaf->at(idx).push_back(f.theta);
                }
            };
        };

        double const vsw = 400.;
        Date date(d[0], d[1], d[2], d[3], d[4]);
        int iopt = (kTS89Model==external ? ioptparmod[0] : 1);
        double const* parmod = (kTS89Model==external && kTSNone==external ? NULL : &ioptparmod.front());
        TSFieldModel fm(date, vsw, internal, iopt, parmod, external);

        _ c(fm, !isCartesianGrid);
        c.setNThreads(n_threads);
        c.setIonoR(ionoR);
        c.setDs(ds);
        c.setD(Point(dx, dy, 1.));
        c.setNPhi(n_phi);
        c.setNTheta(n_theta);
        c.shouldKeepContour = shouldKeepContour;

        c.x0 = &x0;
        c.y0 = &y0;
        c.pa0 = &pa0;
        c.Ls = &Ls;
        c.K = &K;
        c.XorRc = &XorRc;
        c.YorPhic = &YorPhic;
        c.Phif = &Phif;
        c.Thetaf = &Thetaf;

        c.start();

        Phi0[0] = c.Phi0();

#ifdef DEBUG
        fprintf(stderr, "%%%% DEBUG:%s:%d: Loop block exit.\n", __FUNCTION__, __LINE__);
#endif
    }

    //
    // Push result
    //
    long firstDim = 0;
    for (long idx=0; idx<M; idx++) {
        if (XorRc[idx].size() > firstDim) {
            firstDim = XorRc[idx].size();
        }
    }

    IDL_VectorCopy(Ls, argv[0]);
    IDL_VectorCopy(K, argv[1]);
    IDL_VectorCopy(Phi0, argv[2]);
    if (shouldKeepContour) {
        IDL_ArrayCopy(firstDim, XorRc, argv[3]);
        IDL_ArrayCopy(firstDim, YorPhic, argv[4]);
        IDL_ArrayCopy(firstDim, Phif, argv[5]);
        IDL_ArrayCopy(firstDim, Thetaf, argv[6]);
    }

#ifdef DEBUG
    fprintf(stderr, "%%%% DEBUG:%s:%d: %s Enter.\n", __FILE__, __LINE__, __FUNCTION__);
#endif
}