//
//  UBKFieldLinecxx.cpp
//  UBJDevelopment
//
//  Created by Kyungguk Min on 4/7/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#include "UBKFieldLinecxx.h"
#include "UBKLstarxx.h"
#include "matrix.h"
#include "mex.h"
#include <string>

using namespace UBK;
using namespace std;

//
// Log facility
//
class MsgId {
public:
    string id;

    MsgId() {
        id = __FILE__;
        size_t pos = 0;
        while ( (pos=id.find(".", pos))!=string::npos ) {
            id.replace(pos, 1, "_");
        }
        id += ":AssertFailure";
    };
};
static MsgId msgid;

UBK_INLINE void ASSERT(BOOL test, const char *msg) {
    if (!test) {
        mexErrMsgIdAndTxt(msgid.id.c_str(), msg);
    }
}

//
// Sub iteration
//
class Submain : public FieldLineCoordinator {
public:
    // This is needed because insertion to cell array is not thread-safe.
    Mutex *key;

    // Inputs
    long jdx;
    unsigned long M;
    double const* x0;
    double const* y0;
    double const* z0;
    // Outputs
    mxArray *K;
    mxArray *Bm;
    double *Xmeq;
    double *Ymeq;
    double *Zmeq;
    double *Bmeq;
    double *Xeq;
    double *Yeq;
    double *Zeq;
    double *Beq;
    double *Xfoot;
    double *Yfoot;
    double *Zfoot;
    mxArray *Xfl;
    mxArray *Yfl;
    mxArray *Zfl;

    Submain(FieldModel const& fm) : FieldLineCoordinator(fm) {};

    unsigned long count() const {
        return M;
    };
    Point from(long idx) const {
        return Point(x0[idx], y0[idx], z0[idx]);
    };
    void callback(long idx, FieldLine const& fl) const {
        Bmeq[idx] = fl.magneticEquatorFieldMagnitude();
        Xmeq[idx] = fl.magneticEquator().x;
        Ymeq[idx] = fl.magneticEquator().y;
        Zmeq[idx] = fl.magneticEquator().z;
        Beq[idx] = fl.coordinateEquatorFieldMagnitude();
        Xeq[idx] = fl.coordinateEquator().x;
        Yeq[idx] = fl.coordinateEquator().y;
        Zeq[idx] = fl.coordinateEquator().z;
        Xfoot[idx] = fl.footPoint().x;
        Yfoot[idx] = fl.footPoint().y;
        Zfoot[idx] = fl.footPoint().z;

        long cnt = fl.modifiedInvariants().size();
        long wid = jdx*M + idx;
        if (cnt) {
            double *pxfl, *pyfl, *pzfl, *pK, *pBm;
            {
                Locker l(*key); // Prevent seg. fault while inserting to cell

                // K & Bm
                mxSetCell(K, wid, mxCreateDoubleMatrix(cnt, 1, mxREAL));
                pK = mxGetPr(mxGetCell(K, wid));
                mxSetCell(Bm, wid, mxCreateDoubleMatrix(cnt, 1, mxREAL));
                pBm = mxGetPr(mxGetCell(Bm, wid));

                // FL
                mxSetCell(Xfl, wid, mxCreateDoubleMatrix(cnt, 1, mxREAL));
                pxfl = mxGetPr(mxGetCell(Xfl, wid));
                mxSetCell(Yfl, wid, mxCreateDoubleMatrix(cnt, 1, mxREAL));
                pyfl = mxGetPr(mxGetCell(Yfl, wid));
                mxSetCell(Zfl, wid, mxCreateDoubleMatrix(cnt, 1, mxREAL));
                pzfl = mxGetPr(mxGetCell(Zfl, wid));
            }

            // K
            copy(fl.modifiedInvariants().begin(), fl.modifiedInvariants().end(), pK);

            // Bm
            copy(fl.mirrorMagnitudes().begin(), fl.mirrorMagnitudes().end(), pBm);

            // FL
            for (long it=0; it<cnt; it++) {
                Point const &pt = fl.coordinates().at(it);
                *pxfl++ = pt.x;
                *pyfl++ = pt.y;
                *pzfl++ = pt.z;
            }
        }
    };
};

//
// Main
//
class Main {
    // Share with child workers.
    Mutex key;

    // Inputs
    unsigned long M;
    unsigned long N;
    double const* x0;
    double const* y0;
    double const* z0;
    double const* date;
    GeopackInternalFieldModel internal;
    double const* ioptparmod;
    TSExternalFieldModel external;
    double ionoR;
    double ds;
    long M_threads;
    long N_threads;
    // Outputs
    mxArray *K;
    mxArray *Bm;
    double *Xmeq;
    double *Ymeq;
    double *Zmeq;
    double *Bmeq;
    double *Xeq;
    double *Yeq;
    double *Zeq;
    double *Beq;
    double *Xfoot;
    double *Yfoot;
    double *Zfoot;
    mxArray *Xfl;
    mxArray *Yfl;
    mxArray *Zfl;

public:
    Main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        //
        // Check for nargin and nargout
        //
        ASSERT(16==nlhs && 11==nrhs, "Wrong number of input/output.");

        //
        // (xin [np, nt], yin [np, nt], zin [np, nt], [year, doy, hour, min, sec] (5, nt), ioptparmod [1 or 10, nt], external, internal, ionoR, ds, M_threads, N_threads)
        //
        M = mxGetM(prhs[0]);
        N = mxGetN(prhs[3]);

        x0 = mxGetPr(prhs[0]);
        y0 = mxGetPr(prhs[1]);
        z0 = mxGetPr(prhs[2]);
        date = mxGetPr(prhs[3]);
        ioptparmod = mxGetPr(prhs[4]);
        external = round( mxGetScalar(prhs[5]) );
        internal = round( mxGetScalar(prhs[6]) );
        ionoR = mxGetScalar(prhs[7]);
        ds = mxGetScalar(prhs[8]);
        M_threads = round( mxGetScalar(prhs[9]) );
        N_threads = round( mxGetScalar(prhs[10]) );

        //
        // Validity
        //
        ASSERT(5==mxGetM(prhs[3]), "size(data,1) != 5.");
        ASSERT((kGeopackDipoleField==internal) ||
               (kGeopackIGRFField==internal), "Invalid internal magnetic field component.");
        ASSERT((kTSNone==external) ||
               (kTS89Model==external) ||
               (kTS96Model==external) ||
               (kTS02Model==external) ||
               (kTS05Model==external), "Invalid external magnetic field component.");
        ASSERT((external==kTSNone) ||
               (external==kTS89Model && 1==mxGetM(prhs[4])) ||
               (10==mxGetM(prhs[4])), "Invalid ioptparmod dimension.");
        ASSERT(ionoR >= 1., "ionoR < 1.");
        ASSERT(ds > 0., "ds <= 0.");
        ASSERT(M_threads > 0, "M_threads <= 0.");
        ASSERT(N_threads > 0, "N_threads <= 0.");

        //
        // Output buffer
        //
        {
            K = mxCreateCellMatrix(M, N);
            Bm = mxCreateCellMatrix(M, N);
            mxArray* Xmeq_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* Ymeq_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* Zmeq_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* Bmeq_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* Xeq_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* Yeq_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* Zeq_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* Beq_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* Xfoot_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* Yfoot_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* Zfoot_ = mxCreateDoubleMatrix(M, N, mxREAL);
            Xfl = mxCreateCellMatrix(M, N);
            Yfl = mxCreateCellMatrix(M, N);
            Zfl = mxCreateCellMatrix(M, N);

            Xmeq = mxGetPr(Xmeq_);
            Ymeq = mxGetPr(Ymeq_);
            Zmeq = mxGetPr(Zmeq_);
            Bmeq = mxGetPr(Bmeq_);
            Xeq = mxGetPr(Xeq_);
            Yeq = mxGetPr(Yeq_);
            Zeq = mxGetPr(Zeq_);
            Beq = mxGetPr(Beq_);
            Xfoot = mxGetPr(Xfoot_);
            Yfoot = mxGetPr(Yfoot_);
            Zfoot = mxGetPr(Zfoot_);

            //
            // LHS
            //
            plhs[0] = K;
            plhs[1] = Bm;
            plhs[2] = Xmeq_;
            plhs[3] = Ymeq_;
            plhs[4] = Zmeq_;
            plhs[5] = Bmeq_;
            plhs[6] = Xeq_;
            plhs[7] = Yeq_;
            plhs[8] = Zeq_;
            plhs[9] = Beq_;
            plhs[10] = Xfoot_;
            plhs[11] = Yfoot_;
            plhs[12] = Zfoot_;
            plhs[13] = Xfl;
            plhs[14] = Yfl;
            plhs[15] = Zfl;
        }

        //
        // Log
        //
#ifdef DEBUG
        mexPrintf("%%%% DEBUG:%s:%d:\n", __FUNCTION__, __LINE__);
        mexPrintf("\t[M, N] = [%ld, %ld]\n", M, N);
        mexPrintf("\tinternal = %ld\n", internal);
        mexPrintf("\texternal = %ld\n", external);
        mexPrintf("\tionoR = %f\n", ionoR);
        mexPrintf("\tds = %f\n", ds);
        mexPrintf("\tM_threads = %ld\n", M_threads);
        mexPrintf("\tN_threads = %ld\n", N_threads);
#endif

        //
        // Time iteration
        //
        if ( M*N ) {
#ifdef DEBUG
            mexPrintf("%%%% DEBUG:%s:%d: Outer loop start.\n", __FUNCTION__, __LINE__);
#endif

            ThreadFor<Main &> t(N_threads, N, *this);

#ifdef DEBUG
            mexPrintf("%%%% DEBUG:%s:%d: Outer loop end.\n", __FUNCTION__, __LINE__);
#endif
        }
    };

    void operator()(long jdx) {
        const double vsw = 400.;
        long d_offset = jdx * 5;
        long iopt_offset = jdx * 10;

        //
        // Make TS field model
        //
        int iopt = (kTS89Model==external ? round(ioptparmod[jdx]) : 1);
        double const* parmod = ioptparmod + iopt_offset;
        Date d(date[d_offset + 0], date[d_offset + 1], date[d_offset + 2], date[d_offset + 3], date[d_offset + 4]);
        TSFieldModel fm(d, vsw, internal, iopt, parmod, external);

#ifdef DEBUG
        mexPrintf("%%%% DEBUG:%s:%d: date = [%d, %d, %d, %d, %d].\n", __FUNCTION__, __LINE__, d.year, d.doy, d.hour, d.min, d.sec);
#endif

        //
        // Sub iteration
        //
        {
#ifdef DEBUG
            mexPrintf("%%%% DEBUG:%s:%d: Inner loop start.\n", __FUNCTION__, __LINE__);
#endif

            Submain sub(fm);
            sub.key = &key;

            sub.setNThreads(M_threads);
            sub.setDs(ds);
            sub.setIonoR(ionoR);
            sub.M = M;
            sub.jdx = jdx;

            sub.K = K;
            sub.Bm = Bm;
            sub.Xfl = Xfl;
            sub.Yfl = Yfl;
            sub.Zfl = Zfl;

            long sub_offset = jdx * M;

            sub.x0 = x0 + sub_offset;
            sub.y0 = y0 + sub_offset;
            sub.z0 = z0 + sub_offset;

            sub.Bmeq = Bmeq + sub_offset;
            sub.Xmeq = Xmeq + sub_offset;
            sub.Ymeq = Ymeq + sub_offset;
            sub.Zmeq = Zmeq + sub_offset;
            sub.Beq = Beq + sub_offset;
            sub.Xeq = Xeq + sub_offset;
            sub.Yeq = Yeq + sub_offset;
            sub.Zeq = Zeq + sub_offset;
            sub.Xfoot = Xfoot + sub_offset;
            sub.Yfoot = Yfoot + sub_offset;
            sub.Zfoot = Zfoot + sub_offset;

            sub.start();

#ifdef DEBUG
            mexPrintf("%%%% DEBUG:%s:%d: Inner loop end.\n", __FUNCTION__, __LINE__);
#endif
        }
    };
};

//
// MEX entry
//
void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
#ifdef DEBUG
    mexPrintf("%%%% DEBUG:%s:%d: %s Enter.\n", __FILE__, __LINE__, __FUNCTION__);
#endif

    Main main(nlhs, plhs, nrhs, prhs);

#ifdef DEBUG
    mexPrintf("%%%% DEBUG:%s:%d: %s Exit.\n", __FILE__, __LINE__, __FUNCTION__);
#endif
}
