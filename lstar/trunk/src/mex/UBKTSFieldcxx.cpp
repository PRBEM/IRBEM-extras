//
//  UBKTSFieldcxx.cpp
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

#include "UBKTSFieldcxx.h"
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
class Submain : public BFieldCoordinator {
public:
    // Inputs
    unsigned long M;
    double const* x;
    double const* y;
    double const* z;
    // Outputs
    double *Bx;
    double *By;
    double *Bz;

    Submain(FieldModel const& fm) : BFieldCoordinator(fm) {};

    unsigned long count() const {
        return M;
    };
    Point at(long idx) const {
        return Point(x[idx], y[idx], z[idx]);
    };
    void callback(long idx, Point const& B) const {
        Bx[idx] = B.x;
        By[idx] = B.y;
        Bz[idx] = B.z;
    };
};

//
// Main
//
class Main {
    // Inputs
    unsigned long M;
    unsigned long N;
    double const* x;
    double const* y;
    double const* z;
    double const* date;
    GeopackInternalFieldModel internal;
    double const* ioptparmod;
    TSExternalFieldModel external;
    long co_system;
    long M_threads;
    long N_threads;
    // Outputs
    double *Bx;
    double *By;
    double *Bz;

public:
    Main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        //
        // Check for nargin and nargout
        //
        ASSERT(3==nlhs && 10==nrhs, "Wrong number of input/output.");

        //
        // (xin [np, nt], yin [np, nt], zin [np, nt], [year, doy, hour, min, sec] (5, nt), ioptparmod [1 or 10, nt], external, internal, co_system, M_threads, N_threads)
        //
        M = mxGetM(prhs[0]);
        N = mxGetN(prhs[3]);

        x = mxGetPr(prhs[0]);
        y = mxGetPr(prhs[1]);
        z = mxGetPr(prhs[2]);
        date = mxGetPr(prhs[3]);
        ioptparmod = mxGetPr(prhs[4]);
        external = round( mxGetScalar(prhs[5]) );
        internal = round( mxGetScalar(prhs[6]) );
        co_system = round( mxGetScalar(prhs[7]) );
        M_threads = round( mxGetScalar(prhs[8]) );
        N_threads = round( mxGetScalar(prhs[9]) );

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
        ASSERT((kMagneticFieldSM==co_system || kMagneticFieldGSM==co_system), "Invalid to_co_system.");
        ASSERT(M_threads > 0, "M_threads <= 0.");
        ASSERT(N_threads > 0, "N_threads <= 0.");

        //
        // Output buffer
        //
        {
            mxArray* Bx_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* By_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* Bz_ = mxCreateDoubleMatrix(M, N, mxREAL);

            Bx = mxGetPr(Bx_);
            By = mxGetPr(By_);
            Bz = mxGetPr(Bz_);

            //
            // LHS
            //
            plhs[0] = Bx_;
            plhs[1] = By_;
            plhs[2] = Bz_;
        }

#ifdef DEBUG
        mexPrintf("%%%% DEBUG:%s:%d:\n", __FUNCTION__, __LINE__);
        mexPrintf("\t[M, N] = [%ld, %ld]\n", M, N);
        mexPrintf("\tinternal = %ld\n", internal);
        mexPrintf("\texternal = %ld\n", external);
        mexPrintf("\tco_system = %ld\n", co_system);
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
            sub.setNThreads(M_threads);
            sub.setCoSystem(co_system);
            sub.M = M;

            long sub_offset = jdx * M;

            sub.x = x + sub_offset;
            sub.y = y + sub_offset;
            sub.z = z + sub_offset;
            sub.Bx = Bx + sub_offset;
            sub.By = By + sub_offset;
            sub.Bz = Bz + sub_offset;
            
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
