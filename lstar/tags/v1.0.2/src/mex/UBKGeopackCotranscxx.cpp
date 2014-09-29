//
//  UBKGeopackCotranscxx.cpp
//  UBJDevelopment
//
//  Created by Kyungguk Min on 4/10/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#include "UBKGeopackCotranscxx.h"
#include "UBKLstarxx.h"
#include "matrix.h"
#include "mex.h"
#include <string>

using namespace UBK;
using namespace std;

//
// Log facility
//
static Mutex log_key;

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
    if (!test) {Locker _(log_key);
        mexErrMsgIdAndTxt(msgid.id.c_str(), msg);
    }
}

//
// Sub iteration
//
class Submain : public CotransCoordinator {
public:
    // Inputs
    unsigned long M;
    double const* xin;
    double const* yin;
    double const* zin;
    // Outputs
    double *xout;
    double *yout;
    double *zout;

    Submain(FieldModel const& fm) : CotransCoordinator(fm) {};

    unsigned long count() const {
        return M;
    };
    Point from(long idx) const {
        return Point(xin[idx], yin[idx], zin[idx]);
    };
    void callback(long idx, Point const& pt) const {
        xout[idx] = pt.x;
        yout[idx] = pt.y;
        zout[idx] = pt.z;
    };
};

//
// Main
//
class Main {
    // Inputs
    unsigned long M;
    unsigned long N;
    double const* xin;
    double const* yin;
    double const* zin;
    double const* date;
    long to_co_system;
    long M_threads;
    long N_threads;
    // Outputs
    double *xout;
    double *yout;
    double *zout;
    double *psi;

public:
    Main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        //
        // Check for nargin and nargout
        //
        ASSERT(4==nlhs && 7==nrhs, "Wrong number of input/output.");

        //
        // (xin [np, nt], yin [np, nt], zin [np, nt], [year, doy, hour, min, sec] (5, nt), to_co_system, M_threads, N_threads, shouldUseInterpolatedField)
        //
        M = mxGetM(prhs[0]);
        N = mxGetN(prhs[3]);

        xin = mxGetPr(prhs[0]);
        yin = mxGetPr(prhs[1]);
        zin = mxGetPr(prhs[2]);
        date = mxGetPr(prhs[3]);
        to_co_system = round( mxGetScalar(prhs[4]) );
        M_threads = round( mxGetScalar(prhs[5]) );
        N_threads = round( mxGetScalar(prhs[6]) );

        //
        // Validity
        //
        ASSERT(5==mxGetM(prhs[3]), "size(data,1) != 5.");
        ASSERT((kMagneticFieldSM==to_co_system || kMagneticFieldGSM==to_co_system),
               "Invalid to_co_system.");
        ASSERT(M_threads>0, "M_threads <= 0.");
        ASSERT(N_threads>0, "N_threads <= 0.");

        //
        // Output buffer
        //
        {
            mxArray* xout_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* yout_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* zout_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray *psi_ = mxCreateDoubleMatrix(1, N, mxREAL);

            xout = mxGetPr(xout_);
            yout = mxGetPr(yout_);
            zout = mxGetPr(zout_);
            psi = mxGetPr(psi_);

            //
            // LHS
            //
            plhs[0] = xout_;
            plhs[1] = yout_;
            plhs[2] = zout_;
            plhs[3] = psi_;
        }

#ifdef DEBUG
        {Locker _(log_key);
            mexPrintf("%%%% DEBUG:%s:%d:\n", __FUNCTION__, __LINE__);
            mexPrintf("\t[M, N] = [%ld, %ld]\n", M, N);
            mexPrintf("\tto_co_system = %ld\n", to_co_system);
            mexPrintf("\tM_threads = %ld\n", M_threads);
            mexPrintf("\tN_threads = %ld\n", N_threads);
        }
#endif

        //
        // Time iteration
        //
        if (N) {
#ifdef DEBUG
            {Locker _(log_key);
                mexPrintf("%%%% DEBUG:%s:%d: Outer loop start.\n", __FUNCTION__, __LINE__);
            }
#endif

            ThreadFor<Main &> t(N_threads, N, *this);

#ifdef DEBUG
            {Locker _(log_key);
                mexPrintf("%%%% DEBUG:%s:%d: Outer loop end.\n", __FUNCTION__, __LINE__);
            }
#endif
        }
    };

    void operator()(long jdx) {
        const double vsw = 400.;
        long d_offset = jdx * 5;

        //
        // Make TS field model
        //
        Date d(date[d_offset + 0], date[d_offset + 1], date[d_offset + 2], date[d_offset + 3], date[d_offset + 4]);
        TSFieldModel fm(d, vsw, kGeopackDipoleField, 1, NULL, kTSNone);
        psi[jdx] = fm.dipoleTiltAngleRadian();

#ifdef DEBUG
        {Locker _(log_key);
            mexPrintf("%%%% DEBUG:%s:%d: date = [%d, %d, %d, %d, %d].\n", __FUNCTION__, __LINE__, d.year, d.doy, d.hour, d.min, d.sec);
        }
#endif

        //
        // Sub iteration
        //
        if (M) {
#ifdef DEBUG
            {Locker _(log_key);
                mexPrintf("%%%% DEBUG:%s:%d: Inner loop start.\n", __FUNCTION__, __LINE__);
            }
#endif

            Submain sub(fm);
            sub.setNThreads(M_threads);
            sub.setToCoSystem(to_co_system);
            sub.M = M;

            long sub_offset = jdx * M;

            sub.xin = xin + sub_offset;
            sub.yin = yin + sub_offset;
            sub.zin = zin + sub_offset;
            sub.xout = xout + sub_offset;
            sub.yout = yout + sub_offset;
            sub.zout = zout + sub_offset;

            sub.start();

#ifdef DEBUG
            {Locker _(log_key);
                mexPrintf("%%%% DEBUG:%s:%d: Inner loop end.\n", __FUNCTION__, __LINE__);
            }
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
