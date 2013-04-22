//
//  UBKLstarcxx.cpp
//  UBJDevelopment
//
//  Created by Kyungguk Min on 4/20/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#include "UBKLstarcxx.h"
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
static const MsgId msgid;

UBK_INLINE void ASSERT(BOOL test, const char *msg) {
    if (!test) {
        mexErrMsgIdAndTxt(msgid.id.c_str(), msg);
    }
};


//
// Sub iteration
//
class Submain : public LstarCoordinator {
public:
    Mutex *key; // This is needed because insertion to cell array is not thread-safe.

    long jdx;
    unsigned long M;
    double const* x0;
    double const* y0;
    double const* pa0;
    BOOL shouldKeepContour;
    double *Ls;
    double *K;
    mxArray *XorR;
    mxArray *YorPhi;
    mxArray *Phif;
    mxArray *Thetaf;

    Submain(FieldModel const& fm, BOOL isCylindricalGrid) : LstarCoordinator(fm, isCylindricalGrid) {};

    unsigned long count() const {
        return M;
    };
    Point from(long idx) const {
        return Point(x0[idx], y0[idx], pa0[idx]);
    };
    void callback(long idx, Particle const& ptl) const {
        Ls[idx] = ptl.Lstar();
        K[idx] = ptl.K();

        long cnt = ptl.coordinates().size();
        long wid = jdx*M + idx;
        if (shouldKeepContour && cnt) {
            double *pp, *pt, *px, *py;

            {
                Locker l(*key); // Prevent seg. fault while inserting to cell

                mxSetCell(XorR, wid, mxCreateDoubleMatrix(cnt, 1, mxREAL));
                mxSetCell(YorPhi, wid, mxCreateDoubleMatrix(cnt, 1, mxREAL));
                mxSetCell(Phif, wid, mxCreateDoubleMatrix(cnt, 1, mxREAL));
                mxSetCell(Thetaf, wid, mxCreateDoubleMatrix(cnt, 1, mxREAL));

                pp = mxGetPr(mxGetCell(Phif, wid));
                pt = mxGetPr(mxGetCell(Thetaf, wid));
                px = mxGetPr(mxGetCell(XorR, wid));
                py = mxGetPr(mxGetCell(YorPhi, wid));
            }
            for (long it=0; it<cnt; it++) {
                Point const &c = ptl.coordinates().at(it);
                *px++ = c.x;
                *py++ = c.y;
                Point const &f = ptl.footPoints().at(it);
                *pp++ = f.phi;
                *pt++ = f.theta;
            }
        }
    };
};

//
// Main
//
class Main {
    Mutex key; // Share with child workers.

    unsigned long M;
    unsigned long N;
    double const* x0;
    double const* y0;
    double const* pa0;
    double const* date;
    GeopackInternalFieldModel internal;
    double const* ioptparmod;
    TSExternalFieldModel external;
    double ionoR;
    double ds;
    double dx;
    double dy;
    long n_phi;
    long n_theta;
    BOOL shouldKeepContour;
    BOOL isCartesianGrid;
    long n_threads;
    double *Ls;
    double *K;
    double *Phi0;
    mxArray *XorR;
    mxArray *YorPhi;
    mxArray *Phif;
    mxArray *Thetaf;

public:
    Main(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
        //
        // Check for nargin and nargout
        //
        ASSERT(7==nlhs && 16==nrhs, "Wrong number of input/output.");

        //
        // (xin [np, nt], yin [np, nt], zin [np, nt], [year, doy, hour, min, sec] (5, nt), ioptparmod [1 or 10, nt], external, internal, ionoR, ds, dx, dy, n_phi, n_theta, shouldKeepContour, isCartesian, n_threads)
        //
        M = mxGetM(prhs[0]);
        N = mxGetN(prhs[3]);

        x0 = mxGetPr(prhs[0]);
        y0 = mxGetPr(prhs[1]);
        pa0 = mxGetPr(prhs[2]);
        date = mxGetPr(prhs[3]);
        ioptparmod = mxGetPr(prhs[4]);
        external = round( mxGetScalar(prhs[5]) );
        internal = round( mxGetScalar(prhs[6]) );
        ionoR = mxGetScalar(prhs[7]);
        ds = mxGetScalar(prhs[8]);
        dx = mxGetScalar(prhs[9]);
        dy = mxGetScalar(prhs[10]);
        n_phi = round( mxGetScalar(prhs[11]) );
        n_theta = round( mxGetScalar(prhs[12]) );
        shouldKeepContour = *mxGetLogicals(prhs[13]);
        isCartesianGrid = *mxGetLogicals(prhs[14]);
        n_threads = round( mxGetScalar(prhs[15]) );

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
        ASSERT(dx > 0., "dx <= 0.");
        ASSERT(dy > 0., "dy <= 0.");
        ASSERT(n_phi > 10, "n_phi <= 10.");
        ASSERT(n_theta > 10, "n_theta <= 10.");
        ASSERT(n_threads > 0, "n_threads <= 0.");

        //
        // Output buffer
        //
        {
            mxArray* Ls_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* K_ = mxCreateDoubleMatrix(M, N, mxREAL);
            mxArray* Phi0_ = mxCreateDoubleMatrix(1, N, mxREAL);
            XorR = mxCreateCellMatrix(M, N);
            YorPhi = mxCreateCellMatrix(M, N);
            Phif = mxCreateCellMatrix(M, N);
            Thetaf = mxCreateCellMatrix(M, N);

            Ls = mxGetPr(Ls_);
            K = mxGetPr(K_);
            Phi0 = mxGetPr(Phi0_);

            //
            // LHS
            //
            plhs[0] = Ls_;
            plhs[1] = K_;
            plhs[2] = Phi0_;
            plhs[3] = XorR;
            plhs[4] = YorPhi;
            plhs[5] = Phif;
            plhs[6] = Thetaf;
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
        mexPrintf("\tdx = %f\n", dx);
        mexPrintf("\tdy = %f\n", dy);
        mexPrintf("\tn_phi = %ld\n", n_phi);
        mexPrintf("\tn_theta = %ld\n", n_theta);
        mexPrintf("\tshouldKeepContour = %d\n", !!shouldKeepContour);
        mexPrintf("\tisCartesianGrid = %d\n", !!isCartesianGrid);
        mexPrintf("\tn_threads = %ld\n", n_threads);
#endif

        //
        // Time iteration
        //
        if (N) {
#ifdef DEBUG
            mexPrintf("%%%% DEBUG:%s:%d: Outer loop start.\n", __FUNCTION__, __LINE__);
#endif

            ThreadFor<Main &> t((M<10 ? n_threads : 1), N, *this);

#ifdef DEBUG
            mexPrintf("%%%% DEBUG:%s:%d: Outer loop end.\n", __FUNCTION__, __LINE__);
#endif
        }
    };

    void operator()(long jdx) {
        const double vsw = 400.;
        long offset = jdx * 5;

        //
        // Make TS field model
        //
        int iopt = (kTS89Model==external ? round(ioptparmod[jdx]) : 1);
        double const* parmod = (kTSNone==external || kTS89Model==external ? NULL : ioptparmod + jdx*10);
        Date d(date[offset + 0], date[offset + 1], date[offset + 2], date[offset + 3], date[offset + 4]);
        TSFieldModel fm(d, vsw, internal, iopt, parmod, external);

#ifdef DEBUG
        mexPrintf("%%%% DEBUG:%s:%d: date = [%d, %d, %d, %d, %d].\n", __FUNCTION__, __LINE__, d.year, d.doy, d.hour, d.min, d.sec);
#endif

        //
        // Sub iteration
        //
        if (M) {
#ifdef DEBUG
            mexPrintf("%%%% DEBUG:%s:%d: Inner loop start.\n", __FUNCTION__, __LINE__);
#endif

            Submain sub(fm, !isCartesianGrid);
            sub.key = &key;

            sub.setNThreads(n_threads);
            sub.setIonoR(ionoR);
            sub.setDs(ds);
            sub.setNPhi(n_phi);
            sub.setNTheta(n_theta);
            sub.setD(Point(dx, dy, 1.));
            sub.shouldKeepContour = shouldKeepContour;
            sub.M = M;
            sub.jdx = jdx;

            sub.XorR = XorR;
            sub.YorPhi = YorPhi;
            sub.Phif = Phif;
            sub.Thetaf = Thetaf;

            offset = jdx * M;

            sub.x0 = x0 + offset;
            sub.y0 = y0 + offset;
            sub.pa0 = pa0 + offset;

            sub.Ls = Ls + offset;
            sub.K = K + offset;

            sub.start();

            Phi0[jdx] = sub.Phi0();

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
