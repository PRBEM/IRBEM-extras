//
//  MagneticFlux.cpp
//  UBJDevelopment
//
//  Created by Kyungguk Min on 4/2/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#include "MagneticFlux.h"
#include "FieldModel.h"
#include "Thread.h"
#include <cmath>
#include <climits>
#include <cassert>
#include <cstdio>
#include <cstring>

namespace UBK {
    using namespace std;

    //
    // Global options.
    //
    unsigned long MagneticFlux::_nPhi = 360/5;
    unsigned long MagneticFlux::_nTheta = 2*180;
    unsigned long MagneticFlux::_Nthreads = 8;

    //
    // Test
    //
#ifdef DEBUG
    void MagneticFlux::test()
    {
        double IMFb1 = 0.;
        double compb2 = 4.;
        CompressedDipole fm(IMFb1, compb2);

        MagneticFlux::setNPhi(360/5);
        MagneticFlux mf(&fm);
        printf("F0 = 2pi %f\n", mf.coreFlux()/2./M_PI);

        Point vec[361];
        double const L0 = 6.;
        for (long i=0; i<361; i++) {
            vec[i].r = 1.;
            vec[i].phi = i*M_PI/180.;
            vec[i].theta = asin(sqrt(1./L0));
        }
        double Ls = mf.fluxClosedByFootPoints_count(vec, 361);
        Ls = mf.coreFlux()/Ls;
        printf("L* = %f\n", Ls);

        assert(0);
    }
#endif

    //
    // Constructor
    //
    MagneticFlux::~MagneticFlux()
    {
    }
    MagneticFlux::MagneticFlux (FieldModel const* fm) : _coreFlux(NAN), _BrA(), _dPhi(0.), _dTheta(0.), _cntPhi(0), _cntTheta(0)
    {
        this->evaluateUsingFieldModel(fm);
    }

    //
    // Flux closed by foot points
    //
    UBK_INLINE double THETA (double phi, double a, double phi1, double theta1)
    {
        return ((a)*((phi)-(phi1)) + (theta1));
    }
    UBK_INLINE double D (double t)
    {
        return (t)-(floor(t));
    }
    UBK_INLINE double dF (double const* _BrA, long _cntPhi, long ip, double t)
    {
        return ( _BrA[((long)(t))*(_cntPhi) + ((ip))] +
                D(t)*(_BrA[(((long)(t))+1)*(_cntPhi) + ((ip))] - _BrA[((long)(t))*(_cntPhi) + ((ip))]) );
    }
    double MagneticFlux::fluxClosedByFootPoints_count(const Point *footPoints, long count) const
    {
        double flux = 0.;

        double const* BrA = &_BrA.front();
        double pi = M_PI / _dPhi;
        double sphi = 0.;
        for (long idx=0; idx<count-1; idx++) {
            double phi1 = footPoints[idx].phi / _dPhi;
            double phi2 = footPoints[idx+1].phi / _dPhi;
            double dphi = phi2 - phi1;
            if (fabs(dphi) < 1e-6) {
                continue;
            } else if (dphi < -pi) {
                phi2 += 2.*pi;
                dphi = phi2 - phi1;
            } else if (dphi > pi) {
                phi1 += 2.*pi;
                dphi = phi2 - phi1;
            }
            sphi += dphi;
            double theta1 = footPoints[idx].theta / _dTheta;
            double theta2 = footPoints[idx+1].theta / _dTheta;
            double a = (theta2 - theta1) / (dphi);
            double thAtHalf;
            long ip, np;
            for (ip=phi1+1, np=phi2; ip<=np; ip++) {
                thAtHalf = THETA(.5*(phi1+ip), a, phi1, theta1);
                //flux += dF(ip, thAtHalf) * (ip-phi1);
                flux += dF(BrA, _cntPhi, ip, thAtHalf) * (ip-phi1);
                phi1 = ip;
            }
            thAtHalf = THETA(.5*(phi1+phi2), a, phi1, theta1);
            //flux += dF(ip, thAtHalf) * (phi2-phi1);
            flux += dF(BrA, _cntPhi, ip, thAtHalf) * (phi2-phi1);
        }

        // Treat the magnetic island
        sphi = fabs(sphi);
        if (sphi < 1.*pi/180.) { // If < 1˚, ignore
            flux = NAN;
        } else if (sphi < pi) { // Check for magnetic island
            flux = -fabs(flux);
        } else {
            flux =  fabs(flux);
        }

        return flux;
    }

    //
    // Private heavy lifter
    //
    class Job {
    public:
        vector<double>& _BrA;
        double _dTheta;
        long _cntPhi;
        double _dPhi;
        FieldModel const* _fm;

        Job(vector<double>& BrA, double dTheta, long cntPhi, double dPhi, FieldModel const* fm) : _BrA(BrA), _dTheta(dTheta), _cntPhi(cntPhi), _dPhi(dPhi), _fm(fm) {};

        void operator()(long wid) {
            long j = wid + 2;
            long n;
            Point bVec;
            double theta = (j-1)*_dTheta;
            for (long i=0; i<_cntPhi; i++) {
                n = j*_cntPhi + i;
                double phi = (i-1)*_dPhi;
                _fm->getSphericalFieldInSM_atSphericalPoint(&bVec, Point(1., phi, theta));
                _BrA[n] = fabs(bVec.r);
            }
        };
    };
    void MagneticFlux::evaluateUsingFieldModel (FieldModel const* fm)
    {
        long nPhiLocal = _nPhi;
        long nThetaLocal = _nTheta;

        {
            _dPhi = 2.*M_PI / nPhiLocal;
            _dTheta = M_PI / nThetaLocal;
            _cntPhi = nPhiLocal + 2; // Two ghost cells
            _cntTheta = nThetaLocal+1 + 1; // One ghost cell

            _BrA.resize(_cntPhi*_cntTheta, 0.);

            // Theta-major. Zero theta starts from index 1. Since sin(th)dth at theta = 0, pi is zero, they are excluded.
            ///////////////////////////////////////////////////
            // Parallel evaluation start
            ///////////////////////////////////////////////////
            Job j(_BrA, _dTheta, _cntPhi, _dPhi, fm);
            ThreadFor<Job &> t(_Nthreads, nThetaLocal-1, j);
            /*{
                class _ : public WorkItem {
                public:
                    vector<double>& _BrA;
                    double _dTheta;
                    long _cntPhi;
                    double _dPhi;
                    FieldModel const* _fm;

                    _(vector<double>& BrA, double dTheta, long cntPhi, double dPhi, FieldModel const* fm) : _BrA(BrA), _dTheta(dTheta), _cntPhi(cntPhi), _dPhi(dPhi), _fm(fm) {};

                    void main() {
                        long j = this->wid() + 2;
                        long n;
                        Point bVec;
                        double theta = (j-1)*_dTheta;
                        for (long i=0; i<_cntPhi; i++) {
                            n = j*_cntPhi + i;
                            double phi = (i-1)*_dPhi;
                            _fm->getSphericalFieldInSM_atSphericalPoint(&bVec, Point(1., phi, theta));
                            _BrA[n] = fabs(bVec.r);
                        }
                    };
                    void done() {delete this;};
                };

                //
                // Create jobs
                //
                std::vector<WorkItem *> jobs(nThetaLocal-1);
                for (long idx=0, n=nThetaLocal-1; idx<n; idx++) jobs[idx] = new _(_BrA, _dTheta, _cntPhi, _dPhi, fm);
                WorkQueue q(jobs, jobs.size());
                q.join();
            }*/
            /*{
                for (long itmp=0, n=nThetaLocal-1; itmp<n; itmp++) {
                    long j = itmp + 2;
                    long wid;
                    Point bVec;
                    double theta = (j-1)*_dTheta;
                    for (long i=0; i<_cntPhi; i++) {
                        wid = j*_cntPhi + i;
                        double phi = (i-1)*_dPhi;
                        fm->getSphericalFieldInSM_atSphericalPoint(&bVec, Point(1., phi, theta));
                        _BrA[wid] = fabs(bVec.r);
                    }
                }
            }*/
            ///////////////////////////////////////////////////
            // Parallel evaluation end
            ///////////////////////////////////////////////////
        }

        // Br.A accumulation
        {
            double dPhi_2 = _dPhi*.5;
            double dTheta_2 = _dTheta*.5;
            double increment;
            long wid1, wid0;
            for (long j=2; j<_cntTheta; j++) {
                increment = sin((j-1)*_dTheta)*dPhi_2*dTheta_2;
                for (wid0=j*_cntPhi+1,wid1=wid0+1; wid1<(j+1)*_cntPhi; wid0=wid1++) {
                    // Phi integration
                    _BrA[wid0] = (_BrA[wid0]+_BrA[wid1])*increment;
                    // Theta integration
                    _BrA[wid0-_cntPhi] = _BrA[wid0-2*_cntPhi] + (_BrA[wid0-_cntPhi]+_BrA[wid0]);
                }
                // Fill ghost cells
                _BrA[wid0] = _BrA[wid1-_cntPhi+1];
                _BrA[wid1-_cntPhi] = _BrA[wid0-1];
                _BrA[wid0-_cntPhi] = _BrA[wid1-2*_cntPhi+1];
                _BrA[wid1-2*_cntPhi] = _BrA[wid0-1-_cntPhi];
            }
        }

        // Core flux
        {
            _coreFlux = 0.;
            for (long wid=nThetaLocal*_cntPhi+2; wid<(nThetaLocal+1)*_cntPhi; wid++) {
                _coreFlux += _BrA[wid];
            }
            _coreFlux *= .5;
        }
    }

}