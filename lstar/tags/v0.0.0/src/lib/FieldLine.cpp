//
//  FieldLine.cpp
//  UBJDevelopment
//
//  Created by Kyungguk Min on 4/1/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#include "FieldLine.h"
#include "FieldModel.h"
#include <cmath>
#include <climits>
#include <cstdio>
#include <cassert>
#include <algorithm>

namespace UBK {
    using namespace std;

    //
    // Test
    //
#ifdef DEBUG
    void FieldLine::test()
    {
        double IMFb1 = 0.;
        double compb2 = 4.;
        CompressedDipole fm(IMFb1, compb2);

        FieldLine::setStepSize(.05);
        FieldLine::setRadiusMin(2.);
        Point origin(6., 0., 0.);
        FieldLine fl(origin, &fm);

        printf("%ld\n", fl.mirrorMagnitudes().size());

        assert(0);
    }
#endif

    // Bisection index locator.
    UBK_INLINE long locate (double x, vector<double> const& abscissa) {
        long n = abscissa.size();
        long ju, jm, jl;

        jl = 0;
        ju = n-1;
        while (ju-jl > 1) {
            jm = (ju+jl) >> 1;
            if (x-abscissa[jm]>=0.) {
                jl = jm;
            } else {
                ju = jm;
            }
        }
        return jl;
    }

    /*! Global options
     */
    double FieldLine::_stepSize = 0.05;
    double FieldLine::_radiusMax = 15.;
    double FieldLine::_radiusMin = 1.; // Ionospheric boundary
    long FieldLine::_maxSteps = 10000;

    //
    // Constructor
    //
    FieldLine::FieldLine (Point const origin,
                          FieldModel const* fm,
                          IntegrateODESignature odeSolver) :
    _origin(origin),
    _coordinates(),
    _mirrorMagnitudes(),
    _modifiedInvariants(),
    _coordinateEquator(NAN, NAN),
    _coordinateEquatorFieldMagnitude(NAN),
    _magneticEquator(NAN, NAN, NAN),
    _magneticEquatorFieldMagnitude(NAN),
    _footPoint(NAN, NAN, NAN),
    _valid(NO)
    {
        assert(NULL != fm);
        assert(NULL != odeSolver);

        this->traceCartesianFieldLineUsingFieldBlock(fm, odeSolver);
        this->evaluateInvariant();

        if (!_valid) {
            _coordinates.clear();
            _mirrorMagnitudes.clear();
            _modifiedInvariants.clear();
        }
    }

    //
    // Public methods
    //
    double FieldLine::mirrorMagnitudeAtInvariant(const double K) const
    {
        if (!_valid) {
            return NAN;
        }

        double bm;
        // K validity
        if (K<0. || K-_modifiedInvariants.back()>0.) {
            bm = NAN;
        } else {
            // Linear interpolation
            long j = locate(K, _modifiedInvariants);
            bm = _mirrorMagnitudes[j] + (_mirrorMagnitudes[j+1]-_mirrorMagnitudes[j])
            * ( (K-_modifiedInvariants[j])/(_modifiedInvariants[j+1]-_modifiedInvariants[j]) );
            // Check for bifurcation
            if (bm-_magneticEquatorFieldMagnitude<0.) {
                bm = NAN;
            }
        }
        return bm;
    }

    double FieldLine::invariantAtMirrorMagnitude(const double bm) const
    {
        if (!_valid) {
            return NAN;
        }

        double K;
        // Validity
        if (bm-_mirrorMagnitudes.back()>0.) {
            K = NAN;
        } else {
            // Linear interpolation
            long j;
            if (bm-_mirrorMagnitudes[0]>0.) { // If Bm > B(0)
                // Check for bifurcation
                if (bm-_magneticEquatorFieldMagnitude<0.) {
                    K = NAN;
                } else {
                    j = locate(bm, _mirrorMagnitudes);
                    K = _modifiedInvariants[j] + (_modifiedInvariants[j+1]-_modifiedInvariants[j])
                    * ( (bm-_mirrorMagnitudes[j])/(_mirrorMagnitudes[j+1]-_mirrorMagnitudes[j]) );
                }
            } else { // Bm < B(0),Beq. Return negative (extrapolated) K
                j = 0;
                K = _modifiedInvariants[j] + (_modifiedInvariants[j+1]-_modifiedInvariants[j])
                * ( (bm-_mirrorMagnitudes[j])/(_mirrorMagnitudes[j+1]-_mirrorMagnitudes[j]) );
            }
        }
        
        return K;
    }

    double FieldLine::mirrorMagnitudeOfPitchAngle (double const paRad) const
    {
        if (!_valid) {
            return NAN;
        }

        double sinpa = sin(paRad);
        return _magneticEquatorFieldMagnitude / (sinpa*sinpa);
    }
    double FieldLine::pitchAngleOfmirrorMagnitude (double const bm) const
    {
        if (!_valid) {
            return NAN;
        }

        if (bm<_magneticEquatorFieldMagnitude) {
            return NAN;
        }

        double sinpa2 = _magneticEquatorFieldMagnitude / bm;
        return asin(sqrt(sinpa2));
    }

    //
    // Heavy lifter
    //
    struct _TMP {
        FieldModel const* fm;
        Point coeqtemp;
        double r0;
    };
    void FieldLine::derivatives (double x, double const y[], long count, double *dydx)
    {
        struct _TMP *tmp = (struct _TMP *)(_tmp);
        
        // Calculate field
        Point pt(y[0], y[1], y[2]);
        Point b;
        tmp->fm->getCartesianFieldInSM_atCartesianPoint(&b, pt);// fieldBlock(*(ULPoint const *)y, (ULPoint *)dydx);
        dydx[0] = b.x;
        dydx[1] = b.y;
        dydx[2] = b.z;
        
        // Field magnitude
        double mag = sqrt(dydx[0]*dydx[0] + dydx[1]*dydx[1] + dydx[2]*dydx[2]);
        if (mag<1e-6) {
            fprintf(stderr, "Magnetic field magnitude is too small at (%f, %f, %f).\n", pt.x, pt.y, pt.z);
            mag = NAN;
        }
        
        // Unit vector
        dydx[0] /= mag; // dx
        dydx[2] /= mag; // dy
        dydx[1] /= mag; // dz
    }
    void FieldLine::handleEventAtStep_x_y_count_stop(long istep, double x, const double *y, long count, BOOL *stop)
    {
        struct _TMP *tmp = reinterpret_cast<struct _TMP *>(_tmp);

        //
        *stop |= istep>_maxSteps;
        tmp->r0 = sqrt(y[0]*y[0] + y[1]*y[1] + y[2]*y[2]);
        *stop |= *stop || tmp->r0<1. || tmp->r0>_radiusMax || isnan(tmp->r0);

        Point pt(y[0], y[1], y[2]);

        Point b;
        tmp->fm->getCartesianFieldInSM_atCartesianPoint(&b, pt);
        double bm = sqrt(b.x*b.x + b.y*b.y + b.z*b.z);

        if ((tmp->coeqtemp.z)*(y[2]) <= 0.) {
            double d = tmp->coeqtemp.z / (y[2] - tmp->coeqtemp.z);
            _coordinateEquator.x = tmp->coeqtemp.x + d*(y[0]-tmp->coeqtemp.x);
            _coordinateEquator.y = tmp->coeqtemp.y + d*(y[1]-tmp->coeqtemp.y);
            _coordinateEquatorFieldMagnitude = bm;
        }
        tmp->coeqtemp = pt;

        if (tmp->r0 >= _radiusMin) { // Check for ionospheric boundary
            _mirrorMagnitudes.push_back(bm);
            _coordinates.push_back(tmp->coeqtemp);
        }
    }
    void FieldLine::traceCartesianFieldLineUsingFieldBlock (FieldModel const* fm, IntegrateODESignature odeSolver)
    {
        struct _TMP tmp;
        _tmp = reinterpret_cast<void *>(&tmp);
        tmp.fm = fm;

        // Work spaces
        tmp.r0 = 0.;
        double x0 = 0;
        Point y0 = tmp.coeqtemp;
        long count = 0;

        _coordinates.reserve(1000);
        _mirrorMagnitudes.reserve(1000);

        // Down trace
        {
            tmp.coeqtemp = y0 = _origin;
            x0 = 0.;
            odeSolver((-_stepSize), &x0, y0.xyz, 3, this);
            
            if ((count=_mirrorMagnitudes.size())<5 || tmp.r0>1.) {
                _valid = NO;
                return;
            } else {
                reverse(_mirrorMagnitudes.begin(), _mirrorMagnitudes.end());
                reverse(_coordinates.begin(), _coordinates.end());

                _mirrorMagnitudes.pop_back();
                _coordinates.pop_back();
            }
        }

        // Up trace
        {
            tmp.coeqtemp = y0 = _origin;
            x0 = 0.;
            odeSolver((_stepSize), &x0, y0.xyz, 3, this);

            if (-count+_mirrorMagnitudes.size()<5 || tmp.r0>1.) {
                _valid = NO;
                return;
            }
        }

        // Find foot point
        {
            double r0 = tmp.r0;
            double costh = y0.z / r0;
            double z2 = 1. - (1.-costh*costh)/r0;
            _footPoint.z = (costh<0.?-1.:1.) * sqrt(z2);
            double rhoRatio = sqrt((1.-z2)/(y0.x*y0.x+y0.y*y0.y));
            _footPoint.x = rhoRatio*y0.x;
            _footPoint.y = rhoRatio*y0.y;
        }

        _valid = YES;
    }

    void FieldLine::evaluateInvariant()
    {
        if (!_valid) {
            return;
        }

        // Find extrema
        long dnmin = ULONG_MAX, upmin = ULONG_MAX, eqidx = ULONG_MAX;
        {
            long count = _mirrorMagnitudes.size();
            double btmp = _mirrorMagnitudes.front();
            for (long idx=1; idx<count; idx++) {
                if ( btmp-_mirrorMagnitudes[idx] >= 0 ) {
                    dnmin = idx;
                    btmp = _mirrorMagnitudes[idx];
                } else {
                    break;
                }
            }
            btmp = _mirrorMagnitudes.back();
            for (long idx=count-2; idx>=0; idx--) {
                if ( btmp-_mirrorMagnitudes[idx] >= 0 ) {
                    upmin = idx;
                    btmp = _mirrorMagnitudes[idx];
                } else {
                    break;
                }
            }
            eqidx = dnmin;
            if (upmin-dnmin>1) {
                for (long idx=dnmin+1; idx<=upmin; idx++) {
                    if (btmp-_mirrorMagnitudes[idx]<0) {
                        btmp = _mirrorMagnitudes[idx];
                        eqidx = idx;
                    }
                }
                // Check if Beq is the same as either Bs(dnmin) or Bs(upmin). If so, choose index of MIN(Bs(dnmin), Bs(upmin))
                if (eqidx == dnmin) {
                    eqidx = upmin;
                } else if (eqidx == upmin) {
                    eqidx = dnmin;
                }
            }
            if (!(eqidx<count && upmin<count && dnmin<count) ) {
                fprintf(stderr, "Extrema and equator indexes are greater than count: origin=(%f, %f, %f)\n", _origin.x, _origin.y, _origin.z);
                _valid = NO;
                return;
            }

            // Flip if Bdnmin < Bupmin because K(90˚) won't be 0 otherwise.
            if (_mirrorMagnitudes[dnmin] < _mirrorMagnitudes[upmin]) {
                reverse(_mirrorMagnitudes.begin(), _mirrorMagnitudes.end());
                reverse(_coordinates.begin(), _coordinates.end());

                long mintmp = upmin;
                upmin = count - dnmin-1;
                dnmin = count - mintmp-1;
                eqidx = count - eqidx-1;
            }
            // Equator field
            _magneticEquatorFieldMagnitude = _mirrorMagnitudes[eqidx];
            _magneticEquator = _coordinates[eqidx];
        }

        if (!(upmin>=dnmin && dnmin<=eqidx && upmin>=eqidx) ) {
            fprintf(stderr, "Extrema and equator indexes are invalid: origin=(%f, %f, %f)\n", _origin.x, _origin.y, _origin.z);
            _valid = NO;
            return;
        }

        // Invariant
        {
            double ds_2 = _stepSize*.5;
            double Ku = 0., Kd = 0.;
            double bm, b1, b0;
            double bdnmax = _mirrorMagnitudes.front();

            _modifiedInvariants.reserve(_mirrorMagnitudes.size()/2 + 1);
            _modifiedInvariants.push_back(0.);
            for (long i=upmin+1, count=_mirrorMagnitudes.size(); i<count; i++) {
                bm = _mirrorMagnitudes[i];
                // If bm is exceeding bdnmax, stop
                if (bm-bdnmax > 0.) {
                    break;
                }
                // Upper part
                Ku = 0.;
                Kd = 0.;
                for (long i2=eqidx+1; i2<=i; i2++) {
                    b0 = _mirrorMagnitudes[i2-1];
                    b1 = _mirrorMagnitudes[i2];
                    if (bm>=b1) {
                        if (bm>=b0) { // poleward slope
                            Ku += (sqrt(bm - b1) + sqrt(bm - b0)) * ds_2;
                        } else { // equatorward slope
                            double dx = (bm - b1) / (b0 - b1);
                            Ku += (sqrt(bm - b1)) * ds_2*dx;
                        }
                    }
                }
                // Lower part
                for (long j=eqidx-1; j>0; j--) { // Both ends are not used because it is <rmin.
                    b0 = _mirrorMagnitudes[j+1];
                    b1 = _mirrorMagnitudes[j];
                    if (bm>=b1) {
                        if (bm>=b0) { // poleward slope
                            Kd += (sqrt(bm - b1) + sqrt(bm - b0)) * ds_2;
                        } else { // equatorward slope
                            double dx = (bm - b1) / (b0 - b1);
                            Kd += (sqrt(bm - b1)) * ds_2*dx;
                        }
                    } else if (bm>=b0) { // When bm is btw grid. Otherwise, skip
                        double dx = (bm - b0) / (b1 - b0);
                        Kd += (sqrt(bm - b0)) * ds_2*dx;
                    } else if (j < dnmin) { // poleward slope. bm < b1,b0
                        break;
                    }
                }
                // Add
                Ku = (Ku+Kd);
                _modifiedInvariants.push_back(Ku);
            }
        }

        // Wrap up
        _mirrorMagnitudes.erase(_mirrorMagnitudes.begin(), _mirrorMagnitudes.begin()+upmin);
        _mirrorMagnitudes.resize(_modifiedInvariants.size());

        _coordinates.erase(_coordinates.begin(), _coordinates.begin()+upmin);
        _coordinates.resize(_modifiedInvariants.size());
    }

}