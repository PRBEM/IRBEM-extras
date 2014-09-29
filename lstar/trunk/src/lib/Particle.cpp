//
//  Particle.cpp
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

#include "Particle.h"
#include "MagneticFlux.h"
#include "AffineTransform.h"
#include "FieldModel.h"
#include "TSFieldModel.h"
#include "FieldLineTable.h"
#include "FieldLine.h"
#include "ParticleTemplate.h"
#include <cmath>
#include <climits>
#include <cassert>
#include <cstdio>
#include <stdexcept>

namespace UBK {
    using namespace std;

    //
    // Test
    //
#ifdef DEBUG
    void Particle::test()
    {
        //double IMFb1 = 0.;
        //double compb2 = 4.;
        //CompressedDipole fm(IMFb1, compb2);
        TSFieldModel fm(Date(2001, 1, 1, 0, 0), 400, kGeopackDipoleField, 2, Parmod(3, -10, 3, -10, 0, 0, 0, 0, 0, 0), kTS89Model);

        FieldLine::setStepSize(.05);
        MagneticFlux::setNPhi(360/2);
        MagneticFlux mf(&fm);
        printf("F0 = %f\n", mf.coreFlux());

        AffineTransform tr(Point(.1, M_PI/180., .1), Point(0., 0., 0.));
        Point rl(.9, 0.);
        Point ru(12., M_PI*2.);
        rl = tr.pointConvertedToNormal(rl);
        ru = tr.pointConvertedToNormal(ru);
        FieldLineTableCyl<FieldLineTable, FieldLine> cyltbl(Key(rl.r, rl.phi-1, 0.), Key(ru.r, ru.phi+1, 0.), &fm, &tr);

        AffineTransform tx(Point(.1, .1, .1), Point(0., 0., 0.));
        Point xl(-12., -12.);
        Point xu(12., 12.);
        xl = tx.pointConvertedToNormal(xl);
        xu = tx.pointConvertedToNormal(xu);
        FieldLineTableCart<FieldLineTable, FieldLine> carttbl(Key(xl.x, xl.y, 0.), Key(xu.x, xu.y, 0.), &fm, &tx);

        Point origin(6, 0., 0.);
        double pas[18] = {5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90};
        for (long jdx=17; jdx<18; jdx++) {
            double pa = pas[jdx] * M_PI/180.;
            ParticleCyl<FieldLineTable, FieldLine> cylp(origin, pa, &cyltbl, &mf, &tr, YES, YES);
            ParticleCart<FieldLineTable, FieldLine> cartp(origin, pa, &carttbl, &mf, &tx, YES, YES);
            cylp.evaluateLstar();
            cartp.evaluateLstar();

            printf("cylL* = %f, cartL* = %f, dL* = %f\n", cylp.Lstar(), cartp.Lstar(), fabs(cylp.Lstar()-cartp.Lstar()));
            for (long idx=0; NO&&idx<cylp.coordinates().size(); idx++) {
                printf("%s;...\n", cylp.coordinates()[idx].desc().c_str());
            }
        }

        assert(0);
    }
#endif

    //
    // Constructor
    //
    Particle::Particle (Point const origin,
                        double const pitchAngle,
                        FieldLineTable *table,
                        MagneticFlux const* magneticFlux,
                        AffineTransform const* affineTransform,
                        BOOL const shouldKeepFootPoint,
                        BOOL const shouldKeepContourCoordinate) :
    _origin(origin),
    _pitchAngle(pitchAngle),
    _table(table),
    _magneticFlux(magneticFlux),
    _affineTransform(affineTransform),
    _shouldKeepFootPoints(shouldKeepFootPoint),
    _shouldKeepContourCoordinates(shouldKeepContourCoordinate),
    _closed(NO),
    _orbitType(kOrbitTypeUnknown),
    _Lstar(NAN),
    _K(NAN),
    _coordinates(),
    _footPoints()
    {
        assert(table);
        assert(magneticFlux);
        assert(affineTransform);
    }

    Particle const& Particle::evaluateLstar (Contour const* cPtr)
    {
        if (NULL == cPtr) {
            Contour contour(this);
            this->evaluateLstarWithContour(&contour);
        } else {
            this->evaluateLstarWithContour(cPtr);
        }

        return *this;
    }

    void Particle::getLevel_atNode (double *lvlOut, Key const xy)
    {
        FieldLine const* fl;
        if (NULL == (fl=(*_table)[xy])) {
            *lvlOut = NAN;
            _orbitType = kOrbitTypeOutOfBound;
        } else {
            *lvlOut = fl->mirrorMagnitudeAtInvariant(_K);
            if (!fl->isValid()) {
                _orbitType = kOrbitTypeTraceFailed;
            } else if (this->K()>=fl->modifiedInvariants().back()) { // Precipitation
                _orbitType = kOrbitTypePrecipitation;
            } else if (isnan(*lvlOut)) {
                _orbitType = kOrbitTypeBifurcation;
            }
        }
    }

    ////////////////////////////////////////////////////////////////
    // Heavy lifting
    ////////////////////////////////////////////////////////////////

    //
    // Helper 1
    //
    static OrbitType Table_BmAndKForKey_coordinate_pitchAngle (FieldLineTable *table, double *BmOut, double *KOut, Key const key, Point const coord, double const pa)
    {
        *BmOut = *KOut = NAN;

        FieldLine const* fl;

        if (NULL == (fl=(*table)[key])) return kOrbitTypeOutOfBound;
        double Bmbl = fl->mirrorMagnitudeOfPitchAngle(pa);
        double Kbl = fl->invariantAtMirrorMagnitude(Bmbl);
        if (!fl->isValid()) {
            return kOrbitTypeTraceFailed;
        } else if (Kbl>=fl->modifiedInvariants().back()) { // Precipitation
            return kOrbitTypePrecipitation;
        } else if (isnan(Bmbl)) {
            return kOrbitTypeBifurcation;
        }

        if (NULL == (fl=(*table)[Key(key.x+1, key.y, 0)])) return kOrbitTypeOutOfBound;
        double Bmbr = fl->mirrorMagnitudeOfPitchAngle(pa);
        double Kbr = fl->invariantAtMirrorMagnitude(Bmbr);
        if (!fl->isValid()) {
            return kOrbitTypeTraceFailed;
        } else if (Kbr>=fl->modifiedInvariants().back()) { // Precipitation
            return kOrbitTypePrecipitation;
        } else if (isnan(Bmbr)) {
            return kOrbitTypeBifurcation;
        }

        if (NULL == (fl=(*table)[Key(key.x+1, key.y+1, 0)])) return kOrbitTypeOutOfBound;
        double Bmtr = fl->mirrorMagnitudeOfPitchAngle(pa);
        double Ktr = fl->invariantAtMirrorMagnitude(Bmtr);
        if (!fl->isValid()) {
            return kOrbitTypeTraceFailed;
        } else if (Ktr>=fl->modifiedInvariants().back()) { // Precipitation
            return kOrbitTypePrecipitation;
        } else if (isnan(Bmtr)) {
            return kOrbitTypeBifurcation;
        }

        if (NULL == (fl=(*table)[Key(key.x, key.y+1, 0)])) return kOrbitTypeOutOfBound;
        double Bmtl = fl->mirrorMagnitudeOfPitchAngle(pa);
        double Ktl = fl->invariantAtMirrorMagnitude(Bmtl);
        if (!fl->isValid()) {
            return kOrbitTypeTraceFailed;
        } else if (Ktl>=fl->modifiedInvariants().back()) { // Precipitation
            return kOrbitTypePrecipitation;
        } else if (isnan(Bmtl)) {
            return kOrbitTypeBifurcation;
        }

        double dx = coord.x - key.x;
        double dy = coord.y - key.y;

        //
        // Areal interpolation
        //
        *BmOut = (dx-1.)*(dy-1.)*Bmbl + dy*Bmtl + dx*(Bmbr + (Bmtr-Bmtl-Bmbr)*dy);
        *KOut = (dx-1.)*(dy-1.)*Kbl + dy*Ktl + dx*(Kbr + (Ktr-Ktl-Kbr)*dy);

        return kOrbitTypeClosed;
    }

    //
    // Evaluate contour and L*. It assumes that the real grid space is in cartesian coordinate system.
    //
    void Particle::evaluateLstarWithContour(const Contour *cPtr)
    {
        // Conserved quantities
        Point normal = _affineTransform->pointConvertedToNormal(_origin);
        double Bm = 0.;
        OrbitType ot = Table_BmAndKForKey_coordinate_pitchAngle(_table, &Bm, &_K, Key((int)floor(normal.x), (int)floor(normal.y), 0), normal, _pitchAngle);

        // Contour trace
        ContourTermination closed = cPtr->enumerateCoordinatesStartingAtPoint_usingCallback(normal, this);
        _closed = (closed==kContourClosed);

        if (ot) {
            _orbitType = ot;
        }

        // Evaluate L*
        if (_closed) {
            _Lstar = _magneticFlux->fluxClosedByFootPoints_count(&_footPoints.front(), _footPoints.size());
            _Lstar = _magneticFlux->coreFlux() / _Lstar;

            _orbitType = kOrbitTypeClosed;
        }

        if (_shouldKeepContourCoordinates) {
            _affineTransform->convertPoints_fromNormal_count(&_coordinates.front(), &_coordinates.front(), _coordinates.size());
        } else {
            _coordinates.clear();
        }

        if (!_shouldKeepFootPoints) {
            _footPoints.clear();
        }
    }

    //
    // Helper 2
    //
    static Point Table_footPointForKey_coordinate_stopFlag (FieldLineTable *table, Key key, CellSide const side, Point const coord, BOOL *stop)
    {
        FieldLine const* fl1, *fl2;
        double d = 0.;

        switch (-side) { // Rule is that if -side==kCellBottomSide, "coord" is on the bottom side of the current "key" cell (i.e., the index of the bottom-left corner is key).
            case kCellBottomSide:
            {
                d = coord.x - key.x;
                fl1 = (*table)[key];
                key.x++;
                fl2 = (*table)[key];
            }
                break;
            case kCellRightSide:
            {
                key.x++;
                d = coord.y - key.y;
                fl1 = (*table)[key];
                key.y++;
                fl2 = (*table)[key];
            }
                break;
            case kCellTopSide:
            {
                key.y++;
                d = coord.x - key.x;
                fl1 = (*table)[key];
                key.x++;
                fl2 = (*table)[key];
            }
                break;
            case kCellLeftSide:
            {
                d = coord.y - key.y;
                fl1 = (*table)[key];
                key.y++;
                fl2 = (*table)[key];
            }
                break;
            default:
                throw invalid_argument("Unknown side. This should not reach!!");
                break;
        }

        Point foot;
        foot.x = fl1->footPoint().x + d*(fl2->footPoint().x - fl1->footPoint().x);
        foot.y = fl1->footPoint().y + d*(fl2->footPoint().y - fl1->footPoint().y);
        foot.z = fl1->footPoint().z + d*(fl2->footPoint().z - fl1->footPoint().z);

        //
        // In-place transformation to spherical coordinate
        //
        foot.theta = acos(foot.z);
        foot.phi = atan2(foot.y, foot.x);
        foot.r = 1.;
        if (foot.phi<0.) { // Between 0 and 2pi
            foot.phi += 2.*M_PI;
        }

        *stop = isnan(foot.phi+foot.theta) || (foot.phi<0.) || (foot.phi>2.*M_PI);
        return foot;
    }
    //
    // Contour event handler
    //
    void Particle::handleEventAtIndex_node_cellSide_coordinate_level_stop (long const idx, Key *bottomLeftOut, CellSide const side, Point const coordinate, double const level, BOOL *stop)
    {
        // Foot point interpolation
        Point footPoint = Table_footPointForKey_coordinate_stopFlag(_table, *bottomLeftOut, side, coordinate, stop);
        _coordinates.push_back(coordinate);
        _footPoints.push_back(footPoint);
    }

}