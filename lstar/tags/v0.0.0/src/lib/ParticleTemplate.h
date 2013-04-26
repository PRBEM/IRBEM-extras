//
//  ParticleTemplate.h
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

#ifndef UBJDevelopment_ParticleTemplate_h
#define UBJDevelopment_ParticleTemplate_h

#include "Definitions.h"
#include "FieldLineTable.h"
#include "Particle.h"
#include <cassert>
#include <cmath>

//
// Convenience templates.
//
namespace UBK {

    //!
    //! @name Field line table templates.
    //!
    ///@{
    //!
    //! Concrete cartesian table template.
    //!
    class FieldModel;
    template <class FLTBL, class FL>
    class FieldLineTableCart : public FLTBL {
    private:
        AffineTransform const* _t;
        FieldModel const* _fm;
    public:
        //!
        //! Additionally, initialize with field model and transform objects.
        //!
        FieldLineTableCart (Key const lBound, Key const uBound, FieldModel const* fm, AffineTransform const* t) : FLTBL(lBound,uBound), _t(t), _fm(fm) {
            assert(NULL!=fm && NULL!=t);
        };

        virtual FL const* newFieldLineAtNode (Key const key) const {
            Point origin(key.x, key.y, key.z);
            origin = _t->pointConvertedFromNormal(origin);
            return new FL(origin, _fm);
        };
    };

    //!
    //! Concrete cylindrical table template.
    //! @note Note the wrapping of the azimuthal angle. Keep in mind that cos(or sin)(2pi) != cos(or sin)(0).
    //!
    template <class FLTBL, class FL>
    class FieldLineTableCyl : public FLTBL {
    private:
        AffineTransform const* _t;
        FieldModel const* _fm;
        long _maxPhi;
    public:
        //!
        //! Additionally, initialize with field model and transform objects.
        //!
        FieldLineTableCyl (Key const lBound, Key const uBound, FieldModel const* fm, AffineTransform const* t) : FLTBL(lBound,uBound), _t(t), _fm(fm) {
            assert(NULL!=fm && NULL!=t);
            Point maxBd = _t->pointConvertedToNormal(Point(1., 2*M_PI, 0.)); // .r component is dummy, .phi component is used for wrapping at 2pi
            _maxPhi = floor(maxBd.phi);
        };

        virtual FieldLine const* operator [] (Key const key) {
            //
            // Wrap phi and make sure it is between 0 and 2pi.
            //
            Key wrappedKey(key.r, key.phi%_maxPhi, 0);
            if (wrappedKey.phi>=_maxPhi) {
                wrappedKey.phi -= _maxPhi;
            } else if (wrappedKey.phi<0) {
                wrappedKey.phi += _maxPhi;
            }

            //
            // Reuse parent's implementation.
            //
            return FLTBL::operator [] (wrappedKey);
        }

        virtual FL const* newFieldLineAtNode (Key const key) const {
            Point origin(key.r, key.phi, 0.);
            origin = _t->pointConvertedFromNormal(origin);
            return new FL(Point(origin.r*cos(origin.phi), origin.r*sin(origin.phi), 0.), _fm);
        };
    };
    ///@}

    //!
    //! @name Particle templates.
    //!
    ///@{
    //!
    //! Particle cartesian templates.
    //! This is basically the same as Particle class.
    //!
    template <class FLTBL, class FL>
    class ParticleCart : public Particle {
    public:
        virtual ~ParticleCart() {};
        ParticleCart (Point const origin, double const pitchAngle, FieldLineTableCart<FLTBL, FL> *table, MagneticFlux const* magneticFlux, AffineTransform const* affineTransform, BOOL const shouldKeepFootPoint = NO, BOOL const shouldKeepContourCoordinate = NO) : Particle(origin, pitchAngle, table, magneticFlux, affineTransform, shouldKeepFootPoint, shouldKeepContourCoordinate) {};
    };

    //!
    //! Particle cylindrical templates.
    //! Note also on the wrapping which is necessary for cylindrical coordinates. Otherwise, contour tracing won't stop.
    //! Make sure 0<=origin.phi<2pi.
    //!
    template <class FLTBL, class FL>
    class ParticleCyl : public Particle {
    private:
        long _maxPhi;
    public:
        virtual ~ParticleCyl() {};
        ParticleCyl (Point const origin, double const pitchAngle, FieldLineTableCyl<FLTBL, FL> *table, MagneticFlux const* magneticFlux, AffineTransform const* affineTransform, BOOL const shouldKeepFootPoint = NO, BOOL const shouldKeepContourCoordinate = NO) : Particle(origin, pitchAngle, table, magneticFlux, affineTransform, shouldKeepFootPoint, shouldKeepContourCoordinate) {
            Point maxBd = _affineTransform->pointConvertedToNormal(Point(1., 2*M_PI, 0.)); // .r component is dummy, .phi component is used for wrapping at 2pi
            _maxPhi = floor(maxBd.phi);
        }

        virtual void handleEventAtIndex_node_cellSide_coordinate_level_stop (long const idx, Key *bottomLeftOut, CellSide const side, Point const coordinate, double const level, BOOL *stop) {
            Particle::handleEventAtIndex_node_cellSide_coordinate_level_stop(idx, bottomLeftOut, side, coordinate, level, stop);
            // Wrap around azimuthal angle
            bottomLeftOut->phi %= _maxPhi;
            if (bottomLeftOut->phi>=_maxPhi) {
                bottomLeftOut->phi -= _maxPhi;
            } else if (bottomLeftOut->phi<0) {
                bottomLeftOut->phi += _maxPhi;
            }
        }
    };
    ///@}
    
}

#endif
