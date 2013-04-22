//
//  Particle.h
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

#ifndef UBJDevelopment_Particle_h
#define UBJDevelopment_Particle_h

#include "Definitions.h"
#include "Key.h"
#include "Point.h"
#include "Contour.h"
#include "FieldLineTable.h"
#include <vector>
#include <cmath>

//
// Representation of a particle. Its job is to organize a series of operations to get L*. SM cartesian coordinate (otherwise noted).
//
namespace UBK {

    //
    // Forward declarations.
    //
    class MagneticFlux;
    class AffineTransform;

    //!
    //! Representation of a particle.
    //! Given the initial values of a particle, traces a drift shell and calculates L* if defined.
    //!
    //! Adopts ContourLevel, ContourCallback.
    //!
    //! SM coordinate system otherwise specified.
    //! @sa FieldLineTableCart, FieldLineTableCyl
    //!
    class Particle : public ContourLevel, public ContourCallback {
#ifdef DEBUG
    public:
        static void test();
#endif

    public:
        //
        // Adopted from ContourLevel.
        //
        virtual void getLevel_atNode (double *lvlOut, Key const xy);
        //
        // Adopted from ContourCallback.
        //
        virtual void handleEventAtIndex_node_cellSide_coordinate_level_stop (long const idx, Key *keyInOut, CellSide const side, Point const coordinate, double const level, BOOL *stop);

    protected:
        Point _origin;
        double _pitchAngle;
        FieldLineTable *_table;
        MagneticFlux const* _magneticFlux;
        AffineTransform const* _affineTransform;

        BOOL _shouldKeepContourCoordinates; // Flag to keep contour coordinates. Be caucious for memory if set to YES.
        BOOL _shouldKeepFootPoints; // Flag to keep foot point coordinates of the contour. Be caucious for memory if set to YES.

        BOOL _closed;
        double _Lstar;
        double _K;

        std::vector<Point> _coordinates;
        std::vector<Point> _footPoints;

    public:
        //!
        //! @name Property accessors.
        //!
        ///@{
        //!
        //! Initial location of the particle in xy space, i.e., z=0.
        //!
        Point const& origin () const {return _origin;};

        //!
        //! Initial pitch angle of the particle in radian.
        //! @note This is the pitch angle at the magnetic equator, NOT at origin().
        //!
        double const& pitchAngle () const {return _pitchAngle;};

        //!
        //! Field line table pointer which uses the same coordinate system.
        //!
        FieldLineTable *table () const {return _table;};

        //!
        //! Magnetic flux object that will be used for Lstar.
        //!
        MagneticFlux const* magneticFlux () const {return _magneticFlux;};

        //!
        //! Coordinate tramsformation.
        //!
        AffineTransform const* affineTransform () const {return _affineTransform;};

        //!
        //! Flag for contour connectivity. If it's NO, Lstar property returns NAN.
        //!
        BOOL const& isClosed () const {return _closed;};

        //!
        //! Calculated Lstar. NAN value is returned unless Lstar is found.
        //!
        double const& Lstar () const {return _Lstar;};

        //!
        //! Modified 2nd invariant for this particle.
        //!
        double const& K () const {return _K;};

        //!
        //! Contour coordinates in the same coordinate system. Empty if _shouldKeepContourCoordinates==NO.
        //!
        std::vector<Point> const& coordinates () const {return _coordinates;};

        //!
        //! Foot point coordinates in spherical SM. Empty if _shouldKeepFootPoints==NO.
        //!
        std::vector<Point> const& footPoints () const {return _footPoints;};
        ///@}

    public:
        virtual ~Particle () {};
        //!
        //! Constructor requires initial location, pitch angle, field line table pointer, magnetic flux pointer, affine tramsform.
        //! @param [in] origin Initial location of a particle. z-component is ignored.
        //! @param [in] pitchAngle Initial pitch angle of a particle in radian.
        //! @param [in] table A pointer to FieldLineTable object.
        //! @param [in] magneticFlux A pointer to MagneticFlux object.
        //! @param [in] affineTransform A pointer to AffineTransform object.
        //! @param [in] shouldKeepFootPoint If set to YES, footPoints() returns non-empty vector container filled with (if closed) contour foot points. Default is NO.
        //! @param [in] shouldKeepContourCoordinate If set to YES, coordinates() returns non-empty vector container filled with (if closed) contour coordinates. Default is NO.
        //!
        Particle (Point const origin, double const pitchAngle, FieldLineTable *table, MagneticFlux const* magneticFlux, AffineTransform const* affineTransform, BOOL const shouldKeepFootPoint = NO, BOOL const shouldKeepContourCoordinate = NO);

        //!
        //! Evaluate drift contour and L*.
        //! Caller can supply his/her own Contour object.
        //! Default is Contour object created by Contour() constructor.
        //! @note This is the most time consuming task.
        //!
        virtual Particle const& evaluateLstar (Contour const* cPtr = NULL);

    private:
        //
        // Actural heavy lifting.
        //
        void evaluateLstarWithContour (Contour const* cPtr);
    };

}

#endif
