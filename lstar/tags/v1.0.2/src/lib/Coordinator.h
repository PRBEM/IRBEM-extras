//
//  Coordinator.h
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

#ifndef UBJDevelopment_Coordinator_h
#define UBJDevelopment_Coordinator_h

#include "Definitions.h"
#include "Thread.h"
#include "Point.h"
#include "Key.h"
#include "TSFieldModel.h"
#include "MagneticFlux.h"
#include "AffineTransform.h"
#include "FieldLineTable.h"
#include "Particle.h"
#include "ParticleTemplate.h"
#include <stdexcept>
#include <vector>
#include <cmath>

namespace UBK {

    //!
    //! Coordinate system flag.
    //!
    enum _MagneticFieldCoordinateSystem {
        kMagneticFieldGSM = 0, //!< GSM
        kMagneticFieldSM = 1 //!< SM
    };
    typedef long MagneticFieldCoordinateSystem;

    //!
    //! Coordinater abstract.
    //! Base class of all coordinator classes.
    //!
    //! Data parallelism using ThreadFor.
    //!
    class Coordinator {
    private:
        FieldModel const& _fm;
        unsigned long _n_threads;

    public:
        FieldModel const& fieldModel() const {return _fm;};

    public:
        //!
        //! The number of threads to create.
        //! Default is 8. Must be greater than 0. Otherwise an exception is thrown.
        //!
        void setNThreads(unsigned long n_threads) {_n_threads = n_threads;};

    public:
        virtual ~Coordinator() {};
        Coordinator(FieldModel const& fm) : _fm(fm), _n_threads(8) {};

    public:
        //!
        //! Initiates concurrent evaluation.
        //! Call this if ready.
        //!
        virtual void start();

        //!
        //! Task.
        //! Caller should implement this method.
        //! @param [in] idx Index on which the current thread is working.
        //!
        virtual void operator()(long idx) = 0;

        //!
        //! The number of elements.
        //! Caller should return a valid count.
        //!
        virtual unsigned long count() const = 0;
    };

    //!
    //! Cotrans coordinator abstract.
    //! Parallel execution of GSM<->SM coordinate transformation.
    //!
    class CotransCoordinator : public Coordinator {
#ifdef DEBUG
    public:
        static void test();
#endif

    private:
        MagneticFieldCoordinateSystem _to_co_system;

    public:
        virtual ~CotransCoordinator() {};
        CotransCoordinator(FieldModel const& fm) : Coordinator(fm), _to_co_system(kMagneticFieldSM) {};

    public:
        //!
        //! Set the coordinate system to which the caller wants to convert.
        //! Only transform between GSM and SM is available.
        //!
        void setToCoSystem(MagneticFieldCoordinateSystem to_co_system);

    public:
        void operator()(long idx);

    public:
        //!
        //! A coordinate at the given index.
        //!
        virtual Point from(long idx) const = 0;

        //!
        //! Called upon completion of individual jobs.
        //! @param [in] idx Index on which the thread worked on.
        //! @param [in] pt Converted point.
        //!
        virtual void callback(long idx, Point const& pt) const = 0;
    };

    //!
    //! Model magnetic field coordinator abstract.
    //! Parallel evaluation of magnetic field.
    //! @note Available coordinates are cartesian GSM and SM coordinate systems.
    //!
    class BFieldCoordinator : public Coordinator {
#ifdef DEBUG
    public:
        static void test();
#endif

    private:
        MagneticFieldCoordinateSystem _co_system;

    public:
        virtual ~BFieldCoordinator() {};
        BFieldCoordinator(FieldModel const& fm) : Coordinator(fm), _co_system(kMagneticFieldGSM) {};

    public:
        //!
        //! GSM or SM.
        //!
        void setCoSystem(MagneticFieldCoordinateSystem co_system);

    public:
        void operator()(long idx);

    public:
        //!
        //! Coordinate at which to calculate the magnetic field.
        //!
        virtual Point at(long idx) const = 0;

        //!
        //! Called upon completion of individual calculations.
        //! @param [in] idx Index on which the thread worked on.
        //! @param [in] B Calculated magnetic field.
        //!
        virtual void callback(long idx, Point const& B) const = 0;
    };

    //!
    //! Field line coordinator abstract.
    //! Parallel evaluation of field line parameters.
    //! @note SM cartesian coordinate system otherwise specified.
    //!
    class FieldLineCoordinator : public Coordinator {
#ifdef DEBUG
    public:
        static void test();
#endif

    public:
        virtual ~FieldLineCoordinator() {};
        FieldLineCoordinator(FieldModel const& fm) : Coordinator(fm) {
            this->setIonoR(1.015);
        };

    public:
        //!
        //! Set ionospheric boundary.
        //! Radial distance from the center of the Earth in RE. Default is 1.05 RE (~100 km).
        //!
        void setIonoR(double ionoR) const;

        //!
        //! Field line integration step size in RE.
        //! Default is 0.05 RE.
        //!
        void setDs(double ds) const;

    public:
        void operator()(long idx);

    public:
        //!
        //! Location from which to start tracing.
        //!
        virtual Point from(long idx) const = 0;

        //!
        //! Called upon completion of the job at index, idx.
        //! @param [in] idx Index.
        //! @param [in] fl Reference to FieldLine object.
        //!
        virtual void callback(long idx, FieldLine const& fl) const = 0;
    };

    //!
    //! L* coordinator abstract.
    //! Parallel evaluation of L*.
    //! @note SM cartesian/cylindrical coordinate system otherwise specified.
    //!
    class LstarCoordinator : public Coordinator {
#ifdef DEBUG
    public:
        static void test();
#endif

    private:
        AffineTransform _t;
        MagneticFlux const* _mf;
        FieldLineTable* _ftbl;
        double _Phi0;
        BOOL _isCylindricalGrid;

    public:
        virtual ~LstarCoordinator() {};
        //!
        //! Constructor.
        //! @param [in] fm Reference to FieldModel object.
        //! @param [in] isCylindricalGrid Set to YES to use the cylindrical grid. Cartesian grid is used otherwise.
        //!
        LstarCoordinator(FieldModel const& fm, BOOL isCylindricalGrid) : Coordinator(fm), _t(Point(1., 1., 1.), Point()), _mf(NULL), _ftbl(NULL), _Phi0(NAN), _isCylindricalGrid(isCylindricalGrid) {
            this->setIonoR(1.015);
            this->setDs(.1);
            this->setD(Point(.2, .2, 1.));
        };

    public:
        //!
        //! Set ionospheric radial distance from the center of the Earth in RE.
        //! Default is 1.015 RE.
        //! @note Same effect as calling FieldLine::setRadialMin(ionoR).
        //!
        void setIonoR(double ionoR) const;

        //!
        //! Set the field line integration step size in RE.
        //! Default is 0.05 RE.
        //! @note Same effect as calling FieldLine::setStepSize(ionoR).
        //!
        void setDs(double ds) const;

        //!
        //! Set the number of azimuthal grids for magnetic flux integration.
        //! Default is 360/5.
        //! @note Same effect as calling MagneticFlux::setNPhi(nPhi).
        //!
        void setNPhi(unsigned long nPhi) const;

        //!
        //! Set the number of polar grids for magnetic flux integration.
        //! Default is 180*2.
        //! @note Same effect as calling MagneticFlux::setNTheta(nPhi).
        //!
        void setNTheta(unsigned long nTheta) const;

        //!
        //! Set grid resolution in RE.
        //! Appropriate grid resolution corresponding to the grid type should be set.
        //!
        void setD(Point d);

        //!
        //! Magnetic flux through the area defined by Earth' equator.
        //! Set to non-NaN value (if defined) after start() returns or callback() is called.
        //! @see MagneticFlux::coreFlux
        //!
        double const& Phi0() const {return _Phi0;};

    public:
        void start();
        void operator()(long idx);

    public:
        //!
        //! Initial value of the particle at index.
        //! Return the Point object whose first and second (i.e., x and y or r and phi) components represent the location mapped on the SM coordinate equator along the field line, and whose third represents the initial pitch angle at the MAGNETIC EQUATOR.
        //! For example, to trace a particle at XYZ = (6, 0, 0) and pa0 = 90˚, return Point(6, 0, pi/2), and to trace a particle at RPHITHETA = (6, 180˚, 90˚) and pa0 = 45˚, return Point(6, pi, pi/4).
        //!
        //! Distance in RE and angle in radian.
        //!
        virtual Point from(long idx) const = 0;

        //!
        //! Called upon completion of individual tracing.
        //! @param [in] idx Index.
        //! @param [in] ptl Reference to Particle object.
        //!
        virtual void callback(long idx, Particle const& ptl) const = 0;
    };

}

#endif
