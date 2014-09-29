//
//  FieldLine.h
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

#ifndef UBJDevelopment_FieldLine_h
#define UBJDevelopment_FieldLine_h

#include "Point.h"
#include "ODESolver.h"
#include <vector>

//
// Field line tracer, and Bm and K evaluator. SM coordinate system.
//
namespace UBK {

    class FieldModel; // FieldModel forward declaration. FieldModel object should handle getCartesianFieldInSM_atCartesianPoint call.

    //!
    //! @brief Field line tracer, and Bm and K evaluator.
    //! @details Constructs field line passing at the given point and calculates parameters related to field line geometry (line coordinates, mirror magnitude vs. modified 2nd invariant, equator field magnitude, field line foot point (north)).
    //!
    //! It is considered to be invalid when (1) the tracing is beyond the predefined domain (except crash into Earth), (2) maximum iteration count of integration is passed, (3) total line length is too short and (4) the line is not touching minimum specified radius sphere. When the field line is invalid, valid property is set to NO.
    //!
    //! The caller should initialize with the magnetic field model object.
    //! @note SM coordinate system otherwise specified.
    //!
    class FieldLine : public ODECallback {
#ifdef DEBUG
    public:
        static void test();
#endif

    public:
        //
        // Inherited from ODECallback.
        //
        virtual void derivatives (double x, double const y[], long count, double *dydx);
        virtual void handleEventAtStep_x_y_count_stop (long istep, double x, double const y[], long count, BOOL *stop);

    protected:
        Point _origin;
        std::vector<Point> _coordinates;
        std::vector<double> _mirrorMagnitudes;
        std::vector<double> _modifiedInvariants;
        Point _coordinateEquator;
        double _coordinateEquatorFieldMagnitude;
        Point _magneticEquator;
        double _magneticEquatorFieldMagnitude;
        Point _footPoint;
        BOOL _valid;

    private:
        //
        // Static options.
        //
        static double _stepSize;
        static double _radiusMax;
        static double _radiusMin;
        static long _maxSteps;

    public:
        //!
        //! @name Static options.
        //!
        ///@{
        //@{
        //!
        //! Integration step size in RE.
        //! Default is 0.05 RE.
        //!
        static void setStepSize (double const stepSize) {_stepSize=stepSize;};
        static double stepSize () {return _stepSize;};
        //@}

        //@{
        //!
        //! Maximum outer boundary in RE.
        //! Default is 15 RE.
        //!
        static void setRadiusMax (double const radiusMax) {_radiusMax=radiusMax;};
        static double radiusMax () {return _radiusMax;};
        //@}

        //@{
        //!
        //! Radial distance of the ionospheric boundary in RE.
        //! Default is 1 RE.
        //!
        static void setRadiusMin (double const radiusMin) {_radiusMin=radiusMin;};
        static double radiusMin () {return _radiusMin;};
        //@}

        //@{
        //!
        //! Max integration steps at which tracing will stop.
        //! Default is 10000.
        //!
        static void setMaxSteps (unsigned long const maxSteps) {_maxSteps=maxSteps;};
        static long maxSteps () {return _maxSteps;};
        //@}
        ///@}

        //!
        //! @name Property accessors.
        //!
        ///@{
        //!
        //! Origin from which the field line was traced.
        //! @note Set at initialization.
        //!
        Point const& origin () const {return _origin;};

        //!
        //! Array of line node coordinates in cartesian coordinate system.
        //! The coordinates starts from the magnetic minimum location  to the radiusMin. If the minimum location is on the northern hemisphere, the coordinates are northern portion of the field line, otherwise southern portion.
        //! @note Empty if the field line tracing has been failed.
        //!
        std::vector<Point> const& coordinates () const {return _coordinates;};

        //!
        //! Magnetic mirror strength (Bm).
        //!Starts from the magnetic minimum to radiusMin. If the minimum location is on the northern hemisphere, the Bm is northern portion, otherwise southern portion.
        //! @note Empty if the field line tracing has been failed.
        //!
        std::vector<double> const& mirrorMagnitudes () const {return _mirrorMagnitudes;};

        //!
        //! Modified 2nd invariant (K, Roederer 1970).
        //! @note Empty if the field line tracing has been failed.
        //!
        std::vector<double> const& modifiedInvariants () const {return _modifiedInvariants;};

        //!
        //! (x,y,0) cartesian coordinate that field line passes.
        //!
        Point const& coordinateEquator () const {return _coordinateEquator;};

        //!
        //! Magnetic field strength at coordinate equator.
        //!
        double const& coordinateEquatorFieldMagnitude () const {return _coordinateEquatorFieldMagnitude;};

        //!
        //! Magnetic equator coordinate.
        //! @note NaN if the field line tracing has been failed.
        //!
        Point const& magneticEquator () const {return _magneticEquator;};

        //!
        //! Magnetic field magnitude at magnetic equator.
        //! @note NaN if the field line tracing has been failed.
        //!
        double const& magneticEquatorFieldMagnitude () const {return _magneticEquatorFieldMagnitude;};

        //!
        //! Magnetic foot point coordinate in northern hemisphere. Cartesian point. |footPoint|=1.
        //! @note NaN if the field line tracing has been failed.
        //!
        Point const& footPoint () const {return _footPoint;};

        //!
        //! Validity flag.
        //! NO if the field line tracing has been failed.
        //!
        BOOL isValid() const {return _valid;};
        ///@}

    public:
        virtual ~FieldLine() {};

        //!
        //! Constructor.
        //! Requires initial position and field line model pointer. Caller can supply his/her own ode solver.
        //!
        FieldLine (Point const origin, FieldModel const* fm, IntegrateODESignature odeSolver = rk4IntegrateODE);

        //!
        //! Mirror field magnitude linear intepolator from given K.
        //! @param [in] K Modified invariant.
        //! @return Mirror field magnitude. When bm<B_equator or K>K_at_foot, NAN is returned.
        //!
        virtual double mirrorMagnitudeAtInvariant (double const K) const;

        //!
        //! Modified invariant intepolator from given mirror field magnitude.
        //! @param [in] bm Mirror field magnitude.
        //! @return Modified invariant. When bm<B_at_minimum, the invariant is extrapolated and has negative value. When B_at_minimum<bm<B_equator or bm>B_at_foot, NAN is returned.
        //!
        virtual double invariantAtMirrorMagnitude(double const bm) const;

        //!
        //! Pitch angle to mirror magnitude conversion.
        //!
        virtual double mirrorMagnitudeOfPitchAngle(double const paRad) const;

        //!
        //! Mirror magnitude conversion to pitch angle conversion.
        //!
        virtual double pitchAngleOfmirrorMagnitude(double const bm) const;

    private:
        //
        // Opaque temporary variable for communication btw function calls. Internal use only.
        //
        void *_tmp;

        //
        // Called from the constructor. Does actual heavy lifting.
        //
        void traceCartesianFieldLineUsingFieldBlock (FieldModel const* fm, IntegrateODESignature odeSolver);
        void evaluateInvariant ();
        
    };
    
}

#endif
