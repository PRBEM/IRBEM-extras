//
//  MagneticFlux.h
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

#ifndef UBJDevelopment_MagneticFlux_h
#define UBJDevelopment_MagneticFlux_h

#include "Definitions.h"
#include "Point.h"
#include <vector>

//
// Magnetic flux evaluator defined by a contour. SM spherical coordinate system.
//
namespace UBK {

    //
    // FieldModel forward declaration. The object should handle getSphericalFieldInSM_atSphericalPoint call.
    //
    class FieldModel;

    //!
    //! Magnetic flux evaluator defined by a contour.
    //! @note SM spherical coordinate system otherwise specified.
    //!
    class MagneticFlux {
#ifdef DEBUG
    public:
        static void test();
#endif

    private:
        //
        // Global options
        //
        static unsigned long _nPhi;
        static unsigned long _nTheta;
        static unsigned long _Nthreads;

    private:
        //
        // Member variables
        //
        double _coreFlux;
        std::vector<double> _BrA; //! dPhi table. Refer to Min et al. [2012, 2013].
        double _dPhi; //! phi grid resolution calculated using _nPhi.
        double _dTheta; //! theta grid resolution calculated using _nTheta.
        long _cntPhi; //!
        long _cntTheta; //!

    public:
        //!
        //! @name Global options.
        //!
        ///@{
        //@{
        //!
        //! Number of phi grids.
        //! Default is 360/5.
        //!
        static void setNPhi (unsigned long const nPhi) {_nPhi = nPhi;};
        static long nPhi () {return _nPhi;};
        //@}

        //@{
        //!
        //! Number of theta grids.
        //! Default is 2*180.
        //!
        static void setNTheta (unsigned long const nTheta) {_nTheta = nTheta;};
        static long nTheta () {return _nTheta;};
        //@}

        //@{
        //!
        //! N concurrent executions for initialization.
        //! Default is 8.
        //!
        static void setNthreads (unsigned long Nthreads) {_Nthreads = Nthreads;};
        static unsigned long Nthreads() {return _Nthreads;};
        //@}
        ///@}

        //!
        //! @name Property accessors
        //!
        ///@{
        //!
        //! Calculated Phi0, \f$|\oint A_r d\phi|\f$.
        //!
        double const& coreFlux () const {return _coreFlux;}; // Phi0
        ///@}

    public:
        virtual ~MagneticFlux ();
        //!
        //! Constructor requires field model.
        //!
        MagneticFlux (FieldModel const* fm);

        //!
        //! Flux calculation (Phi).
        //! L* = this->coreFlux()/Phi.
        //! @note A foot point should be within 0<=theta<=pi, 0<=phi<=2pi.
        //! @param [in] Array of foot points.
        //! @param [in] The number of elements of the array.
        //!
        virtual double fluxClosedByFootPoints_count (Point const* footPoints, long count) const;

    private:
        //
        // Initiate preparation.
        //
        void evaluateUsingFieldModel (FieldModel const* fm);
    };
}

#endif
