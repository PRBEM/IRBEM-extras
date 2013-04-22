//
//  FieldModel.h
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

#ifndef UBJDevelopment_FieldModel_h
#define UBJDevelopment_FieldModel_h

#include "Definitions.h"
#include "Point.h"

namespace UBK {

    //!
    //! Field model abstract.
    //!
    class FieldModel {
    public:
        virtual ~FieldModel () {};

        //!
        //! Earth's dipole tilt angle in radian that the model has calculated.
        //!
        virtual double dipoleTiltAngleRadian () const = 0;

        //!
        //! Model field in cartesian GSM coordinate.
        //!
        virtual void getCartesianFieldInGSM_atCartesianPoint (Point *bOut, Point const pt) const = 0;
        //!
        //! Model field in cartesian SM coordinate.
        //!
        virtual void getCartesianFieldInSM_atCartesianPoint (Point *bOut, Point const pt) const = 0;

        //!
        //! Model field in spherical GSM coordinate.
        //!
        virtual void getSphericalFieldInGSM_atSphericalPoint (Point *bOut, Point const pt) const = 0;
        //!
        //! Model field in spherical SM coordinate.
        //!
        virtual void getSphericalFieldInSM_atSphericalPoint (Point *bOut, Point const pt) const = 0;

        //!
        //! GSM->SM cartesian coordinate transform.
        //!
        virtual void convertGSMCoordinate_toSM (Point const gsm, Point *smOut) const = 0;
        //!
        //! SM->GSM cartesian coordinate transform.
        //!
        virtual void convertSMCoordinate_toGSM (Point const sm, Point *gsmOut) const = 0;
    };

#ifdef DEBUG
    //
    // Analytic compressed dipole field.
    // References: Elkington et al. 2003; Kabin et al. 2007
    //
    class CompressedDipole : public FieldModel {
    public:
        static void test ();

    private:
        double _IMFb1, _compb2, _B0;

    public:
        virtual ~CompressedDipole () {};
        //
        // IMFb1 and compb2 correspond to b1 and b2 of the compressed dipole model, respectively. Refer to the references for the meaning.
        // B0: Earth's dipole moment.
        //
        CompressedDipole (double const IMFb1, double const compb2, double const B0 = 31200.);

        //
        // Always returns 0.
        //
        virtual double dipoleTiltAngleRadian () const {return 0.;};

        virtual void getCartesianFieldInGSM_atCartesianPoint (Point *bOut, Point const pt) const;
        virtual void getCartesianFieldInSM_atCartesianPoint (Point *bOut, Point const pt) const;

        virtual void getSphericalFieldInGSM_atSphericalPoint (Point *bOut, Point const pt) const;
        virtual void getSphericalFieldInSM_atSphericalPoint (Point *bOut, Point const pt) const;

        virtual void convertGSMCoordinate_toSM (Point const gsm, Point *smOut) const;
        virtual void convertSMCoordinate_toGSM (Point const sm, Point *gsmOut) const;

    private:
        double cospsi() const;
        double sinpsi() const;
    };
#endif

}

#endif