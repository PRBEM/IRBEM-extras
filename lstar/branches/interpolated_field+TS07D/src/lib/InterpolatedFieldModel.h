//
//  InterpolatedFieldModel.h
//  UBJDevelopment
//
//  Created by KYUNGGUK MIN on 7/10/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

#ifndef UBJDevelopment_InterpolatedFieldModel_h
#define UBJDevelopment_InterpolatedFieldModel_h

#include "FieldModel.h"
#include "Key.h"
#include <vector>
#include <cmath>

namespace UBK {

    namespace __helper {
        struct HashTable;
        struct F;
    }

    //!
    //! Interpolation order flag.
    //!
    enum _PolyInterpolationOrder {
        k1st = 1, //!< Linear.
        k2nd //!< 2nd order.
    };
    typedef long PolyInterpolationOrder;

    //!
    //! Field model abstract.
    //!
    class InterpolatedFieldModel : public FieldModel {
#ifdef DEBUG
    public:
        static void test ();
#endif

    public:
        //!
        //! Field vector insertion failure exception.
        //! Critical error.
        //!
        struct hash_table_insertion_failure_exception {};
        //!
        //! Thrown when the input interpolation order is invalid.
        //!
        struct wrong_interpolation_order_exception {};

    public:
        //
        // Grid spacing and size to determine the table size.
        //
        static const double dx;
        static const double dy;
        static const double dz;
        static const signed ixmin;
        static const signed iymin;
        static const signed izmin;
        static const signed ixmax;
        static const signed iymax;
        static const signed izmax;

    private:
        PolyInterpolationOrder _order;
        __helper::F const* _f_ptr;
        std::vector<__helper::HashTable*> _table;
        FieldModel const& _field_model;

    public:
        //!
        //! Get interpolation order of the object.
        //!
        PolyInterpolationOrder interpolationOrder() const {return _order;};
        //!
        //! Get the underlying field model.
        //!
        FieldModel const& fieldModel() const {return _field_model;};

    public:
        virtual ~InterpolatedFieldModel ();
        InterpolatedFieldModel(FieldModel const& field_model, PolyInterpolationOrder order = k2nd);

        //!
        //! Indexed access to the tabulated field vector.
        //! Caller should calculate the proper index based on the grid spacing and table size.
        //!
        Point const& operator[](Key const& key) const;

        //!
        //! Earth's dipole tilt angle in radian that the model has calculated.
        //!
        virtual double dipoleTiltAngleRadian () const {return _field_model.dipoleTiltAngleRadian();};
        //!
        //! Interpolated model field in cartesian GSM coordinate.
        //!
        virtual void getCartesianFieldInGSM_atCartesianPoint (Point *bOut, Point pt) const;
        //!
        //! Interpolated model field in cartesian SM coordinate.
        //!
        virtual void getCartesianFieldInSM_atCartesianPoint (Point *bOut, Point const pt) const;

        //!
        //! Interpolated model field in spherical GSM coordinate.
        //!
        virtual void getSphericalFieldInGSM_atSphericalPoint (Point *bOut, Point const pt) const;
        //!
        //! Interpolated model field in spherical SM coordinate.
        //!
        virtual void getSphericalFieldInSM_atSphericalPoint (Point *bOut, Point const pt) const;

        //!
        //! GSM->SM cartesian coordinate transform.
        //!
        virtual void convertGSMCoordinate_toSM (Point const gsm, Point *smOut) const;
        //!
        //! SM->GSM cartesian coordinate transform.
        //!
        virtual void convertSMCoordinate_toGSM (Point const sm, Point *gsmOut) const;

    private:
        double cospsi() const {
            return std::cos(this->dipoleTiltAngleRadian());
        };
        double sinpsi() const {
            return std::sin(this->dipoleTiltAngleRadian());
        };
    };

}

#endif
