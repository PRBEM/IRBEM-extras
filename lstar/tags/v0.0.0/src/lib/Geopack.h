//
//  Geopack.h
//  UBJDevelopment
//
//  Created by Kyungguk Min on 4/5/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#ifndef UBJDevelopment_Geopack_h
#define UBJDevelopment_Geopack_h

#include "TSFieldComponent.h"

//
// Geopack class.
//
namespace UBK {

    //
    // Opaque object. Internal use only.
    //
    typedef void* GeopackContextRef;
    typedef void* GeopackFieldHelper;

    //!
    //! Internal model flag.
    //!
    enum _GeopackInternalFieldModel {
        kGeopackDipoleField = 0, //!< Dipole.
        kGeopackIGRFField //!< IGRF.
    };
    //!
    //! Internal model flag type.
    //!
    typedef long GeopackInternalFieldModel;

    //!
    //! Wrapper to GEOPACK module.
    //!
    //! Reference: http://geo.phys.spbu.ru/~tsyganenko/modeling.html
    //!
    class Geopack : public TSFieldComponent {
#ifdef DEBUG
    public:
        static void test ();
#endif

    private:
        int _year;
        int _doy;
        int _hour;
        int _min;
        int _sec;
        Point _vgse;
        GeopackContextRef _ctx;
        GeopackFieldHelper _helper;

    public:
        //!
        //! @name Accessors.
        //!
        ///@{
        //!
        //! Year.
        //!
        int const& year () const {return _year;};

        //!
        //! Day of year.
        //!
        int const& doy () const {return _doy;};

        //!
        //! Hour.
        //!
        int const& hour () const {return _hour;};

        //!
        //! Minute.
        //!
        int const& min () const {return _min;};

        //!
        //! Second.
        //!
        int const& sec () const {return _sec;};

        //!
        //! Solar wind velocity to determine the GSW coordinate system.
        //!
        Point const& vgse () const {return _vgse;};

        //@{
        //!
        //! Dipole tilt angle in radian.
        //!
        //! @note User can arbitrary set this tilt angle (for experimental purpose only).
        //!
        void setPsi (double const psi);
        double const& psi () const;
        //@}

        //!
        //! sin(psi()).
        //!
        double const& sinpsi () const;
        //!
        //! cos(psi()).
        //!
        double const& cospsi () const;
        ///@}

    public:
        ~Geopack ();
        //@{
        //!
        //! Constructor.
        //!
        Geopack (int year, int doy, int hour, int min, int sec, Point vgse, GeopackInternalFieldModel flag = kGeopackIGRFField);
        Geopack (int year, int doy, int hour, int min, int sec, double vsw, GeopackInternalFieldModel flag = kGeopackIGRFField);
        //@}

        //
        // Copy semantic
        //
        Geopack (Geopack const& g);
        Geopack &operator = (Geopack const& g);

        //!
        //! Update the context.
        //! Automatically called on instantiation.
        //!
        //! @note Call if date changes.
        //!
        void recaluateWithDate_vGSE (int year, int doy, int hour, int min, int sec, Point vgse);

        //
        // Adopted from TSFieldComponent.
        //
        void getFieldInGSW_atPoint (Point *bOut, Point const ptgsw) const;
    };

}

#endif
