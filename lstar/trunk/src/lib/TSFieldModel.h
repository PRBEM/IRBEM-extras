//
//  TSFieldModel.h
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

#ifndef UBJDevelopment_TSFieldModel_h
#define UBJDevelopment_TSFieldModel_h

#include "Definitions.h"
#include "FieldModel.h"
#include "Geopack.h"
#include "TSExternalField.h"

//
// Tsyganenko field model with geopack.
// References: Tsyganenko 1987, 1989, 1995, 2002a, 2002b; Tsyganenko and Sitnov 2005;
//
namespace UBK {

    //!
    //! Convenient date representation.
    //!
    struct Date {
        int year, doy, hour, min, sec;
        Date (int y, int d, int h, int m, int s) : year(y), doy(d), hour(h), min(m), sec(s) {};
    };
    
    //!
    //! Parameters for external field.
    //! Refer to the references.
    //!
    union Parmod {
        double parmod[10];
        struct {
            double Pdyn, DST, ByIMF, BzIMF, W1, W2, W3, W4, W5, W6;
        };
        Parmod (double pdyn, double dst, double byimf, double bzimf, double w1, double w2, double w3, double w4, double w5, double w6) : Pdyn(pdyn), DST(dst), ByIMF(byimf), BzIMF(bzimf), W1(w1), W2(w2), W3(w3), W4(w4), W5(w5), W6(w6) {};
    };

    //!
    //! Tsyganenko field model with geopack.
    //!
    //! References: Tsyganenko 1987, 1989, 1995, 2002a, 2002b; Tsyganenko and Sitnov 2005;
    //!
    class TSFieldModel : public FieldModel {
#ifdef DEBUG
    public:
        static void test();
#endif

    private:
        Geopack const* _geopack;
        TSFieldComponent const* _external;

    private:
        //
        // Prevent copying
        //
        TSFieldModel(TSFieldModel const& fm);
        TSFieldModel& operator = (TSFieldModel const& fm);

    public:
        ~TSFieldModel ();
        //@{
        //!
        //! Contructor requires Date, solar wind speed, internal model flag, iopt (for T89) or parmod[10] (for others), and external model flag.
        //!
        TSFieldModel (Date date, double vsw, GeopackInternalFieldModel intModel, int iopt, double const parmod[], TSExternalFieldModel extModel);
        TSFieldModel (Date date, double vsw, GeopackInternalFieldModel intModel, int iopt, Parmod parmod, TSExternalFieldModel extModel);
        //@}

        //
        // Adopted from FieldModel
        //
        double dipoleTiltAngleRadian () const;

        void getCartesianFieldInGSM_atCartesianPoint (Point *bOut, Point const pt) const;
        void getCartesianFieldInSM_atCartesianPoint (Point *bOut, Point const pt) const;

        void getSphericalFieldInGSM_atSphericalPoint (Point *bOut, Point const pt) const;
        void getSphericalFieldInSM_atSphericalPoint (Point *bOut, Point const pt) const;

        void convertGSMCoordinate_toSM (Point const gsm, Point *smOut) const;
        void convertSMCoordinate_toGSM (Point const sm, Point *gsmOut) const;
    };

}

#endif

