//
//  TSFieldModel.cpp
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

#include "TSFieldModel.h"
#include "Geopack.h"
#include "T89.h"
#include "T96.h"
#include "T02.h"
#include "TS05.h"
#include <cmath>
#include <cassert>
#include <cstdio>
#include <stdexcept>
#include <algorithm>

namespace UBK {
    using namespace std;

    //
    // Test
    //
#ifdef DEBUG
    void TSFieldModel::test()
    {
        Date date(2007, 1, 12, 00, 00);
        double vsw = 400;
        int iopt = 1;
        double parmod[10] = {2, -10, 3, -10, };

        {
            TSFieldModel fm(date, vsw, kGeopackDipoleField, iopt, parmod, kTS96Model);

            Point pt(6.);
            Point b;
            fm.getCartesianFieldInGSM_atCartesianPoint(&b, pt);
            printf("B(%s) = %s\n", pt.desc().c_str(), b.desc().c_str());
        }

        assert(0);
    }
#endif

    //
    // None external
    //
    class TSNoneExtern : public TSFieldComponent {
    public:
        void getFieldInGSW_atPoint (Point *bOut, Point const ptgsw) const {
            fill(bOut->xyz, bOut->xyz+3, 0.);
        };
    };

    TSFieldModel::~TSFieldModel ()
    {
        delete _geopack;
        delete _external;
    }
    TSFieldModel::TSFieldModel (Date date, double vsw, GeopackInternalFieldModel intModel, int iopt, double const parmod[], TSExternalFieldModel extModel)
    {
        _geopack = new Geopack(date.year, date.doy, date.hour, date.min, date.sec, vsw, intModel);
        switch (extModel) {
            case kTSNone:
                _external = new TSNoneExtern();
                break;
            case kTS89Model:
                _external = new T89(_geopack, iopt);
                break;
            case kTS96Model:
                _external = new T96(_geopack, parmod);
                break;
            case kTS02Model:
                _external = new T02(_geopack, parmod);
                break;
            case kTS05Model:
                _external = new TS05(_geopack, parmod);
                break;
            default:
                throw invalid_argument("Invalid external model name.");
                break;
        }
    }
    TSFieldModel::TSFieldModel (Date date, double vsw, GeopackInternalFieldModel intModel, int iopt, Parmod parmod, TSExternalFieldModel extModel)
    {
        _geopack = new Geopack(date.year, date.doy, date.hour, date.min, date.sec, vsw, intModel);
        switch (extModel) {
            case kTSNone:
                _external = new TSNoneExtern();
                break;
            case kTS89Model:
                _external = new T89(_geopack, iopt);
                break;
            case kTS96Model:
                _external = new T96(_geopack, parmod.parmod);
                break;
            case kTS02Model:
                _external = new T02(_geopack, parmod.parmod);
                break;
            case kTS05Model:
                _external = new TS05(_geopack, parmod.parmod);
                break;
            default:
                throw invalid_argument("Invalid external model name.");
                break;
        }
    }

    /*! Earth's dipole tilt angle that the model has calculated. */
    double TSFieldModel::dipoleTiltAngleRadian () const
    {
        return _geopack->psi();
    }

    /*! Model field calculation in cartesian coordinate.
     */
    void TSFieldModel::getCartesianFieldInGSM_atCartesianPoint (Point *bOut, Point const pt) const
    {
        _geopack->getFieldInGSW_atPoint(bOut, pt);
        Point h;
        _external->getFieldInGSW_atPoint(&h, pt);

        bOut->x += h.x;
        bOut->y += h.y;
        bOut->z += h.z;
    }
    void TSFieldModel::getCartesianFieldInSM_atCartesianPoint (Point *bOut, Point const pt) const
    {
        {
            Point gsm;
            this->convertSMCoordinate_toGSM(pt, &gsm);
            this->getCartesianFieldInGSM_atCartesianPoint(bOut, gsm);
        }
        this->convertGSMCoordinate_toSM(*bOut, bOut);
    }

    /*! Model field calculation in spherical coordinate.
     */
    void TSFieldModel::getSphericalFieldInGSM_atSphericalPoint (Point *bOut, Point const pt) const
    {
        Point ptCart;
        double costh, sinth, cosph, sinph;
        ptCart.z = pt.r * (costh=cos(pt.theta));
        ptCart.y = pt.r * (sinth=sin(pt.theta)); // Temp for rho
        ptCart.x = ptCart.y * (cosph=cos(pt.phi));
        ptCart.y = ptCart.y * (sinph=sin(pt.phi));

        Point bCart;
        this->getCartesianFieldInGSM_atCartesianPoint(&bCart, ptCart);

        bOut->r = bCart.x*sinth*cosph + bCart.y*sinth*sinph + bCart.z*costh;
        bOut->theta = bCart.x*costh*cosph + bCart.y*costh*sinph - bCart.z*sinth;
        bOut->phi = bCart.y*cosph - bCart.x*sinph;
    }
    void TSFieldModel::getSphericalFieldInSM_atSphericalPoint (Point *bOut, Point const pt) const
    {
        Point ptCart;
        double costh, sinth, cosph, sinph;
        ptCart.z = pt.r * (costh=cos(pt.theta));
        ptCart.y = pt.r * (sinth=sin(pt.theta)); // Temp for rho
        ptCart.x = ptCart.y * (cosph=cos(pt.phi));
        ptCart.y = ptCart.y * (sinph=sin(pt.phi));

        Point bCart;
        this->getCartesianFieldInSM_atCartesianPoint(&bCart, ptCart);

        bOut->r = bCart.x*sinth*cosph + bCart.y*sinth*sinph + bCart.z*costh;
        bOut->theta = bCart.x*costh*cosph + bCart.y*costh*sinph - bCart.z*sinth;
        bOut->phi = bCart.y*cosph - bCart.x*sinph;
    }

    /*! GSM<->SM coordinate transform
     */
    void TSFieldModel::convertGSMCoordinate_toSM (Point const gsm, Point *smOut) const
    {
        smOut->x = gsm.x*_geopack->cospsi() - gsm.z*_geopack->sinpsi();
        smOut->y = gsm.y;
        smOut->z = gsm.x*_geopack->sinpsi() + gsm.z*_geopack->cospsi();
    }
    void TSFieldModel::convertSMCoordinate_toGSM (Point const sm, Point *gsmOut) const
    {
        gsmOut->x = sm.x*_geopack->cospsi() + sm.z*_geopack->sinpsi();
        gsmOut->y = sm.y;
        gsmOut->z = sm.z*_geopack->cospsi() - sm.x*_geopack->sinpsi();
    }

}