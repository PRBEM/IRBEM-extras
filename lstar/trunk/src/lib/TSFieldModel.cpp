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
#include "TS07.h"
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
        // Hard coded model parameter for specific time moment. THIS IS AN EXAMPLE!!
        double const A[101] = {1.00000,-17.2196,-43.8456,-73.3322,-97.5089,-119.160,44.1685,-28.3264,-10.0919,-6.92003,-9.75857,15.2687,-34.5620,-4.25799,-24.7348,-2.83235,-0.496528,-24.5300,-24.1017,17.9590,18.7953,10.5520,11.1135,-15.2804,-0.598420,4.85943,-10.3045,-6.15506,13.5167,-3.53762,-5.13419,5.62481,11.6524,4.45178,14.5666,-19.2603,-11.6760,5.27181,16.2248,-4.03140,-10.1715,-9.40554,-8.63108,19.2596,8.45173,5.43871,-7.09198,-5.21992,5.06120,8.87807,15.7564,46.7266,-48.0707,-8.33984,1.25046,-23.5229,0.877495,-31.7131,4.87073,-8.23912,1.57636,-33.3436,-21.1656,-3.80171,-14.3754,6.42494,-2.94579,-18.7112,37.4036,10.8904,-1.33068,6.50348,1.51085,-6.63102,4.40149,-2.77839,-3.40469,-4.87805,17.3867,-7.80144,0.872749,9.48005,5.42796,-2.30417,17.9994,0.485395,-8.85228,1.79156,-9.94661,0.396798E-01,7.30606,0.421574,0.575977E-01,-0.246429,0.678498E-01,2.63071,9.15388,32.1001,0.820630,1.71374,0.141481E-01};

        Date date(2007, 1, 1, 1, 1);
        double vsw = 400;
        int iopt = 2;
        Parmod parmod(1., A);

        {
            TSFieldModel fm(date, vsw, kGeopackDipoleField, iopt, parmod, kTS07Model);

            Point pt(-6., 0., 0.);
            Point b;
            fm.getCartesianFieldInSM_atCartesianPoint(&b, pt);
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
            case kTS07Model:
                _external = new TS07(_geopack, parmod);
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
            case kTS07Model:
                _external = new TS07(_geopack, parmod.parmod);
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