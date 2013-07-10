//
//  FieldModel.cpp
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

#include "FieldModel.h"
#include <cmath>
#include <iostream>
#include <cassert>

namespace UBK {

    //
    // Analytic compressed dipole field
    // References: Elkington et al. 2003; Kabin et al. 2007
    //
    
#ifdef DEBUG
    //
    // Test
    //
    void CompressedDipole::test()
    {
        double IMFb1 = 0.;
        double compb2 = 4.;
        FieldModel *fm = new CompressedDipole(IMFb1, compb2);

        Point pt(6., 0., M_PI_2);
        Point b;
        fm->getSphericalFieldInSM_atSphericalPoint(&b, pt);
        std::cout << "Bx=" << b.x << ", By=" << b.y << ", Bz=" << b.z << std::endl;

        delete fm;

        assert(0);
    }

    //
    // Constructor
    //
    CompressedDipole::CompressedDipole (double const IMFb1, double const compb2, double const B0) : _IMFb1(IMFb1), _compb2(compb2), _B0(B0)
    {
        assert(IMFb1 >= 0.);
        assert(compb2 >= 0.);
        assert(B0 >= 0.);
    }

    /*! Model field calculation in cartesian coordinate.
     */
    void CompressedDipole::getCartesianFieldInGSM_atCartesianPoint(Point *bOut, Point const pt) const
    {
        {
            Point sm;
            this->convertGSMCoordinate_toSM(pt, &sm);
            this->getCartesianFieldInSM_atCartesianPoint(bOut, sm);
        }
        this->convertSMCoordinate_toGSM(*bOut, bOut);
    }
    void CompressedDipole::getCartesianFieldInSM_atCartesianPoint(Point *bOut, Point const pt) const
    {
        double rho2 = pt.x*pt.x + pt.y*pt.y;
        double rho = sqrt(rho2);
        double r2 = rho2 + pt.z*pt.z;
        double i_r2 = 1./r2;
        double i_r = sqrt(i_r2);
        double i_r3 = i_r2*i_r;
        double costh = pt.z * i_r;
        double sinth = rho * i_r;
        double sinph, cosph;
        if (rho<1e-10) {
            sinph = 0.;
            cosph = 1.;
        } else {
            sinph = pt.y / rho;
            cosph = pt.x / rho;
        }

        bOut->r = - (2.*_B0*i_r3 - _IMFb1*(1. + _compb2*cosph)) * costh;
        bOut->theta = - (_B0*i_r3 + _IMFb1*(1. + _compb2*cosph)) * sinth;

        bOut->y = bOut->r*sinth + bOut->theta*costh; // Temp for Brho
        bOut->z = bOut->r*costh - bOut->theta*sinth;
        bOut->x = bOut->y*cosph;
        bOut->y = bOut->y*sinph;
    }

    /*! Model field calculation in spherical coordinate.
     */
    void CompressedDipole::getSphericalFieldInGSM_atSphericalPoint(Point *bOut, Point const pt) const
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
    void CompressedDipole::getSphericalFieldInSM_atSphericalPoint(Point *bOut, Point const pt) const
    {
        double r3 = pt.r*pt.r*pt.r;
        double i_r3 = 1./r3;
        double costh = cos(pt.theta);
        double sinth = sin(pt.theta);
        double cosph = cos(pt.phi);

        bOut->r = - (2.*_B0*i_r3 - _IMFb1*(1. + _compb2*cosph)) * costh;
        bOut->theta = - (_B0*i_r3 + _IMFb1*(1. + _compb2*cosph)) * sinth;
        bOut->phi = 0.;
    }

    /*! GSM<->SM coordinate transform
     */
    void CompressedDipole::convertGSMCoordinate_toSM(Point const gsm, Point *smOut) const
    {
        smOut->x = gsm.x*this->cospsi() - gsm.z*this->sinpsi();
        smOut->y = gsm.y;
        smOut->z = gsm.x*this->sinpsi() + gsm.z*this->cospsi();
    }
    void CompressedDipole::convertSMCoordinate_toGSM(Point const sm, Point *gsmOut) const
    {
        gsmOut->x = sm.x*this->cospsi() + sm.z*this->sinpsi();
        gsmOut->y = sm.y;
        gsmOut->z = sm.z*this->cospsi() - sm.x*this->sinpsi();
    }

    double CompressedDipole::cospsi() const
    {
        return cos(this->dipoleTiltAngleRadian());
    }
    double CompressedDipole::sinpsi() const
    {
        return sin(this->dipoleTiltAngleRadian());
    }
#endif

}