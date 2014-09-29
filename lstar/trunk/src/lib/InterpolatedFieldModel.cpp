//
//  InterpolatedFieldModel.cpp
//  UBJDevelopment
//
//  Created by KYUNGGUK MIN on 7/10/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

#include "InterpolatedFieldModel.h"
#include "Thread.h"
#include <map>
#include <stdexcept>
#include <cassert>
#include <iostream>
#include "TSFieldModel.h"

using namespace std;

namespace UBK {

#ifdef DEBUG
    void InterpolatedFieldModel::test()
    {
        Date date(2007, 1, 1, 1, 1);
        double vsw = 400;
        int iopt = 2;
        Parmod parmod(1., 10., -5, 1, 0, 0, 0, 0, 0, 0);

        TSFieldModel fm(date, vsw, kGeopackDipoleField, iopt, parmod, kTS05Model);
        InterpolatedFieldModel ifm(fm);

        Point gsm(6.2, 0., -9.);
        Point ba, bi;
        ifm.getCartesianFieldInGSM_atCartesianPoint(&bi, gsm);
        fm.getCartesianFieldInGSM_atCartesianPoint(&ba, gsm);
        cout << "Ba={" << ba.x << ", " << ba.y << ", " << ba.z << "}, ";
        cout << "Bi={" << bi.x << ", " << bi.y << ", " << bi.z << "}, " << endl;

        assert(0);
    }
#endif

    ////////////////////////////////////////////////////////////
    // InterpolatedFieldModel static constants
    ////////////////////////////////////////////////////////////
    const double InterpolatedFieldModel::dx = .1;
    const double InterpolatedFieldModel::dy = .1;
    const double InterpolatedFieldModel::dz = .1;
    const signed InterpolatedFieldModel::ixmin = -150;
    const signed InterpolatedFieldModel::iymin = -150;
    const signed InterpolatedFieldModel::izmin = -100;
    const signed InterpolatedFieldModel::ixmax = 150;
    const signed InterpolatedFieldModel::iymax = 150;
    const signed InterpolatedFieldModel::izmax = 100;

    ////////////////////////////////////////////////////////////
    // Interpolation helpers
    ////////////////////////////////////////////////////////////
    namespace __helper {

        static const unsigned ND = 3; // 3D

        typedef double Tx;
        typedef Point Ty;
        typedef Point T;

        struct F {
            InterpolatedFieldModel const& _field_model;
            virtual ~F() {};
            F(InterpolatedFieldModel const& field_model) : _field_model(field_model) {};
            virtual T operator()(T const& pt) const {
                Key key(round(pt.x),
                        round(pt.y),
                        round(pt.z));
                return _field_model[key];
            };
        };

        namespace __1st {
            inline Ty poly_interp (Tx const dx[], Ty y[]) {
                Point yi;
                yi.x =
                (dx[1] * y[0].x * (-1.)) +
                (dx[0] * y[1].x * ( 1.));
                yi.y =
                (dx[1] * y[0].y * (-1.)) +
                (dx[0] * y[1].y * ( 1.));
                yi.z =
                (dx[1] * y[0].z * (-1.)) +
                (dx[0] * y[1].z * ( 1.));
                return yi;
            }

            inline T poly3D_interpZ (T const& pt, F const& f) {
                Tx z1 = round(pt.z);
                Tx z0 = z1 - 1.;

                T Y[ND];
                Y[0] = f(T(pt.x,pt.y,z0));
                Y[1] = f(T(pt.x,pt.y,z1));

                Tx dz[ND];
                dz[0] = pt.z-z0;
                dz[1] = pt.z-z1;
                return poly_interp(dz, Y);
            }

            inline T poly3D_interpY (T const& pt, F const& f) {
                Tx y1 = round(pt.y);
                Tx y0 = y1 - 1.;

                T Y[ND];
                Y[0] = poly3D_interpZ(T(pt.x,y0,pt.z), f);
                Y[1] = poly3D_interpZ(T(pt.x,y1,pt.z), f);

                Tx dy[ND];
                dy[0] = pt.y-y0;
                dy[1] = pt.y-y1;
                return poly_interp(dy, Y);
            }

            inline T poly3D_interp (T const& pt, F const& f) {
                Tx x1 = round(pt.x);
                Tx x0 = x1 - 1.;

                T Y[ND];
                Y[0] = poly3D_interpY(T(x0,pt.y,pt.z), f);
                Y[1] = poly3D_interpY(T(x1,pt.y,pt.z), f);

                Tx dx[ND];
                dx[0] = pt.x-x0;
                dx[1] = pt.x-x1;
                return poly_interp(dx, Y);
            }
        }

        namespace __2nd {
            inline Ty poly_interp (Tx const dx[], Ty y[]) {
                Point yi;
                yi.x =
                (dx[1]*dx[2] * y[0].x * (0.5)) +
                (dx[0]*dx[2] * y[1].x * (-1.)) +
                (dx[0]*dx[1] * y[2].x * (0.5));
                yi.y =
                (dx[1]*dx[2] * y[0].y * (0.5)) +
                (dx[0]*dx[2] * y[1].y * (-1.)) +
                (dx[0]*dx[1] * y[2].y * (0.5));
                yi.z =
                (dx[1]*dx[2] * y[0].z * (0.5)) +
                (dx[0]*dx[2] * y[1].z * (-1.)) +
                (dx[0]*dx[1] * y[2].z * (0.5));
                return yi;
            }

            inline T poly3D_interpZ (T const& pt, F const& f) {
                Tx z1 = round(pt.z);
                Tx z0 = z1 - 1.;
                Tx z2 = z1 + 1.;

                T Y[ND];
                Y[0] = f(T(pt.x,pt.y,z0));
                Y[1] = f(T(pt.x,pt.y,z1));
                Y[2] = f(T(pt.x,pt.y,z2));

                Tx dz[ND];
                dz[0] = pt.z-z0;
                dz[1] = pt.z-z1;
                dz[2] = pt.z-z2;
                return poly_interp(dz, Y);
            }

            inline T poly3D_interpY (T const& pt, F const& f) {
                Tx y1 = round(pt.y);
                Tx y0 = y1 - 1.;
                Tx y2 = y1 + 1.;

                T Y[ND];
                Y[0] = poly3D_interpZ(T(pt.x,y0,pt.z), f);
                Y[1] = poly3D_interpZ(T(pt.x,y1,pt.z), f);
                Y[2] = poly3D_interpZ(T(pt.x,y2,pt.z), f);

                Tx dy[ND];
                dy[0] = pt.y-y0;
                dy[1] = pt.y-y1;
                dy[2] = pt.y-y2;
                return poly_interp(dy, Y);
            }

            inline T poly3D_interp (T const& pt, F const& f) {
                Tx x1 = round(pt.x);
                Tx x0 = x1 - 1.;
                Tx x2 = x1 + 1.;

                T Y[ND];
                Y[0] = poly3D_interpY(T(x0,pt.y,pt.z), f);
                Y[1] = poly3D_interpY(T(x1,pt.y,pt.z), f);
                Y[2] = poly3D_interpY(T(x2,pt.y,pt.z), f);
                
                Tx dx[ND];
                dx[0] = pt.x-x0;
                dx[1] = pt.x-x1;
                dx[2] = pt.x-x2;
                return poly_interp(dx, Y);
            }
        }
    }

    ////////////////////////////////////////////////////////////
    // Table proxy
    ////////////////////////////////////////////////////////////
    namespace __helper {
        struct HashTable {
            mutable Mutex _key;
            mutable map<int, Point> _hash_table;
            FieldModel const& _field_model;

            Point const& operator[](Key const& key) const {Locker _(_key);
                map<int, Point>::const_iterator it;
                int ind = key.x;
                if (_hash_table.end() == (it=_hash_table.find(ind))) {
                    Point pt(key.x*InterpolatedFieldModel::dx, key.y*InterpolatedFieldModel::dy, key.z*InterpolatedFieldModel::dz);
                    _field_model.getCartesianFieldInGSM_atCartesianPoint(&pt, pt);
                    std::pair<map<int, Point>::iterator, bool> result = _hash_table.insert(std::make_pair(ind, pt));
                    if (!result.second) {
                        throw InterpolatedFieldModel::hash_table_insertion_failure_exception();
                    }
                    it = result.first;
                }
                return it->second;
            };

            ~HashTable() {};

            HashTable(FieldModel const& model) : _field_model(model) {};
        };
    }

    ////////////////////////////////////////////////////////////
    // InterpolatedFieldModel implementation
    ////////////////////////////////////////////////////////////
    InterpolatedFieldModel::~InterpolatedFieldModel ()
    {
        for (long idx=0, cnt=_table.size(); idx<cnt; idx++) {
            delete _table.at(idx);
        }
        delete _f_ptr;
    }

    InterpolatedFieldModel::InterpolatedFieldModel(FieldModel const& field_model, PolyInterpolationOrder order) : _field_model(field_model), _order(order)
    {
        if (k1st!=_order && k2nd!=_order) {
            throw wrong_interpolation_order_exception();
        }

        for (long idz=InterpolatedFieldModel::izmin; idz<=InterpolatedFieldModel::izmax; idz++) {
            for (long idy=InterpolatedFieldModel::iymin; idy<=InterpolatedFieldModel::iymax; idy++) {
                _table.push_back(new __helper::HashTable(_field_model));
            }
        }

        _f_ptr = new __helper::F(*this);
    }

    Point const& InterpolatedFieldModel::operator[](Key const& key) const
    {
        if (
            (InterpolatedFieldModel::ixmin>key.x || InterpolatedFieldModel::ixmax<key.x) ||
            (InterpolatedFieldModel::iymin>key.y || InterpolatedFieldModel::iymax<key.y) ||
            (InterpolatedFieldModel::izmin>key.z || InterpolatedFieldModel::izmax<key.z)) {
            throw out_of_range("Key.x out of range.");
        }

        static const long ny = (InterpolatedFieldModel::iymax-InterpolatedFieldModel::iymin+1);
        long ind = (key.z-InterpolatedFieldModel::izmin)*ny + (key.y-InterpolatedFieldModel::iymin);
        return _table.at(ind)->operator[](key);
    }

    void InterpolatedFieldModel::getCartesianFieldInGSM_atCartesianPoint (Point *bOut, Point pt) const
    {
        try {
            pt.x /= InterpolatedFieldModel::dx;
            pt.y /= InterpolatedFieldModel::dy;
            pt.z /= InterpolatedFieldModel::dz;
            if (k1st == _order) {
                *bOut = __helper::__1st::poly3D_interp(pt, *_f_ptr);
            } else {
                *bOut = __helper::__2nd::poly3D_interp(pt, *_f_ptr);
            }
        } catch (out_of_range&) {
            bOut->x=NAN, bOut->y=NAN, bOut->z=NAN;
        }
    }

    void InterpolatedFieldModel::getCartesianFieldInSM_atCartesianPoint (Point *bOut, Point const pt) const
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
    void InterpolatedFieldModel::getSphericalFieldInGSM_atSphericalPoint (Point *bOut, Point const pt) const
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
    void InterpolatedFieldModel::getSphericalFieldInSM_atSphericalPoint (Point *bOut, Point const pt) const
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
    void InterpolatedFieldModel::convertGSMCoordinate_toSM (Point const gsm, Point *smOut) const
    {
        smOut->x = gsm.x*this->cospsi() - gsm.z*this->sinpsi();
        smOut->y = gsm.y;
        smOut->z = gsm.x*this->sinpsi() + gsm.z*this->cospsi();
    }
    void InterpolatedFieldModel::convertSMCoordinate_toGSM (Point const sm, Point *gsmOut) const
    {
        gsmOut->x = sm.x*this->cospsi() + sm.z*this->sinpsi();
        gsmOut->y = sm.y;
        gsmOut->z = sm.z*this->cospsi() - sm.x*this->sinpsi();
    }

}