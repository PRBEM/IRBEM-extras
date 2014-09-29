//
//  Coordinator.cpp
//  UBJDevelopment
//
//  Created by Kyungguk Min on 4/20/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

#include "Coordinator.h"
#include "UBKLstarxx.h"

namespace UBK {
    using namespace std;

    //
    // Coordinator
    //
    void Coordinator::start()
    {
        ThreadFor<Coordinator &> t(_n_threads, this->count(), *this);
    }

    //
    // Cotrans coordinator
    //
#ifdef DEBUG
    void CotransCoordinator::test()
    {
        Date date(2007, 1, 12, 00, 00);
        double vsw = 400;
        int iopt = 1;
        double parmod[10] = {2, -10, 3, -10, };
        GeopackInternalFieldModel internal = kGeopackDipoleField;
        TSExternalFieldModel external = kTSNone;

        TSFieldModel fm(date, vsw, internal, iopt, parmod, external);
        {
            class _ : public CotransCoordinator {
            public:
                _(FieldModel const& fm) : CotransCoordinator(fm) {};

                unsigned long count() const {return 1;};
                Point from(long idx) const {return Point(0., 6.);};
                void callback(long idx, Point const& pt) const {
                    fprintf(stderr, "idx = %ld, pt = %s\n", idx, pt.desc().c_str());
                };
            };

            _ f(fm);
            f.setToCoSystem(kMagneticFieldSM);
            f.start();
        }

        assert(0);
    }
#endif

    void CotransCoordinator::operator()(long idx)
    {
        Point to;

        switch (_to_co_system) {
            case kMagneticFieldGSM:
                this->fieldModel().convertSMCoordinate_toGSM(this->from(idx), &to);
                break;
            case kMagneticFieldSM:
                this->fieldModel().convertGSMCoordinate_toSM(this->from(idx), &to);
                break;
            default:
                throw std::invalid_argument("Invalid to_co_system parameter.");
                break;
        }

        this->callback(idx, to);
    }

    void CotransCoordinator::setToCoSystem(MagneticFieldCoordinateSystem to_co_system) {
        if (kMagneticFieldGSM!=to_co_system && kMagneticFieldSM!=to_co_system) {
            throw std::invalid_argument("to_co_system should be either kMagneticFieldGSM or kMagneticFieldSM.");
        }
        _to_co_system = to_co_system;
    }

    //
    // Model magnetic field coordinator.
    //
#ifdef DEBUG
    void BFieldCoordinator::test()
    {
        Date date(2001, 1, 1, 0, 0);
        double vsw = 400;
        int iopt = 1;
        double parmod[10] = {2, -10, 3, -10, };
        GeopackInternalFieldModel internal = kGeopackIGRFField;
        TSExternalFieldModel external = kTS05Model;

        TSFieldModel fm(date, vsw, internal, iopt, parmod, external);
        {
            class _ : public BFieldCoordinator {
            public:
                _(FieldModel const& fm) : BFieldCoordinator(fm) {};

                unsigned long count() const {return 1;};
                Point at(long idx) const {return Point(6.);};
                void callback(long idx, Point const& B) const {
                    fprintf(stderr, "idx = %ld, pt = %s\n", idx, B.desc().c_str());
                };
            };

            _ f(fm);
            f.setCoSystem(kMagneticFieldSM);
            f.start();
        }
        
        assert(0);
    }
#endif

    void BFieldCoordinator::operator()(long idx)
    {
        Point B;

        switch (_co_system) {
            case kMagneticFieldGSM:
                this->fieldModel().getCartesianFieldInGSM_atCartesianPoint(&B, this->at(idx));
                break;
            case kMagneticFieldSM:
                this->fieldModel().getCartesianFieldInSM_atCartesianPoint(&B, this->at(idx));
                break;
            default:
                throw std::invalid_argument("Invalid co_system parameter.");
                break;
        }

        this->callback(idx, B);
    }

    void BFieldCoordinator::setCoSystem(MagneticFieldCoordinateSystem co_system) {
        if (kMagneticFieldGSM!=co_system && kMagneticFieldSM!=co_system) {
            throw std::invalid_argument("co_system should be either kMagneticFieldGSM or kMagneticFieldSM.");
        }
        _co_system = co_system;
    }

    //
    // Field line coordinator abstract.
    //
#ifdef DEBUG
    void FieldLineCoordinator::test()
    {
        Date date(2001, 1, 1, 0, 0);
        double vsw = 400;
        int iopt = 1;
        double parmod[10] = {2, -10, 3, -10, 0, 0, 0, 0, 0, 0};
        GeopackInternalFieldModel internal = kGeopackIGRFField;
        TSExternalFieldModel external = kTS05Model;

        TSFieldModel fm_tmp(date, vsw, internal, iopt, parmod, external);
        InterpolatedFieldModel fm(fm_tmp);
        {
            class _ : public FieldLineCoordinator {
            public:
                _(FieldModel const& fm) : FieldLineCoordinator(fm) {};

                unsigned long count() const {return 1;};
                Point from(long idx) const {return Point(-6.);};
                void callback(long idx, FieldLine const& fl) const {
                    fprintf(stderr, "idx = %ld, XYZmeq = %s, Bmeq = %f, Rflend = %s\n", idx, fl.magneticEquator().desc().c_str(), fl.magneticEquatorFieldMagnitude(), fl.coordinates().back().desc().c_str());
                    fprintf(stderr, "%ld, %ld\n", fl.modifiedInvariants().size(), fl.coordinates().size());
                };
            };

            _ f(fm);
            f.setDs(.05);
            f.setIonoR(1.);
            f.start();
        }

        assert(0);
    }
#endif

    void FieldLineCoordinator::setIonoR(double ionoR) const {
        if (ionoR < 1.) {
            throw std::invalid_argument("ionoR < 1 RE.");
        }
        FieldLine::setRadiusMin(ionoR);
    }
    void FieldLineCoordinator::setDs(double ds) const {
        if (ds <= 0.) {
            throw std::invalid_argument("ds <= 0 RE.");
        }
        FieldLine::setStepSize(ds);
    }

    void FieldLineCoordinator::operator()(long idx) {
        FieldLine fl(this->from(idx), &this->fieldModel());
        this->callback(idx, fl);
    }

    //
    // Cartesian L* coordinator abstract.
    //
#ifdef DEBUG
    void LstarCoordinator::test()
    {
        Date date(2001, 1, 1, 0, 0);
        double vsw = 400;
        int iopt = 2;
        double parmod[10] = {2, -10, 3, -10, };
        GeopackInternalFieldModel internal = kGeopackDipoleField;
        TSExternalFieldModel external = kTSNone;

        TSFieldModel fm_tmp(date, vsw, internal, iopt, parmod, external);
        InterpolatedFieldModel fm(fm_tmp, k2nd);
        {
            class _ : public LstarCoordinator {
            public:
                _(FieldModel const& fm) : LstarCoordinator(fm, YES) {};

                unsigned long count() const {return 18;};
                Point from(long idx) const {return Point(6.1, 2.5 * M_PI/180., (idx/(this->count()-1.)*80+10.) * M_PI/180.);};
                void callback(long idx, Particle const& ptl) const {
                    fprintf(stderr, "idx = %ld, L* = %f, K = %f\n", idx, ptl.Lstar(), ptl.K());
                };
            };

            _ f(fm);
            f.setNPhi(360/5);
            //f.setD(Point(.1, M_PI*2./MagneticFlux::nPhi()));
            f.setD(Point(.2, 2.*M_PI/MagneticFlux::nPhi()));
            f.setDs(.1);
            f.start();
        }
        
        assert(0);
    }
#endif

    void LstarCoordinator::setIonoR(double ionoR) const {
        if (ionoR < 1.) {
            throw std::invalid_argument("ionoR < 1 RE.");
        }
        FieldLine::setRadiusMin(ionoR);
    }
    void LstarCoordinator::setDs(double ds) const {
        if (ds <= 0.) {
            throw std::invalid_argument("ds <= 0 RE.");
        }
        FieldLine::setStepSize(ds);
    }
    void LstarCoordinator::setNPhi(unsigned long nPhi) const {
        if (nPhi<10) {
            throw std::invalid_argument("nPhi < 10.");
        }
        MagneticFlux::setNPhi(nPhi);
    }
    void LstarCoordinator::setNTheta(unsigned long nTheta) const {
        if (nTheta<10) {
            throw std::invalid_argument("nTheta < 10.");
        }
        MagneticFlux::setNTheta(nTheta);
    }
    void LstarCoordinator::setD(Point d) {
        if (d.x<=0. || d.y<=0.) {
            throw std::invalid_argument("d <= 0 RE.");
        }
        d.z = 1.;
        _t = AffineTransform(d, Point());
    }

    void LstarCoordinator::start() {
        //
        // Limit
        //
        Key l, u;
        if (_isCylindricalGrid) {
            Point lb = _t.pointConvertedToNormal( Point(.9, 0., 0.) );
            Point ub = _t.pointConvertedToNormal( Point(15., M_PI*2., 0.) );

            l = Key(lb.r, lb.phi-2, 0);
            u = Key(ub.r, ub.phi+2, 0);
        } else {
            Point lb = _t.pointConvertedToNormal( Point(-15., -15., 0.) );
            Point ub = _t.pointConvertedToNormal( Point( 15.,  15., 0.) );
            using namespace std;
            l = Key(round(lb.x), round(lb.y), 0);
            u = Key(round(ub.x), round(ub.y), 0);
        }

        //
        // Instantiate magnetic flux
        //
        UBK::MagneticFlux mf(&this->fieldModel());
        _mf = &mf;
        _Phi0 = mf.coreFlux();

        //
        // Instantiate field line table
        //
        if (_isCylindricalGrid) {
            _ftbl = new FieldLineTableCyl<FieldLineTable, FieldLine>(l, u, &this->fieldModel(), &_t);
        } else {
            _ftbl = new FieldLineTableCart<FieldLineTable, FieldLine>(l, u, &this->fieldModel(), &_t);
        }

        Coordinator::start();

        delete _ftbl;
    }
    void LstarCoordinator::operator()(long idx) {
        Point origin = this->from(idx);
        Particle *ptl;
        if (_isCylindricalGrid) {
            ptl = new ParticleCyl<FieldLineTable, FieldLine>(origin, origin.theta, dynamic_cast< FieldLineTableCyl<FieldLineTable, FieldLine>* >(_ftbl), _mf, &_t, YES, YES);
        } else {
            ptl = new ParticleCart<FieldLineTable, FieldLine>(origin, origin.theta, dynamic_cast< FieldLineTableCart<FieldLineTable, FieldLine>* >(_ftbl), _mf, &_t, YES, YES);
        }

        ptl->evaluateLstar();
        this->callback(idx, *ptl);

        delete ptl;
    }

}