//
//  main.cpp
//  UBKLstarxx
//
//  Created by Kyungguk Min on 3/30/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#include "UBKLstarxx.h"
#include <cstdio>
#include <unistd.h>

//
// Example 1. Multi-threading
//
static void multi_threading();

//
// Example 2. High-level interfaces
//
static void high_level();

//
// Example 3. Low-level interface
//
static void low_level();

int main(int argc, const char * argv[])
{
    using namespace std;
    using namespace UBK;

    fprintf(stdout, "%%%% Multi-threading.\n");
    multi_threading();

    fprintf(stdout, "%%%% High level interfaces.\n");
    high_level();

    fprintf(stdout, "%%%% Low level interfaces.\n");
    low_level();

    return 0;
}

/////////////////////////////////////////////////////////////
//
// Implementation
//
/////////////////////////////////////////////////////////////

// 1. Multi-threading
class ATask { // Some task
public:
    void operator()() {
        fprintf(stdout, "A task started.\n");
        sleep(2);
        fprintf(stdout, "A task ended.\n");
    };
};
class Sum { // Loop task
public:
    long sum;
    UBK::Mutex key;

    Sum() : sum(0) {};
    void operator()(long idx) {
        UBK::Locker l(key); // Race condition
        sum += idx+1;
    };
};
static void multi_threading()
{
    using namespace std;
    using namespace UBK;

    // 1.1 Excute a task concurrently
    {
        // Make Thread object and submit the task
        ATask aTask;
        Thread<ATask> t(aTask);
        //
        // Do something...
        //

        // Join the thread
        t.join();
    }

    // 1.2 Parallel for loop
    {
        Sum s;
        long Nthreads = 4;
        long count = 10;
        ThreadFor<Sum &> t(Nthreads, count, s);
        fprintf(stdout, "1 + 2 + ... 10 = %ld.\n", s.sum);
    }
}

// 2. High-level
// 2.1 Extend LstarCoordinator abstract class
class Lstar : public UBK::LstarCoordinator {
public:
    double *Ls;

    Lstar(UBK::FieldModel const& fm, BOOL isCylindricalGrid) : LstarCoordinator(fm,isCylindricalGrid) {};
    unsigned long count() const { // Number of iterations
        return 17;
    };
    UBK::Point from(long idx) const {
        double pa = 80.*idx/(this->count()-1.) + 10;
        return UBK::Point(6,M_PI,pa * M_PI/180.);
    };
    void callback(long idx, UBK::Particle const& ptl) const { // Called by worker thread when its job is done
        Ls[idx] = ptl.Lstar();
        UBK::Point p0 = this->from(idx);
        p0.phi *= 180./M_PI;
        p0.theta *= 180./M_PI;
        fprintf(stdout, "idx = %ld, Initial = %s, L* = %f\n", idx, p0.desc().c_str(), ptl.Lstar());
    };
};
static void high_level()
{
    using namespace std;
    using namespace UBK;

    // Field model
    double const vsw = 400.;
    Date d(2007, 1, 12, 0, 0);
    GeopackInternalFieldModel const internal = kGeopackDipoleField;
    TSExternalFieldModel const external = kTS89Model;
    int iopt = 3;
    TSFieldModel fm(d, vsw, internal, iopt, NULL, external);

    // Instantiate Lstar
    BOOL isCylindrical = YES;
    Lstar lstar(fm, isCylindrical);
    // Some settings ...

    // Start
    double *Ls = new double[lstar.count()];
    lstar.Ls = Ls;
    lstar.start();
    delete Ls;
}

// 3. Low-level
static void low_level()
{
    using namespace std;
    using namespace UBK;

    // Field model. Users can supply their own model by inheriting FieldModel abstract class.
    TSFieldModel fm(Date(2001, 1, 1, 0, 0), 400, kGeopackDipoleField, 3, Parmod(3, -10, 3, -10, 0, 0, 0, 0, 0, 0), kTS89Model);

    // Global options
    FieldLine::setStepSize(.05);
    MagneticFlux::setNPhi(360/2);
    MagneticFlux mf(&fm);
    printf("F0 = %f\n", mf.coreFlux());

    // Grid resolution.
    double dr = .2;
    double dphi = 2.*M_PI/MagneticFlux::nPhi();
    AffineTransform tr(Point(dr, dphi, 1.), Point(0., 0., 0.));

    // Index space bound. .9<=R<=12RE, and phi is of course from 0 to 2pi.
    Point rl(.9, 0.); // Lower bound.
    Point ru(12., M_PI*2.); // Upper bound.
    // Translate to index space.
    rl = tr.pointConvertedToNormal(rl);
    ru = tr.pointConvertedToNormal(ru);
    // Field line table for cylindrical grid.
    FieldLineTableCyl<FieldLineTable, FieldLine> cyltbl(Key(rl.r, rl.phi-1, 0.), Key(ru.r, ru.phi+1, 0.), &fm, &tr);

    // Example for cartesian grid
//    AffineTransform tx(Point(.1, .1, .1), Point(0., 0., 0.));
//    Point xl(-12., -12.);
//    Point xu(12., 12.);
//    xl = tx.pointConvertedToNormal(xl);
//    xu = tx.pointConvertedToNormal(xu);
//    FieldLineTableCart<FieldLineTable, FieldLine> carttbl(Key(xl.x, xl.y, 0.), Key(xu.x, xu.y, 0.), &fm, &tx);

    // Initial values.
    Point origin(6, M_PI, 0.);
    double pas[17] = {10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90};
    // Calculation
    for (long jdx=0; jdx<17; jdx++) {
        double pa = pas[jdx] * M_PI/180.;
        ParticleCyl<FieldLineTable, FieldLine> cylp(origin, pa, &cyltbl, &mf, &tr, YES, YES);
        //ParticleCart<FieldLineTable, FieldLine> cartp(origin, pa, &carttbl, &mf, &tx, YES, YES);
        cylp.evaluateLstar();
        //cartp.evaluateLstar();

        fprintf(stdout, "idx = %ld, Initial = [%f, %f, %f], L* = %f\n", jdx, origin.r, origin.phi*180/M_PI, pas[jdx], cylp.Lstar());
    }
}