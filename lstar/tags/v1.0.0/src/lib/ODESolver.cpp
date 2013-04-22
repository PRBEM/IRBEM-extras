//
//  ODESolver.cpp
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

#include "ODESolver.h"
#include <cmath>
#include <climits>
#include <cassert>

namespace UBK {
    using namespace std;

    /*! C function version of RK4 solver.
     */
    void rk4IntegrateODE (double delta, double *xInOut, double *yInOut, long count, ODECallback *odeCallbacks)
    {
        assert(odeCallbacks);

        // Inputs
        double hh = delta*.5;
        double h6 = delta/6.;
        // Working space
        double *wspace = new double[4*count];// (double *)malloc(4*count*sizeof(double));
        double *dydx = wspace;
        double *dym = wspace + count;
        double *dyt = wspace + 2*count;
        double *yt = wspace + 3*count;

        BOOL stop = NO;
        long istep = 0;
        // Zeroth
        odeCallbacks->handleEventAtStep_x_y_count_stop(istep, *xInOut, yInOut, count, &stop);// handler(istep, *xInOut, yInOut, count, &stop, handlerContext);
        // Iteration
        for (istep=1; !stop&&istep<LONG_MAX; istep++) {
            double xh = *xInOut + hh;
            odeCallbacks->derivatives(*xInOut, yInOut, count, dydx);// derivs(*xInOut, yInOut, count, dydx, derivsContext);
            for (long i=0; i<count; i++) {
                yt[i] = yInOut[i] + hh*dydx[i];
            }
            odeCallbacks->derivatives(xh, yt, count, dyt);// derivs(xh, yt, count, dyt, derivsContext);
            for (long i=0; i<count; i++) {
                yt[i] = yInOut[i] + hh*dyt[i];
            }
            odeCallbacks->derivatives(xh, yt, count, dym);// derivs(xh, yt, count, dym, derivsContext);
            for (long i=0; i<count; i++) {
                yt[i] = yInOut[i] + delta*dym[i];
                dym[i] += dyt[i];
            }
            odeCallbacks->derivatives((*xInOut+=delta), yt, count, dyt);// derivs((*xInOut+=delta), yt, count, dyt, derivsContext);
            for (long i=0; i<count; i++) {
                yInOut[i] += h6*(dydx[i] + dyt[i] + 2.*dym[i]);
            }
            // Call user block
            odeCallbacks->handleEventAtStep_x_y_count_stop(istep, *xInOut, yInOut, count, &stop);// handler(istep, *xInOut, yInOut, count, &stop, handlerContext);
        }
        
        delete [] wspace;// free(wspace);
    }

}