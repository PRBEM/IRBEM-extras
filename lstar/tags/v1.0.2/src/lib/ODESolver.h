//
//  ODESolver.h
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

#ifndef UBJDevelopment_ODESolver_h
#define UBJDevelopment_ODESolver_h

#include "Definitions.h"

//! @file ODESolver.h
//! ODE solver and callback definitions.
//! Here fixed step sized RK4 scheme is implemented and used for field line tracing. One can implement a different solver that has the same function signature.
//!

namespace UBK {

    //!
    //! ODE callback abstract.
    //! The callbacks are defined within an abstract class. The ODE solver calls "derivatives" for derivatives and "handleEventAtStep_x_y_count_stop" at each integration step (including zeroth step).
    //! @sa FieldLine
    //!
    class ODECallback {
    public:
        virtual ~ODECallback() {};
        //!
        //! Derivative of a function, i.e., f(x, y) = dy/dx.
        //! Here the count is the array size of y and dydx.
        //!
        virtual void derivatives (double x, double const y[], long count, double *dydx/*!<[out] dy/dx. */) = 0;
        //!
        //! Integration event handler.
        //! Is called at every step of integration.
        //! @param [in] istep Current step index. Starts from 0 which is just initial condition.
        //! @param [out] stop It is important to note that a caller decide when to stop.
        //!
        virtual void handleEventAtStep_x_y_count_stop (long istep, double x, double const y[], long count, BOOL *stop) = 0;
    };

    //!
    //! ODE solver signature.
    //!
    typedef void (*IntegrateODESignature) (double delta, double *xInOut, double *yInOut, long count, ODECallback *odeCallbacks);

    //!
    //! @brief C function version of RK4 solver.
    //! @details Integrates ODE using RK4 scheme. Once it is set, the step size, delta, cannot be adjusted during integeration. Initial values are passed by pointer and replaced by last values on exit. The derivatives provided by callback object provides f(x,y) = dy/dx. At each step, handler member function of the callback object is called, passing intermediate values. It is caller's responsibility to continue or stop by setting necessary boolean to *stop. The integration continues until LONG_MAX reaches.
    //! @note The handler function is called at zeroth step. The 'istep' variable in the handler block will be set to 0 at this step.
    //!
    extern void rk4IntegrateODE (double delta, double *xInOut, double *yInOut, long count, ODECallback *odeCallbacks);

}

#endif
