//
//  Point.h
//  UBJDevelopment
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

#ifndef UBJDevelopment_Point_h
#define UBJDevelopment_Point_h

#include "Definitions.h"
#include <string>

namespace UBK {

    //!
    //! 3 element vector representation.
    //!
    union Point {
        double xyz[3];
        struct {double x, y, z;};
        struct {double r, phi, theta;};

        Point (double x1=0., double x2=0., double x3=0.) : x(x1), y(x2), z(x3) {};
        std::string const desc () const;
    };

}

#endif
