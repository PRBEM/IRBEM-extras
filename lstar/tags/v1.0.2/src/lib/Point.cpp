//
//  Point.cpp
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

#include "Point.h"
#include <sstream>
#include <string>

namespace UBK {

    std::string const Point::desc () const
    {
        std::ostringstream os;
        os << "[" << this->x << ", " << this->y << ", " << this->z << "]";
        return os.str();
    }

}