//
//  Key.cpp
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

#include "Key.h"
#include <sstream>

namespace UBK {

    std::string const Key::desc () const
    {
        std::ostringstream os;
        os << "[" << this->x << ", " << this->y << ", " << this->z << "]";
        return os.str();
    }

}