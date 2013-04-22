//
//  Key.h
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

#ifndef UBJDevelopment_Key_h
#define UBJDevelopment_Key_h

#include "Definitions.h"
#include <string>

namespace UBK {

    //!
    //! Index vector and hash key.
    //!
    typedef union Key {
        long key; //! Unique key for the first two elements.
        struct {int i1, i2, i3;};
        struct {int x, y, z;};
        struct {int r, phi, theta;};

        Key (int i=0, int j=0, int k=0) : x(i), y(j), z(k) {};
        std::string const desc () const;
    } Key;
    
}

#endif
