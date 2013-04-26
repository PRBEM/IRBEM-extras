//
//  T02.h
//  UBJDevelopment
//
//  Created by Kyungguk Min on 4/6/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#ifndef UBJDevelopment_T02_h
#define UBJDevelopment_T02_h

#include "TSExternalField.h"

namespace UBK {

    class Geopack;

    //!
    //! T02 field model.
    //!
    //! References: Tsyganenko 2002a and b;
    //!
    class T02 : public TSExternalField {
    public:
        T02 (Geopack const* geopack, double const parmod[]);

        void getFieldInGSW_atPoint (Point *bOut, Point const ptgsw) const;
    };
    
}

#endif
