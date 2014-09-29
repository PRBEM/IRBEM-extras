//
//  T96.h
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

#ifndef UBJDevelopment_T96_h
#define UBJDevelopment_T96_h

#include "TSExternalField.h"

namespace UBK {

    class Geopack;

    //!
    //! T96 field model.
    //!
    //! References: Tsyganenko 1999;
    //!
    class T96 : public TSExternalField {
    public:
        T96 (Geopack const* geopack, double const parmod[]);

        void getFieldInGSW_atPoint (Point *bOut, Point const ptgsw) const;
    };
    
}

#endif
