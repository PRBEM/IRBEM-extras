//
//  TS05.h
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

#ifndef UBJDevelopment_TS05_h
#define UBJDevelopment_TS05_h

#include "TSExternalField.h"

namespace UBK {

    class Geopack;

    //!
    //! TS05 field model.
    //!
    //! References: Tsyganenko and Sitnov 2005;
    //!
    class TS05 : public TSExternalField {
    public:
        TS05 (Geopack const* geopack, double const parmod[]);

        void getFieldInGSW_atPoint (Point *bOut, Point const ptgsw) const;
    };
    
}

#endif
