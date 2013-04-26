//
//  T89.h
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

#ifndef UBJDevelopment_T89_h
#define UBJDevelopment_T89_h

#include "TSExternalField.h"

namespace UBK {

    class Geopack;

    //!
    //! T89 field model.
    //!
    //! References: Tsyganenko 1987, 1989;
    //!
    class T89 : public TSExternalField {
    public:
        T89 (Geopack const* geopack, int iopt);

        void getFieldInGSW_atPoint (Point *bOut, Point const ptgsw) const;
    };

}

#endif
