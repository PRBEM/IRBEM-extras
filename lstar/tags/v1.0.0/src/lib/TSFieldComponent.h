//
//  TSFieldComponent.h
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

#ifndef UBJDevelopment_TSFieldComponent_h
#define UBJDevelopment_TSFieldComponent_h

#include "Definitions.h"
#include "Point.h"

namespace UBK {

    //!
    //! Intermediate TS field component abstract.
    //!
    class TSFieldComponent {
    public:
        virtual ~TSFieldComponent() {};

        //!
        //! Calculates field vector at a given point in GSW.
        //!
        virtual void getFieldInGSW_atPoint (Point *bOut, Point const ptgsw) const = 0;
    };

}

#endif
