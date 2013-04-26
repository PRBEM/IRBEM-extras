//
//  TSExternalField.cpp
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

#include "TSExternalField.h"
#include "Geopack.h"
#include <cassert>

namespace UBK {

    TSExternalField::TSExternalField (Geopack const* geopack, int iopt, double const parmod[]) : _geopack(geopack), _iopt(iopt)
    {
        assert(NULL != geopack);

        this->setParmod(parmod);
    }

}