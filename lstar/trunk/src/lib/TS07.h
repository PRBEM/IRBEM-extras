//
//  TS07.h
//  UBJDevelopment
//
//  Created by KYUNGGUK MIN on 7/8/13.
//  Copyright (c) 2013 Kyungguk.com. All rights reserved.
//

//
// $Author$
// $LastChangedDate$
// $Revision$
// $Id$
//

#ifndef UBJDevelopment_TS07_h
#define UBJDevelopment_TS07_h

#include "TSExternalField.h"

namespace UBK {

    class Geopack;

    //!
    //! TS07 field model.
    //!
    //! References: Tsyganenko and Sitnov 2007;
    //!
    class TS07 : public TSExternalField {
    private:
        double _parmod[102];
    public:
        virtual double const* parmod () const {return _parmod;};
        virtual void setParmod (double const parmod[]) {if (parmod) std::copy(parmod, parmod+102, _parmod);};

    public:
        TS07 (Geopack const* geopack, double const parmod[]);

        void getFieldInGSW_atPoint (Point *bOut, Point const ptgsw) const;
    };
    
}

#endif
