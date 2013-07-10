//
//  TSExternalField.h
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

#ifndef UBJDevelopment_TSExternalField_h
#define UBJDevelopment_TSExternalField_h

#include "TSFieldComponent.h"
#include <algorithm>

//
// Abstract class for the external field.
//
namespace UBK {

    //
    // Geopack forward declaration.
    //
    class Geopack;

    //!
    //! External field flag.
    //!
    enum _TSExternalFieldModel {
        kTSNone = 0, //!< No external component.
        kTS89Model, //!< T89 model.
        kTS96Model, //!< T96 model.
        kTS02Model, //!< T02 model.
        kTS05Model, //!< TS05 model.
        kTS07Model //!< TS07 model.
    };
    //!
    //! External field flag type.
    //!
    typedef long TSExternalFieldModel;

    //!
    //! External field abstract class.
    //!
    class TSExternalField : public TSFieldComponent {
    private:
        int _iopt;
        double _parmod[10];
        Geopack const* _geopack;

    public:
        //!
        //! @name Property.
        //!
        ///@{
        //@{
        //!
        //! IOPT parameter.
        //!
        int const& iopt () const {return _iopt;};
        void setIopt (int iopt) {_iopt = iopt;};
        //@}

        //@{
        //!
        //! PARMOD parameter.
        //!
        virtual double const* parmod () const {return _parmod;};
        virtual void setParmod (double const parmod[]) {if (parmod) std::copy(parmod, parmod+10, _parmod);};
        //@}

        //!
        //! Pointer to Geopack object.
        //!
        Geopack const* geopack () const {return _geopack;};
        ///@}

    protected:
        TSExternalField (Geopack const* geopack, int iopt, double const parmod[]);
    };
    
}

#endif
